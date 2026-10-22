from dataclasses import FrozenInstanceError
from copy import copy
from pathlib import Path
import re
import shlex
import shutil
import subprocess

import pandas as pd
import pytest
from pyomo import environ as pe
from pyomo.core.expr.visitor import identify_variables
from pyomo.repn.standard_repn import generate_standard_repn

from hatchet import config
from hatchet.bin.HATCHet import (
    execute,
    load_fixed_matrix_entries,
    makeBaseCMD,
    parsing_arguments,
    runningDiploid,
    runningTetraploid,
)
from hatchet.utils.commands import commands
from hatchet.utils.fixed_matrix_editor import main as fixed_matrix_editor
from hatchet.utils.fixed_matrix import (
    FixedMatrixFormatError,
    parse_fixed_matrix,
    resolve_fixed_matrix,
    validate_fixed_matrix,
    verify_fixed_matrix_solution,
)
from hatchet.utils.solve.ilp_subset import ILPSubset, ILPSubsetSplit
from hatchet.utils.solve import _finalize_solution, solver_available
from hatchet.utils.solve.utils import canonicalize_clone_order


FORMAT_MARKER = '# HATCHET_FIXED_ENTRIES 1'
DATA_FOLDER = Path(__file__).parent / 'data'
BINNING_SCENARIOS = (
    pytest.param(
        DATA_FOLDER / 'fw' / 'bulk.seg',
        DATA_FOLDER / 'fw' / 'fixed-matrix.tsv',
        id='fixed-width',
    ),
    pytest.param(
        DATA_FOLDER / 'vw' / 'bbc' / 'bulk.seg',
        DATA_FOLDER / 'vw' / 'fixed-matrix.tsv',
        id='variable-width',
    ),
)


def _write_fixed_matrix(tmp_path, contents):
    path = tmp_path / 'fixed.tsv'
    path.write_text(contents)
    return path


def test_parse_fixed_matrix_preserves_sparse_entries_and_source_lines(tmp_path):
    path = _write_fixed_matrix(
        tmp_path,
        f"""{FORMAT_MARKER}
# An optional comment.
# clones 3
matrix\trow\tcolumn\tvalue
U\tclone2\tTumor sample\t0.25
A\tcluster-01\tclone1\t2
B\tcluster-01\tnormal\t1
""",
    )

    fixed = parse_fixed_matrix(path)

    assert fixed.clones == 3
    assert fixed.source == str(path)
    assert tuple(entry.coordinate for entry in fixed.entries) == (
        ('U', 'clone2', 'Tumor sample'),
        ('A', 'cluster-01', 'clone1'),
        ('B', 'cluster-01', 'normal'),
    )
    assert tuple(entry.value for entry in fixed.entries) == (0.25, 2, 1)
    assert tuple(entry.line_number for entry in fixed.entries) == (5, 6, 7)
    assert fixed.entries_for('A') == (fixed.entries[1],)

    with pytest.raises(FrozenInstanceError):
        fixed.clones = 4
    with pytest.raises(FrozenInstanceError):
        fixed.entries[0].value = 0.5


def test_parse_fixed_matrix_accepts_utf8_bom_blank_lines_and_arbitrary_row_order(tmp_path):
    path = tmp_path / 'fixed-bom.tsv'
    path.write_bytes(
        (
            f'{FORMAT_MARKER}\n'
            '# clones 3\n'
            '\n'
            '# Entries intentionally appear in B, U, A order.\n'
            'matrix\trow\tcolumn\tvalue\n'
            'B\tcluster1\tclone2\t1\n'
            '\n'
            'U\tclone1\tsample1\t0.3\n'
            'A\tcluster1\tclone1\t2\n'
        ).encode('utf-8-sig')
    )

    fixed = parse_fixed_matrix(path)

    assert tuple(entry.matrix for entry in fixed.entries) == ('B', 'U', 'A')
    assert tuple(entry.value for entry in fixed.entries) == (1, 0.3, 2)


def test_parse_fixed_matrix_rejects_invalid_utf8(tmp_path):
    path = tmp_path / 'invalid-encoding.tsv'
    path.write_bytes(b'\xff\xfe\x00')

    with pytest.raises(FixedMatrixFormatError, match='must use UTF-8 encoding'):
        parse_fixed_matrix(path)


@pytest.mark.parametrize(
    ('contents', 'message'),
    [
        ('', 'expected format marker'),
        ('# HATCHET_FIXED_ENTRIES 2\n', 'expected format marker'),
        (f'{FORMAT_MARKER}\nmatrix\trow\tcolumn\tvalue\n', 'clone-count metadata must appear before'),
        (f'{FORMAT_MARKER}\n# clones 1\n', 'clone count must be at least 2'),
        (f'{FORMAT_MARKER}\n# clones 2\n# clones 2\n', 'duplicate clone-count metadata'),
        (f'{FORMAT_MARKER}\n# clones 2\nwrong\theader\n', 'expected tab-separated header'),
        (f'{FORMAT_MARKER}\n# clones 2\nmatrix\trow\tcolumn\tvalue\n', 'at least one fixed entry'),
        (
            f'{FORMAT_MARKER}\n# clones 2\nmatrix\trow\tcolumn\tvalue\nA\tcluster1\tclone1\n',
            'expected 4 tab-separated fields',
        ),
        (
            f'{FORMAT_MARKER}\n# clones 2\nmatrix\trow\tcolumn\tvalue\nC\tcluster1\tclone1\t2\n',
            'unknown matrix',
        ),
        (
            f'{FORMAT_MARKER}\n# clones 2\nmatrix\trow\tcolumn\tvalue\nA\t\tclone1\t2\n',
            'row identifier must not be empty',
        ),
        (
            f'{FORMAT_MARKER}\n# clones 2\nmatrix\trow\tcolumn\tvalue\nA\tcluster1\tclone0\t2\n',
            'invalid clone label',
        ),
        (
            f'{FORMAT_MARKER}\n# clones 2\nmatrix\trow\tcolumn\tvalue\nA\tcluster1\tclone2\t2\n',
            'outside the declared clone count',
        ),
        (
            f'{FORMAT_MARKER}\n# clones 2\nmatrix\trow\tcolumn\tvalue\nA\tcluster1\tclone1\t1.5\n',
            'must be an integer',
        ),
        (
            f'{FORMAT_MARKER}\n# clones 2\nmatrix\trow\tcolumn\tvalue\nB\tcluster1\tclone1\t-1\n',
            'must be nonnegative',
        ),
        (
            f'{FORMAT_MARKER}\n# clones 2\nmatrix\trow\tcolumn\tvalue\nA\tcluster1\tnormal\t2\n',
            'normal clone must equal 1',
        ),
        (
            f'{FORMAT_MARKER}\n# clones 2\nmatrix\trow\tcolumn\tvalue\nU\tclone1\tsample1\tnan\n',
            'must be finite',
        ),
        (
            f'{FORMAT_MARKER}\n# clones 2\nmatrix\trow\tcolumn\tvalue\nU\tclone1\tsample1\t1.1\n',
            'must be in [0, 1]',
        ),
        (
            f'{FORMAT_MARKER}\n# clones 2\nmatrix\trow\tcolumn\tvalue\n'
            'A\tcluster1\tclone1\t2\nA\tcluster1\tclone1\t2\n',
            'duplicate fixed entry',
        ),
    ],
)
def test_parse_fixed_matrix_rejects_invalid_files(tmp_path, contents, message):
    path = _write_fixed_matrix(tmp_path, contents)

    with pytest.raises(FixedMatrixFormatError, match=re.escape(message)) as exc_info:
        parse_fixed_matrix(path)

    assert str(path) in str(exc_info.value)


def _parse_and_resolve(tmp_path, rows, clones=3, cluster_ids=('cluster1', 'cluster2'), sample_ids=('sample1',)):
    path = _write_fixed_matrix(
        tmp_path,
        f'{FORMAT_MARKER}\n# clones {clones}\nmatrix\trow\tcolumn\tvalue\n{rows}',
    )
    return resolve_fixed_matrix(parse_fixed_matrix(path), cluster_ids, sample_ids)


def test_resolve_fixed_matrix_uses_exact_textual_seg_identifiers(tmp_path):
    fixed = _parse_and_resolve(
        tmp_path,
        'A\t1\tclone1\t2\nU\tclone2\tTumor sample\t0.25\n',
        cluster_ids=(1, 2),
        sample_ids=('Tumor sample',),
    )

    assert tuple((entry.matrix, entry.row_index, entry.column_index) for entry in fixed.entries) == (
        ('A', 0, 1),
        ('U', 2, 0),
    )
    with pytest.raises(FrozenInstanceError):
        fixed.entries[0].row_index = 1


@pytest.mark.parametrize(
    ('rows', 'message'),
    [
        ('A\tmissing\tclone1\t2\n', 'unknown SEG cluster identifier'),
        ('U\tclone1\tmissing\t0.5\n', 'unknown SEG sample identifier'),
    ],
)
def test_resolve_fixed_matrix_rejects_unknown_seg_identifiers(tmp_path, rows, message):
    path = _write_fixed_matrix(
        tmp_path,
        f'{FORMAT_MARKER}\n# clones 3\nmatrix\trow\tcolumn\tvalue\n{rows}',
    )

    with pytest.raises(FixedMatrixFormatError, match=message):
        resolve_fixed_matrix(parse_fixed_matrix(path), ('cluster1',), ('sample1',))


def test_validate_fixed_matrix_accepts_compatible_partial_entries(tmp_path):
    _parse_and_resolve(
        tmp_path,
        'A\tcluster1\tclone1\t2\n' 'B\tcluster1\tclone1\t1\n' 'U\tnormal\tsample1\t0.2\n' 'U\tclone1\tsample1\t0.3\n',
    )


@pytest.mark.parametrize('fixed_proportion', (0.0, 0.029, 0.03))
def test_validate_fixed_matrix_accepts_minprop_boundaries(tmp_path, fixed_proportion):
    fixed = _parse_and_resolve(tmp_path, f'U\tclone1\tsample1\t{fixed_proportion}\n')

    assert validate_fixed_matrix(fixed, cn_max=6, mu=0.03) is fixed

    assert (
        validate_fixed_matrix(
            fixed,
            cn_max=6,
            mu=0.03,
            copy_numbers={'cluster2': (1, 1)},
            purities=None,
            ampdel=False,
            weights={'cluster1': 1, 'cluster2': 1},
        )
        is fixed
    )


@pytest.mark.parametrize(
    ('rows', 'kwargs', 'message'),
    [
        ('A\tcluster1\tclone1\t7\n', {}, 'exceeds copy-number bound'),
        (
            'A\tcluster1\tclone1\t4\nB\tcluster1\tclone1\t3\n',
            {},
            'fixed A+B total',
        ),
        (
            'A\tcluster1\tclone1\t2\n',
            {'copy_numbers': {'cluster1': (1, 1)}},
            'conflicts with clonal copy number',
        ),
        ('U\tclone1\tsample1\t0.01\n', {}, 'between 0 and minprop'),
        (
            'U\tnormal\tsample1\t0.7\nU\tclone1\tsample1\t0.4\n',
            {},
            'sum to 1.1, which exceeds 1',
        ),
        (
            'U\tnormal\tsample1\t0.4\nU\tclone1\tsample1\t0.5\n',
            {'clones': 2},
            'fully fixed U entries',
        ),
        (
            'U\tnormal\tsample1\t0.99\nU\tclone2\tsample1\t0\n',
            {},
            'remaining U proportion',
        ),
        (
            'U\tnormal\tsample1\t0.4\n',
            {'purities': [0.5]},
            'conflicts with purity',
        ),
        (
            'A\tcluster1\tclone1\t1\nA\tcluster1\tclone2\t3\n',
            {'ampdel': True, 'base': 2},
            'force both amplification and deletion',
        ),
        (
            'A\tcluster1\tclone1\t0\nB\tcluster1\tclone1\t0\n',
            {'weights': {'cluster1': 1, 'cluster2': 1}},
            'cannot have fixed tumor copy number',
        ),
        (
            'A\tcluster1\tclone1\t1\nB\tcluster1\tclone1\t0\n' 'A\tcluster1\tclone2\t2\nB\tcluster1\tclone2\t0\n',
            {'cnstates': 1},
            'exceeding cnstates',
        ),
    ],
)
def test_validate_fixed_matrix_rejects_local_solver_conflicts(tmp_path, rows, kwargs, message):
    validation_kwargs = dict(kwargs)
    clones = validation_kwargs.pop('clones', 3)
    fixed = _parse_and_resolve(tmp_path, rows, clones=clones)
    options = {
        'cn_max': 6,
        'mu': 0.03,
        'copy_numbers': {},
        'purities': None,
        'ampdel': False,
        'base': 1,
        'cnstates': None,
        'weights': None,
    }
    options.update(validation_kwargs)

    with pytest.raises(FixedMatrixFormatError, match=re.escape(message)) as exc_info:
        validate_fixed_matrix(fixed, **options)

    assert str(fixed.source) in str(exc_info.value)


def test_validate_fixed_matrix_accepts_any_fully_fixed_clone_order(tmp_path):
    # The solver no longer imposes a clone order, so a heavier clone1 than clone2 is valid
    fixed = _parse_and_resolve(
        tmp_path,
        'A\tcluster1\tclone1\t2\nB\tcluster1\tclone1\t1\n'
        'A\tcluster2\tclone1\t2\nB\tcluster2\tclone1\t1\n'
        'A\tcluster1\tclone2\t1\nB\tcluster1\tclone2\t1\n'
        'A\tcluster2\tclone2\t1\nB\tcluster2\tclone2\t1\n',
        clones=3,
    )

    assert (
        validate_fixed_matrix(
            fixed,
            cn_max=6,
            mu=0.03,
            copy_numbers={},
            purities=None,
            ampdel=False,
            base=1,
            cnstates=None,
            weights=None,
        )
        is fixed
    )
    assert fixed.addresses_tumor_clones


def _compute_cn_input(tmp_path):
    prefix = tmp_path / 'bulk'
    prefix.with_suffix('.seg').touch()
    prefix.with_suffix('.bbc').touch()
    return str(prefix)


def test_compute_cn_fixed_matrix_argument_requires_exact_clone_count(tmp_path, monkeypatch):
    monkeypatch.setattr(config.compute_cn, 'solver', 'gurobi')
    fixed_path = _write_fixed_matrix(
        tmp_path,
        f'{FORMAT_MARKER}\n# clones 3\nmatrix\trow\tcolumn\tvalue\nA\tcluster1\tclone1\t2\n',
    )
    input_prefix = _compute_cn_input(tmp_path)

    parsed = parsing_arguments(['-i', input_prefix, '-x', str(tmp_path), '-n', '3', '--fixed-matrix', str(fixed_path)])
    assert parsed['fixedmatrix'] == str(fixed_path)
    assert parsed['ln'] == parsed['un'] == 3

    with pytest.raises(ValueError, match='requires one exact number of clones'):
        parsing_arguments(['-i', input_prefix, '-x', str(tmp_path), '-n', '2,3', '--fixed-matrix', str(fixed_path)])


def test_compute_cn_fixed_matrix_uses_ini_default_and_checks_path(tmp_path, monkeypatch):
    monkeypatch.setattr(config.compute_cn, 'solver', 'gurobi')
    input_prefix = _compute_cn_input(tmp_path)
    missing_path = tmp_path / 'missing.tsv'
    monkeypatch.setattr(config.compute_cn, 'fixedmatrix', str(missing_path))

    with pytest.raises(ValueError, match='Fixed-matrix file .* not found'):
        parsing_arguments(['-i', input_prefix, '-x', str(tmp_path), '-n', '3'])


def test_compute_cn_fixed_matrix_cli_ini_and_environment_precedence(tmp_path, monkeypatch):
    monkeypatch.setattr(config.compute_cn, 'solver', 'gurobi')
    input_prefix = _compute_cn_input(tmp_path)
    ini_path = _write_fixed_matrix(
        tmp_path,
        f'{FORMAT_MARKER}\n# clones 3\nmatrix\trow\tcolumn\tvalue\nA\tini\tclone1\t2\n',
    )
    environment_path = tmp_path / 'environment.tsv'
    environment_path.write_text(f'{FORMAT_MARKER}\n# clones 3\nmatrix\trow\tcolumn\tvalue\nA\tenvironment\tclone1\t2\n')
    cli_path = tmp_path / 'cli.tsv'
    cli_path.write_text(f'{FORMAT_MARKER}\n# clones 3\nmatrix\trow\tcolumn\tvalue\nA\tcli\tclone1\t2\n')
    monkeypatch.setattr(config.compute_cn, 'fixedmatrix', str(ini_path))

    parsed_ini = parsing_arguments(['-i', input_prefix, '-x', str(tmp_path), '-n', '3'])
    assert parsed_ini['fixedmatrix'] == str(ini_path)

    monkeypatch.setenv('HATCHET_COMPUTE_CN_FIXEDMATRIX', str(environment_path))
    parsed_environment = parsing_arguments(['-i', input_prefix, '-x', str(tmp_path), '-n', '3'])
    assert parsed_environment['fixedmatrix'] == str(environment_path)

    parsed_cli = parsing_arguments(
        ['-i', input_prefix, '-x', str(tmp_path), '-n', '3', '--fixed-matrix', str(cli_path)]
    )
    assert parsed_cli['fixedmatrix'] == str(cli_path)


def test_compute_cn_max_homdel_size_cli_config_and_environment(tmp_path, monkeypatch):
    monkeypatch.setattr(config.compute_cn, 'solver', 'gurobi')
    input_prefix = _compute_cn_input(tmp_path)
    monkeypatch.setattr(config.compute_cn, 'maxhomdelsize', 0.01)

    parsed_config = parsing_arguments(['-i', input_prefix, '-x', str(tmp_path), '-n', '3'])
    assert parsed_config['max_homdel_size'] == pytest.approx(0.01)

    monkeypatch.setenv('HATCHET_COMPUTE_CN_MAXHOMDELSIZE', '0.02')
    parsed_environment = parsing_arguments(['-i', input_prefix, '-x', str(tmp_path), '-n', '3'])
    assert parsed_environment['max_homdel_size'] == pytest.approx(0.02)

    parsed_cli = parsing_arguments(['-i', input_prefix, '-x', str(tmp_path), '-n', '3', '--max-homdel-size', '0.03'])
    assert parsed_cli['max_homdel_size'] == pytest.approx(0.03)

    for invalid in ('-0.01', '1.01', 'nan', 'inf'):
        with pytest.raises(ValueError, match='homozygous-deletion cluster-size fraction'):
            parsing_arguments(['-i', input_prefix, '-x', str(tmp_path), '-n', '3', '--max-homdel-size', invalid])


def test_load_fixed_matrix_entries_checks_metadata_and_seg_identifiers(tmp_path):
    fixed_path = _write_fixed_matrix(
        tmp_path,
        f'{FORMAT_MARKER}\n# clones 3\nmatrix\trow\tcolumn\tvalue\nA\tcluster1\tclone1\t2\n',
    )
    args = {'fixedmatrix': str(fixed_path), 'ln': 3}

    fixed = load_fixed_matrix_entries(args, ('cluster1',), ('sample1',))
    assert fixed.clones == 3

    args['ln'] = 2
    with pytest.raises(ValueError, match='declares 3 clones'):
        load_fixed_matrix_entries(args, ('cluster1',), ('sample1',))

    args['ln'] = 3
    with pytest.raises(FixedMatrixFormatError, match='unknown SEG cluster identifier'):
        load_fixed_matrix_entries(args, ('other-cluster',), ('sample1',))


def test_load_fixed_matrix_entries_preserves_empty_option():
    assert load_fixed_matrix_entries({'fixedmatrix': None}, (), ()) is None


def test_cpp_command_forwards_fixed_matrix_path():
    command = makeBaseCMD(
        {
            'solver': 'solve',
            'input': 'bulk',
            'ampdel': False,
            'd': None,
            'j': None,
            'p': None,
            'u': None,
            'm': None,
            's': None,
            'f': None,
            'r': None,
            'M': None,
            'fixedmatrix': '/tmp/fixed entries.tsv',
            'max_homdel_size': 0.02,
        },
        6,
    )
    tokens = shlex.split(command)
    fixed_matrix_index = tokens.index('-F')
    assert tokens.count('-F') == 1
    assert tokens[fixed_matrix_index + 1] == '/tmp/fixed entries.tsv'
    homdel_index = tokens.index('--max-homdel-size')
    assert tokens[homdel_index + 1] == '0.02'


def test_cpp_execute_parses_final_objective_at_quiet_verbosity(monkeypatch, tmp_path):
    class CompletedProcess:
        returncode = 0

        def communicate(self):
            return None, 'Final objective:\t12.5\n'

    monkeypatch.setattr(config.compute_cn, 'solver', 'cpp')
    monkeypatch.setattr('hatchet.bin.HATCHet.subprocess.Popen', lambda *args, **kwargs: CompletedProcess())

    objective = execute(
        {'p': 1, 'v': 0, 'fixed': object()},
        'solve bulk',
        2,
        str(tmp_path / 'result'),
    )

    assert objective == 12.5


@pytest.mark.parametrize(('seg_path', 'fixed_path'), BINNING_SCENARIOS)
def test_fixed_matrix_fixture_matches_binning_scenario(seg_path, fixed_path):
    segments = pd.read_csv(seg_path, sep='\t')
    cluster_ids = tuple(segments['#ID'].drop_duplicates())
    sample_ids = tuple(segments['SAMPLE'].drop_duplicates())
    weights = segments.groupby('#ID', sort=False)['#BINS'].first()

    fixed = parse_fixed_matrix(fixed_path)
    resolved = resolve_fixed_matrix(fixed, cluster_ids, sample_ids)

    assert (
        validate_fixed_matrix(
            resolved,
            cn_max=6,
            mu=0.03,
            ampdel=False,
            weights=weights,
        )
        is resolved
    )
    assert resolved.entries


@pytest.fixture(scope='module')
def cpp_fixed_matrix_reader(tmp_path_factory):
    compiler = shutil.which('g++')
    if compiler is None:
        pytest.skip('g++ is required for the native fixed-matrix parser test')
    repository = Path(__file__).parents[1]
    executable = tmp_path_factory.mktemp('fixed-matrix-cpp') / 'reader'
    subprocess.run(
        [
            compiler,
            '-std=c++11',
            '-I' + str(repository / 'src'),
            str(repository / 'tests' / 'fixed_matrix_cpp_test.cpp'),
            str(repository / 'src' / 'fixed_matrix.cpp'),
            str(repository / 'src' / 'sample_clonal.cpp'),
            str(repository / 'src' / 'utils.cpp'),
            '-pthread',
            '-o',
            str(executable),
        ],
        check=True,
    )
    return executable


def test_cpp_fixed_matrix_parser_and_operations(tmp_path, cpp_fixed_matrix_reader):
    fixed_path = _write_fixed_matrix(
        tmp_path,
        f'{FORMAT_MARKER}\n# clones 3\nmatrix\trow\tcolumn\tvalue\n'
        'A\tcluster1\tclone1\t2\nB\tcluster2\tclone2\t1\n'
        'U\tnormal\tsample1\t0.2\nU\tclone1\tsample1\t0.3\n',
    )
    subprocess.run([str(cpp_fixed_matrix_reader), str(fixed_path)], check=True)


def test_cpp_canonical_clone_order_sorts_tumor_clones_stably(cpp_fixed_matrix_reader):
    subprocess.run([str(cpp_fixed_matrix_reader), '--canonicalize'], check=True)


def test_cpp_fixed_matrix_parser_accepts_utf8_bom_blank_lines_and_arbitrary_order(
    tmp_path,
    cpp_fixed_matrix_reader,
):
    fixed_path = tmp_path / 'native-bom.tsv'
    fixed_path.write_bytes(
        (
            f'{FORMAT_MARKER}\n'
            '# clones 3\n'
            '\n'
            'matrix\trow\tcolumn\tvalue\n'
            'U\tclone1\tsample1\t0.3\n'
            'B\tcluster2\tclone2\t1\n'
            '\n'
            'A\tcluster1\tclone1\t2\n'
            'U\tnormal\tsample1\t0.2\n'
        ).encode('utf-8-sig')
    )

    subprocess.run([str(cpp_fixed_matrix_reader), str(fixed_path)], check=True)


def test_cpp_fixed_matrix_parser_rejects_invalid_clone_labels(tmp_path, cpp_fixed_matrix_reader):
    fixed_path = _write_fixed_matrix(
        tmp_path,
        f'{FORMAT_MARKER}\n# clones 3\nmatrix\trow\tcolumn\tvalue\nA\tcluster1\tclone+1\t2\n',
    )
    result = subprocess.run(
        [str(cpp_fixed_matrix_reader), str(fixed_path)],
        text=True,
        capture_output=True,
    )
    assert result.returncode == 1
    assert 'invalid clone label' in result.stderr


@pytest.mark.parametrize(
    ('rows', 'scenario', 'message'),
    [
        ('A\tcluster1\tclone1\t7\n', 'bound', 'exceeds the copy-number bound'),
        ('A\tcluster1\tclone1\t2\n', 'clonal', 'conflicts with a clonal copy number'),
        ('U\tclone1\tsample1\t0.01\n', 'minprop', 'between 0 and the minimum clone proportion'),
        (
            'U\tnormal\tsample1\t0.7\nU\tclone1\tsample1\t0.4\n',
            'mixture',
            'sum to more than 1',
        ),
        (
            'A\tcluster1\tclone1\t1\nA\tcluster1\tclone2\t3\n',
            'ampdel',
            'force both amplification and deletion',
        ),
        (
            'A\tcluster1\tclone1\t1\nB\tcluster1\tclone1\t0\n' 'A\tcluster1\tclone2\t2\nB\tcluster1\tclone2\t0\n',
            'cnstates',
            'exceed the maximum number of copy-number states',
        ),
        (
            'A\tcluster1\tclone1\t0\nB\tcluster1\tclone1\t0\n',
            'large-zero',
            'large cluster cannot have fixed tumor copy number',
        ),
    ],
)
def test_cpp_fixed_matrix_validation_matches_documented_constraints(
    tmp_path,
    cpp_fixed_matrix_reader,
    rows,
    scenario,
    message,
):
    fixed_path = _write_fixed_matrix(
        tmp_path,
        f'{FORMAT_MARKER}\n# clones 3\nmatrix\trow\tcolumn\tvalue\n{rows}',
    )

    result = subprocess.run(
        [str(cpp_fixed_matrix_reader), str(fixed_path), '3', scenario],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 1
    assert message in result.stderr


def test_cpp_fixed_matrix_validation_skips_ampdel_exempt_clusters(tmp_path, cpp_fixed_matrix_reader):
    fixed_path = _write_fixed_matrix(
        tmp_path,
        f'{FORMAT_MARKER}\n# clones 3\nmatrix\trow\tcolumn\tvalue\n' 'A\tcluster1\tclone1\t1\nA\tcluster1\tclone2\t3\n',
    )

    subprocess.run([str(cpp_fixed_matrix_reader), str(fixed_path), '3', 'ampdel-exempt'], check=True)


def test_cpp_fixed_matrix_homdel_validation_uses_custom_inclusive_cutoff(
    tmp_path,
    cpp_fixed_matrix_reader,
):
    fixed_path = _write_fixed_matrix(
        tmp_path,
        f'{FORMAT_MARKER}\n# clones 3\nmatrix\trow\tcolumn\tvalue\n' 'A\tcluster1\tclone1\t0\nB\tcluster1\tclone1\t0\n',
    )

    subprocess.run([str(cpp_fixed_matrix_reader), str(fixed_path), '3', 'large-zero', '0.5'], check=True)
    result = subprocess.run(
        [str(cpp_fixed_matrix_reader), str(fixed_path), '3', 'large-zero', str(1 / 3)],
        text=True,
        capture_output=True,
    )
    assert result.returncode == 1
    assert 'cannot have fixed tumor copy number' in result.stderr


def test_cpp_fixed_matrix_homdel_validation_ignores_zero_weight_clusters(
    tmp_path,
    cpp_fixed_matrix_reader,
):
    fixed_path = _write_fixed_matrix(
        tmp_path,
        f'{FORMAT_MARKER}\n# clones 3\nmatrix\trow\tcolumn\tvalue\n' 'A\tcluster1\tclone1\t0\nB\tcluster1\tclone1\t0\n',
    )

    subprocess.run(
        [str(cpp_fixed_matrix_reader), str(fixed_path), '3', 'zero-weight', '0'],
        check=True,
    )


@pytest.mark.parametrize(('_seg_path', 'fixed_path'), BINNING_SCENARIOS)
def test_cpp_fixed_matrix_parser_accepts_binning_scenario_fixture(
    _seg_path,
    fixed_path,
    cpp_fixed_matrix_reader,
):
    subprocess.run([str(cpp_fixed_matrix_reader), str(fixed_path), '2'], check=True)


def test_fixed_matrix_editor_is_registered_and_packaged(monkeypatch, capsys):
    opened = []
    monkeypatch.setattr(
        'hatchet.utils.fixed_matrix_editor.webbrowser.open', lambda *args, **kwargs: opened.append(args)
    )

    editor_path = fixed_matrix_editor(['--no-open'])

    assert 'fixed-matrix-editor' in commands
    assert editor_path.is_file()
    assert editor_path.name == 'fixed_matrix_editor.html'
    assert capsys.readouterr().out.strip() == str(editor_path)
    assert opened == []


def test_fixed_matrix_editor_opens_packaged_file(monkeypatch):
    opened = []
    monkeypatch.setattr(
        'hatchet.utils.fixed_matrix_editor.webbrowser.open', lambda *args, **kwargs: opened.append(args)
    )

    editor_path = fixed_matrix_editor([])

    assert opened == [(editor_path.as_uri(),)]


def test_fixed_matrix_editor_asset_is_self_contained_and_feature_complete():
    editor_path = fixed_matrix_editor(['--no-open'])
    source = editor_path.read_text()

    assert 'id="seg-file"' in source
    assert 'id="constraints-file"' in source
    assert 'data-matrix="A"' in source
    assert 'data-matrix="B"' in source
    assert 'data-matrix="U"' in source
    assert 'function parseSeg' in source
    assert 'function parseConstraintFile' in source
    assert 'function canonicalTsv' in source
    assert "elements.grid.addEventListener('paste'" in source
    assert 'https://' not in source
    assert 'http://' not in source


def test_compute_cn_documentation_covers_fixed_matrix_contract():
    documentation = (Path(__file__).parents[1] / 'docs' / 'source' / 'doc_compute_cn.md').read_text()

    assert '# HATCHET_FIXED_ENTRIES 1' in documentation
    assert 'hatchet fixed-matrix-editor' in documentation
    assert '--fixed-matrix' in documentation
    assert 'fixedmatrix = /path/to/patient.fixed-matrix.tsv' in documentation
    assert 'exact total number of clones' in documentation
    assert 'HATCHET_COMPUTE_CN_FIXEDMATRIX' in documentation
    assert 'same file is applied to both model-selection branches' in documentation


def test_compute_cn_documentation_covers_max_homdel_size_contract():
    documentation = (Path(__file__).parents[1] / 'docs' / 'source' / 'doc_compute_cn.md').read_text()

    assert '--max-homdel-size' in documentation
    assert 'maxhomdelsize = 0.005' in documentation
    assert 'HATCHET_COMPUTE_CN_MAXHOMDELSIZE' in documentation
    assert 'cutoff is inclusive' in documentation
    assert 'fraction in `[0, 1]`' in documentation


def _build_ilp(n=4, fixed_matrix=None, weights=None, max_homdel_size=0.005):
    cluster_ids = pd.Index(['cluster1', 'cluster2'], name='#ID')
    sample_ids = ['sample1']
    f_a = pd.DataFrame([[1.0], [2.0]], index=cluster_ids, columns=sample_ids)
    f_b = pd.DataFrame([[1.0], [1.0]], index=cluster_ids, columns=sample_ids)
    if weights is None:
        weights = pd.Series([1.0, 2.0], index=cluster_ids)

    ilp = ILPSubset(
        n=n,
        cn_max=6,
        d=-1,
        mu=0.03,
        ampdel=False,
        copy_numbers={},
        f_a=f_a,
        f_b=f_b,
        w=weights,
        purities=None,
        fixed_matrix=fixed_matrix,
        max_homdel_size=max_homdel_size,
    )
    ilp.create_model()
    return ilp


def _fixed_entry_values(model):
    """Copy-number and proportion variables pinned by bounds, excluding the normal clone's copy numbers."""
    assert len(model.fixed_value_conflicts) == 0
    fixed = {}
    for variable in model.component_data_objects(pe.Var):
        prefix, *indices = variable.name.split('_')
        if prefix not in ('cA', 'cB', 'u') or (prefix != 'u' and indices[1] == '1'):
            continue
        if variable.lb is not None and variable.lb == variable.ub:
            fixed[variable.name] = float(variable.lb)
    return fixed


def test_python_models_enforce_and_preserve_sparse_fixed_entries(tmp_path):
    fixed = _parse_and_resolve(
        tmp_path,
        'A\tcluster1\tclone1\t2\nB\tcluster2\tclone2\t1\nU\tclone1\tsample1\t0.4\n',
    )
    ilp = _build_ilp(n=3, fixed_matrix=fixed)

    assert _fixed_entry_values(ilp.model) == {'cA_1_2': 2.0, 'cB_2_3': 1.0, 'u_2_1': 0.4}

    carch = copy(ilp)
    carch.fix_u([[0.2], [0.4], [0.4]])
    carch.create_model()
    assert carch.fixed_matrix is fixed
    assert carch._fixed_u[1][0] == 0.4
    assert _fixed_entry_values(carch.model) == {'cA_1_2': 2.0, 'cB_2_3': 1.0}

    uarch = copy(ilp)
    uarch.fix_c([[1, 0, 0], [1, 0, 0]], [[1, 0, 0], [1, 0, 0]])
    uarch.create_model()
    assert uarch.fixed_matrix is fixed
    assert uarch._fixed_cA[0][1] == 2
    assert uarch._fixed_cB[1][2] == 1
    assert _fixed_entry_values(uarch.model) == {'u_2_1': 0.4}


def _build_split_ilp(fixed_matrix, weights=None, max_homdel_size=0.005):
    cluster_ids = pd.Index(['cluster1', 'cluster2'], name='#ID')
    sample_ids = ['sample1']
    f_a = pd.DataFrame([[1.0], [2.0]], index=cluster_ids, columns=sample_ids)
    f_b = pd.DataFrame([[1.0], [1.0]], index=cluster_ids, columns=sample_ids)
    bins_a = {
        cluster_id: pd.DataFrame([[f_a.loc[cluster_id, 'sample1']]], columns=sample_ids) for cluster_id in cluster_ids
    }
    bins_b = {
        cluster_id: pd.DataFrame([[f_b.loc[cluster_id, 'sample1']]], columns=sample_ids) for cluster_id in cluster_ids
    }
    lengths = {cluster_id: [1.0] for cluster_id in cluster_ids}
    ilp = ILPSubsetSplit(
        n=3,
        cn_max=6,
        d=-1,
        mu=0.03,
        ampdel=False,
        copy_numbers={},
        f_a=f_a,
        f_b=f_b,
        binsA=bins_a,
        binsB=bins_b,
        lengths=lengths,
        fixed_matrix=fixed_matrix,
        w=weights,
        max_homdel_size=max_homdel_size,
    )
    ilp.create_model()
    return ilp


def test_binwise_models_enforce_and_preserve_sparse_fixed_entries(tmp_path):
    fixed = _parse_and_resolve(
        tmp_path,
        'A\tcluster1\tclone1\t2\nB\tcluster2\tclone2\t1\nU\tclone1\tsample1\t0.4\n',
    )
    ilp = _build_split_ilp(fixed)
    assert _fixed_entry_values(ilp.model) == {'cA_1_2': 2.0, 'cB_2_3': 1.0, 'u_2_1': 0.4}

    carch = copy(ilp)
    carch.fix_u([[0.2], [0.4], [0.4]])
    carch.create_model()
    assert _fixed_entry_values(carch.model) == {'cA_1_2': 2.0, 'cB_2_3': 1.0}

    uarch = copy(ilp)
    uarch.fix_c([[1, 0, 0], [1, 0, 0]], [[1, 0, 0], [1, 0, 0]])
    uarch.create_model()
    assert _fixed_entry_values(uarch.model) == {'u_2_1': 0.4}


def _constraint_variable_names(constraints):
    return {
        frozenset(variable.name for variable in generate_standard_repn(constraint.body).linear_vars)
        for constraint in constraints.values()
    }


def _build_clonal_ilp(fixed_matrix=None):
    cluster_ids = pd.Index(['cluster1', 'cluster2'], name='#ID')
    f_a = pd.DataFrame([[1.0], [2.0]], index=cluster_ids, columns=['sample1'])
    f_b = pd.DataFrame([[1.0], [1.0]], index=cluster_ids, columns=['sample1'])
    ilp = ILPSubset(
        n=3,
        cn_max=6,
        d=-1,
        mu=0.03,
        ampdel=False,
        copy_numbers={'cluster1': (1, 1)},
        f_a=f_a,
        f_b=f_b,
        w=pd.Series([1.0, 2.0], index=cluster_ids),
        purities=None,
        fixed_matrix=fixed_matrix,
    )
    ilp.create_model()
    return ilp


def test_python_models_fix_clonal_clusters_with_bounds_and_omit_implied_rows():
    ilp = _build_clonal_ilp()

    assert _fixed_entry_values(ilp.model) == {'cA_1_2': 1.0, 'cB_1_2': 1.0, 'cA_1_3': 1.0, 'cB_1_3': 1.0}
    free_tumor_cells = {frozenset(('cA_2_2', 'cB_2_2')), frozenset(('cA_2_3', 'cB_2_3'))}
    assert _constraint_variable_names(ilp.model.homdel_constraints) == free_tumor_cells
    total_rows = {
        names
        for names in _constraint_variable_names(ilp.model.constraints)
        if all(name.startswith(('cA_', 'cB_')) for name in names)
    }
    assert total_rows == free_tumor_cells


def test_python_models_enforce_conflicting_fixed_values_as_equalities(tmp_path):
    fixed = _parse_and_resolve(tmp_path, 'A\tcluster1\tclone1\t2\n')
    ilp = _build_clonal_ilp(fixed)

    assert (ilp.model.cA_1_2.lb, ilp.model.cA_1_2.ub) == (1, 1)
    conflicts = [
        (constraint.body.name, float(constraint.lower()), float(constraint.upper()))
        for constraint in ilp.model.fixed_value_conflicts.values()
    ]
    assert conflicts == [('cA_1_2', 2.0, 2.0)]


@pytest.mark.parametrize('split', (False, True))
def test_python_models_apply_custom_homdel_size_in_full_and_carch_modes(split):
    weights = pd.Series([0.25, 0.75], index=['cluster1', 'cluster2'])
    if split:
        ilp = _build_split_ilp(None, weights=weights, max_homdel_size=0.75)
    else:
        ilp = _build_ilp(n=3, weights=weights, max_homdel_size=0.75)

    expected = {
        frozenset(('cA_2_2', 'cB_2_2')),
        frozenset(('cA_2_3', 'cB_2_3')),
    }
    assert _constraint_variable_names(ilp.model.homdel_constraints) == expected

    carch = copy(ilp)
    carch.fix_u([[0.5], [0.5], [0.0]])
    carch.create_model()
    assert carch.max_homdel_size == pytest.approx(0.75)
    assert _constraint_variable_names(carch.model.homdel_constraints) == expected


def test_fixed_matrix_homdel_validation_uses_custom_inclusive_cutoff(tmp_path):
    fixed = _parse_and_resolve(tmp_path, 'A\tcluster1\tclone1\t0\nB\tcluster1\tclone1\t0\n')
    options = {
        'cn_max': 6,
        'mu': 0.03,
        'ampdel': False,
        'weights': {'cluster1': 1, 'cluster2': 99},
    }

    assert validate_fixed_matrix(fixed, max_homdel_size=0.02, **options) is fixed
    with pytest.raises(FixedMatrixFormatError, match='cannot have fixed tumor copy number'):
        validate_fixed_matrix(fixed, max_homdel_size=0.01, **options)


def test_fixed_matrix_homdel_validation_ignores_zero_weight_clusters(tmp_path):
    fixed = _parse_and_resolve(tmp_path, 'A\tcluster1\tclone1\t0\nB\tcluster1\tclone1\t0\n')

    assert (
        validate_fixed_matrix(
            fixed,
            cn_max=6,
            mu=0.03,
            ampdel=False,
            weights={'cluster1': 0, 'cluster2': 1},
            max_homdel_size=0,
        )
        is fixed
    )


def test_fixed_entries_condition_coordinate_descent_seeds_and_hot_starts(tmp_path):
    fixed = _parse_and_resolve(
        tmp_path,
        'A\tcluster1\tclone1\t2\nB\tcluster2\tclone2\t1\n' 'U\tnormal\tsample1\t0.2\nU\tclone1\tsample1\t0.3\n',
    )
    ilp = _build_ilp(n=3, fixed_matrix=fixed)

    seed = ilp.build_random_u(random_seed=7)
    assert seed[0][0] == 0.2
    assert seed[1][0] == 0.3
    assert seed[2][0] == pytest.approx(0.5)
    assert sum(seed[:, 0]) == pytest.approx(1)

    start_a, start_b = ilp.first_hot_start()
    assert start_a[0][1] == 2
    assert start_b[1][2] == 1

    ilp.hot_start(start_a, start_b)
    assert ilp.cA[0][1].value == 2
    assert ilp.cB[1][2].value == 1

    ilp.hot_start(start_a, start_b, [[0.2], [0.3], [0.5]])
    assert [ilp.u[clone_index][0].value for clone_index in range(3)] == [0.2, 0.3, 0.5]
    assert ilp.model.x_2_1.value == 1
    assert (ilp.model.bitcA_1_1_2.value, ilp.model.bitcA_2_1_2.value) == (0, 1)


def test_hot_start_keeps_clone_columns_in_place():
    ilp = _build_ilp(n=4)
    start_a = [[1, 3, 2, 1], [1, 3, 2, 1]]
    start_b = [[1, 1, 1, 0], [1, 1, 1, 0]]

    ilp.hot_start(start_a, start_b)

    assert [[ilp.cA[m][n].value for n in range(4)] for m in range(2)] == start_a
    assert [[ilp.cB[m][n].value for n in range(4)] for m in range(2)] == start_b


def test_fixed_subproblem_matrices_reject_invalid_shapes_and_sums():
    ilp = _build_ilp(n=3)

    with pytest.raises(ValueError, match='Fixed U must have shape'):
        ilp.fix_u([[1]])
    with pytest.raises(ValueError, match='sums to'):
        ilp.fix_u([[0.2], [0.3], [0.4]])
    with pytest.raises(ValueError, match='Fixed A must have shape'):
        ilp.fix_c([[1, 1, 1]], [[1, 1, 1], [1, 1, 1]])


def test_verify_fixed_matrix_solution_checks_all_requested_cells(tmp_path):
    fixed = _parse_and_resolve(
        tmp_path,
        'A\tcluster1\tclone1\t2\nB\tcluster2\tclone2\t1\nU\tclone1\tsample1\t0.4\n',
    )
    cA = [[1, 2, 0], [1, 0, 0]]
    cB = [[1, 0, 0], [1, 0, 1]]
    u = [[0.2], [0.4005], [0.3995]]

    verify_fixed_matrix_solution(fixed, cA, cB, u)

    cA[0][1] = 3
    with pytest.raises(RuntimeError, match='solver returned 3.*expected 2'):
        verify_fixed_matrix_solution(fixed, cA, cB, u)


def test_fixed_matrix_solve_preserves_constraints_with_positive_residual(tmp_path):
    if not solver_available('gurobi'):
        pytest.skip('a working Gurobi Pyomo solver is unavailable')
    fixed = _parse_and_resolve(
        tmp_path,
        'A\tcluster1\tclone1\t3\nU\tnormal\tsample1\t0.5\nU\tclone1\tsample1\t0.5\n',
        clones=2,
    )
    ilp = _build_ilp(n=2, fixed_matrix=fixed)

    objective, c_a, c_b, u, _, _ = ilp.run(solver_type='gurobi', timelimit=30)

    assert objective > 0
    assert c_a[0][1] == 3
    assert u[0][0] == pytest.approx(0.5)
    assert u[1][0] == pytest.approx(0.5)
    verify_fixed_matrix_solution(fixed, c_a, c_b, u)


def test_partially_fixed_matrix_can_pass_local_validation_but_make_solver_infeasible(tmp_path):
    if not solver_available('gurobi'):
        pytest.skip('a working Gurobi Pyomo solver is unavailable')
    # clone2 has zero copies everywhere, so it cannot receive the fixed positive proportion
    fixed = _parse_and_resolve(
        tmp_path,
        'A\tcluster1\tclone2\t0\nB\tcluster1\tclone2\t0\n'
        'A\tcluster2\tclone2\t0\nB\tcluster2\tclone2\t0\n'
        'U\tclone2\tsample1\t0.5\n',
        clones=3,
    )
    assert (
        validate_fixed_matrix(
            fixed,
            cn_max=6,
            mu=0.03,
            ampdel=False,
            weights={'cluster1': 1, 'cluster2': 1},
            max_homdel_size=1,
        )
        is fixed
    )
    ilp = _build_ilp(n=3, fixed_matrix=fixed, max_homdel_size=1)

    assert ilp.run(solver_type='gurobi', timelimit=30) is None


def test_unforced_model_selection_routes_same_fixed_matrix_to_both_branches(tmp_path, monkeypatch):
    input_prefix = DATA_FOLDER / 'fw' / 'bulk'
    fixed_path = DATA_FOLDER / 'fw' / 'fixed-matrix.tsv'
    monkeypatch.setattr(config.compute_cn, 'solver', 'gurobi')
    seen = []

    def capture_diploid(neutral, args, clonal=None):
        seen.append(('diploid', args['fixed']))
        return [(2, 0.0, str(tmp_path / 'diploid'))]

    def capture_tetraploid(clonal, scale, size, args):
        seen.append(('tetraploid', args['fixed']))
        return [(2, 0.0, str(tmp_path / 'tetraploid'))]

    monkeypatch.setattr('hatchet.bin.HATCHet.runningDiploid', capture_diploid)
    monkeypatch.setattr(
        'hatchet.bin.HATCHet.findClonalClusters', lambda **kwargs: ({'1': (2, 2), '2': (2, 1)}, ('1', '2'))
    )
    monkeypatch.setattr('hatchet.bin.HATCHet.runningTetraploid', capture_tetraploid)
    monkeypatch.setattr('hatchet.bin.HATCHet.select', lambda **kwargs: None)

    from hatchet.bin.HATCHet import main as compute_cn

    compute_cn(
        args=[
            '-i',
            str(input_prefix),
            '-x',
            str(tmp_path),
            '-n',
            '2',
            '--fixed-matrix',
            str(fixed_path),
            '-v',
            '0',
        ]
    )

    assert [branch for branch, _ in seen] == ['diploid', 'tetraploid']
    assert seen[0][1] is seen[1][1]


@pytest.mark.parametrize(
    ('runner', 'expected_bound'),
    ((runningDiploid, 6), (runningTetraploid, 12)),
)
def test_fixed_matrix_uses_active_ploidy_copy_bound(tmp_path, monkeypatch, runner, expected_bound):
    fixed = object()
    args = {
        'sample_clonal_entries': None,
        'solver': 'solve',
        'input': 'bulk',
        'ampdel': False,
        'd': None,
        'j': None,
        'p': 1,
        'u': None,
        'm': None,
        's': None,
        'f': None,
        'r': None,
        'M': 1,
        'fixedmatrix': '/tmp/fixed.tsv',
        'fixed': fixed,
        'purities': None,
        'eD': 6,
        'eT': 12,
        'ln': 2,
        'un': 2,
        'v': 0,
        'x': str(tmp_path),
    }
    calls = []
    monkeypatch.setattr(
        'hatchet.bin.HATCHet.execute',
        lambda passed_args, basecmd, n, outprefix: calls.append((dict(passed_args), basecmd)) or 0,
    )

    if runner is runningDiploid:
        runner(neutral='1', args=args)
    else:
        runner(clonal={'1': (2, 2), '2': (2, 1)}, scale=('1', '2'), size={}, args=args)

    passed_args, basecmd = calls[0]
    assert passed_args['fixed'] is fixed
    assert passed_args['e'] == expected_bound
    assert f'-e {expected_bound}' in basecmd
    assert '-F /tmp/fixed.tsv' in basecmd


def test_unforced_model_selection_reports_branch_specific_fixed_matrix_bound_conflict(tmp_path, monkeypatch):
    if not solver_available('gurobi'):
        pytest.skip('a working Gurobi Pyomo solver is unavailable')
    fixed_path = _write_fixed_matrix(
        tmp_path,
        f'{FORMAT_MARKER}\n# clones 2\nmatrix\trow\tcolumn\tvalue\nA\t2\tclone1\t7\n',
    )
    monkeypatch.setattr(config.compute_cn, 'solver', 'gurobi')
    from hatchet.bin.HATCHet import main as compute_cn

    with pytest.raises(FixedMatrixFormatError, match='copy-number bound 6'):
        compute_cn(
            args=[
                '-i',
                str(DATA_FOLDER / 'fw' / 'bulk'),
                '-x',
                str(tmp_path),
                '-n',
                '2',
                '-p',
                '1',
                '-j',
                '1',
                '-v',
                '0',
                '--mode',
                '1',
                '-s',
                '30',
                '-eD',
                '8',
                '-eT',
                '6',
                '--fixed-matrix',
                str(fixed_path),
            ]
        )


def _tumor_copy_number_footprint(constraint):
    clusters, clones = set(), set()
    for variable in identify_variables(constraint.body):
        prefix, *indices = variable.name.split('_')
        if prefix in ('cA', 'cB') and int(indices[1]) > 1:
            clusters.add(indices[0])
            clones.add(indices[1])
    return clusters, clones


def test_python_models_have_no_clone_order_constraints():
    ilp = _build_ilp(n=4)
    carch = copy(ilp)
    carch.fix_u([[0.4], [0.3], [0.2], [0.1]])
    carch.create_model()

    for model in (ilp.model, carch.model):
        for constraint in model.component_data_objects(pe.Constraint):
            clusters, clones = _tumor_copy_number_footprint(constraint)
            assert len(clusters) < 2 or len(clones) < 2, constraint.name


UNSORTED_SOLUTION = (
    [[1, 3, 1, 1], [1, 2, 1, 1]],
    [[1, 1, 0, 0], [1, 1, 1, 1]],
    [[0.1], [0.5], [0.3], [0.1]],
)
SORTED_SOLUTION = (
    [[1, 1, 1, 3], [1, 1, 1, 2]],
    [[1, 0, 0, 1], [1, 1, 1, 1]],
    [[0.1], [0.3], [0.1], [0.5]],
)


def test_canonicalize_clone_order_sorts_tumor_clones_stably():
    # Burdens: clone1 = 10, clone2 = 5, clone3 = 5; the clone2/clone3 tie is kept in input order
    assert canonicalize_clone_order(*UNSORTED_SOLUTION) == SORTED_SOLUTION


def test_finalize_solution_sorts_clones_unless_the_fixed_matrix_addresses_them(tmp_path):
    result = (0.0, *UNSORTED_SOLUTION, ('cluster1', 'cluster2'), ('sample1',))

    assert _finalize_solution(result, None)[1:4] == SORTED_SOLUTION

    normal_only = _parse_and_resolve(tmp_path, 'U\tnormal\tsample1\t0.1\n', clones=4)
    assert not normal_only.addresses_tumor_clones
    assert _finalize_solution(result, normal_only)[1:4] == SORTED_SOLUTION

    addressed = _parse_and_resolve(tmp_path, 'A\tcluster1\tclone1\t3\n', clones=4)
    assert addressed.addresses_tumor_clones
    assert _finalize_solution(result, addressed)[1:4] == UNSORTED_SOLUTION


def test_cpp_models_have_no_clone_order_constraints():
    source = (Path(__file__).parents[1] / 'src' / 'ilp-min.cpp').read_text()
    assert 'buildSymmetryBreaking' not in source
    assert 'sort_indexes' not in source
