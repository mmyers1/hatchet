from dataclasses import FrozenInstanceError
from copy import copy
from pathlib import Path
import re
import shlex
import shutil
import subprocess

import pandas as pd
import pytest
from pyomo.repn.standard_repn import generate_standard_repn

from hatchet import config
from hatchet.bin.HATCHet import execute, load_fixed_matrix_entries, makeBaseCMD, parsing_arguments
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
    fixed = _parse_and_resolve(
        tmp_path,
        'A\tcluster1\tclone1\t2\n' 'B\tcluster1\tclone1\t1\n' 'U\tnormal\tsample1\t0.2\n' 'U\tclone1\tsample1\t0.3\n',
    )

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
        (
            'A\tcluster1\tclone1\t2\nB\tcluster1\tclone1\t1\n'
            'A\tcluster2\tclone1\t2\nB\tcluster2\tclone1\t1\n'
            'A\tcluster1\tclone2\t1\nB\tcluster1\tclone2\t1\n'
            'A\tcluster2\tclone2\t1\nB\tcluster2\tclone2\t1\n',
            {},
            'fully fixed clone1 burden',
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
        },
        6,
    )
    tokens = shlex.split(command)
    fixed_matrix_index = tokens.index('-F')
    assert tokens.count('-F') == 1
    assert tokens[fixed_matrix_index + 1] == '/tmp/fixed entries.tsv'


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


def _build_ilp(n=4, fixed_matrix=None):
    cluster_ids = pd.Index(['cluster1', 'cluster2'], name='#ID')
    sample_ids = ['sample1']
    f_a = pd.DataFrame([[1.0], [2.0]], index=cluster_ids, columns=sample_ids)
    f_b = pd.DataFrame([[1.0], [1.0]], index=cluster_ids, columns=sample_ids)
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
    )
    ilp.create_model()
    return ilp


def _fixed_equalities(model):
    return {constraint.body.name: float(constraint.lower()) for constraint in model.fixed_matrix_constraints.values()}


def test_python_models_enforce_and_preserve_sparse_fixed_entries(tmp_path):
    fixed = _parse_and_resolve(
        tmp_path,
        'A\tcluster1\tclone1\t2\nB\tcluster2\tclone2\t1\nU\tclone1\tsample1\t0.4\n',
    )
    ilp = _build_ilp(n=3, fixed_matrix=fixed)

    assert _fixed_equalities(ilp.model) == {'cA_1_2': 2.0, 'cB_2_3': 1.0, 'u_2_1': 0.4}

    carch = copy(ilp)
    carch.fix_u([[0.2], [0.4], [0.4]])
    carch.create_model()
    assert carch.fixed_matrix is fixed
    assert carch._fixed_u[1][0] == 0.4
    assert _fixed_equalities(carch.model) == {'cA_1_2': 2.0, 'cB_2_3': 1.0}

    uarch = copy(ilp)
    uarch.fix_c([[1, 0, 0], [1, 0, 0]], [[1, 0, 0], [1, 0, 0]])
    uarch.create_model()
    assert uarch.fixed_matrix is fixed
    assert uarch._fixed_cA[0][1] == 2
    assert uarch._fixed_cB[1][2] == 1
    assert _fixed_equalities(uarch.model) == {'u_2_1': 0.4}


def _build_split_ilp(fixed_matrix):
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
    )
    ilp.create_model()
    return ilp


def test_binwise_models_enforce_and_preserve_sparse_fixed_entries(tmp_path):
    fixed = _parse_and_resolve(
        tmp_path,
        'A\tcluster1\tclone1\t2\nB\tcluster2\tclone2\t1\nU\tclone1\tsample1\t0.4\n',
    )
    ilp = _build_split_ilp(fixed)
    assert _fixed_equalities(ilp.model) == {'cA_1_2': 2.0, 'cB_2_3': 1.0, 'u_2_1': 0.4}

    carch = copy(ilp)
    carch.fix_u([[0.2], [0.4], [0.4]])
    carch.create_model()
    assert _fixed_equalities(carch.model) == {'cA_1_2': 2.0, 'cB_2_3': 1.0}

    uarch = copy(ilp)
    uarch.fix_c([[1, 0, 0], [1, 0, 0]], [[1, 0, 0], [1, 0, 0]])
    uarch.create_model()
    assert _fixed_equalities(uarch.model) == {'u_2_1': 0.4}


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


def _linear_coefficients(constraint):
    representation = generate_standard_repn(constraint.body)
    return dict(zip((variable.name for variable in representation.linear_vars), representation.linear_coefs))


def test_python_symmetry_breaking_orders_adjacent_tumor_clones():
    ilp = _build_ilp()
    constraint_count = len(ilp.model.constraints)
    symmetry_constraints = [ilp.model.constraints[index] for index in range(constraint_count - 1, constraint_count + 1)]

    for tumor_clone, constraint in enumerate(symmetry_constraints, start=1):
        coefficients = _linear_coefficients(constraint)
        for cluster_index in range(ilp.m):
            weight = ilp.symmCoeff(cluster_index)
            left_clone = tumor_clone + 1
            right_clone = tumor_clone + 2
            assert coefficients[f'cA_{cluster_index + 1}_{left_clone}'] == weight
            assert coefficients[f'cB_{cluster_index + 1}_{left_clone}'] == weight
            assert coefficients[f'cA_{cluster_index + 1}_{right_clone}'] == -weight
            assert coefficients[f'cB_{cluster_index + 1}_{right_clone}'] == -weight


def test_cpp_symmetry_breaking_orders_adjacent_tumor_clones():
    source = (Path(__file__).parents[1] / 'src' / 'ilp-min.cpp').read_text()
    assert 'sumI1 += _cA[s][i+1] * symmCoeff(s) + _cB[s][i+1] * symmCoeff(s);' in source
