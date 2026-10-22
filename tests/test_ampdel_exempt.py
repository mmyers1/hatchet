from copy import copy
from pathlib import Path
import os
import shlex
import subprocess

import pandas as pd
import pytest

from hatchet import config
from hatchet.bin.HATCHet import main as compute_cn_main
from hatchet.bin.HATCHet import makeBaseCMD, parsing_arguments
from hatchet.utils.fixed_matrix import (
    FixedMatrixFormatError,
    parse_fixed_matrix,
    resolve_fixed_matrix,
    validate_fixed_matrix,
)
from hatchet.utils.sample_clonal import (
    SampleClonalFormatError,
    parse_sample_clonal_inline,
    resolve_sample_clonal,
    validate_sample_clonal_constraints,
)
from hatchet.utils.solve import solve, solver_available
from hatchet.utils.solve.ilp_subset import ILPSubset, ILPSubsetSplit


DATA_FOLDER = Path(__file__).parent / 'data'
REPOSITORY = Path(__file__).parents[1]
FIXED_MATRIX_MARKER = '# HATCHET_FIXED_ENTRIES 1'
CLUSTER_IDS = pd.Index(['cluster1', 'cluster2'], name='#ID')


def _build_ilp(model_class=ILPSubset, n=3, ampdel_exempt=None, copy_numbers=None, fixed_matrix=None, f_a=None):
    sample_ids = ['sample1', 'sample2']
    if f_a is None:
        f_a = pd.DataFrame([[1.0, 1.0], [2.0, 2.0]], index=CLUSTER_IDS, columns=sample_ids)
    f_b = pd.DataFrame([[1.0, 1.0], [1.0, 1.0]], index=CLUSTER_IDS, columns=sample_ids)
    options = dict(
        n=n,
        cn_max=6,
        d=-1,
        mu=0.03,
        ampdel=True,
        copy_numbers=copy_numbers or {},
        f_a=f_a,
        f_b=f_b,
        w=pd.Series([1.0, 2.0], index=CLUSTER_IDS),
        purities=None,
        fixed_matrix=fixed_matrix,
        base=1,
        max_homdel_size=1,
        ampdel_exempt=ampdel_exempt,
    )
    if model_class is ILPSubsetSplit:
        options['binsA'] = {c: f_a.loc[[c]].reset_index(drop=True) for c in CLUSTER_IDS}
        options['binsB'] = {c: f_b.loc[[c]].reset_index(drop=True) for c in CLUSTER_IDS}
        options['lengths'] = {c: [1.0] for c in CLUSTER_IDS}
    return model_class(**options)


def _resolve_fixed_matrix(tmp_path, rows, clones=3):
    path = tmp_path / 'fixed.tsv'
    path.write_text(f'{FIXED_MATRIX_MARKER}\n# clones {clones}\nmatrix\trow\tcolumn\tvalue\n{rows}')
    return resolve_fixed_matrix(parse_fixed_matrix(path), tuple(CLUSTER_IDS), ('sample1', 'sample2'))


@pytest.mark.parametrize('model_class', (ILPSubset, ILPSubsetSplit))
def test_exempt_cluster_has_no_ampdel_variables_or_constraints(model_class):
    n = 3
    enforced = _build_ilp(model_class, n=n)
    enforced.create_model()
    exempt = _build_ilp(model_class, n=n, ampdel_exempt={'cluster2'})
    exempt.create_model()

    assert set(enforced._adA) == {0, 1}
    assert set(exempt._adA) == set(exempt._adB) == {0}
    assert exempt.model.find_component('adA_2') is None
    assert enforced.model.nconstraints() - exempt.model.nconstraints() == 4 * (n - 1)
    assert enforced.model.nvariables() - exempt.model.nvariables() == 2

    carch = copy(exempt)
    assert carch.ampdel_exempt == frozenset({'cluster2'})
    carch.fix_u([[0.2, 0.2], [0.4, 0.4], [0.4, 0.4]])
    carch.create_model()
    assert set(carch._adA) == {0}


def test_exempting_a_clonal_cluster_changes_nothing():
    clonal = {'cluster1': (1, 1)}
    without = _build_ilp(copy_numbers=clonal)
    without.create_model()
    with_exempt = _build_ilp(copy_numbers=clonal, ampdel_exempt={'cluster1'})
    with_exempt.create_model()

    assert set(without._adA) == set(with_exempt._adA) == {1}
    assert without.model.nconstraints() == with_exempt.model.nconstraints()


def test_first_hot_start_does_not_clamp_exempt_clusters():
    # Clone 1 starts from sample2 (A=3, amplified); clone 2 starts from sample1 (A=0, deleted)
    f_a = pd.DataFrame([[0.0, 3.0], [2.0, 2.0]], index=CLUSTER_IDS, columns=['sample1', 'sample2'])

    enforced_a, _ = _build_ilp(f_a=f_a).first_hot_start()
    exempt_a, _ = _build_ilp(f_a=f_a, ampdel_exempt={'cluster1'}).first_hot_start()

    assert list(enforced_a[0]) == [1, 3, 1]
    assert list(exempt_a[0]) == [1, 3, 0]


def test_validate_fixed_matrix_skips_exempt_clusters(tmp_path):
    fixed = _resolve_fixed_matrix(tmp_path, 'A\tcluster1\tclone1\t1\nA\tcluster1\tclone2\t3\n')
    options = {'cn_max': 6, 'mu': 0.03, 'ampdel': True, 'base': 2}

    with pytest.raises(FixedMatrixFormatError, match='force both amplification and deletion'):
        validate_fixed_matrix(fixed, **options)
    with pytest.raises(FixedMatrixFormatError, match='force both amplification and deletion'):
        validate_fixed_matrix(fixed, ampdel_exempt={'cluster2'}, **options)
    assert validate_fixed_matrix(fixed, ampdel_exempt={'cluster1'}, **options) is fixed


def test_validate_sample_clonal_skips_exempt_clusters():
    resolved = resolve_sample_clonal(
        parse_sample_clonal_inline('TumorA=c1:0:1,c2:2:1;TumorB=c1:2:1,c3:1:1'),
        ('c1', 'c2', 'c3'),
        ('TumorA', 'TumorB'),
    )
    options = {
        'cn_max': 6,
        'base': 1,
        'ampdel': True,
        'weights': pd.Series({'c1': 1.0, 'c2': 1.0, 'c3': 1.0}),
    }

    with pytest.raises(SampleClonalFormatError, match='force both amplification and deletion'):
        validate_sample_clonal_constraints(resolved, [0.5, 0.5], **options)
    validate_sample_clonal_constraints(resolved, [0.5, 0.5], ampdel_exempt={'c1'}, **options)


def test_solve_resolves_exempt_labels_to_seg_identifiers(monkeypatch):
    monkeypatch.setattr('hatchet.utils.solve.solver_available', lambda solver=None: True)
    captured = {}

    class StopAtCoordinateDescent(Exception):
        pass

    def fake_coordinate_descent(**kwargs):
        captured.update(kwargs)
        raise StopAtCoordinateDescent()

    monkeypatch.setattr('hatchet.utils.solve.CoordinateDescent', fake_coordinate_descent)
    files = {'bbc_file': DATA_FOLDER / 'fw' / 'bulk.bbc', 'seg_file': DATA_FOLDER / 'fw' / 'bulk.seg'}

    with pytest.raises(StopAtCoordinateDescent):
        solve(clonal=None, n=3, ampdel_exempt=['3', '7'], **files)
    # The seg file's #ID column is numeric, so string labels must resolve to the integer identifiers
    assert captured['ampdel_exempt'] == frozenset({3, 7})

    with pytest.raises(ValueError, match='not found in the seg file: 99'):
        solve(clonal=None, n=3, ampdel_exempt=['3', '99'], **files)
    with pytest.raises(ValueError, match='require the amp-del constraint'):
        solve(clonal=None, n=3, ampdel=False, ampdel_exempt=['3'], **files)


def _compute_cn_input(tmp_path):
    prefix = tmp_path / 'bulk'
    prefix.with_suffix('.seg').touch()
    prefix.with_suffix('.bbc').touch()
    return str(prefix)


def test_compute_cn_parses_ampdel_exempt_from_cli_and_config(tmp_path, monkeypatch):
    monkeypatch.setattr(config.compute_cn, 'solver', 'gurobi')
    base_args = ['-i', _compute_cn_input(tmp_path), '-x', str(tmp_path), '-n', '3']

    assert parsing_arguments(base_args)['ampdel_exempt'] is None
    assert parsing_arguments(base_args + ['--ampdel-exempt', '3,7, 3,'])['ampdel_exempt'] == ['3', '7']

    # The ini reader converts a single numeric value to int
    monkeypatch.setattr(config.compute_cn, 'ampdelexempt', 5)
    assert parsing_arguments(base_args)['ampdel_exempt'] == ['5']

    with pytest.raises(ValueError, match='cannot be used with --noampdel'):
        parsing_arguments(base_args + ['-f'])
    with pytest.raises(ValueError, match='No cluster IDs'):
        parsing_arguments(base_args + ['--ampdel-exempt', ' , '])


def test_compute_cn_rejects_exempt_clusters_missing_from_seg(tmp_path, monkeypatch):
    monkeypatch.setattr(config.compute_cn, 'solver', 'gurobi')
    args = ['-i', str(DATA_FOLDER / 'fw' / 'bulk'), '-x', str(tmp_path), '-n', '3', '--ampdel-exempt', '3,99']

    with pytest.raises(ValueError, match='not found in SEG file: 99'):
        compute_cn_main(args)


@pytest.mark.parametrize(
    ('ampdel', 'ampdel_exempt', 'expected'),
    [(True, ['3', '7'], '3,7'), (True, None, None), (False, ['3'], None)],
)
def test_cpp_command_forwards_ampdel_exempt_only_with_ampdel(ampdel, ampdel_exempt, expected):
    args = {key: None for key in ('d', 'j', 'p', 'u', 'm', 's', 'f', 'r', 'M', 'fixedmatrix')}
    args.update({'solver': 'solve', 'input': 'bulk', 'ampdel': ampdel, 'ampdel_exempt': ampdel_exempt})

    tokens = shlex.split(makeBaseCMD(args, 6))

    if expected is None:
        assert '--ampdel-exempt' not in tokens
    else:
        assert tokens[tokens.index('--ampdel-exempt') + 1] == expected
        assert '-f' in tokens


def test_solver_allows_amplification_and_deletion_only_in_exempt_cluster(tmp_path):
    if not solver_available('gurobi'):
        pytest.skip('a working Gurobi Pyomo solver is unavailable')
    fixed = _resolve_fixed_matrix(tmp_path, 'A\tcluster1\tclone1\t3\nA\tcluster1\tclone2\t0\n')

    enforced = _build_ilp(fixed_matrix=fixed)
    enforced.create_model()
    assert enforced.run(solver_type='gurobi', timelimit=30) is None

    exempt = _build_ilp(fixed_matrix=fixed, ampdel_exempt={'cluster1'})
    exempt.create_model()
    _, c_a, _, _, _, _ = exempt.run(solver_type='gurobi', timelimit=30)
    assert c_a[0][1:] == [3, 0]


@pytest.fixture(scope='module')
def native_solver():
    configured = os.getenv('HATCHET_TEST_CPP_SOLVER')
    candidates = [Path(configured)] if configured else []
    candidates += [REPOSITORY / 'build-native' / 'solve', REPOSITORY / 'src' / 'hatchet' / 'solve']
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            result = subprocess.run([str(candidate), '-h'], text=True, capture_output=True)
            if '--ampdel-exempt' in result.stdout:
                return candidate
    pytest.skip('build the native solver with --ampdel-exempt support or set HATCHET_TEST_CPP_SOLVER')


def _write_native_inputs(tmp_path):
    """Cluster 1 is the clonal diploid cluster; cluster 2 is forced to be amplified in clone1 and deleted in clone2."""
    prefix = tmp_path / 'native-input'
    prefix.with_suffix('.seg').write_text(
        '#ID\tSAMPLE\t#BINS\tRD\t#SNPS\tCOV\tALPHA\tBETA\tBAF\n'
        '1\tTumor\t1\t0.6666666666666666\t20\t30\t10\t10\t0.5\n'
        '2\tTumor\t1\t0.8333333333333334\t20\t30\t12\t8\t0.4\n'
    )
    prefix.with_suffix('.bbc').write_text(
        '#CHR\tSTART\tEND\tSAMPLE\tRD\t#SNPS\tCOV\tALPHA\tBETA\tBAF\tCLUSTER\n'
        'chr1\t0\t100\tTumor\t0.6666666666666666\t20\t30\t10\t10\t0.5\t1\n'
        'chr1\t100\t200\tTumor\t0.8333333333333334\t20\t30\t12\t8\t0.4\t2\n'
    )
    fixed = tmp_path / 'native-fixed.tsv'
    fixed.write_text(
        f'{FIXED_MATRIX_MARKER}\n# clones 3\nmatrix\trow\tcolumn\tvalue\nA\t2\tclone1\t3\nA\t2\tclone2\t0\n'
    )
    return prefix, fixed


def _native_command(solver, prefix, fixed, tmp_path, *extra, mode=1):
    return [
        str(solver),
        str(prefix),
        '-n', '3', '-c', '1:1:1', '-e', '6', '-j', '1', '-p', '2', '-i', '2', '-u', '0.1',
        '-s', '30', '-r', '7', '-M', str(mode), '-v', '0', '-F', str(fixed),
        '-o', str(tmp_path / f'native-out-m{mode}'),
        *extra,
    ]  # fmt: skip


@pytest.mark.parametrize(
    ('extra', 'message'),
    [
        (['-f'], 'force both amplification and deletion'),
        (['-f', '--ampdel-exempt', '99'], 'amp-del exempt cluster with ID 99 has not been found'),
        (['--ampdel-exempt', '2'], 'require the amp-del constraint -f'),
    ],
)
def test_native_solver_rejects_invalid_ampdel_exemptions_before_solving(tmp_path, native_solver, extra, message):
    prefix, fixed = _write_native_inputs(tmp_path)

    command = _native_command(native_solver, prefix, fixed, tmp_path, *extra)
    result = subprocess.run(command, text=True, capture_output=True)

    assert result.returncode != 0
    assert message in result.stderr + result.stdout


@pytest.mark.parametrize('mode', (0, 1, 2), ids=('cd-plus-ilp', 'ilp', 'cd'))
def test_native_solver_allows_amplification_and_deletion_in_exempt_cluster(tmp_path, native_solver, mode):
    if not solver_available('cpp'):
        pytest.skip('GRB_LICENSE_FILE is not set; load the Gurobi module to run native solves')
    prefix, fixed = _write_native_inputs(tmp_path)

    subprocess.run(
        _native_command(native_solver, prefix, fixed, tmp_path, '-f', '--ampdel-exempt', '2', mode=mode), check=True
    )

    bbc = pd.read_table(tmp_path / f'native-out-m{mode}.bbc.ucn.tsv')
    cluster2 = bbc[bbc['CLUSTER'] == 2].iloc[0]
    tumor_a = sorted(int(cluster2[f'cn_clone{i}'].split('|')[0]) for i in (1, 2))
    assert tumor_a == [0, 3]


def test_documentation_and_ini_cover_ampdel_exempt():
    documentation = (REPOSITORY / 'docs' / 'source' / 'doc_compute_cn.md').read_text()
    assert '`--ampdel-exempt`' in documentation
    assert '`ampdelexempt`' in documentation
    for ini in (REPOSITORY / 'src' / 'hatchet' / 'hatchet.ini', REPOSITORY / 'script' / 'hatchet.ini'):
        assert 'ampdelexempt =' in ini.read_text()
