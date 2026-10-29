from copy import copy
from pathlib import Path
import os
import re
import shlex
import shutil
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
CLUSTER_IDS = pd.Index(['c1', 'c2', 'c3'], name='#ID')
SAMPLE_IDS = ['s1', 's2', 's3']
# Tumor clone t makes up 90% of sample t, next to 10% normal cells. Cluster c2 needs three tumor states; c1 and c3
# are clonal.
TRUE_A = {'c1': (1, 1, 1), 'c2': (2, 1, 1), 'c3': (2, 2, 2)}
TRUE_B = {'c1': (1, 1, 1), 'c2': (1, 1, 0), 'c3': (1, 1, 1)}
MAJOR_U = [[0.1, 0.1, 0.1], [0.9, 0.0, 0.0], [0.0, 0.9, 0.0], [0.0, 0.0, 0.9]]


def _fractions(true_a=TRUE_A, true_b=TRUE_B):
    """Fractional copy numbers of the mixtures in MAJOR_U (tumor clone t = entry t of a cluster's tuple)."""

    def mix(true):
        return [[0.1 + 0.9 * true[c][sample] for sample in range(len(SAMPLE_IDS))] for c in CLUSTER_IDS]

    f_a = pd.DataFrame(mix(true_a), index=CLUSTER_IDS, columns=SAMPLE_IDS)
    f_b = pd.DataFrame(mix(true_b), index=CLUSTER_IDS, columns=SAMPLE_IDS)
    return f_a, f_b


def _build_ilp(
    model_class=ILPSubset, n=4, d=-1, cluster_cnstates=None, copy_numbers=None, f_a=None, f_b=None, fixed_matrix=None
):
    if f_a is None:
        f_a, f_b = _fractions()
    options = dict(
        n=n,
        cn_max=6,
        d=d,
        mu=0.03,
        ampdel=False,
        copy_numbers=copy_numbers or {},
        f_a=f_a,
        f_b=f_b,
        w=pd.Series(1.0, index=f_a.index),
        purities=None,
        fixed_matrix=fixed_matrix,
        base=1,
        max_homdel_size=1,
        cluster_cnstates=cluster_cnstates,
    )
    if model_class is ILPSubsetSplit:
        options['binsA'] = {c: f_a.loc[[c]].reset_index(drop=True) for c in f_a.index}
        options['binsB'] = {c: f_b.loc[[c]].reset_index(drop=True) for c in f_a.index}
        options['lengths'] = {c: [1.0] for c in f_a.index}
    return model_class(**options)


def _slot_counts(model):
    """Number of slot-state variables per 1-based cluster index."""
    counts = {}
    for variable in model.component_objects():
        match = re.fullmatch(r'sA_(\d+)_(\d+)', variable.name)
        if match:
            counts[int(match.group(1))] = counts.get(int(match.group(1)), 0) + 1
    return counts


def _component_names(model, prefix):
    return [variable.name for variable in model.component_objects() if variable.name.startswith(prefix)]


def _tumor_states(c_a, c_b, cluster_index):
    return {(c_a[cluster_index][t], c_b[cluster_index][t]) for t in range(1, len(c_a[cluster_index]))}


@pytest.mark.parametrize('model_class', (ILPSubset, ILPSubsetSplit))
def test_state_limit_variables_exist_only_for_limited_non_clonal_clusters(model_class):
    ilp = _build_ilp(model_class, d=1, cluster_cnstates={'c2': None}, copy_numbers={'c1': (1, 1)})
    ilp.create_model()

    # c1 is clonal and c2 is exempt, so only c3 gets a slot
    assert _slot_counts(ilp.model) == {3: 1}
    # restricted growth: tumor clone t may only use slots l < t, so a single slot gives one z per tumor clone
    assert sorted(_component_names(ilp.model, 'z_')) == ['z_3_2_1', 'z_3_3_1', 'z_3_4_1']


def test_cluster_limits_override_the_global_limit():
    ilp = _build_ilp(d=1, cluster_cnstates={'c2': 2})
    ilp.create_model()

    assert _slot_counts(ilp.model) == {1: 1, 2: 2, 3: 1}
    # clone 2 (tumor clone 1) only uses slot 1; tumor clones 2 and 3 may use slots 1 and 2
    assert sorted(_component_names(ilp.model, 'z_2_')) == ['z_2_2_1', 'z_2_3_1', 'z_2_3_2', 'z_2_4_1', 'z_2_4_2']


def test_cluster_limits_apply_without_a_global_limit():
    ilp = _build_ilp(d=-1, cluster_cnstates={'c2': 1})
    ilp.create_model()

    assert _slot_counts(ilp.model) == {2: 1}


def test_limit_of_at_least_the_tumor_clone_count_adds_nothing():
    unlimited = _build_ilp(d=-1)
    unlimited.create_model()
    vacuous = _build_ilp(d=3)  # n=4 has 3 tumor clones
    vacuous.create_model()

    assert _slot_counts(vacuous.model) == {}
    assert vacuous.model.nconstraints() == unlimited.model.nconstraints()
    assert vacuous.model.nvariables() == unlimited.model.nvariables()


def test_carch_models_with_a_limit_have_no_bit_variables():
    carch = copy(_build_ilp(d=1, cluster_cnstates={'c2': 2}))
    assert carch.cluster_cnstates == {'c2': 2}
    carch.fix_u(MAJOR_U)
    carch.create_model()

    assert _component_names(carch.model, 'bitc') == []
    assert _slot_counts(carch.model) == {1: 1, 2: 2, 3: 1}


def test_first_hot_start_has_at_most_the_limited_number_of_states():
    # d < min(n, k) used to read past the end of the target array (IndexError). c1 and c3 vary across samples, so
    # the all-samples start would give them three states.
    f_a, f_b = _fractions(dict(TRUE_A, c1=(1, 2, 3), c3=(3, 2, 1)))
    hc_a, hc_b = _build_ilp(d=1, cluster_cnstates={'c2': None}, f_a=f_a, f_b=f_b).first_hot_start()

    for cluster_index in (0, 2):
        assert len(_tumor_states(hc_a, hc_b, cluster_index)) == 1
    # an unlimited cluster keeps the all-samples start: clone i starts from sample i % min(n, k)
    assert list(hc_a[1]) == [1, 1, 1, 2]


def test_first_hot_start_allows_a_limit_above_the_sample_count():
    f_a = pd.DataFrame([[1.0, 2.0]], index=pd.Index(['c1'], name='#ID'), columns=['s1', 's2'])
    f_b = pd.DataFrame([[1.0, 1.0]], index=f_a.index, columns=f_a.columns)

    hc_a, hc_b = _build_ilp(n=5, d=3, f_a=f_a, f_b=f_b).first_hot_start()

    assert len(_tumor_states(hc_a, hc_b, 0)) <= 3


# Pyomo MIP solver for the solving tests; e.g. HATCHET_TEST_MIP_SOLVER=cbc where the Gurobi licence is unreachable
MIP_SOLVER = os.getenv('HATCHET_TEST_MIP_SOLVER', 'gurobi')


def _run_carch(ilp):
    carch = copy(ilp)
    carch.fix_u(MAJOR_U)
    carch.create_model()
    return carch.run(solver_type=MIP_SOLVER, timelimit=30)


def _require_mip_solver():
    if not solver_available(MIP_SOLVER):
        pytest.skip(f'the Pyomo solver {MIP_SOLVER} is unavailable')


def test_limit_of_two_lets_the_first_two_tumor_clones_differ():
    _require_mip_solver()
    # c2 has states (2,1), (1,1), (1,1): two states, with tumor clones 1 and 2 different. The old prefix-sum
    # symmetry rows forced tumor clones 1 and 2 to share a state in every cluster.
    true_a = dict(TRUE_A, c2=(2, 1, 1))
    true_b = dict(TRUE_B, c2=(1, 1, 1))
    f_a, f_b = _fractions(true_a, true_b)

    obj, c_a, c_b, _, _, _ = _run_carch(_build_ilp(d=2, f_a=f_a, f_b=f_b))

    assert obj == pytest.approx(0.0, abs=1e-6)
    assert [c_a[1][t] for t in (1, 2, 3)] == [2, 1, 1]


def test_exempt_and_limited_clusters_in_coordinate_descent_steps():
    _require_mip_solver()
    limited_obj, c_a, c_b, _, _, _ = _run_carch(_build_ilp(d=1))
    assert limited_obj > 0.1
    assert all(len(_tumor_states(c_a, c_b, cluster_index)) == 1 for cluster_index in range(3))

    exempt_obj, c_a, c_b, _, _, _ = _run_carch(_build_ilp(d=1, cluster_cnstates={'c2': None}))
    assert exempt_obj == pytest.approx(0.0, abs=1e-6)
    assert _tumor_states(c_a, c_b, 1) == {(2, 1), (1, 1), (1, 0)}

    two_obj, c_a, c_b, _, _, _ = _run_carch(_build_ilp(d=1, cluster_cnstates={'c2': 2}))
    assert len(_tumor_states(c_a, c_b, 1)) == 2
    assert 0.0 < two_obj < limited_obj


def test_exempt_cluster_in_the_full_ilp():
    _require_mip_solver()
    ilp = _build_ilp(d=1, cluster_cnstates={'c2': None})
    ilp.create_model()

    obj, c_a, c_b, _, _, _ = ilp.run(solver_type=MIP_SOLVER, timelimit=60)

    assert obj == pytest.approx(0.0, abs=1e-4)
    assert len(_tumor_states(c_a, c_b, 1)) == 3
    assert len(_tumor_states(c_a, c_b, 0)) == len(_tumor_states(c_a, c_b, 2)) == 1


def _resolve_fixed_matrix(tmp_path, rows, clones=4):
    path = tmp_path / 'fixed.tsv'
    path.write_text(f'{FIXED_MATRIX_MARKER}\n# clones {clones}\nmatrix\trow\tcolumn\tvalue\n{rows}')
    return resolve_fixed_matrix(parse_fixed_matrix(path), tuple(CLUSTER_IDS), tuple(SAMPLE_IDS))


def test_validate_fixed_matrix_applies_each_cluster_limit(tmp_path):
    fixed = _resolve_fixed_matrix(tmp_path, 'A\tc1\tclone1\t1\nB\tc1\tclone1\t0\nA\tc1\tclone2\t2\nB\tc1\tclone2\t0\n')
    options = {'cn_max': 6, 'mu': 0.03, 'ampdel': False, 'base': 1}

    with pytest.raises(FixedMatrixFormatError, match='exceeding cnstates 1'):
        validate_fixed_matrix(fixed, cnstates=1, **options)
    with pytest.raises(FixedMatrixFormatError, match='exceeding cnstates 1'):
        validate_fixed_matrix(fixed, cnstates=1, cluster_cnstates={'c2': None}, **options)
    with pytest.raises(FixedMatrixFormatError, match='exceeding cnstates 1'):
        validate_fixed_matrix(fixed, cluster_cnstates={'c1': 1}, **options)
    assert validate_fixed_matrix(fixed, cnstates=1, cluster_cnstates={'c1': None}, **options) is fixed
    assert validate_fixed_matrix(fixed, cnstates=1, cluster_cnstates={'c1': 2}, **options) is fixed


def test_validate_sample_clonal_applies_each_cluster_limit():
    resolved = resolve_sample_clonal(
        parse_sample_clonal_inline('TumorA=c1:0:1,c2:2:1;TumorB=c1:2:1,c3:1:1'),
        ('c1', 'c2', 'c3'),
        ('TumorA', 'TumorB'),
    )
    options = {'cn_max': 6, 'base': 1, 'ampdel': False, 'weights': pd.Series({'c1': 1.0, 'c2': 1.0, 'c3': 1.0})}

    with pytest.raises(SampleClonalFormatError, match='exceeding cnstates 1'):
        validate_sample_clonal_constraints(resolved, [0.5, 0.5], cnstates=1, **options)
    with pytest.raises(SampleClonalFormatError, match='exceeding cnstates 1'):
        validate_sample_clonal_constraints(resolved, [0.5, 0.5], cluster_cnstates={'c1': 1}, **options)
    validate_sample_clonal_constraints(resolved, [0.5, 0.5], cnstates=1, cluster_cnstates={'c1': None}, **options)


def test_solve_resolves_cluster_limit_labels_to_seg_identifiers(monkeypatch):
    monkeypatch.setattr('hatchet.utils.solve.solver_available', lambda solver=None: True)
    captured = {}

    class StopAtCoordinateDescent(Exception):
        pass

    def fake_coordinate_descent(**kwargs):
        captured.update(kwargs)
        raise StopAtCoordinateDescent()

    monkeypatch.setattr('hatchet.utils.solve.CoordinateDescent', fake_coordinate_descent)
    files = {'bbc_file': DATA_FOLDER / 'fw' / 'bulk.bbc', 'seg_file': DATA_FOLDER / 'fw' / 'bulk.seg'}

    # The seg file's #ID column is numeric, so string labels must resolve to the integer identifiers
    with pytest.raises(StopAtCoordinateDescent):
        solve(clonal=None, n=3, d=1, cnstates_exempt=['3', '7'], **files)
    assert captured['d'] == 1
    assert captured['cluster_cnstates'] == {3: None, 7: None}

    with pytest.raises(StopAtCoordinateDescent):
        solve(clonal=None, n=3, cluster_cnstates={'3': 2, '7': 4}, **files)
    assert captured['cluster_cnstates'] == {3: 2, 7: 4}

    with pytest.raises(ValueError, match='not found in the seg file: 99'):
        solve(clonal=None, n=3, d=1, cnstates_exempt=['3', '99'], **files)
    with pytest.raises(ValueError, match='not found in the seg file: 99'):
        solve(clonal=None, n=3, cluster_cnstates={'99': 2}, **files)
    with pytest.raises(ValueError, match='require a global --cnstates'):
        solve(clonal=None, n=3, cnstates_exempt=['3'], **files)
    with pytest.raises(ValueError, match='only one of'):
        solve(clonal=None, n=3, d=1, cnstates_exempt=['3'], cluster_cnstates={'7': 2}, **files)
    with pytest.raises(ValueError, match='positive integer'):
        solve(clonal=None, n=3, cluster_cnstates={'3': 0}, **files)


def _compute_cn_input(tmp_path):
    prefix = tmp_path / 'bulk'
    prefix.with_suffix('.seg').touch()
    prefix.with_suffix('.bbc').touch()
    return str(prefix)


def test_compute_cn_parses_cluster_limits_from_cli_and_config(tmp_path, monkeypatch):
    monkeypatch.setattr(config.compute_cn, 'solver', 'gurobi')
    base_args = ['-i', _compute_cn_input(tmp_path), '-x', str(tmp_path), '-n', '3']

    parsed = parsing_arguments(base_args)
    assert parsed['cnstates_exempt'] is None and parsed['cluster_cnstates'] is None
    assert parsing_arguments(base_args + ['-d', '2', '--cnstates-exempt', '3,7, 3,'])['cnstates_exempt'] == ['3', '7']
    assert parsing_arguments(base_args + ['--cluster-cnstates', '3:2, 7:4'])['cluster_cnstates'] == {'3': 2, '7': 4}

    # The ini reader converts a single numeric value to int
    monkeypatch.setattr(config.compute_cn, 'cnstatesexempt', 5)
    assert parsing_arguments(base_args + ['-d', '2'])['cnstates_exempt'] == ['5']
    with pytest.raises(ValueError, match='only one of'):
        parsing_arguments(base_args + ['-d', '2', '--cluster-cnstates', '3:2'])


@pytest.mark.parametrize(
    ('extra', 'message'),
    [
        (['--cnstates-exempt', '3'], 'requires a global -d/--cnstates'),
        (['-d', '2', '--cnstates-exempt', ' , '], 'No cluster IDs'),
        (['--cluster-cnstates', '3'], 'ID:INT'),
        (['--cluster-cnstates', '3:x'], 'ID:INT'),
        (['--cluster-cnstates', '3:0'], 'positive integer'),
        (['--cluster-cnstates', '3:2,3:4'], 'more than once'),
        (['-d', '2', '--cnstates-exempt', '3', '--cluster-cnstates', '7:2'], 'only one of'),
    ],
)
def test_compute_cn_rejects_invalid_cluster_limits(tmp_path, monkeypatch, extra, message):
    monkeypatch.setattr(config.compute_cn, 'solver', 'gurobi')
    with pytest.raises(ValueError, match=message):
        parsing_arguments(['-i', _compute_cn_input(tmp_path), '-x', str(tmp_path), '-n', '3'] + extra)


@pytest.mark.parametrize(
    'extra',
    [['-d', '1', '--cnstates-exempt', '3,99'], ['--cluster-cnstates', '3:1,99:2']],
)
def test_compute_cn_rejects_cluster_limits_missing_from_seg(tmp_path, monkeypatch, extra):
    monkeypatch.setattr(config.compute_cn, 'solver', 'gurobi')
    args = ['-i', str(DATA_FOLDER / 'fw' / 'bulk'), '-x', str(tmp_path), '-n', '3'] + extra

    with pytest.raises(ValueError, match='not found in SEG file: 99'):
        compute_cn_main(args)


@pytest.mark.parametrize(
    ('d', 'cnstates_exempt', 'cluster_cnstates', 'flag', 'expected'),
    [
        (2, ['3', '7'], None, '--cnstates-exempt', '3,7'),
        (None, None, {'3': 2, '7': 4}, '--cluster-cnstates', '3:2,7:4'),
        (2, None, None, '--cnstates-exempt', None),
    ],
)
def test_cpp_command_forwards_cluster_limits(d, cnstates_exempt, cluster_cnstates, flag, expected):
    args = {key: None for key in ('j', 'p', 'u', 'm', 's', 'f', 'r', 'M', 'fixedmatrix')}
    args.update(
        {
            'solver': 'solve',
            'input': 'bulk',
            'ampdel': False,
            'd': d,
            'cnstates_exempt': cnstates_exempt,
            'cluster_cnstates': cluster_cnstates,
        }
    )

    tokens = shlex.split(makeBaseCMD(args, 6))

    if expected is None:
        assert flag not in tokens
    else:
        assert tokens[tokens.index(flag) + 1] == expected


def test_documentation_and_ini_cover_cluster_limits():
    documentation = (REPOSITORY / 'docs' / 'source' / 'doc_compute_cn.md').read_text()
    assert '`--cnstates-exempt`' in documentation
    assert '`--cluster-cnstates`' in documentation
    assert 'This option is deprecated' not in documentation
    ini = (REPOSITORY / 'src' / 'hatchet' / 'hatchet.ini').read_text()
    assert re.search(r'^cnstatesexempt =', ini, re.MULTILINE)
    assert re.search(r'^clustercnstates =', ini, re.MULTILINE)


@pytest.fixture(scope='module')
def cpp_validators(tmp_path_factory):
    """The fixed-matrix and sample-clonal C++ validator drivers of test_fixed_matrix.py and test_sample_clonal.py."""
    compiler = shutil.which('g++')
    if compiler is None:
        pytest.skip('g++ is required for the native validator tests')
    built = {}
    for name, sources in (
        ('fixed_matrix', ('fixed_matrix.cpp', 'sample_clonal.cpp', 'utils.cpp')),
        ('sample_clonal', ('sample_clonal.cpp', 'utils.cpp')),
    ):
        executable = tmp_path_factory.mktemp(f'{name}-cpp') / 'reader'
        subprocess.run(
            [compiler, '-std=c++11', '-I' + str(REPOSITORY / 'src'), str(REPOSITORY / 'tests' / f'{name}_cpp_test.cpp')]
            + [str(REPOSITORY / 'src' / source) for source in sources]
            + ['-pthread', '-o', str(executable)],
            check=True,
        )
        built[name] = executable
    return built


@pytest.mark.parametrize(('scenario', 'returncode'), [('cnstates', 1), ('cnstates-exempt', 0)])
def test_cpp_fixed_matrix_validation_applies_each_cluster_limit(tmp_path, cpp_validators, scenario, returncode):
    path = tmp_path / 'fixed.tsv'
    path.write_text(
        f'{FIXED_MATRIX_MARKER}\n# clones 3\nmatrix\trow\tcolumn\tvalue\n'
        'A\tcluster1\tclone1\t1\nB\tcluster1\tclone1\t0\nA\tcluster1\tclone2\t2\nB\tcluster1\tclone2\t0\n'
    )

    command = [str(cpp_validators['fixed_matrix']), str(path), '3', scenario]
    result = subprocess.run(command, text=True, capture_output=True)

    assert result.returncode == returncode, result.stderr
    if returncode:
        assert 'exceed the maximum number of copy-number states' in result.stderr


@pytest.mark.parametrize(('scenario', 'returncode'), [('cnstates', 1), ('cnstates-exempt', 0)])
def test_cpp_sample_clonal_validation_applies_each_cluster_limit(tmp_path, cpp_validators, scenario, returncode):
    path = tmp_path / 'sample-clonal.tsv'
    path.write_text(
        '# HATCHET_SAMPLE_CLONAL 1\nsample\tcluster\tcn_a\tcn_b\n'
        'Tumor\tc1\t1\t1\nTumor\tc2\t2\t1\nTumor2\tc1\t3\t1\nTumor2\tc2\t2\t1\n'
    )

    result = subprocess.run(
        [str(cpp_validators['sample_clonal']), str(path), scenario, '0.005'], text=True, capture_output=True
    )

    assert result.returncode == returncode, result.stderr
    if returncode:
        assert 'require more copy-number states than allowed' in result.stderr


# ---- native solver


@pytest.fixture(scope='module')
def native_solver():
    configured = os.getenv('HATCHET_TEST_CPP_SOLVER')
    candidates = [Path(configured)] if configured else []
    candidates += [REPOSITORY / 'build-native' / 'solve', REPOSITORY / 'src' / 'hatchet' / 'solve']
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            result = subprocess.run([str(candidate), '-h'], text=True, capture_output=True)
            if '--cnstates-exempt' in result.stdout:
                return candidate
    pytest.skip('build the native solver with --cnstates-exempt support or set HATCHET_TEST_CPP_SOLVER')


def _require_native_license():
    if not solver_available('cpp'):
        pytest.skip('GRB_LICENSE_FILE is not set; load the Gurobi module to run native solves')


def _write_native_inputs(tmp_path):
    """Cluster 1 is the clonal diploid cluster, cluster 2 holds (2,1), (1,1) and (1,0) in pure samples Tumor1-3,
    and cluster 3 is a clonal (2,1)."""
    rows = {
        1: [(1.0, 0.5)] * 3,
        2: [(1.5, 1 / 3), (1.0, 0.5), (0.5, 0.0)],
        3: [(1.5, 1 / 3)] * 3,
    }
    prefix = tmp_path / 'native-input'
    seg = ['#ID\tSAMPLE\t#BINS\tRD\t#SNPS\tCOV\tALPHA\tBETA\tBAF']
    bbc = ['#CHR\tSTART\tEND\tSAMPLE\tRD\t#SNPS\tCOV\tALPHA\tBETA\tBAF\tCLUSTER']
    for cluster, values in rows.items():
        for sample, (rd, baf) in enumerate(values, start=1):
            seg.append(f'{cluster}\tTumor{sample}\t1\t{rd}\t20\t30\t10\t10\t{baf}')
            start = 100 * (cluster - 1)
            bbc.append(f'chr1\t{start}\t{start + 100}\tTumor{sample}\t{rd}\t20\t30\t10\t10\t{baf}\t{cluster}')
    prefix.with_suffix('.seg').write_text('\n'.join(seg) + '\n')
    prefix.with_suffix('.bbc').write_text('\n'.join(bbc) + '\n')
    return prefix


def _native_command(solver, prefix, tmp_path, *extra, mode=1, name='out'):
    return [
        str(solver),
        str(prefix),
        '-n',
        '4',
        '-c',
        '1:1:1',
        '-e',
        '6',
        '-j',
        '1',
        '-p',
        '2',
        '-i',
        '2',
        '-u',
        '0.1',
        '-s',
        '60',
        '-r',
        '7',
        '-M',
        str(mode),
        '-v',
        '1',
        '-o',
        str(tmp_path / f'native-{name}-m{mode}'),
        *extra,
    ]


def _native_states(tmp_path, mode, name='out'):
    """Distinct tumor (A, B) states per cluster in the native solver's output."""
    bbc = pd.read_table(tmp_path / f'native-{name}-m{mode}.bbc.ucn.tsv')
    clones = [column for column in bbc.columns if column.startswith('cn_clone')]
    first = bbc.groupby('CLUSTER')[clones].first()
    return {cluster: set(row.values) for cluster, row in first.iterrows()}


def _native_objective(result):
    values = re.findall(r'Final objective:\s*([0-9.eE+-]+)', result.stdout + result.stderr)
    assert values, result.stdout + result.stderr
    return float(values[-1])


@pytest.mark.parametrize(
    ('extra', 'message'),
    [
        (['--cnstates-exempt', '2'], 'require a global -d'),
        (['-d', '1', '--cnstates-exempt', '99'], 'cluster with ID 99 has not been found'),
        (['--cluster-cnstates', '99:2'], 'cluster with ID 99 has not been found'),
        (['-d', '1', '--cnstates-exempt', '2', '--cluster-cnstates', '3:2'], 'only one of'),
        (['--cluster-cnstates', '2:0'], 'positive integer'),
        (['--cluster-cnstates', '2'], 'ID:INT'),
        (['--cluster-cnstates', '2:2,2:3'], 'more than once'),
        (['-d', '1', '--cnstates-formulation', 'other'], 'pairwise or slot'),
    ],
)
def test_native_solver_rejects_invalid_cluster_limits(tmp_path, native_solver, extra, message):
    prefix = _write_native_inputs(tmp_path)

    result = subprocess.run(_native_command(native_solver, prefix, tmp_path, *extra), text=True, capture_output=True)

    assert result.returncode != 0
    assert message in result.stderr + result.stdout


@pytest.mark.parametrize('mode', (0, 1, 2), ids=('cd-plus-ilp', 'ilp', 'cd'))
def test_native_solver_enforces_the_global_limit(tmp_path, native_solver, mode):
    _require_native_license()
    prefix = _write_native_inputs(tmp_path)

    subprocess.run(_native_command(native_solver, prefix, tmp_path, '-d', '1', mode=mode), check=True)

    assert all(len(states) == 1 for states in _native_states(tmp_path, mode).values())


def test_native_solver_frees_exempt_and_bounds_per_cluster_limits(tmp_path, native_solver):
    _require_native_license()
    prefix = _write_native_inputs(tmp_path)

    subprocess.run(
        _native_command(native_solver, prefix, tmp_path, '-d', '1', '--cnstates-exempt', '2', name='exempt'),
        check=True,
    )
    states = _native_states(tmp_path, 1, name='exempt')
    assert states[2] == {'2|1', '1|1', '1|0'}
    assert len(states[3]) == 1

    subprocess.run(
        _native_command(native_solver, prefix, tmp_path, '--cluster-cnstates', '2:2', name='map'),
        check=True,
    )
    assert len(_native_states(tmp_path, 1, name='map')[2]) == 2


def test_native_pairwise_and_slot_formulations_reach_the_same_optimum(tmp_path, native_solver):
    _require_native_license()
    prefix = _write_native_inputs(tmp_path)

    objectives = {}
    for formulation in ('pairwise', 'slot'):
        result = subprocess.run(
            _native_command(
                native_solver, prefix, tmp_path, '-d', '2', '--cnstates-formulation', formulation, name=formulation
            ),
            text=True,
            capture_output=True,
            check=True,
        )
        objectives[formulation] = _native_objective(result)

    assert objectives['slot'] > 0.01  # cluster 2 needs three states
    assert objectives['pairwise'] == pytest.approx(objectives['slot'], abs=1e-4)


def test_native_coordinate_descent_reports_model_sizes(tmp_path, native_solver):
    _require_native_license()
    prefix = _write_native_inputs(tmp_path)
    command = _native_command(native_solver, prefix, tmp_path, '-d', '1', mode=2)
    command[command.index('-v') + 1] = '2'

    subprocess.run(command, check=True, cwd=tmp_path, capture_output=True)

    steps = pd.read_csv(tmp_path / 'coordinate_descence.tsv', sep='\t', header=None)
    assert steps.shape[1] == 9
    carch = steps[steps[2] == 'C']
    assert (carch[6] > 0).all() and (carch[7] > 0).all() and (carch[8] > 0).all()
