import os
import subprocess
import pytest
import pandas as pd
from pandas.testing import assert_frame_equal

import hatchet
from hatchet.bin.HATCHet import main as main
from hatchet import config
from hatchet.utils.fixed_matrix import parse_fixed_matrix
from hatchet.utils.solve import solver_available

this_dir = os.path.dirname(__file__)
DATA_FOLDER = os.path.join(this_dir, 'data', 'fw')
SOLVE = os.getenv('HATCHET_TEST_CPP_SOLVER', os.path.join(os.path.dirname(hatchet.__file__), 'solve'))
FIXED_MATRIX_SCENARIOS = (
    pytest.param(
        os.path.join(this_dir, 'data', 'fw', 'bulk'),
        os.path.join(this_dir, 'data', 'fw', 'fixed-matrix.tsv'),
        id='fixed-width',
    ),
    pytest.param(
        os.path.join(this_dir, 'data', 'vw', 'bbc', 'bulk'),
        os.path.join(this_dir, 'data', 'vw', 'fixed-matrix.tsv'),
        id='variable-width',
    ),
)
FIXED_MATRIX_SYMMETRY_SCENARIOS = (
    pytest.param(
        os.path.join(this_dir, 'data', 'fw', 'bulk'),
        os.path.join(this_dir, 'data', 'fw', 'fixed-matrix-n3.tsv'),
        id='fixed-width',
    ),
    pytest.param(
        os.path.join(this_dir, 'data', 'vw', 'bbc', 'bulk'),
        os.path.join(this_dir, 'data', 'vw', 'fixed-matrix-n3.tsv'),
        id='variable-width',
    ),
)
SOLVE_MODES = (
    pytest.param(0, id='cd-plus-ilp'),
    pytest.param(1, id='ilp'),
    pytest.param(2, id='cd'),
)


@pytest.fixture(scope='module')
def output_folder():
    out = os.path.join(this_dir, 'out', 'solver')
    os.makedirs(out, exist_ok=True)
    return out


@pytest.fixture(scope='session')
def gurobi_runtime_available():
    environment = None
    model = None
    try:
        import gurobipy as gp

        environment = gp.Env(empty=True)
        environment.setParam('OutputFlag', 0)
        environment.start()
        model = gp.Model(env=environment)
        variable = model.addVar(lb=0, ub=1)
        model.setObjective(variable, gp.GRB.MINIMIZE)
        model.optimize()
        return model.Status == gp.GRB.OPTIMAL
    except Exception:
        return False
    finally:
        if model is not None:
            model.dispose()
        if environment is not None:
            environment.dispose()


def require_solver(solver, gurobi_runtime_available):
    if not gurobi_runtime_available:
        pytest.skip('a working Gurobi runtime license is unavailable')
    if solver == 'cpp' and not os.path.isfile(SOLVE):
        pytest.skip('compiled solver is unavailable')
    if solver != 'cpp' and not solver_available(solver):
        pytest.skip(f'{solver} solver is unavailable')


def assert_fixed_matrix_entries(result_path, fixed_path):
    result = pd.read_csv(result_path, sep='\t')
    fixed = parse_fixed_matrix(fixed_path)

    for entry in fixed.entries:
        if entry.matrix in ('A', 'B'):
            rows = result[result['CLUSTER'].astype(str) == entry.row]
            column = 'cn_normal' if entry.column == 'normal' else f'cn_{entry.column}'
            allele_index = 0 if entry.matrix == 'A' else 1
            observed = rows[column].astype(str).str.split('|').str[allele_index].astype(int)
            assert not observed.empty
            assert (observed == entry.value).all()
        else:
            rows = result[result['SAMPLE'].astype(str) == entry.column]
            column = 'u_normal' if entry.row == 'normal' else f'u_{entry.row}'
            assert not rows.empty
            assert (rows[column] - entry.value).abs().max() <= 0.001


def assert_solve_output(result_path, input_prefix, fixed_path, clones):
    result = pd.read_csv(result_path, sep='\t')
    input_columns = pd.read_csv(input_prefix + '.bbc', sep='\t', nrows=0).columns.tolist()
    expected_extra = []
    for clone in range(clones):
        label = 'normal' if clone == 0 else f'clone{clone}'
        expected_extra.extend((f'cn_{label}', f'u_{label}'))

    assert result.columns.tolist() == input_columns + expected_extra
    assert (result['cn_normal'] == '1|1').all()
    clone_labels = ['normal'] + [f'clone{i}' for i in range(1, clones)]
    proportion_columns = [f'u_{label}' for label in clone_labels]
    proportions = result.groupby('SAMPLE', sort=False)[proportion_columns].first()
    assert ((proportions.sum(axis=1) - 1).abs() <= 0.001).all()
    assert_fixed_matrix_entries(result_path, fixed_path)


def assert_canonical_clone_order(result_path, seg_path, clones):
    result = pd.read_csv(result_path, sep='\t')
    cluster_order = pd.read_csv(seg_path, sep='\t')['#ID'].astype(str).drop_duplicates().tolist()
    burdens = []
    for clone in range(1, clones):
        burden = 0
        for cluster_index, cluster_id in enumerate(cluster_order):
            states = result.loc[result['CLUSTER'].astype(str) == cluster_id, f'cn_clone{clone}'].unique()
            assert len(states) == 1
            allele_a, allele_b = (int(value) for value in str(states[0]).split('|'))
            burden += (cluster_index + 1) * (allele_a + allele_b)
        burdens.append(burden)
    assert burdens == sorted(burdens)


def run_fixed_matrix_solve(output_folder, input_prefix, fixed_path, solver, mode):
    old_value = config.compute_cn.solver
    config.compute_cn.solver = solver
    try:
        args = [SOLVE] if solver == 'cpp' else []
        args.extend(
            [
                '-x',
                str(output_folder),
                '-i',
                input_prefix,
                '-n{}'.format(parse_fixed_matrix(fixed_path).clones),
                '-p',
                '1',
                '-v',
                '0',
                '-u',
                '0.03',
                '--mode',
                str(mode),
                '--maxiterations',
                '2',
                '-s',
                '60',
                '-r',
                '6700',
                '-j',
                '1',
                '-eD',
                '6',
                '-td',
                '0.1',
                '--diploid',
                '--fixed-matrix',
                fixed_path,
            ]
        )
        main(args=args)
    finally:
        config.compute_cn.solver = old_value

    return output_folder / 'best.bbc.ucn'


@pytest.mark.parametrize(('input_prefix', 'fixed_path'), FIXED_MATRIX_SCENARIOS)
@pytest.mark.parametrize('solver', ('cpp', 'gurobi'))
@pytest.mark.parametrize('mode', SOLVE_MODES)
def test_fixed_matrix_solve_binning_scenario(
    tmp_path,
    input_prefix,
    fixed_path,
    solver,
    mode,
    gurobi_runtime_available,
):
    require_solver(solver, gurobi_runtime_available)
    result_path = run_fixed_matrix_solve(tmp_path, input_prefix, fixed_path, solver, mode)
    assert_solve_output(result_path, input_prefix, fixed_path, clones=2)


@pytest.mark.parametrize(('input_prefix', 'fixed_path'), FIXED_MATRIX_SYMMETRY_SCENARIOS)
@pytest.mark.parametrize('solver', ('cpp', 'gurobi'))
def test_fixed_matrix_solve_preserves_canonical_clone_order(
    tmp_path,
    input_prefix,
    fixed_path,
    solver,
    gurobi_runtime_available,
):
    require_solver(solver, gurobi_runtime_available)
    result_path = run_fixed_matrix_solve(tmp_path, input_prefix, fixed_path, solver, mode=1)
    assert_solve_output(result_path, input_prefix, fixed_path, clones=3)
    assert_canonical_clone_order(result_path, input_prefix + '.seg', clones=3)


@pytest.mark.skipif(os.getenv('GRB_LICENSE_FILE') is None, reason='No license for Gurobi found')
def test_solve_binary():
    cmd = [SOLVE, os.path.join(DATA_FOLDER, 'bulk')]
    # -M 1 for ILP solve mode
    cmd.extend([_x for _x in '-f -e 6 -p 400 -u 0.03 -r 6700 -M 1 -v 2 -c 5:1:1 -n 2'.split()])

    print(' '.join(cmd))

    p = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    p.communicate()
    assert p.returncode == 0


@pytest.mark.skipif(os.getenv('GRB_LICENSE_FILE') is None, reason='No license for Gurobi found')
def test_solve_command_cpp(output_folder):
    old_value = config.compute_cn.solver
    config.compute_cn.solver = 'cpp'

    main(
        args=[
            SOLVE,
            '-x',
            os.path.join(output_folder),
            '-i',
            os.path.join(this_dir, 'data/fw/bulk'),
            '-n2',
            '-p',
            '400',
            '-v',
            '3',
            '-u',
            '0.03',
            '--mode',
            '1',
            '-r',
            '6700',
            '-j',
            '8',
            '-eD',
            '6',
            '-eT',
            '12',
            '-g',
            '0.35',
            '-l',
            '0.6',
        ]
    )

    df1 = pd.read_csv(os.path.join(output_folder, 'best.bbc.ucn'), sep='\t')
    df2 = pd.read_csv(os.path.join(this_dir, 'data', 'fw', 'best.bbc.ucn'), sep='\t')
    assert_frame_equal(df1, df2)

    config.compute_cn.solver = old_value


@pytest.mark.skipif(
    not solver_available('gurobipy'),
    reason='gurobipy solver not available for pyomo',
)
def test_solve_command_gurobipy(output_folder):
    old_value = config.compute_cn.solver
    config.compute_cn.solver = 'gurobipy'

    main(
        args=[
            '-x',
            os.path.join(output_folder),
            '-i',
            os.path.join(this_dir, 'data/fw/bulk'),
            '-n2',
            '-p',
            '400',
            '-v',
            '3',
            '-u',
            '0.03',
            '--mode',
            '1',
            '-r',
            '6700',
            '-j',
            '8',
            '-eD',
            '6',
            '-eT',
            '12',
            '-g',
            '0.35',
            '-l',
            '0.6',
        ]
    )

    df1 = pd.read_csv(os.path.join(output_folder, 'best.bbc.ucn'), sep='\t')
    df2 = pd.read_csv(os.path.join(this_dir, 'data', 'fw', 'best.bbc.ucn'), sep='\t')
    assert_frame_equal(df1, df2)

    config.compute_cn.solver = old_value


@pytest.mark.skipif(
    not solver_available('gurobi'),
    reason='gurobi solver not available for pyomo',
)
def test_solve_command_gurobi(output_folder):
    old_value = config.compute_cn.solver
    config.compute_cn.solver = 'gurobi'

    main(
        args=[
            '-x',
            os.path.join(output_folder),
            '-i',
            os.path.join(this_dir, 'data/fw/bulk'),
            '-n2',
            '-p',
            '400',
            '-v',
            '3',
            '-u',
            '0.03',
            '--mode',
            '1',
            '-r',
            '6700',
            '-j',
            '8',
            '-eD',
            '6',
            '-eT',
            '12',
            '-g',
            '0.35',
            '-l',
            '0.6',
        ]
    )

    df1 = pd.read_csv(os.path.join(output_folder, 'best.bbc.ucn'), sep='\t')
    df2 = pd.read_csv(os.path.join(this_dir, 'data', 'fw', 'best.bbc.ucn'), sep='\t')
    assert_frame_equal(df1, df2)

    config.compute_cn.solver = old_value


@pytest.mark.skipif(not solver_available('cbc'), reason='cbc solver not available for pyomo')
def test_solve_command_cbc(output_folder):
    old_value = config.compute_cn.solver
    config.compute_cn.solver = 'cbc'

    main(
        args=[
            '-x',
            os.path.join(output_folder),
            '-i',
            os.path.join(this_dir, 'data/fw/bulk'),
            '-n2',
            '-p',
            '400',
            '-v',
            '3',
            '-u',
            '0.03',
            '--mode',
            '1',
            '-r',
            '6700',
            '-j',
            '8',
            '-eD',
            '6',
            '-eT',
            '12',
            '-g',
            '0.35',
            '-l',
            '0.6',
        ]
    )

    df1 = pd.read_csv(os.path.join(output_folder, 'best.bbc.ucn'), sep='\t')
    df2 = pd.read_csv(os.path.join(this_dir, 'data', 'fw', 'best.bbc.ucn'), sep='\t')
    assert_frame_equal(df1, df2)

    config.compute_cn.solver = old_value
