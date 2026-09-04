from dataclasses import FrozenInstanceError
from pathlib import Path
import os
import re
import shutil
import subprocess

import pandas as pd
import pytest
from pyomo.environ import value

from hatchet import config
from hatchet.bin.HATCHet import (
    execute,
    execute_python,
    load_sample_clonal_entries,
    main as compute_cn,
    parsing_arguments,
    runningDiploid,
    runningTetraploid,
)
from hatchet.utils.solve import solve
from hatchet.utils.solve.ilp_subset import ILPSubset, ILPSubsetSplit
from hatchet.utils.fixed_matrix import (
    FixedMatrixFormatError,
    parse_fixed_matrix,
    resolve_fixed_matrix,
    validate_fixed_matrix,
)
from hatchet.utils.sample_clonal import (
    FORMAT_MARKER,
    SampleClonalFormatError,
    derive_sample_clonal_calibration,
    parse_sample_clonal_file,
    parse_sample_clonal_inline,
    resolve_sample_clonal,
    validate_sample_clonal_constraints,
    verify_sample_clonal_solution,
)


def _write_sample_clonal(tmp_path, contents):
    path = tmp_path / 'sample-clonal.tsv'
    path.write_text(contents)
    return path


def test_inline_and_file_specs_resolve_identically(tmp_path):
    path = _write_sample_clonal(
        tmp_path,
        f"""{FORMAT_MARKER}
# comment
sample\tcluster\tcn_a\tcn_b
Tumor 1\tcluster-a\t1\t1
Tumor 1\tcluster-b\t2\t1
Tumor 2\tcluster-b\t3\t1
Tumor 2\tcluster-c\t1\t0
""",
    )
    from_file = parse_sample_clonal_file(path)
    inline = parse_sample_clonal_inline('Tumor 1=cluster-a:1:1,cluster-b:2:1;Tumor 2=cluster-b:3:1,cluster-c:1:0')

    cluster_ids = ('cluster-a', 'cluster-b', 'cluster-c')
    sample_ids = ('Tumor 1', 'Tumor 2')
    resolved_file = resolve_sample_clonal(from_file, cluster_ids, sample_ids)
    resolved_inline = resolve_sample_clonal(inline, cluster_ids, sample_ids)

    values = lambda resolved: tuple(
        (entry.sample_index, entry.cluster_index, entry.cn_a, entry.cn_b) for entry in resolved.entries
    )
    assert values(resolved_file) == values(resolved_inline)
    with pytest.raises(FrozenInstanceError):
        from_file.entries[0].cn_a = 2


@pytest.mark.parametrize(
    ('contents', 'message'),
    [
        ('', 'expected format marker'),
        (f'{FORMAT_MARKER}\nwrong\theader\n', 'expected tab-separated header'),
        (f'{FORMAT_MARKER}\nsample\tcluster\tcn_a\tcn_b\n', 'at least one entry'),
        (
            f'{FORMAT_MARKER}\nsample\tcluster\tcn_a\tcn_b\ns1\tc1\t1\n',
            'expected 4 tab-separated fields',
        ),
        (
            f'{FORMAT_MARKER}\nsample\tcluster\tcn_a\tcn_b\ns1\tc1\t1.5\t1\ns1\tc2\t2\t1\n',
            'must be an integer',
        ),
        (
            f'{FORMAT_MARKER}\nsample\tcluster\tcn_a\tcn_b\ns1\tc1\t-1\t1\ns1\tc2\t2\t1\n',
            'must be nonnegative',
        ),
        (
            f'{FORMAT_MARKER}\nsample\tcluster\tcn_a\tcn_b\ns1\tc1\t1\t1\n',
            'must have exactly two entries',
        ),
        (
            f'{FORMAT_MARKER}\nsample\tcluster\tcn_a\tcn_b\n' 's1\tc1\t1\t1\ns1\tc1\t2\t1\n',
            'duplicate sample-clonal entry',
        ),
        (
            f'{FORMAT_MARKER}\nsample\tcluster\tcn_a\tcn_b\n' 's1\tc1\t1\t1\ns1\tc2\t2\t0\n',
            'distinct total copy-number states',
        ),
    ],
)
def test_file_parser_rejects_invalid_specs(tmp_path, contents, message):
    path = _write_sample_clonal(tmp_path, contents)

    with pytest.raises(SampleClonalFormatError, match=re.escape(message)) as exc_info:
        parse_sample_clonal_file(path)

    assert str(path) in str(exc_info.value)


@pytest.mark.parametrize(
    ('value', 'message'),
    [
        ('', 'empty sample specification'),
        ('sample:c1:1:1,c2:2:1', 'expected SAMPLE='),
        ('sample=c1:1:1', 'exactly two entries'),
        ('sample=c1:1:1,c2:2', 'expected each state'),
        ('sample=c1:1:1,c1:2:1', 'duplicate sample-clonal entry'),
        ('sample=c1:1:1,c2:2:0', 'distinct total copy-number states'),
    ],
)
def test_inline_parser_rejects_invalid_specs(value, message):
    with pytest.raises(SampleClonalFormatError, match=re.escape(message)):
        parse_sample_clonal_inline(value)


def test_resolution_requires_exact_complete_seg_identifiers():
    parsed = parse_sample_clonal_inline('Tumor=1:1:1,2:2:1')

    resolved = resolve_sample_clonal(parsed, (1, 2), ('Tumor',))
    assert tuple((entry.cluster_index, entry.sample_index) for entry in resolved.entries) == ((0, 0), (1, 0))

    with pytest.raises(SampleClonalFormatError, match='unknown SEG cluster identifier'):
        resolve_sample_clonal(parsed, (1, 3), ('Tumor',))
    with pytest.raises(SampleClonalFormatError, match='unknown SEG sample identifier'):
        resolve_sample_clonal(parsed, (1, 2), ('Other',))
    with pytest.raises(SampleClonalFormatError, match='missing sample-clonal entries'):
        resolve_sample_clonal(parsed, (1, 2), ('Tumor', 'Other'))


def _derive(value, rdr):
    parsed = parse_sample_clonal_inline(value)
    resolved = resolve_sample_clonal(parsed, rdr.index, rdr.columns)
    return derive_sample_clonal_calibration(resolved, rdr)


def test_calibration_uses_general_theorem_and_is_pair_order_invariant():
    rdr = pd.DataFrame({'Tumor': [2.0 / 3.0, 2.5 / 3.0]}, index=('c1', 'c2'))

    forward = _derive('Tumor=c1:1:1,c2:2:1', rdr)
    reverse = _derive('Tumor=c2:2:1,c1:1:1', rdr)

    assert forward.gamma_by_sample() == pytest.approx({'Tumor': 3.0})
    assert forward.purity_by_sample() == pytest.approx({'Tumor': 0.5})
    assert reverse.gamma == pytest.approx(forward.gamma)
    assert reverse.purity == pytest.approx(forward.purity)


def test_calibration_supports_two_nonbaseline_states():
    rdr = pd.DataFrame({'Tumor': [2.6 / 2.4, 3.2 / 2.4]}, index=('c1', 'c2'))

    calibration = _derive('Tumor=c1:2:1,c2:3:1', rdr)

    assert calibration.gamma == pytest.approx((2.4,))
    assert calibration.purity == pytest.approx((0.6,))


def test_calibration_rejects_singular_or_invalid_results():
    singular_rdr = pd.DataFrame({'Tumor': [1.0, 2.0]}, index=('c1', 'c2'))
    with pytest.raises(SampleClonalFormatError, match='singular'):
        _derive('Tumor=c1:2:1,c2:3:1', singular_rdr)

    invalid_rdr = pd.DataFrame({'Tumor': [1.0, 2.0]}, index=('c1', 'c2'))
    with pytest.raises(SampleClonalFormatError, match='invalid purity'):
        _derive('Tumor=c1:1:1,c2:2:1', invalid_rdr)


def _compute_cn_input(tmp_path):
    prefix = tmp_path / 'bulk'
    prefix.with_suffix('.seg').touch()
    prefix.with_suffix('.bbc').touch()
    return str(prefix)


def test_compute_cn_accepts_inline_or_file_sample_clonal(tmp_path, monkeypatch):
    monkeypatch.setattr(config.compute_cn, 'solver', 'gurobi')
    input_prefix = _compute_cn_input(tmp_path)
    sample_clonal_path = _write_sample_clonal(
        tmp_path,
        f'{FORMAT_MARKER}\nsample\tcluster\tcn_a\tcn_b\nTumor\t1\t1\t1\nTumor\t2\t2\t1\n',
    )

    inline = parsing_arguments(
        ['-i', input_prefix, '-x', str(tmp_path), '--diploid', '--sample-clonal', 'Tumor=1:1:1,2:2:1']
    )
    from_file = parsing_arguments(
        ['-i', input_prefix, '-x', str(tmp_path), '--tetraploid', '--sample-clonal-file', str(sample_clonal_path)]
    )

    assert inline['sampleclonal'] == 'Tumor=1:1:1,2:2:1'
    assert inline['sampleclonalfile'] is None
    assert from_file['sampleclonal'] is None
    assert from_file['sampleclonalfile'] == str(sample_clonal_path)


def test_compute_cn_uses_sample_clonal_ini_defaults(tmp_path, monkeypatch):
    monkeypatch.setattr(config.compute_cn, 'solver', 'gurobi')
    monkeypatch.setattr(config.compute_cn, 'sampleclonal', 'Tumor=1:1:1,2:2:1')
    monkeypatch.setattr(config.compute_cn, 'diploid', True)

    parsed = parsing_arguments(['-i', _compute_cn_input(tmp_path), '-x', str(tmp_path)])

    assert parsed['sampleclonal'] == 'Tumor=1:1:1,2:2:1'


@pytest.mark.parametrize(
    ('extra_args', 'message'),
    [
        (['--sample-clonal', 'Tumor=1:1:1,2:2:1'], 'requires either --diploid or --tetraploid'),
        (
            ['--diploid', '--tetraploid', '--sample-clonal', 'Tumor=1:1:1,2:2:1'],
            'Diploid and tetraploid cannot be forced at the same time',
        ),
        (
            ['--diploid', '--sample-clonal', 'Tumor=1:1:1,2:2:1', '--clonal', '1:1:1'],
            'cannot be used with --clonal',
        ),
        (
            ['--diploid', '--sample-clonal', 'Tumor=1:1:1,2:2:1', '--sample-clonal-file', 'missing.tsv'],
            'Sample-clonal file .* not found',
        ),
    ],
)
def test_compute_cn_rejects_invalid_sample_clonal_option_combinations(tmp_path, monkeypatch, extra_args, message):
    monkeypatch.setattr(config.compute_cn, 'solver', 'gurobi')

    with pytest.raises(ValueError, match=message):
        parsing_arguments(['-i', _compute_cn_input(tmp_path), '-x', str(tmp_path)] + extra_args)


def test_compute_cn_rejects_both_sample_clonal_input_forms(tmp_path, monkeypatch):
    monkeypatch.setattr(config.compute_cn, 'solver', 'gurobi')
    sample_clonal_path = _write_sample_clonal(
        tmp_path,
        f'{FORMAT_MARKER}\nsample\tcluster\tcn_a\tcn_b\nTumor\t1\t1\t1\nTumor\t2\t2\t1\n',
    )

    with pytest.raises(ValueError, match='cannot be used together'):
        parsing_arguments(
            [
                '-i',
                _compute_cn_input(tmp_path),
                '-x',
                str(tmp_path),
                '--diploid',
                '--sample-clonal',
                'Tumor=1:1:1,2:2:1',
                '--sample-clonal-file',
                str(sample_clonal_path),
            ]
        )


def test_load_sample_clonal_entries_resolves_both_forms(tmp_path):
    inline_args = {
        'sampleclonal': 'Tumor=1:1:1,2:2:1',
        'sampleclonalfile': None,
    }
    resolved_inline = load_sample_clonal_entries(inline_args, (1, 2), ('Tumor',))
    assert tuple(entry.cluster_index for entry in resolved_inline.entries) == (0, 1)

    path = _write_sample_clonal(
        tmp_path,
        f'{FORMAT_MARKER}\nsample\tcluster\tcn_a\tcn_b\nTumor\t1\t1\t1\nTumor\t2\t2\t1\n',
    )
    file_args = {'sampleclonal': None, 'sampleclonalfile': str(path)}
    resolved_file = load_sample_clonal_entries(file_args, (1, 2), ('Tumor',))
    assert tuple(entry.cluster_index for entry in resolved_file.entries) == (0, 1)

    assert load_sample_clonal_entries({'sampleclonal': None, 'sampleclonalfile': None}, (), ()) is None


def _write_solver_inputs(tmp_path):
    prefix = tmp_path / 'solver-input'
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
    return prefix


def _write_conflicting_solver_inputs(tmp_path):
    prefix = tmp_path / 'conflicting-input'
    prefix.with_suffix('.seg').write_text(
        '#ID\tSAMPLE\t#BINS\tRD\t#SNPS\tCOV\tALPHA\tBETA\tBAF\n'
        '1\tTumorA\t1\t0.6666666666666666\t20\t30\t10\t10\t0.5\n'
        '1\tTumorB\t1\t0.8333333333333334\t20\t30\t12\t8\t0.4\n'
        '2\tTumorA\t1\t0.8333333333333334\t20\t30\t12\t8\t0.4\n'
        '2\tTumorB\t1\t0.6666666666666666\t20\t30\t10\t10\t0.5\n'
        '3\tTumorA\t1\t0.6666666666666666\t20\t30\t10\t10\t0.5\n'
        '3\tTumorB\t1\t0.6666666666666666\t20\t30\t10\t10\t0.5\n'
    )
    prefix.with_suffix('.bbc').write_text(
        '#CHR\tSTART\tEND\tSAMPLE\tRD\t#SNPS\tCOV\tALPHA\tBETA\tBAF\tCLUSTER\n'
        'chr1\t0\t100\tTumorA\t0.6666666666666666\t20\t30\t10\t10\t0.5\t1\n'
        'chr1\t0\t100\tTumorB\t0.8333333333333334\t20\t30\t12\t8\t0.4\t1\n'
        'chr1\t100\t200\tTumorA\t0.8333333333333334\t20\t30\t12\t8\t0.4\t2\n'
        'chr1\t100\t200\tTumorB\t0.6666666666666666\t20\t30\t10\t10\t0.5\t2\n'
        'chr1\t200\t300\tTumorA\t0.6666666666666666\t20\t30\t10\t10\t0.5\t3\n'
        'chr1\t200\t300\tTumorB\t0.6666666666666666\t20\t30\t10\t10\t0.5\t3\n'
    )
    return prefix


def test_python_solver_uses_sample_clonal_scaling_and_explicit_base(tmp_path, monkeypatch):
    captured = {}

    class CapturingILP:
        def __init__(self, n, cn_max, **kwargs):
            captured.update(n=n, cn_max=cn_max, **kwargs)

        def create_model(self, pprint=False):
            captured['pprint'] = pprint

        def run(self, solver_type, timelimit=None):
            captured['solver_type'] = solver_type
            return None

    monkeypatch.setattr('hatchet.utils.solve.solver_available', lambda solver: True)
    monkeypatch.setattr('hatchet.utils.solve.ILPSubset', CapturingILP)
    prefix = _write_solver_inputs(tmp_path)
    sample_clonal = parse_sample_clonal_inline('Tumor=1:1:1,2:2:1')

    solve(
        clonal=None,
        bbc_file=str(prefix.with_suffix('.bbc')),
        seg_file=str(prefix.with_suffix('.seg')),
        n=2,
        solver='gurobi',
        solve_mode='ilp',
        cn_max=6,
        mu=0.1,
        base=1,
        sample_clonal=sample_clonal,
    )

    assert captured['copy_numbers'] == {}
    assert captured['purities'] == pytest.approx([0.5])
    assert captured['base'] == 1
    assert captured['sample_clonal'].sample_ids == ('Tumor',)
    assert captured['f_a'].loc[1, 'Tumor'] == pytest.approx(1.0)
    assert captured['f_b'].loc[1, 'Tumor'] == pytest.approx(1.0)
    assert captured['f_a'].loc[2, 'Tumor'] == pytest.approx(1.5)
    assert captured['f_b'].loc[2, 'Tumor'] == pytest.approx(1.0)


def test_python_solver_checks_explicit_base_and_purity_assertion(tmp_path, monkeypatch):
    monkeypatch.setattr('hatchet.utils.solve.solver_available', lambda solver: True)
    prefix = _write_solver_inputs(tmp_path)
    sample_clonal = parse_sample_clonal_inline('Tumor=1:1:1,2:2:1')
    common = {
        'clonal': None,
        'bbc_file': str(prefix.with_suffix('.bbc')),
        'seg_file': str(prefix.with_suffix('.seg')),
        'n': 2,
        'solver': 'gurobi',
        'solve_mode': 'ilp',
        'cn_max': 6,
        'sample_clonal': sample_clonal,
    }

    with pytest.raises(ValueError, match='explicit factorization base'):
        solve(**common)
    with pytest.raises(ValueError, match='does not match sample-clonal purity'):
        solve(**common, base=1, purities=[0.8])
    conflicting = dict(common, clonal='1:1:1')
    with pytest.raises(ValueError, match='cannot be combined with global clonal'):
        solve(**conflicting, base=1)


@pytest.mark.parametrize('provided_purity', (0.5, 0.45, 0.55))
def test_python_solver_accepts_matching_sample_clonal_purity_within_tolerance(
    tmp_path,
    monkeypatch,
    provided_purity,
):
    captured = {}

    class CapturingILP:
        def __init__(self, n, cn_max, **kwargs):
            captured.update(kwargs)

        def create_model(self, pprint=False):
            pass

        def run(self, solver_type, timelimit=None):
            return None

    monkeypatch.setattr('hatchet.utils.solve.solver_available', lambda solver: True)
    monkeypatch.setattr('hatchet.utils.solve.ILPSubset', CapturingILP)
    prefix = _write_solver_inputs(tmp_path)

    solve(
        clonal=None,
        bbc_file=str(prefix.with_suffix('.bbc')),
        seg_file=str(prefix.with_suffix('.seg')),
        n=2,
        solver='gurobi',
        solve_mode='ilp',
        cn_max=6,
        base=1,
        purities=[provided_purity],
        sample_clonal=parse_sample_clonal_inline('Tumor=1:1:1,2:2:1'),
    )

    assert captured['purities'] == pytest.approx([0.5])


@pytest.mark.parametrize('provided_purity', (0.449, 0.551))
def test_python_solver_rejects_sample_clonal_purity_outside_tolerance(tmp_path, monkeypatch, provided_purity):
    monkeypatch.setattr('hatchet.utils.solve.solver_available', lambda solver: True)
    prefix = _write_solver_inputs(tmp_path)

    with pytest.raises(ValueError, match='does not match sample-clonal purity'):
        solve(
            clonal=None,
            bbc_file=str(prefix.with_suffix('.bbc')),
            seg_file=str(prefix.with_suffix('.seg')),
            n=2,
            solver='gurobi',
            solve_mode='ilp',
            cn_max=6,
            base=1,
            purities=[provided_purity],
            sample_clonal=parse_sample_clonal_inline('Tumor=1:1:1,2:2:1'),
        )


@pytest.mark.parametrize('solve_mode', ('ilp', 'cd', 'both'))
def test_python_gurobi_solves_enforced_sample_clonal_fixture(tmp_path, solve_mode):
    from hatchet.utils.solve import solver_available

    if not solver_available('gurobi'):
        pytest.skip('a working Gurobi Pyomo solver is unavailable')
    prefix = _write_solver_inputs(tmp_path)

    result = solve(
        clonal=None,
        bbc_file=str(prefix.with_suffix('.bbc')),
        seg_file=str(prefix.with_suffix('.seg')),
        n=2,
        solver='gurobi',
        solve_mode=solve_mode,
        cn_max=6,
        mu=0.1,
        base=1,
        sample_clonal=parse_sample_clonal_inline('Tumor=1:1:1,2:2:1'),
        n_seed=1,
        n_worker=1,
        random_seed=7,
        max_iters=2,
        timelimit=30,
    )

    objective, c_a, c_b, u, _, sample_ids = result
    assert objective == pytest.approx(0.0, abs=1e-6)
    assert c_a == [[1, 1], [1, 2]]
    assert c_b == [[1, 1], [1, 1]]
    assert sample_ids.tolist() == ['Tumor']
    assert [row[0] for row in u] == pytest.approx([0.5, 0.5])


def _assert_public_sample_clonal_outputs(output_folder):
    bbc_path = output_folder / 'best.bbc.ucn'
    seg_path = output_folder / 'best.seg.ucn'
    assert bbc_path.is_file()
    assert seg_path.is_file()

    bbc = pd.read_csv(bbc_path, sep='\t')
    assert set(bbc['cn_normal']) == {'1|1'}
    assert set(bbc.loc[bbc['CLUSTER'] == 1, 'cn_clone1']) == {'1|1'}
    assert set(bbc.loc[bbc['CLUSTER'] == 2, 'cn_clone1']) == {'2|1'}
    assert (bbc['u_normal'] - 0.5).abs().max() <= 0.001
    assert (bbc['u_clone1'] - 0.5).abs().max() <= 0.001

    seg = pd.read_csv(seg_path, sep='\t')
    assert set(seg['cn_normal']) == {'1|1'}
    assert set(seg['cn_clone1']) == {'1|1', '2|1'}
    assert (seg['u_normal'] - 0.5).abs().max() <= 0.001
    assert (seg['u_clone1'] - 0.5).abs().max() <= 0.001


def _run_public_sample_clonal_solve(tmp_path, monkeypatch, solver, mode, input_form, cpp_solver=None):
    input_prefix = _write_solver_inputs(tmp_path)
    output_folder = tmp_path / 'output'
    output_folder.mkdir()
    monkeypatch.setattr(config.compute_cn, 'solver', solver)

    args = []
    if solver == 'cpp':
        args.append(str(cpp_solver))
    args.extend(
        [
            '-i',
            str(input_prefix),
            '-x',
            str(output_folder),
            '-n',
            '2',
            '-p',
            '1',
            '-v',
            '0',
            '-u',
            '0.1',
            '--mode',
            str(mode),
            '--maxiterations',
            '2',
            '-s',
            '30',
            '-r',
            '7',
            '-j',
            '1',
            '-eD',
            '6',
            '--diploid',
        ]
    )
    if input_form == 'inline':
        args.extend(['--sample-clonal', 'Tumor=1:1:1,2:2:1'])
    else:
        sample_clonal_path = _write_sample_clonal(
            tmp_path,
            f'{FORMAT_MARKER}\nsample\tcluster\tcn_a\tcn_b\n' 'Tumor\t1\t1\t1\nTumor\t2\t2\t1\n',
        )
        args.extend(['--sample-clonal-file', str(sample_clonal_path)])

    compute_cn(args=args)
    _assert_public_sample_clonal_outputs(output_folder)


@pytest.mark.parametrize('mode', (0, 1, 2), ids=('cd-plus-ilp', 'ilp', 'cd'))
@pytest.mark.parametrize('input_form', ('inline', 'file'))
def test_compute_cn_solves_with_sample_clonal_parameter_python(tmp_path, monkeypatch, mode, input_form):
    from hatchet.utils.solve import solver_available

    if not solver_available('gurobi'):
        pytest.skip('a working Gurobi Pyomo solver is unavailable')
    _run_public_sample_clonal_solve(tmp_path, monkeypatch, 'gurobi', mode, input_form)


def _supports_native_sample_clonal(solver):
    result = subprocess.run([str(solver), '-h'], text=True, capture_output=True)
    return result.returncode == 0 and '-S SAMPLE-CLONAL' in result.stdout


@pytest.fixture(scope='session')
def sample_clonal_cpp_solver(tmp_path_factory):
    configured = os.getenv('HATCHET_TEST_CPP_SOLVER')
    candidates = [Path(configured)] if configured else []
    packaged = Path(__file__).parents[1] / 'src' / 'hatchet' / 'solve'
    candidates.append(packaged)
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK) and _supports_native_sample_clonal(candidate):
            return candidate

    gurobi_home = os.getenv('GUROBI_HOME')
    cmake = shutil.which('cmake')
    if not gurobi_home or cmake is None:
        pytest.skip('load the Gurobi module or set HATCHET_TEST_CPP_SOLVER for native sample-clonal solve tests')

    gurobi_root = Path(gurobi_home)
    if gurobi_root.name in ('linux64', 'mac64'):
        gurobi_root = gurobi_root.parent
    repository = Path(__file__).parents[1]
    build_dir = tmp_path_factory.mktemp('sample-clonal-native-build')
    subprocess.run(
        [
            cmake,
            '-S',
            str(repository),
            '-B',
            str(build_dir),
            '-DCMAKE_BUILD_TYPE=Release',
            f'-DGUROBI_HOME={gurobi_root}',
        ],
        check=True,
    )
    subprocess.run([cmake, '--build', str(build_dir), '--config', 'Release'], check=True)
    solver = build_dir / 'solve'
    assert _supports_native_sample_clonal(solver)
    return solver


@pytest.mark.parametrize('mode', (0, 1, 2), ids=('cd-plus-ilp', 'ilp', 'cd'))
@pytest.mark.parametrize('input_form', ('inline', 'file'))
def test_compute_cn_solves_with_sample_clonal_parameter_cpp(
    tmp_path,
    monkeypatch,
    mode,
    input_form,
    sample_clonal_cpp_solver,
):
    _run_public_sample_clonal_solve(
        tmp_path,
        monkeypatch,
        'cpp',
        mode,
        input_form,
        cpp_solver=sample_clonal_cpp_solver,
    )


def _run_conflicting_sample_clonal_solve(tmp_path, monkeypatch, solver, clone_spec, cpp_solver=None):
    input_prefix = _write_conflicting_solver_inputs(tmp_path)
    output_folder = tmp_path / 'output'
    output_folder.mkdir()
    monkeypatch.setattr(config.compute_cn, 'solver', solver)
    args = [str(cpp_solver)] if solver == 'cpp' else []
    args.extend(
        [
            '-i',
            str(input_prefix),
            '-x',
            str(output_folder),
            '-n',
            clone_spec,
            '-p',
            '1',
            '-v',
            '0',
            '-u',
            '0.1',
            '--mode',
            '1',
            '-s',
            '30',
            '-r',
            '7',
            '-j',
            '1',
            '-eD',
            '6',
            '--diploid',
            '--sample-clonal',
            'TumorA=1:1:1,2:2:1;TumorB=1:2:1,3:1:1',
        ]
    )
    compute_cn(args=args)
    return output_folder


@pytest.mark.parametrize('solver', ('gurobi', 'cpp'))
def test_compute_cn_separates_samples_with_conflicting_sample_clonal_states(
    tmp_path,
    monkeypatch,
    solver,
    sample_clonal_cpp_solver,
):
    if solver == 'gurobi':
        from hatchet.utils.solve import solver_available

        if not solver_available('gurobi'):
            pytest.skip('a working Gurobi Pyomo solver is unavailable')
    output_folder = _run_conflicting_sample_clonal_solve(
        tmp_path,
        monkeypatch,
        solver,
        '3',
        cpp_solver=sample_clonal_cpp_solver,
    )
    result = pd.read_csv(output_folder / 'best.bbc.ucn', sep='\t')
    proportions = result.groupby('SAMPLE', sort=False)[['u_clone1', 'u_clone2']].first()
    assert (result.groupby('SAMPLE')['u_normal'].first() - 0.5).abs().max() <= 0.001

    positive_clone = {}
    for sample in ('TumorA', 'TumorB'):
        present = [clone for clone in (1, 2) if proportions.loc[sample, f'u_clone{clone}'] > 0.001]
        assert len(present) == 1
        positive_clone[sample] = present[0]
    assert positive_clone['TumorA'] != positive_clone['TumorB']

    states = result.groupby('CLUSTER', sort=False).first()
    clone_a = positive_clone['TumorA']
    clone_b = positive_clone['TumorB']
    assert states.loc[1, f'cn_clone{clone_a}'] == '1|1'
    assert states.loc[2, f'cn_clone{clone_a}'] == '2|1'
    assert states.loc[1, f'cn_clone{clone_b}'] == '2|1'
    assert states.loc[3, f'cn_clone{clone_b}'] == '1|1'


def test_compute_cn_sample_clonal_supports_clone_count_range(tmp_path, monkeypatch):
    from hatchet.utils.solve import solver_available

    if not solver_available('gurobi'):
        pytest.skip('a working Gurobi Pyomo solver is unavailable')
    input_prefix = _write_solver_inputs(tmp_path)
    output_folder = tmp_path / 'output'
    output_folder.mkdir()
    monkeypatch.setattr(config.compute_cn, 'solver', 'gurobi')

    compute_cn(
        args=[
            '-i',
            str(input_prefix),
            '-x',
            str(output_folder),
            '-n',
            '2,3',
            '-p',
            '1',
            '-v',
            '0',
            '--mode',
            '1',
            '-s',
            '30',
            '-eD',
            '6',
            '--diploid',
            '--sample-clonal',
            'Tumor=1:1:1,2:2:1',
        ]
    )

    assert (output_folder / 'results.diploid.n2.bbc.ucn.tsv').is_file()
    assert (output_folder / 'results.diploid.n3.bbc.ucn.tsv').is_file()
    assert (output_folder / 'best.bbc.ucn').is_file()


def test_compute_cn_matching_purity_is_only_a_check_and_derived_purity_is_used(tmp_path, monkeypatch):
    from hatchet.utils.solve import solver_available

    if not solver_available('gurobi'):
        pytest.skip('a working Gurobi Pyomo solver is unavailable')
    input_prefix = _write_solver_inputs(tmp_path)
    output_folder = tmp_path / 'output'
    output_folder.mkdir()
    monkeypatch.setattr(config.compute_cn, 'solver', 'gurobi')

    compute_cn(
        args=[
            '-i',
            str(input_prefix),
            '-x',
            str(output_folder),
            '-n',
            '2',
            '--mode',
            '1',
            '-eD',
            '6',
            '--diploid',
            '--purities',
            '0.55',
            '--sample-clonal',
            'Tumor=1:1:1,2:2:1',
        ]
    )

    result = pd.read_csv(output_folder / 'best.bbc.ucn', sep='\t')
    assert (result['u_normal'] - 0.5).abs().max() <= 0.001
    assert (result['u_clone1'] - 0.5).abs().max() <= 0.001


def test_compute_cn_sample_clonal_reports_insufficient_clone_count(tmp_path, monkeypatch):
    from hatchet.utils.solve import solver_available

    if not solver_available('gurobi'):
        pytest.skip('a working Gurobi Pyomo solver is unavailable')
    monkeypatch.setattr(config.compute_cn, 'solver', 'gurobi')

    with pytest.raises(RuntimeError, match='(?i)infeasible|feasible solution|no solution'):
        _run_conflicting_sample_clonal_solve(tmp_path, monkeypatch, 'gurobi', '2')


def _runner_args(tmp_path, sample_clonal):
    return {
        'sample_clonal_entries': sample_clonal,
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
        'M': 2,
        'fixedmatrix': None,
        'fixed': None,
        'purities': None,
        'eD': 6,
        'eT': 12,
        'ln': 2,
        'un': 2,
        'v': 0,
        'x': str(tmp_path),
    }


@pytest.mark.parametrize(('runner', 'expected_base'), ((runningDiploid, 1), (runningTetraploid, 2)))
def test_forced_ploidy_runners_bypass_global_clonal_input(tmp_path, monkeypatch, runner, expected_base):
    sample_clonal = parse_sample_clonal_inline('Tumor=1:1:1,2:2:1')
    args = _runner_args(tmp_path, sample_clonal)
    calls = []
    monkeypatch.setattr(
        'hatchet.bin.HATCHet.execute',
        lambda passed_args, basecmd, n, outprefix: calls.append((dict(passed_args), basecmd, n, outprefix)) or 0,
    )

    if runner is runningDiploid:
        runner(neutral=None, args=args)
    else:
        runner(clonal=None, scale=None, size={}, args=args)

    passed_args, basecmd, clone_count, _ = calls[0]
    assert passed_args['base'] == expected_base
    assert passed_args['c'] is None
    assert ' -c ' not in basecmd
    assert clone_count == 2


def test_execute_python_passes_sample_clonal_payload_and_base(tmp_path, monkeypatch):
    sample_clonal = parse_sample_clonal_inline('Tumor=1:1:1,2:2:1')
    args = _runner_args(tmp_path, sample_clonal)
    args.update({'c': None, 'base': 2, 'e': 12, 'bbc': 'bulk.bbc', 'seg': 'bulk.seg', 'binwise': False})
    captured = {}

    def capture_solve(**kwargs):
        captured.update(kwargs)
        return 0.0, [[1, 1]], [[1, 1]], [[0.5], [0.5]], (1,), ('Tumor',)

    monkeypatch.setattr('hatchet.bin.HATCHet.solve', capture_solve)
    monkeypatch.setattr('hatchet.bin.HATCHet.segmentation', lambda *args, **kwargs: None)

    execute_python('gurobi', args, 2, str(tmp_path / 'result'))

    assert captured['clonal'] is None
    assert captured['base'] == 2
    assert captured['sample_clonal'] is sample_clonal


def test_execute_cpp_uses_normalized_temporary_sample_clonal_file(tmp_path, monkeypatch):
    captured = {}

    class CompletedProcess:
        returncode = 0

        def communicate(self):
            return None, 'Final objective:\t4.5\n'

    def capture_popen(command, **kwargs):
        sample_path = Path(command[command.index('-S') + 1])
        captured['command'] = command
        captured['sample_path'] = sample_path
        captured['contents'] = sample_path.read_text()
        return CompletedProcess()

    monkeypatch.setattr(config.compute_cn, 'solver', 'cpp')
    monkeypatch.setattr('hatchet.bin.HATCHet.subprocess.Popen', capture_popen)
    sample_clonal = resolve_sample_clonal(
        parse_sample_clonal_inline('Tumor=1:1:1,2:2:1'),
        (1, 2),
        ('Tumor',),
    )
    args = _runner_args(tmp_path, sample_clonal)
    args['base'] = 2

    objective = execute(args, 'solve bulk', 2, str(tmp_path / 'result'))

    assert objective == 4.5
    assert captured['command'][captured['command'].index('-b') + 1] == '2'
    assert captured['contents'] == (
        f'{FORMAT_MARKER}\nsample\tcluster\tcn_a\tcn_b\n' 'Tumor\t1\t1\t1\nTumor\t2\t2\t1\n'
    )
    assert not captured['sample_path'].exists()


def _sample_clonal_ilp(mu=0.1, split=False):
    cluster_ids = ('c1', 'c2')
    sample_ids = ('Tumor',)
    sample_clonal = resolve_sample_clonal(
        parse_sample_clonal_inline('Tumor=c1:1:1,c2:2:1'),
        cluster_ids,
        sample_ids,
    )
    f_a = pd.DataFrame({'Tumor': [1.2, 1.5]}, index=cluster_ids)
    f_b = pd.DataFrame({'Tumor': [0.8, 1.0]}, index=cluster_ids)
    common = {
        'n': 3,
        'cn_max': 6,
        'd': -1,
        'mu': mu,
        'ampdel': False,
        'copy_numbers': {},
        'f_a': f_a,
        'f_b': f_b,
        'purities': [0.5],
        'base': 1,
        'sample_clonal': sample_clonal,
    }
    if split:
        common.update(
            binsA={cluster: pd.DataFrame({'Tumor': [f_a.loc[cluster, 'Tumor']]}) for cluster in cluster_ids},
            binsB={cluster: pd.DataFrame({'Tumor': [f_b.loc[cluster, 'Tumor']]}) for cluster in cluster_ids},
            lengths={cluster: [1.0] for cluster in cluster_ids},
        )
        return ILPSubsetSplit(**common)
    return ILPSubset(w=pd.Series({'c1': 1.0, 'c2': 1.0}), **common)


def _constraint_satisfied(constraint, tolerance=1e-9):
    body = value(constraint.body)
    if constraint.lower is not None and body < value(constraint.lower) - tolerance:
        return False
    if constraint.upper is not None and body > value(constraint.upper) + tolerance:
        return False
    return True


@pytest.mark.parametrize('split', (False, True))
def test_full_model_conditionally_enforces_sample_clonal_states(split):
    ilp = _sample_clonal_ilp(split=split)
    ilp.create_model()

    assert len(ilp.model.sample_clonal_constraints) == 16
    ilp.model.x_2_1.value = 1
    ilp.model.x_3_1.value = 0
    for cluster_index, (cn_a, cn_b) in enumerate(((1, 1), (2, 1)), start=1):
        getattr(ilp.model, f'cA_{cluster_index}_2').value = cn_a
        getattr(ilp.model, f'cB_{cluster_index}_2').value = cn_b
        getattr(ilp.model, f'cA_{cluster_index}_3').value = 0
        getattr(ilp.model, f'cB_{cluster_index}_3').value = 0
    assert all(_constraint_satisfied(constraint) for constraint in ilp.model.sample_clonal_constraints.values())

    ilp.model.cA_1_2.value = 2
    assert not all(_constraint_satisfied(constraint) for constraint in ilp.model.sample_clonal_constraints.values())


def test_full_model_creates_presence_variables_when_minprop_is_zero():
    ilp = _sample_clonal_ilp(mu=0)
    ilp.create_model()

    assert hasattr(ilp.model, 'x_2_1')
    assert len(ilp.model.sample_clonal_constraints) == 16


def test_carch_fixes_states_only_for_clones_present_in_fixed_u():
    ilp = _sample_clonal_ilp()
    ilp.fix_u([[0.5], [0.5], [0.0]])
    ilp.create_model()

    assert len(ilp.model.sample_clonal_constraints) == 4
    constrained_variables = {constraint.body.name for constraint in ilp.model.sample_clonal_constraints.values()}
    assert constrained_variables == {'cA_1_2', 'cB_1_2', 'cA_2_2', 'cB_2_2'}


def test_uarch_excludes_clones_with_incompatible_fixed_states():
    ilp = _sample_clonal_ilp()
    ilp.fix_c(
        [[1, 1, 2], [1, 2, 2]],
        [[1, 1, 1], [1, 1, 1]],
    )
    ilp.create_model()

    assert len(ilp.model.sample_clonal_constraints) == 1
    constraint = next(iter(ilp.model.sample_clonal_constraints.values()))
    assert constraint.body.name == 'u_3_1'
    assert value(constraint.lower) == value(constraint.upper) == 0


def _conflicting_sample_clonal_ilp(n):
    cluster_ids = ('c1', 'c2', 'c3')
    sample_ids = ('TumorA', 'TumorB')
    sample_clonal = resolve_sample_clonal(
        parse_sample_clonal_inline('TumorA=c1:1:1,c2:2:1;TumorB=c1:2:1,c3:1:1'),
        cluster_ids,
        sample_ids,
    )
    f_a = pd.DataFrame(
        {'TumorA': [1.2, 1.5, 1.0], 'TumorB': [1.5, 1.0, 1.2]},
        index=cluster_ids,
    )
    f_b = pd.DataFrame(
        {'TumorA': [0.8, 1.0, 1.0], 'TumorB': [1.0, 1.0, 0.8]},
        index=cluster_ids,
    )
    return ILPSubset(
        n=n,
        cn_max=6,
        d=-1,
        mu=0.1,
        ampdel=False,
        copy_numbers={},
        f_a=f_a,
        f_b=f_b,
        w=pd.Series({cluster: 1.0 for cluster in cluster_ids}),
        purities=[0.5, 0.6],
        base=1,
        sample_clonal=sample_clonal,
    )


def test_sample_clonal_seed_fixes_purity_and_separates_incompatible_samples():
    ilp = _conflicting_sample_clonal_ilp(n=3)

    seed = ilp.build_random_u(random_seed=7)

    assert seed[0] == pytest.approx([0.5, 0.4])
    assert seed.sum(axis=0) == pytest.approx([1.0, 1.0])
    assert seed[1:, :].sum(axis=0) == pytest.approx([0.5, 0.6])
    assert not any(seed[clone_index, 0] > 0 and seed[clone_index, 1] > 0 for clone_index in (1, 2))


def test_sample_clonal_seed_rejects_insufficient_tumor_clones():
    ilp = _conflicting_sample_clonal_ilp(n=2)

    with pytest.raises(ValueError, match='No feasible tumor-clone assignment'):
        ilp.build_random_u(random_seed=7)


@pytest.mark.parametrize(
    ('inline', 'options', 'message'),
    [
        (
            'Tumor=c1:1:1,c2:2:1',
            {'cn_max': 2, 'base': 1, 'ampdel': False},
            'exceeds copy-number bound',
        ),
        (
            'TumorA=c1:0:1,c2:2:1;TumorB=c1:2:1,c3:1:1',
            {'cn_max': 6, 'base': 1, 'ampdel': True},
            'force both amplification and deletion',
        ),
        (
            'TumorA=c1:1:1,c2:2:1;TumorB=c1:2:1,c3:1:1',
            {'cn_max': 6, 'base': 1, 'ampdel': False, 'cnstates': 1},
            'exceeding cnstates',
        ),
        (
            'Tumor=c1:0:0,c2:1:0',
            {'cn_max': 6, 'base': 1, 'ampdel': False},
            'large cluster',
        ),
    ],
)
def test_validate_sample_clonal_constraints_rejects_infeasible_requirements(inline, options, message):
    parsed = parse_sample_clonal_inline(inline)
    samples = tuple(entry.sample for entry in parsed.entries[::2])
    clusters = tuple(dict.fromkeys(entry.cluster for entry in parsed.entries))
    resolved = resolve_sample_clonal(parsed, clusters, samples)

    with pytest.raises(SampleClonalFormatError, match=message):
        validate_sample_clonal_constraints(
            resolved,
            [0.5] * len(samples),
            weights=pd.Series({cluster: 1.0 for cluster in clusters}),
            **options,
        )


def _resolved_fixed_matrix(tmp_path, rows, sample_clonal, clones=3):
    path = tmp_path / 'fixed.tsv'
    path.write_text('# HATCHET_FIXED_ENTRIES 1\n' f'# clones {clones}\n' 'matrix\trow\tcolumn\tvalue\n' f'{rows}')
    fixed = parse_fixed_matrix(path)
    return resolve_fixed_matrix(fixed, sample_clonal.cluster_ids, sample_clonal.sample_ids)


def test_fixed_matrix_rejects_fixed_positive_incompatible_clone(tmp_path):
    sample_clonal = _conflicting_sample_clonal_ilp(n=3).sample_clonal
    fixed = _resolved_fixed_matrix(
        tmp_path,
        'A\tc1\tclone1\t2\nU\tclone1\tTumorA\t0.5\n',
        sample_clonal,
    )

    with pytest.raises(FixedMatrixFormatError, match='fixed-positive clone1 is incompatible'):
        validate_fixed_matrix(
            fixed,
            cn_max=6,
            mu=0.1,
            copy_numbers={},
            purities=[0.5, 0.6],
            sample_clonal=sample_clonal,
            ampdel=False,
        )


def _write_compatible_fixed_matrix(tmp_path, normal_proportion=0.5):
    path = tmp_path / 'sample-clonal-fixed.tsv'
    path.write_text(
        '# HATCHET_FIXED_ENTRIES 1\n'
        '# clones 2\n'
        'matrix\trow\tcolumn\tvalue\n'
        'A\t1\tclone1\t1\n'
        'B\t1\tclone1\t1\n'
        'A\t2\tclone1\t2\n'
        'B\t2\tclone1\t1\n'
        f'U\tnormal\tTumor\t{normal_proportion}\n'
        'U\tclone1\tTumor\t0.5\n'
    )
    return path


@pytest.mark.parametrize('solver', ('gurobi', 'cpp'))
def test_compute_cn_solves_with_sample_clonal_and_fixed_matrix(
    tmp_path,
    monkeypatch,
    solver,
    sample_clonal_cpp_solver,
):
    if solver == 'gurobi':
        from hatchet.utils.solve import solver_available

        if not solver_available('gurobi'):
            pytest.skip('a working Gurobi Pyomo solver is unavailable')
    input_prefix = _write_solver_inputs(tmp_path)
    output_folder = tmp_path / 'output'
    output_folder.mkdir()
    fixed_path = _write_compatible_fixed_matrix(tmp_path)
    monkeypatch.setattr(config.compute_cn, 'solver', solver)
    args = [str(sample_clonal_cpp_solver)] if solver == 'cpp' else []
    args.extend(
        [
            '-i',
            str(input_prefix),
            '-x',
            str(output_folder),
            '-n',
            '2',
            '-p',
            '1',
            '-v',
            '0',
            '--mode',
            '1',
            '-s',
            '30',
            '-eD',
            '6',
            '--diploid',
            '--sample-clonal',
            'Tumor=1:1:1,2:2:1',
            '--fixed-matrix',
            str(fixed_path),
        ]
    )

    compute_cn(args=args)

    _assert_public_sample_clonal_outputs(output_folder)


def test_compute_cn_rejects_fixed_normal_proportion_conflicting_with_sample_clonal_purity(tmp_path, monkeypatch):
    from hatchet.utils.solve import solver_available

    if not solver_available('gurobi'):
        pytest.skip('a working Gurobi Pyomo solver is unavailable')
    input_prefix = _write_solver_inputs(tmp_path)
    output_folder = tmp_path / 'output'
    output_folder.mkdir()
    fixed_path = tmp_path / 'conflicting-normal.tsv'
    fixed_path.write_text(
        '# HATCHET_FIXED_ENTRIES 1\n' '# clones 2\n' 'matrix\trow\tcolumn\tvalue\n' 'U\tnormal\tTumor\t0.4\n'
    )
    monkeypatch.setattr(config.compute_cn, 'solver', 'gurobi')

    with pytest.raises(FixedMatrixFormatError, match='conflicts with purity'):
        compute_cn(
            args=[
                '-i',
                str(input_prefix),
                '-x',
                str(output_folder),
                '-n',
                '2',
                '--mode',
                '1',
                '-eD',
                '6',
                '--diploid',
                '--sample-clonal',
                'Tumor=1:1:1,2:2:1',
                '--fixed-matrix',
                str(fixed_path),
            ]
        )


@pytest.mark.parametrize(
    ('ploidy_flag', 'bound_flag'),
    (('--diploid', '--diploidcmax'), ('--tetraploid', '--tetraploidcmax')),
)
def test_compute_cn_applies_active_ploidy_copy_bound_to_sample_clonal_states(
    tmp_path,
    monkeypatch,
    ploidy_flag,
    bound_flag,
):
    from hatchet.utils.solve import solver_available

    if not solver_available('gurobi'):
        pytest.skip('a working Gurobi Pyomo solver is unavailable')
    input_prefix = _write_solver_inputs(tmp_path)
    output_folder = tmp_path / 'output'
    output_folder.mkdir()
    monkeypatch.setattr(config.compute_cn, 'solver', 'gurobi')

    with pytest.raises(SampleClonalFormatError, match='exceeds copy-number bound 2'):
        compute_cn(
            args=[
                '-i',
                str(input_prefix),
                '-x',
                str(output_folder),
                '-n',
                '2',
                '--mode',
                '1',
                ploidy_flag,
                bound_flag,
                '2',
                '--sample-clonal',
                'Tumor=1:1:1,2:2:1',
            ]
        )


def test_fixed_matrix_rejects_shared_clone_fixed_positive_in_incompatible_samples(tmp_path):
    sample_clonal = _conflicting_sample_clonal_ilp(n=3).sample_clonal
    fixed = _resolved_fixed_matrix(
        tmp_path,
        'U\tclone1\tTumorA\t0.5\nU\tclone1\tTumorB\t0.6\n',
        sample_clonal,
    )

    with pytest.raises(FixedMatrixFormatError, match='assigned to samples with incompatible'):
        validate_fixed_matrix(
            fixed,
            cn_max=6,
            mu=0.1,
            copy_numbers={},
            purities=[0.5, 0.6],
            sample_clonal=sample_clonal,
            ampdel=False,
        )


def test_fixed_matrix_rejects_when_no_clone_can_receive_remaining_purity(tmp_path):
    sample_clonal = _conflicting_sample_clonal_ilp(n=3).sample_clonal
    fixed = _resolved_fixed_matrix(
        tmp_path,
        'U\tclone1\tTumorA\t0\nU\tclone2\tTumorA\t0\n',
        sample_clonal,
    )

    with pytest.raises(FixedMatrixFormatError, match='no tumor clone is eligible'):
        validate_fixed_matrix(
            fixed,
            cn_max=6,
            mu=0.1,
            copy_numbers={},
            purities=[0.5, 0.6],
            sample_clonal=sample_clonal,
            ampdel=False,
        )


def test_verify_sample_clonal_solution_checks_purity_and_present_clone_states():
    sample_clonal = _sample_clonal_ilp().sample_clonal
    c_a = [[1, 1, 0], [1, 2, 0]]
    c_b = [[1, 1, 0], [1, 1, 0]]
    u = [[0.5], [0.5], [0.0]]

    verify_sample_clonal_solution(sample_clonal, [0.5], c_a, c_b, u)

    with pytest.raises(RuntimeError, match='violates A'):
        verify_sample_clonal_solution(sample_clonal, [0.5], [[1, 2, 0], [1, 2, 0]], c_b, u)
    with pytest.raises(RuntimeError, match='normal proportion'):
        verify_sample_clonal_solution(sample_clonal, [0.5], c_a, c_b, [[0.4], [0.6], [0.0]])


@pytest.fixture(scope='module')
def cpp_sample_clonal_reader(tmp_path_factory):
    compiler = shutil.which('g++')
    if compiler is None:
        pytest.skip('g++ is required for the native sample-clonal parser test')
    repository = Path(__file__).parents[1]
    executable = tmp_path_factory.mktemp('sample-clonal-cpp') / 'reader'
    subprocess.run(
        [
            compiler,
            '-std=c++11',
            '-I' + str(repository / 'src'),
            str(repository / 'tests' / 'sample_clonal_cpp_test.cpp'),
            str(repository / 'src' / 'sample_clonal.cpp'),
            str(repository / 'src' / 'utils.cpp'),
            '-pthread',
            '-o',
            str(executable),
        ],
        check=True,
    )
    return executable


def test_cpp_sample_clonal_parser_and_calibration_match_python(tmp_path, cpp_sample_clonal_reader):
    path = _write_sample_clonal(
        tmp_path,
        f'{FORMAT_MARKER}\nsample\tcluster\tcn_a\tcn_b\nTumor\tc1\t1\t1\nTumor\tc2\t2\t1\n',
    )

    subprocess.run([str(cpp_sample_clonal_reader), str(path)], check=True)


@pytest.mark.parametrize(
    ('contents', 'message'),
    [
        ('wrong\n', 'expected format marker'),
        (
            f'{FORMAT_MARKER}\nsample\tcluster\tcn_a\tcn_b\nTumor\tmissing\t1\t1\nTumor\tc2\t2\t1\n',
            'unknown SEG cluster identifier',
        ),
        (
            f'{FORMAT_MARKER}\nsample\tcluster\tcn_a\tcn_b\nTumor\tc1\t1\t1\n',
            'exactly two sample-clonal entries',
        ),
    ],
)
def test_cpp_sample_clonal_parser_rejects_invalid_files(tmp_path, cpp_sample_clonal_reader, contents, message):
    path = _write_sample_clonal(tmp_path, contents)

    result = subprocess.run([str(cpp_sample_clonal_reader), str(path)], text=True, capture_output=True)

    assert result.returncode == 1
    assert message in result.stderr


def test_compute_cn_documentation_covers_sample_clonal_contract():
    documentation = (Path(__file__).parents[1] / 'docs' / 'source' / 'doc_compute_cn.md').read_text()

    assert '--sample-clonal' in documentation
    assert '--sample-clonal-file' in documentation
    assert '# HATCHET_SAMPLE_CLONAL 1' in documentation
    assert 'HATCHET_COMPUTE_CN_SAMPLECLONAL' in documentation
    assert 'HATCHET_COMPUTE_CN_SAMPLECLONALFILE' in documentation
    assert 'Exactly one of `--diploid` or `--tetraploid` is required' in documentation
    assert 'cannot be combined with global `--clonal`' in documentation
