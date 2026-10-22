import os
import numpy as np
import pandas as pd
from pyomo import environ as pe
from collections import OrderedDict

from hatchet.utils.solve.ilp_subset import ILPSubset, ILPSubsetSplit
from hatchet.utils.solve.cd import CoordinateDescent, CoordinateDescentSplit
from hatchet.utils.fixed_matrix import resolve_fixed_matrix, validate_fixed_matrix, verify_fixed_matrix_solution
from hatchet.utils.sample_clonal import (
    PURITY_TOLERANCE,
    derive_sample_clonal_calibration,
    reindex_sample_clonal,
    validate_sample_clonal_constraints,
    verify_sample_clonal_solution,
)
from hatchet.utils.solve.utils import canonicalize_clone_order, parse_clonal, scale_rdr
from hatchet import config
import hatchet.utils.Supporting as sp


def solver_available(solver=None):
    solver = solver or config.compute_cn.solver
    if solver == 'cpp':
        return os.getenv('GRB_LICENSE_FILE') is not None
    elif solver == 'gurobipy':
        return pe.SolverFactory('gurobi', solver_io='python').available(exception_flag=False)
    return pe.SolverFactory(solver).available(exception_flag=False)


def _finalize_solution(result, fixed_matrix, sample_clonal=None, purities=None):
    """Verify user constraints and report interchangeable tumor clones in canonical order."""
    if result is not None:
        obj, cA, cB, u, cluster_ids, sample_ids = result
        if fixed_matrix is not None:
            verify_fixed_matrix_solution(fixed_matrix, cA, cB, u)
        if sample_clonal is not None:
            verify_sample_clonal_solution(sample_clonal, purities, cA, cB, u)
        # The model has no clone-order constraint; clones addressed by the fixed matrix keep their labels
        if fixed_matrix is None or not fixed_matrix.addresses_tumor_clones:
            cA, cB, u = canonicalize_clone_order(cA, cB, u)
        result = obj, cA, cB, u, cluster_ids, sample_ids
    return result


def solve(
    clonal,
    bbc_file,
    seg_file,
    n,
    solver='gurobi',
    solve_mode='cd',
    d=-1,
    cn_max=-1,
    mu=0.01,
    diploid_threshold=0.1,
    ampdel=True,
    n_seed=400,
    n_worker=8,
    random_seed=None,
    max_iters=None,
    timelimit=None,
    binwise=False,
    purities=None,
    fixed_matrix=None,
    base=None,
    sample_clonal=None,
    max_homdel_size=0.005,
    ampdel_exempt=None,
):

    if not np.isfinite(max_homdel_size) or not 0 <= max_homdel_size <= 1:
        raise ValueError('The maximum homozygous-deletion cluster-size fraction must be in [0, 1]')

    assert solve_mode in ('ilp', 'cd', 'both'), 'Unrecognized solve_mode'
    assert solver_available(solver), f'Solver {solver} not available or not licensed'

    if max_iters is None:
        max_iters = 10

    df = pd.read_csv(seg_file, sep='\t')
    df = df.sort_values(['#ID', 'SAMPLE'])

    bbc = pd.read_table(bbc_file)
    bbc = bbc.sort_values(by=['#CHR', 'START', 'SAMPLE'])

    # sanity-check
    sample_ids = np.sort(df['SAMPLE'].unique())
    for _cluster_id, _df in df.groupby('#ID'):
        _sample_ids = _df['SAMPLE'].values
        if not np.all(_sample_ids == sample_ids):
            raise ValueError(f'Sample IDs for cluster {_cluster_id} do not match {sample_ids}')

    rdr = df.pivot(index='#ID', columns='SAMPLE', values='RD')
    baf = df.pivot(index='#ID', columns='SAMPLE', values='BAF')

    if ampdel_exempt:
        if not ampdel:
            raise ValueError('Amp-del exempt clusters require the amp-del constraint to be active')
        # Match labels as strings so that IDs given on the command line resolve to the seg file's #ID dtype
        known = {str(cluster_id): cluster_id for cluster_id in rdr.index}
        unknown = [label for label in map(str, ampdel_exempt) if label not in known]
        if unknown:
            raise ValueError(f'Amp-del exempt clusters not found in the seg file: {", ".join(unknown)}')
        ampdel_exempt = frozenset(known[label] for label in map(str, ampdel_exempt))

    bins = {}  # cluster_id => no. of bins
    for cluster_id, _df in df.groupby('#ID'):
        if len(_df['#BINS'].unique()) != 1:
            raise ValueError(f'Bin sizes for cluster {cluster_id} across tumor samples are not identical!')
        bins[cluster_id] = _df.iloc[0]['#BINS']
    bins = pd.Series(bins)

    weights = 100 * bins / sum(bins)

    if sample_clonal is not None:
        if clonal is not None:
            raise ValueError('Sample-clonal input cannot be combined with global clonal copy numbers')
        if base not in (1, 2):
            raise ValueError('Sample-clonal input requires an explicit factorization base of 1 or 2')

        resolved_sample_clonal = reindex_sample_clonal(sample_clonal, rdr.index, rdr.columns)
        calibration = derive_sample_clonal_calibration(resolved_sample_clonal, rdr)
        derived_purities = list(calibration.purity)
        if purities is not None:
            if len(purities) != len(derived_purities):
                raise ValueError(f'Purity count {len(purities)} does not match sample count {len(derived_purities)}')
            for sample_id, provided, derived in zip(rdr.columns, purities, derived_purities):
                difference = abs(provided - derived)
                if difference > PURITY_TOLERANCE and not np.isclose(
                    difference,
                    PURITY_TOLERANCE,
                    rtol=0.0,
                    atol=1e-12,
                ):
                    raise ValueError(
                        f'Purity {provided} for sample {sample_id!r} does not match '
                        f'sample-clonal purity {derived} within tolerance {PURITY_TOLERANCE}'
                    )

        purities = derived_purities
        gamma = pd.Series(calibration.gamma, index=calibration.sample_ids, dtype=float)
        copy_numbers = {}
    else:
        resolved_sample_clonal = None
        if clonal is None:
            _candidate_cluster_ids = np.all(rdr > 0.5 - diploid_threshold, axis=1)
            if not np.any(_candidate_cluster_ids):
                raise RuntimeError(f'Unable to determine cluster with diploid RDR threshold {diploid_threshold}')
            dipoid_cluster_id = (_candidate_cluster_ids * weights).idxmax()
            copy_numbers = {dipoid_cluster_id: (1, 1)}
        else:
            copy_numbers = parse_clonal(clonal)
        gamma = scale_rdr(rdr, copy_numbers)
        if base is None:
            base = min(2, len(copy_numbers))

    sp.log(
        msg='Computed scaling factor gamma = \n' + str(gamma) + '\n',
        level='INFO',
    )
    if resolved_sample_clonal is not None:
        sp.log(
            msg='Computed sample-clonal purities = \n' + str(pd.Series(purities, index=rdr.columns)) + '\n',
            level='INFO',
        )
    rdr = rdr * gamma
    f_b = rdr * baf
    f_a = rdr - f_b

    if cn_max is None or cn_max <= 0:
        cn_max = max(1, int(np.ceil(np.nanmax(rdr.values))))
        if resolved_sample_clonal is not None:
            cn_max = max(cn_max, max(entry.total_copy_number for entry in resolved_sample_clonal.entries))

    if resolved_sample_clonal is not None:
        validate_sample_clonal_constraints(
            resolved_sample_clonal,
            purities,
            cn_max=cn_max,
            base=base,
            ampdel=ampdel,
            cnstates=d,
            weights=weights,
            max_homdel_size=max_homdel_size,
            ampdel_exempt=ampdel_exempt,
        )

    resolved_fixed_matrix = None
    if fixed_matrix is not None:
        if fixed_matrix.clones != n:
            raise ValueError(
                f'Fixed-matrix file declares {fixed_matrix.clones} clones, but the solver was configured for {n}'
            )
        resolved_fixed_matrix = resolve_fixed_matrix(fixed_matrix, f_a.index, f_a.columns)
        validate_fixed_matrix(
            resolved_fixed_matrix,
            cn_max=cn_max,
            mu=mu,
            copy_numbers=copy_numbers,
            purities=purities,
            sample_clonal=resolved_sample_clonal,
            ampdel=ampdel,
            base=base,
            cnstates=d,
            weights=weights,
            max_homdel_size=max_homdel_size,
            ampdel_exempt=ampdel_exempt,
        )

    if not binwise:
        if solve_mode == 'ilp':
            ilp = ILPSubset(
                n,
                cn_max,
                d=d,
                mu=mu,
                ampdel=ampdel,
                copy_numbers=copy_numbers,
                f_a=f_a,
                f_b=f_b,
                w=weights,
                purities=purities,
                fixed_matrix=resolved_fixed_matrix,
                base=base,
                sample_clonal=resolved_sample_clonal,
                max_homdel_size=max_homdel_size,
                ampdel_exempt=ampdel_exempt,
            )
            ilp.create_model(pprint=True)
            return _finalize_solution(
                ilp.run(solver_type=solver, timelimit=timelimit),
                resolved_fixed_matrix,
                resolved_sample_clonal,
                purities,
            )
        elif solve_mode == 'cd':
            cd = CoordinateDescent(
                f_a=f_a,
                f_b=f_b,
                n=n,
                mu=mu,
                d=d,
                cn_max=cn_max,
                w=weights,
                ampdel=ampdel,
                cn=copy_numbers,
                purities=purities,
                fixed_matrix=resolved_fixed_matrix,
                base=base,
                sample_clonal=resolved_sample_clonal,
                max_homdel_size=max_homdel_size,
                ampdel_exempt=ampdel_exempt,
            )
            return _finalize_solution(
                cd.run(
                    solver_type=solver,
                    max_iters=max_iters,
                    n_seed=n_seed,
                    j=n_worker,
                    random_seed=random_seed,
                    timelimit=timelimit,
                ),
                resolved_fixed_matrix,
                resolved_sample_clonal,
                purities,
            )
        else:
            cd = CoordinateDescent(
                f_a=f_a,
                f_b=f_b,
                n=n,
                mu=mu,
                d=d,
                cn_max=cn_max,
                w=weights,
                ampdel=ampdel,
                cn=copy_numbers,
                purities=purities,
                fixed_matrix=resolved_fixed_matrix,
                base=base,
                sample_clonal=resolved_sample_clonal,
                max_homdel_size=max_homdel_size,
                ampdel_exempt=ampdel_exempt,
            )
            cd_result = cd.run(
                solver_type=solver,
                max_iters=max_iters,
                n_seed=n_seed,
                j=n_worker,
                random_seed=random_seed,
                timelimit=timelimit,
            )
            _finalize_solution(cd_result, resolved_fixed_matrix, resolved_sample_clonal, purities)
            _, cA, cB, u, _, _ = cd_result

            ilp = ILPSubset(
                n,
                cn_max,
                d=d,
                mu=mu,
                ampdel=ampdel,
                copy_numbers=copy_numbers,
                f_a=f_a,
                f_b=f_b,
                w=weights,
                purities=purities,
                fixed_matrix=resolved_fixed_matrix,
                base=base,
                sample_clonal=resolved_sample_clonal,
                max_homdel_size=max_homdel_size,
                ampdel_exempt=ampdel_exempt,
            )
            ilp.create_model()
            ilp.hot_start(cA, cB, u)
            return _finalize_solution(
                ilp.run(solver_type=solver, timelimit=timelimit),
                resolved_fixed_matrix,
                resolved_sample_clonal,
                purities,
            )

    else:
        bins = OrderedDict()  # cluster_id => RDR for cluster
        for cluster_id, _df in df.groupby('#ID'):
            if len(_df['#BINS'].unique()) != 1:
                raise ValueError(f'Bin sizes for cluster {cluster_id} across tumor samples are not identical!')
            my_bbc = bbc[bbc.CLUSTER == cluster_id]
            if _df.iloc[0]['#BINS'] * len(_df) != len(my_bbc):   # seg file should have 1 row per sample
                raise ValueError(f'BBC and SEG files describe inconsisitent # bins for cluster {cluster_id}!')
            bins[cluster_id] = my_bbc.RD

        binned_length = (bbc[bbc.SAMPLE == bbc.iloc[0].SAMPLE].END - bbc[bbc.SAMPLE == bbc.iloc[0].SAMPLE].START).sum()

        bins_rdr = {
            k: my_df.pivot(index=['#CHR', 'START'], columns='SAMPLE', values='RD')
            for k, my_df in bbc.groupby('CLUSTER')
        }
        bins_baf = {
            k: my_df.pivot(index=['#CHR', 'START'], columns='SAMPLE', values='BAF')
            for k, my_df in bbc.groupby('CLUSTER')
        }
        bins_length = {
            k: (
                (
                    my_df.pivot(index=['#CHR', 'START'], columns='SAMPLE', values='END')
                    - my_df.pivot(
                        index=['#CHR', 'START'],
                        columns='SAMPLE',
                        values='START',
                    )
                ).values[:, 0]
                * 100
            )
            / binned_length
            for k, my_df in bbc.groupby('CLUSTER')
        }

        binsA = {}
        binsB = {}
        for k in bins_rdr.keys():
            bins_rdr[k] = bins_rdr[k] * gamma
            binsB[k] = bins_rdr[k] * bins_baf[k]
            binsA[k] = bins_rdr[k] - binsB[k]

        if solve_mode == 'ilp':
            ilp = ILPSubsetSplit(
                n,
                cn_max,
                d=d,
                mu=mu,
                ampdel=ampdel,
                copy_numbers=copy_numbers,
                f_a=f_a,
                f_b=f_b,
                binsA=binsA,
                binsB=binsB,
                lengths=bins_length,
                purities=purities,
                fixed_matrix=resolved_fixed_matrix,
                base=base,
                sample_clonal=resolved_sample_clonal,
                w=weights,
                max_homdel_size=max_homdel_size,
                ampdel_exempt=ampdel_exempt,
            )
            ilp.create_model(pprint=True)
            return _finalize_solution(
                ilp.run(solver_type=solver, timelimit=timelimit),
                resolved_fixed_matrix,
                resolved_sample_clonal,
                purities,
            )
        elif solve_mode == 'cd':
            cd = CoordinateDescentSplit(
                f_a=f_a,
                f_b=f_b,
                n=n,
                mu=mu,
                d=d,
                cn_max=cn_max,
                ampdel=ampdel,
                cn=copy_numbers,
                binsA=binsA,
                binsB=binsB,
                lengths=bins_length,
                purities=purities,
                fixed_matrix=resolved_fixed_matrix,
                base=base,
                sample_clonal=resolved_sample_clonal,
                w=weights,
                max_homdel_size=max_homdel_size,
                ampdel_exempt=ampdel_exempt,
            )
            return _finalize_solution(
                cd.run(
                    solver_type=solver,
                    max_iters=max_iters,
                    n_seed=n_seed,
                    j=n_worker,
                    random_seed=random_seed,
                    timelimit=timelimit,
                ),
                resolved_fixed_matrix,
                resolved_sample_clonal,
                purities,
            )
        else:
            cd = CoordinateDescentSplit(
                f_a=f_a,
                f_b=f_b,
                n=n,
                mu=mu,
                d=d,
                cn_max=cn_max,
                ampdel=ampdel,
                cn=copy_numbers,
                binsA=binsA,
                binsB=binsB,
                lengths=bins_length,
                purities=purities,
                fixed_matrix=resolved_fixed_matrix,
                base=base,
                sample_clonal=resolved_sample_clonal,
                w=weights,
                max_homdel_size=max_homdel_size,
                ampdel_exempt=ampdel_exempt,
            )
            cd_result = cd.run(
                solver_type=solver,
                max_iters=max_iters,
                n_seed=n_seed,
                j=n_worker,
                random_seed=random_seed,
                timelimit=timelimit,
            )
            _finalize_solution(cd_result, resolved_fixed_matrix, resolved_sample_clonal, purities)
            _, cA, cB, u, _, _ = cd_result

            ilp = ILPSubsetSplit(
                n,
                cn_max,
                d=d,
                mu=mu,
                ampdel=ampdel,
                copy_numbers=copy_numbers,
                f_a=f_a,
                f_b=f_b,
                binsA=binsA,
                binsB=binsB,
                lengths=bins_length,
                purities=purities,
                fixed_matrix=resolved_fixed_matrix,
                base=base,
                sample_clonal=resolved_sample_clonal,
                w=weights,
                max_homdel_size=max_homdel_size,
                ampdel_exempt=ampdel_exempt,
            )
            ilp.create_model()
            ilp.hot_start(cA, cB, u)
            return _finalize_solution(
                ilp.run(solver_type=solver, timelimit=timelimit),
                resolved_fixed_matrix,
                resolved_sample_clonal,
                purities,
            )
