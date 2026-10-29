"""Summarize the runs of bench_cnstates.sh as Markdown tables, one section per scenario and solving mode.

usage: summarize_cnstates.py <WORK>

Reads <WORK>/<scenario>/input/truth/ (from simulate.py) and, for every run <WORK>/<scenario>/mode<M>/<config>/:
time.json (from measure.py), log.err (compute-cn's log), coordinate_descence.tsv (the native solver's per-step times
and, since the per-cluster limits, model sizes), results.*.bbc.ucn.tsv and best.bbc.ucn.

Accuracy of the chosen solution (best.bbc.ucn) against the truth, independent of clone labels:
- frac MAE: mean |inferred - true| fractional A plus B copy number per cluster and sample, weighted by cluster bins
- purity MAE: mean |inferred - true| tumor purity over samples
- over / under: genome fraction of clusters with more / fewer distinct tumor states (over clones present in at least
  one sample) than the truth; hot: mean inferred / true states of the hot clusters
- exact CN: genome fraction with the true allele-specific state, averaged over tumor clones matched to true clones by
  their proportions across samples (Hungarian matching); unmatched clones are reported as extra or missing
"""
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
import hashlib

ANSI = re.compile(r'\x1b\[[0-9;]*m')


sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '2026-09-29-compute-cn-efficiency'))
CONFIG_ORDER = [
    'free', 'free_old', 'd', 'd_exempt', 'd_exempt_partial', 'map_oracle', 'd_exempt_pairwise', 'd1_exempt', 'd3_exempt'
]  # fmt: skip
PRESENT = 1e-3  # a clone with a larger proportion in some sample counts as present


def load_truth(input_dir):
    tdir = os.path.join(input_dir, 'truth')
    states = pd.read_csv(os.path.join(tdir, 'cluster_states.tsv'), sep='\t').set_index('cluster')
    A = pd.read_csv(os.path.join(tdir, 'A.tsv'), sep='\t', index_col=0)
    B = pd.read_csv(os.path.join(tdir, 'B.tsv'), sep='\t', index_col=0)
    U = pd.read_csv(os.path.join(tdir, 'U.tsv'), sep='\t', index_col=0)
    params = json.load(open(os.path.join(tdir, 'params.json')))
    return {'states': states, 'A': A, 'B': B, 'U': U, 'params': params}


def step_stats(path):
    """C-step count, mean and total Gurobi runtime, and mean model size (native solver only)."""
    tsv = os.path.join(path, 'coordinate_descence.tsv')
    if not os.path.exists(tsv) or os.path.getsize(tsv) == 0:
        return None
    steps = pd.read_csv(tsv, sep='\t', header=None)
    carch = steps[steps[2] == 'C']
    stats = {'count': len(carch), 'mean_s': carch[5].mean(), 'total_s': carch[5].sum()}
    if steps.shape[1] >= 9:
        stats.update(rows=carch[6].mean(), discrete=carch[7].mean(), nonzeros=carch[8].mean())
    return stats


def solution(best):
    """Per cluster: tumor (A, B) per clone; per sample: proportions of normal and tumor clones."""
    result = pd.read_csv(best, sep='\t')
    clones = sorted(
        (column for column in result.columns if re.fullmatch(r'cn_clone\d+', column)), key=lambda c: int(c[8:])
    )
    first = result.drop_duplicates(['CLUSTER', 'SAMPLE'])
    states = first.drop_duplicates('CLUSTER').set_index('CLUSTER')
    A = pd.DataFrame({c: states[c].map(lambda v: int(str(v).split('|')[0])) for c in clones})
    B = pd.DataFrame({c: states[c].map(lambda v: int(str(v).split('|')[1])) for c in clones})
    U = first.drop_duplicates('SAMPLE').set_index('SAMPLE')[['u_normal'] + [f'u_{c[3:]}' for c in clones]].T
    U.index = ['normal'] + clones
    return A, B, U.astype(float)


def md5(path):
    with open(path, 'rb') as handle:
        return hashlib.md5(handle.read()).hexdigest()


def chosen(path):
    """The results file that compute-cn copied to best.bbc.ucn, e.g. 'diploid.n4'."""
    best = os.path.join(path, 'best.bbc.ucn')
    if not os.path.exists(best):
        return 'none'
    digest = md5(best)
    for result in sorted(glob.glob(os.path.join(path, 'results.*.bbc.ucn.tsv'))):
        if md5(result) == digest:
            return os.path.basename(result)[len('results.') : -len('.bbc.ucn.tsv')]
    return '?'


def objectives(log_err):
    """(ploidy, clones) -> objective, from compute-cn's '# Best objective found with N clones' lines."""
    found, ploidy = {}, None
    if not os.path.exists(log_err):
        return found
    with open(log_err, errors='replace') as log:
        for line in log:
            line = ANSI.sub('', line)
            match = re.search(r'Running (diploid|tetraploid) with (\d+) clones', line)
            if match:
                ploidy = match.group(1)
                continue
            match = re.search(r'Best objective found with (\d+) clones: (\S+)', line)
            if match and ploidy is not None:
                found[(ploidy, int(match.group(1)))] = float(match.group(2))
    return found


def accuracy(best, truth):
    A, B, U = solution(best)
    tA, tB, tU = truth['A'], truth['B'], truth['U']
    samples = list(tU.columns)
    U = U.reindex(columns=samples)
    clusters = truth['states'].index
    weights = truth['states']['bins'].reindex(clusters).astype(float)
    weights /= weights.sum()
    A, B = A.reindex(clusters), B.reindex(clusters)

    def fractions(a, b, u):
        fa = u.loc['normal'].values[None, :] + a.values @ u.iloc[1:].values
        fb = u.loc['normal'].values[None, :] + b.values @ u.iloc[1:].values
        return fa, fb

    fa, fb = fractions(A, B, U)
    ta, tb = fractions(tA, tB, tU)
    frac_mae = float((weights.values[:, None] * (np.abs(fa - ta) + np.abs(fb - tb))).sum() / len(samples))
    purity_mae = float(np.abs((1 - U.loc['normal']) - (1 - tU.loc['normal'])).mean())

    present = [c for c in A.columns if U.loc[c].max() > PRESENT]
    inferred_states = pd.Series(
        [len(set(zip(A.loc[c, present], B.loc[c, present]))) if present else 0 for c in clusters], index=clusters
    )
    true_states = truth['states']['true_states']
    hot = truth['states']['kind'] == 'hot'
    over = float(weights[inferred_states > true_states].sum())
    under = float(weights[inferred_states < true_states].sum())
    hot_cell = f'{inferred_states[hot].mean():.2f} / {true_states[hot].mean():.2f}' if hot.any() else 'n/a'

    true_clones = list(tU.index[1:])
    cost = np.array([[np.abs(U.loc[c].values - tU.loc[t].values).sum() for t in true_clones] for c in present])
    exact, extra, missing = float('nan'), len(present), len(true_clones)
    if len(present):
        rows, cols = linear_sum_assignment(cost)
        matched = [
            float(weights[(A[present[r]] == tA[true_clones[t]]) & (B[present[r]] == tB[true_clones[t]])].sum())
            for r, t in zip(rows, cols)
        ]
        exact = float(np.mean(matched))
        extra, missing = len(present) - len(rows), len(true_clones) - len(rows)
    return {
        'tumor_clones': len(present),
        'frac_mae': frac_mae,
        'purity_mae': purity_mae,
        'over': over,
        'under': under,
        'hot': hot_cell,
        'exact': exact,
        'extra': extra,
        'missing': missing,
    }


def config_rank(config):
    return CONFIG_ORDER.index(config) if config in CONFIG_ORDER else len(CONFIG_ORDER)


def runs(work):
    found = {}
    for path in glob.glob(os.path.join(work, '*', 'mode*', '*')):
        if not os.path.isdir(path):
            continue
        mode_dir = os.path.dirname(path)
        scenario = os.path.basename(os.path.dirname(mode_dir))
        mode = os.path.basename(mode_dir)[len('mode') :]
        found.setdefault((scenario, mode), []).append(path)
    for key in found:
        found[key].sort(key=lambda path: config_rank(os.path.basename(path)))
    return dict(sorted(found.items()))


def fmt(value, spec, missing='n/a'):
    return missing if value is None or (isinstance(value, float) and np.isnan(value)) else format(value, spec)


def main(work):
    found = runs(work)
    if not found:
        sys.exit(f'No runs under {work}/<scenario>/mode*/')
    print(f'# cnstates benchmark: {work}\n')
    for (scenario, mode), paths in found.items():
        truth = load_truth(os.path.join(work, scenario, 'input'))
        params = truth['params']
        print(
            f'## {scenario}, mode {mode}\n\n'
            f'{params["clusters"]} clusters ({int((truth["states"]["kind"] == "hot").sum())} hot, '
            f'{int((truth["states"]["kind"] == "background").sum())} background), {params["samples"]} samples, '
            f'{params["tumor_clones"]} true tumor clones, {"WGD" if params["wgd"] else "diploid"}, D = {params["D"]}.\n'
        )
        print(
            '| config | exit | wall | peak RSS | C steps | mean C step | mean C rows / discrete | chosen '
            '| tumor clones | frac MAE | purity MAE | over | under | hot states | exact CN | extra / missing |'
        )
        print('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
        objective_rows = {}
        for path in paths:
            config = os.path.basename(path)
            timing_path = os.path.join(path, 'time.json')
            timing = json.load(open(timing_path)) if os.path.exists(timing_path) else {}
            steps = step_stats(path)
            best = os.path.join(path, 'best.bbc.ucn')
            acc = accuracy(best, truth) if os.path.exists(best) else {}
            size = 'n/a'
            if steps and 'rows' in steps:
                size = f'{steps["rows"]:.0f} / {steps["discrete"]:.0f}'
            print(
                f'| {config} | {timing.get("returncode", "?")} | {fmt(timing.get("wall_s"), ".0f")} s '
                f'| {fmt(timing.get("peak_tree_rss_mb", float("nan")) / 1024, ".2f")} GB '
                f'| {steps["count"] if steps else "n/a"} | {fmt(steps["mean_s"] * 1000 if steps else None, ".1f")} ms '
                f'| {size} | {chosen(path)} | {acc.get("tumor_clones", "n/a")} '
                f'| {fmt(acc.get("frac_mae"), ".4f")} | {fmt(acc.get("purity_mae"), ".3f")} '
                f'| {fmt(acc.get("over"), ".3f")} | {fmt(acc.get("under"), ".3f")} | {acc.get("hot", "n/a")} '
                f'| {fmt(acc.get("exact"), ".3f")} | {acc.get("extra", "n/a")} / {acc.get("missing", "n/a")} |'
            )
            objective_rows[config] = objectives(os.path.join(path, 'log.err'))

        keys = sorted({key for values in objective_rows.values() for key in values})
        if keys:
            configs = list(objective_rows)
            print('\nObjective per ploidy and clone count (compute-cn counts the normal clone):\n')
            print('| ploidy | clones | ' + ' | '.join(configs) + ' |')
            print('|---|---|' + '---|' * len(configs))
            for ploidy, clones in keys:
                cells = [fmt(objective_rows[c].get((ploidy, clones)), '.4f') for c in configs]
                print(f'| {ploidy} | {clones} | ' + ' | '.join(cells) + ' |')
        print()


if __name__ == '__main__':
    main(sys.argv[1])
