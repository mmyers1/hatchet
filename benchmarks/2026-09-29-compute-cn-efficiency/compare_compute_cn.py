"""Summarize the runs of bench_compute_cn.sh as Markdown tables.

usage: compare_compute_cn.py <WORK> <SOLVER>

Reads <WORK>/mode<M>/<version>/ for every mode and version: time.json (from measure.py), log.err (compute-cn's log),
coordinate_descence.tsv (per-step times written by the native solver), results.*.bbc.ucn.tsv and best.bbc.ucn.
"""
import glob
import hashlib
import json
import os
import re
import sys

import pandas as pd

ANSI = re.compile(r'\x1b\[[0-9;]*m')
VERSION_ORDER = {'old': 0, 'new': 1}


def runs(work):
    found = []
    for path in glob.glob(os.path.join(work, 'mode*', '*')):
        mode = os.path.basename(os.path.dirname(path))[len('mode') :]
        version = os.path.basename(path)
        if os.path.isdir(path):
            found.append((mode, version, path))
    return sorted(found, key=lambda run: (run[0], VERSION_ORDER.get(run[1], 2), run[1]))


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


def step_times(path):
    """Count, total and mean Gurobi runtime of the C and U steps of coordinate descent (native solver only)."""
    tsv = os.path.join(path, 'coordinate_descence.tsv')
    if not os.path.exists(tsv) or os.path.getsize(tsv) == 0:
        return {}
    # newer native solvers append model-size columns; only the first six are used here
    steps = pd.read_csv(tsv, sep='\t', header=None).iloc[:, :6]
    steps.columns = ['seed', 'iteration', 'step', 'obj', 'gap', 'runtime']
    stats = {}
    for step, group in steps.groupby('step'):
        stats[step] = (len(group), group['runtime'].sum(), group['runtime'].mean())
    return stats


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


def clone_order(path, seg, solver):
    """Whether the tumor clones in best.bbc.ucn have nondecreasing sum_s (s + 1)(A + B), with the solver's cluster
    order (the native solver: first appearance in the SEG file; Pyomo: sorted cluster IDs)."""
    best = os.path.join(path, 'best.bbc.ucn')
    if not os.path.exists(best):
        return 'n/a'
    result = pd.read_csv(best, sep='\t')
    clones = [column for column in result.columns if re.fullmatch(r'cn_clone\d+', column)]
    if len(clones) < 2:
        return 'n/a (one tumor clone)'
    clusters = pd.read_csv(seg, sep='\t')['#ID'].drop_duplicates().tolist()
    if solver != 'cpp':
        clusters = sorted(clusters)
    states = result.drop_duplicates('CLUSTER').set_index('CLUSTER')
    burdens = []
    for clone in clones:
        burden = 0
        for index, cluster in enumerate(clusters):
            if cluster in states.index:
                allele_a, allele_b = (int(value) for value in str(states.loc[cluster, clone]).split('|'))
                burden += (index + 1) * (allele_a + allele_b)
        burdens.append(burden)
    label = 'sorted' if burdens == sorted(burdens) else 'not sorted'
    return f'{label} {burdens}'


def main(work, solver):
    found = runs(work)
    if not found:
        sys.exit(f'No runs under {work}/mode*/')
    seg = None
    first_time = os.path.join(found[0][2], 'time.json')
    if os.path.exists(first_time):
        cmd = json.load(open(first_time))['cmd']
        seg = cmd[cmd.index('-i') + 1] + '.seg'

    print(f'# compute-cn benchmark: {work} (solver {solver})\n')
    print(
        '| mode | version | exit | wall | CPU | peak RSS | C steps (n, total, mean) | U steps (n, total, mean) '
        '| chosen | clone order of best |'
    )
    print('|---|---|---|---|---|---|---|---|---|---|')
    for mode, version, path in found:
        timing_path = os.path.join(path, 'time.json')
        timing = json.load(open(timing_path)) if os.path.exists(timing_path) else {}
        steps = step_times(path)

        def step(name):
            if name not in steps:
                return 'n/a'
            count, total, mean = steps[name]
            return f'{count}, {total:.1f} s, {mean * 1000:.1f} ms'

        order = clone_order(path, seg, solver) if seg else 'n/a'
        print(
            f"| {mode} | {version} | {timing.get('returncode', '?')} | {timing.get('wall_s', float('nan')):.0f} s "
            f"| {timing.get('cpu_s', float('nan')):.0f} s "
            f"| {timing.get('peak_tree_rss_mb', float('nan')) / 1024:.2f} GB "
            f"| {step('C')} | {step('U')} | {chosen(path)} | {order} |"
        )

    print('\nObjective per ploidy and clone count (lower is better):\n')
    by_mode = {}
    for mode, version, path in found:
        by_mode.setdefault(mode, {})[version] = objectives(os.path.join(path, 'log.err'))
    versions = sorted({version for _, version, _ in found}, key=lambda version: VERSION_ORDER.get(version, 2))
    print('| mode | ploidy | clones | ' + ' | '.join(versions) + ' |')
    print('|---|---|---|' + '---|' * len(versions))
    for mode in sorted(by_mode):
        keys = sorted({key for values in by_mode[mode].values() for key in values})
        for ploidy, clones in keys:
            values = [by_mode[mode].get(version, {}).get((ploidy, clones)) for version in versions]
            cells = ['n/a' if value is None else f'{value:.6f}' for value in values]
            print(f'| {mode} | {ploidy} | {clones} | ' + ' | '.join(cells) + ' |')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else 'cpp')
