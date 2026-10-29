"""Simulate compute-cn inputs with a known clone tree, allele-specific copy numbers and clone proportions.

usage: simulate.py <scenarios.json> <scenario> <outdir> [--smoke] [--self-check]

Writes <outdir>/bulk.bbc and <outdir>/bulk.seg (compute-cn's -i <outdir>/bulk), the ground truth in <outdir>/truth/
and the benchmark runs in <outdir>/runs.tsv and <outdir>/common_args.txt (read by bench_cnstates.sh). Every random
draw comes from one numpy Generator seeded by the scenario's 'seed', so a scenario is reproduced exactly.

Model (all parameters in scenarios.json):
- Tumor clones form a random tree (clone 0 is the trunk; clone t > 0 descends from a clone < t).
- Every cluster has a clonal state drawn from a fixed distribution (doubled under WGD). Cluster 1 is the largest and
  copy-neutral, cluster 2 the second largest and clonal but not neutral, so compute-cn's heuristics find them.
- 'Background' clusters get up to background_states - 1 subclonal events and 'hot' clusters hot_states - 1 events.
  An event adds +-1 copy of one allele in a subtree of the clone tree, keeping 0 <= B <= A and 1 <= A + B <= cmax.
  Hot clusters are resampled until they have more than D distinct tumor states, so a limit of D binds on them.
- Each sample has a purity in 'purity' and one dominant tumor clone (every clone dominates at least one sample when
  there are enough samples); tumor proportions below minprop are set to 0.
- Fractional copy numbers become a read-depth ratio (normalized by the sample's mean ploidy) and a B-allele frequency.
  Gaussian noise is added per cluster (rd_noise relative, baf_noise absolute; none for the neutral cluster) and per
  bin (bin_rd_noise, bin_baf_noise). The SEG values are the means of their bins.
"""
import argparse
import json
import math
import os
import shlex

import numpy as np
import pandas as pd

DIPLOID_STATES = {
    (1, 1): 0.40,
    (2, 1): 0.15,
    (1, 0): 0.15,
    (2, 0): 0.07,
    (2, 2): 0.08,
    (3, 1): 0.06,
    (3, 2): 0.03,
    (3, 0): 0.03,
    (4, 1): 0.03,
}
TOTAL_ALLELES = 1200  # ALPHA + BETA per bin
SNPS = 40
COVERAGE = 30.0


def load_scenario(path, name, smoke):
    with open(path) as handle:
        spec = json.load(handle)
    if name not in spec['scenarios']:
        raise SystemExit(f'Unknown scenario {name}; choose from {", ".join(spec["scenarios"])}')
    cfg = dict(spec['defaults'])
    cfg.update(spec['scenarios'][name])
    if smoke:
        cfg.update(spec['smoke'])
    cfg.setdefault('D', cfg['background_states'])
    cfg['name'] = name
    return cfg


def clone_tree(rng, n):
    """Parent of every tumor clone (-1 for the trunk) and the subtree (clone and descendants) of each clone."""
    parent = [-1] + [int(rng.integers(0, t)) for t in range(1, n)]
    subtree = [{t} for t in range(n)]
    for t in range(n - 1, 0, -1):
        subtree[parent[t]] |= subtree[t]
    return parent, subtree


def valid(state, cmax):
    a, b = state
    return 0 <= b <= a and 1 <= a + b <= cmax


def apply_event(rng, states, clones, cmax):
    """Add +-1 copy of one allele to every clone in `clones`; False when no change keeps all of them valid."""
    deltas = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    for index in rng.permutation(len(deltas)):
        da, db = deltas[index]
        changed = [(states[t][0] + da, states[t][1] + db) for t in clones]
        if all(valid(state, cmax) for state in changed):
            for t, state in zip(clones, changed):
                states[t] = state
            return True
    return False


def clonal_distribution(cfg):
    scale = 2 if cfg['wgd'] else 1
    scaled = {(a * scale, b * scale): w for (a, b), w in DIPLOID_STATES.items()}
    states = {state: w for state, w in scaled.items() if valid(state, cfg['cmax_truth'])}
    total = sum(states.values())
    return list(states), np.array([w / total for w in states.values()])


def simulate_states(cfg, rng, subtree):
    """A[m][n], B[m][n] tumor copy numbers and the cluster kinds ('clonal', 'background', 'hot')."""
    m, n, cmax = cfg['clusters'], cfg['tumor_clones'], cfg['cmax_truth']
    choices, weights = clonal_distribution(cfg)
    scale = 2 if cfg['wgd'] else 1
    neutral, second = (scale, scale), (2 * scale, scale)

    kinds = ['clonal'] * m
    others = rng.permutation(np.arange(2, m))
    n_hot = int(round(cfg['frac_hot'] * m))
    if cfg['frac_hot'] > 0 and n > 1:
        n_hot = max(n_hot, 1)
    n_background = int(round(cfg['frac_background'] * m)) if cfg['background_states'] > 1 and n > 1 else 0
    if n_hot + n_background > len(others):
        raise SystemExit('frac_hot + frac_background leave no room for the two clonal reference clusters')
    for index in others[:n_hot]:
        kinds[index] = 'hot'
    for index in others[n_hot : n_hot + n_background]:
        kinds[index] = 'background'

    A = np.zeros((m, n), dtype=int)
    B = np.zeros((m, n), dtype=int)
    for c in range(m):
        if c == 0:
            base = neutral
        elif c == 1:
            base = second
        else:
            base = choices[rng.choice(len(choices), p=weights)]
        for attempt in range(100):
            states = [base] * n
            if kinds[c] == 'background':
                events = int(rng.integers(1, cfg['background_states']))
            elif kinds[c] == 'hot':
                events = min(cfg['hot_states'] - 1, n - 1)
            else:
                events = 0
            for node in rng.choice(np.arange(1, n), size=min(events, n - 1), replace=False) if events else []:
                apply_event(rng, states, sorted(subtree[node]), cmax)
            distinct = len(set(states))
            if kinds[c] != 'hot' or distinct > cfg['D']:
                break
        else:
            raise SystemExit(
                f'could not give hot cluster {c + 1} more than D={cfg["D"]} states; lower D or raise hot_states'
            )
        A[c] = [state[0] for state in states]
        B[c] = [state[1] for state in states]
    return A, B, kinds


def simulate_proportions(cfg, rng):
    """U[n + 1][k]: row 0 is the normal clone, rows 1..n the tumor clones."""
    n, k = cfg['tumor_clones'], cfg['samples']
    low, high = cfg['purity']
    dominant = rng.permutation(np.arange(k)) % n  # each clone dominates at least one sample when k >= n
    U = np.zeros((n + 1, k))
    for p in range(k):
        purity = rng.uniform(low, high)
        weights = rng.dirichlet(np.full(n, cfg['dirichlet_alpha']))
        weights[dominant[p]] += cfg['dominance']
        tumor = purity * weights / weights.sum()
        while True:
            small = (tumor > 0) & (tumor < cfg['minprop'])
            if not small.any():
                break
            tumor[small] = 0.0
            tumor = purity * tumor / tumor.sum()
        U[0, p] = 1.0 - purity
        U[1:, p] = tumor
    return U


def cluster_sizes(cfg, rng):
    m, total = cfg['clusters'], cfg['bins']
    rest = rng.lognormal(0.0, cfg['size_sigma'], size=m - 2)
    others = (1 - cfg['neutral_share'] - cfg['second_share']) * rest / rest.sum()
    shares = np.concatenate([[cfg['neutral_share'], cfg['second_share']], others])
    raw = shares * (total - cfg['min_bins'] * m)
    sizes = np.floor(raw).astype(int) + cfg['min_bins']
    remainder = total - sizes.sum()
    sizes[np.argsort(raw - np.floor(raw))[::-1][:remainder]] += 1
    return sizes


def layout(cfg, rng, sizes):
    """(cluster index, chromosome, start) of every bin: each cluster is split into 1-3 segments placed at random."""
    segments = []
    for c, size in enumerate(sizes):
        pieces = int(min(rng.integers(1, 4), size))
        cuts = np.sort(rng.choice(np.arange(1, size), size=pieces - 1, replace=False)) if pieces > 1 else []
        for length in np.diff(np.concatenate([[0], cuts, [size]])):
            segments.append((c, int(length)))
    order = rng.permutation(len(segments))
    clusters = np.concatenate([[segments[i][0]] * segments[i][1] for i in order])
    total = len(clusters)
    chromosome = 1 + (np.arange(total) * 22) // total
    start = np.zeros(total, dtype=int)
    for chrom in range(1, 23):
        where = np.flatnonzero(chromosome == chrom)
        start[where] = np.arange(len(where)) * cfg['bin_size']
    return clusters, chromosome, start


def simulate(cfg):
    rng = np.random.default_rng(cfg['seed'])
    n, k, _ = cfg['tumor_clones'], cfg['samples'], cfg['clusters']
    parent, subtree = clone_tree(rng, n)
    A, B, kinds = simulate_states(cfg, rng, subtree)
    U = simulate_proportions(cfg, rng)
    sizes = cluster_sizes(cfg, rng)

    # fractional copy numbers (normal cells are diploid (1, 1) also under WGD)
    FA = U[0][None, :] + A @ U[1:]
    FB = U[0][None, :] + B @ U[1:]
    F = FA + FB
    ploidy = (sizes[:, None] * F).sum(axis=0) / sizes.sum()
    rd = F / ploidy[None, :]
    baf = FB / F
    cluster_rd = rd * (1 + rng.normal(0.0, cfg['rd_noise'], size=rd.shape))
    cluster_baf = baf + rng.normal(0.0, cfg['baf_noise'], size=baf.shape)
    cluster_rd[0], cluster_baf[0] = rd[0], baf[0]  # keep the neutral cluster balanced for compute-cn's -td check

    clusters, chromosome, start = layout(cfg, rng, sizes)
    bin_rd = cluster_rd[clusters] * (1 + rng.normal(0.0, cfg['bin_rd_noise'], size=(len(clusters), k)))
    bin_baf = np.clip(cluster_baf[clusters] + rng.normal(0.0, cfg['bin_baf_noise'], size=(len(clusters), k)), 0.0, 1.0)
    bin_rd = np.clip(bin_rd, 0.0, None)

    truth = {
        'parent': parent,
        'A': A,
        'B': B,
        'U': U,
        'kinds': kinds,
        'sizes': sizes,
        'rd': rd,
        'baf': baf,
        'ploidy': ploidy,
        'F': F,
        'FB': FB,
    }
    return truth, (clusters, chromosome, start, bin_rd, bin_baf)


def write_inputs(outdir, cfg, truth, bins):
    clusters, chromosome, start, bin_rd, bin_baf = bins
    samples = [f'Sample{p + 1}' for p in range(cfg['samples'])]
    rows = []
    for p, sample in enumerate(samples):
        beta = np.rint(bin_baf[:, p] * TOTAL_ALLELES).astype(int)
        rows.append(
            pd.DataFrame(
                {
                    '#CHR': [f'chr{c}' for c in chromosome],
                    'START': start,
                    'END': start + cfg['bin_size'],
                    'SAMPLE': sample,
                    'RD': bin_rd[:, p],
                    '#SNPS': SNPS,
                    'COV': COVERAGE,
                    'ALPHA': TOTAL_ALLELES - beta,
                    'BETA': beta,
                    'BAF': bin_baf[:, p],
                    'CLUSTER': clusters + 1,
                }
            )
        )
    bbc = pd.concat(rows).sort_values(
        ['#CHR', 'START', 'SAMPLE'], key=lambda column: column.map(_chrom_key) if column.name == '#CHR' else column
    )
    bbc.to_csv(os.path.join(outdir, 'bulk.bbc'), sep='\t', index=False)

    seg = (
        bbc.groupby(['CLUSTER', 'SAMPLE'])
        .agg(
            **{
                '#BINS': ('RD', 'size'),
                'RD': ('RD', 'mean'),
                '#SNPS': ('#SNPS', 'sum'),
                'COV': ('COV', 'mean'),
                'ALPHA': ('ALPHA', 'sum'),
                'BETA': ('BETA', 'sum'),
                'BAF': ('BAF', 'mean'),
            }
        )
        .reset_index()
        .rename(columns={'CLUSTER': '#ID'})
    )
    seg[['#ID', 'SAMPLE', '#BINS', 'RD', '#SNPS', 'COV', 'ALPHA', 'BETA', 'BAF']].to_csv(
        os.path.join(outdir, 'bulk.seg'), sep='\t', index=False
    )
    return samples


def _chrom_key(value):
    return int(str(value)[3:])


def write_truth(outdir, cfg, truth, samples):
    tdir = os.path.join(outdir, 'truth')
    os.makedirs(tdir, exist_ok=True)
    m, n = truth['A'].shape
    cluster_ids = np.arange(1, m + 1)
    clones = [f'clone{t + 1}' for t in range(n)]
    for matrix in ('A', 'B'):
        frame = pd.DataFrame(truth[matrix], index=pd.Index(cluster_ids, name='cluster'), columns=clones)
        frame.to_csv(os.path.join(tdir, f'{matrix}.tsv'), sep='\t')
    frame = pd.DataFrame(truth['U'], index=pd.Index(['normal'] + clones, name='clone'), columns=samples)
    frame.to_csv(os.path.join(tdir, 'U.tsv'), sep='\t')
    states = [len(set(zip(truth['A'][c], truth['B'][c]))) for c in range(m)]
    pd.DataFrame(
        {'cluster': cluster_ids, 'bins': truth['sizes'], 'kind': truth['kinds'], 'true_states': states}
    ).to_csv(os.path.join(tdir, 'cluster_states.tsv'), sep='\t', index=False)
    hot = [int(c) for c, kind in zip(cluster_ids, truth['kinds']) if kind == 'hot']
    with open(os.path.join(tdir, 'hot_clusters.txt'), 'w') as handle:
        handle.write(','.join(map(str, hot)) + '\n')
    with open(os.path.join(tdir, 'params.json'), 'w') as handle:
        json.dump(dict(cfg, parent=truth['parent']), handle, indent=2)
    return hot, states


def write_runs(outdir, cfg, hot, states):
    """runs.tsv: config, version (new, old or pairwise) and the extra compute-cn arguments of every run."""
    D = cfg['D']
    hot_ids = ','.join(map(str, hot))
    half = ','.join(map(str, hot[: max(1, math.ceil(len(hot) / 2))]))
    oracle = ','.join(f'{c + 1}:{count}' for c, count in enumerate(states))
    available = {
        'free': ('new', []),
        'free_old': ('old', []),
        'd': ('new', ['-d', str(D)]),
        'd_exempt': ('new', ['-d', str(D), '--cnstates-exempt', hot_ids]),
        'd_exempt_partial': ('new', ['-d', str(D), '--cnstates-exempt', half]),
        'map_oracle': ('new', ['--cluster-cnstates', oracle]),
        'd_exempt_pairwise': ('pairwise', ['-d', str(D), '--cnstates-exempt', hot_ids]),
        'd1_exempt': ('new', ['-d', '1', '--cnstates-exempt', hot_ids]),
        'd3_exempt': ('new', ['-d', '3', '--cnstates-exempt', hot_ids]),
    }
    with open(os.path.join(outdir, 'runs.tsv'), 'w') as handle:
        for config in cfg['configs']:
            if config not in available:
                raise SystemExit(f'Unknown configuration {config}; choose from {", ".join(available)}')
            version, args = available[config]
            if '--cnstates-exempt' in args and not hot:
                continue  # nothing to exempt in this scenario
            handle.write(f'{config}\t{version}\t{shlex.join(args)}\n')
    with open(os.path.join(outdir, 'common_args.txt'), 'w') as handle:
        handle.write(shlex.join(['-n', cfg['clones']] + cfg['compute_cn_args']) + '\n')
    with open(os.path.join(outdir, 'modes.txt'), 'w') as handle:
        timelimit = '' if cfg['timelimit'] is None else str(cfg['timelimit'])
        handle.write(' '.join(map(str, cfg['modes'])) + '\n' + timelimit + '\n')


def self_check(cfg):
    """Without noise, the SEG values must reproduce the truth exactly; the truth must satisfy the model's rules."""
    quiet = dict(cfg, rd_noise=0.0, baf_noise=0.0, bin_rd_noise=0.0, bin_baf_noise=0.0)
    truth, bins = simulate(quiet)
    clusters, _, _, bin_rd, bin_baf = bins
    for c in range(cfg['clusters']):
        assert np.allclose(bin_rd[clusters == c].mean(axis=0) * truth['ploidy'], truth['F'][c]), c
        assert np.allclose(bin_baf[clusters == c].mean(axis=0), truth['FB'][c] / truth['F'][c]), c
    U = truth['U']
    assert np.allclose(U.sum(axis=0), 1.0)
    assert ((U[1:] == 0) | (U[1:] >= cfg['minprop'] - 1e-12)).all()
    if cfg['samples'] >= cfg['tumor_clones']:
        assert (U[1:].max(axis=1) > 0).all(), 'a tumor clone is absent from all samples'
    for c, kind in enumerate(truth['kinds']):
        distinct = len(set(zip(truth['A'][c], truth['B'][c])))
        assert all(valid(state, cfg['cmax_truth']) for state in zip(truth['A'][c], truth['B'][c])), c
        if kind == 'clonal':
            assert distinct == 1, c
        elif kind == 'background':
            assert distinct <= cfg['background_states'], c
        else:
            assert cfg['D'] < distinct <= cfg['hot_states'], c
    assert truth['sizes'].sum() == cfg['bins'] and truth['sizes'].argmax() == 0
    print(
        f'self-check OK: {cfg["name"]} ({cfg["clusters"]} clusters, {cfg["samples"]} samples, '
        f'{cfg["tumor_clones"]} tumor clones)'
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('scenarios')
    parser.add_argument('scenario')
    parser.add_argument('outdir', nargs='?')
    parser.add_argument('--smoke', action='store_true', help='shrink the scenario for a quick end-to-end check')
    parser.add_argument('--self-check', action='store_true', help='verify the simulator without writing files')
    args = parser.parse_args()
    cfg = load_scenario(args.scenarios, args.scenario, args.smoke)
    if args.self_check:
        self_check(cfg)
        return
    if args.outdir is None:
        parser.error('outdir is required unless --self-check is given')
    os.makedirs(args.outdir, exist_ok=True)
    truth, bins = simulate(cfg)
    samples = write_inputs(args.outdir, cfg, truth, bins)
    hot, states = write_truth(args.outdir, cfg, truth, samples)
    write_runs(args.outdir, cfg, hot, states)
    print(
        f'{cfg["name"]}: {cfg["clusters"]} clusters ({len(hot)} hot, '
        f'{sum(kind == "background" for kind in truth["kinds"])} background), {cfg["samples"]} samples, '
        f'{cfg["tumor_clones"]} tumor clones, {cfg["bins"]} bins -> {args.outdir}'
    )


if __name__ == '__main__':
    main()
