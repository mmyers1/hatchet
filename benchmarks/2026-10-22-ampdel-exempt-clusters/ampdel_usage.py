"""Summarize per-cluster amp+del usage in compute-cn outputs (base = 1 copy per allele, diploid)."""
import sys
import pandas as pd

bench, base = sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 1
configs = ['default', 'exempt1', 'exempt16', 'noampdel']
profiles = {}
for config in configs:
    bbc = pd.read_table(f'{bench}/{config}/results.diploid.n4.bbc.ucn.tsv')
    clone_cols = [c for c in bbc.columns if c.startswith('cn_clone')]
    weight = bbc.groupby('CLUSTER').size() / len(bbc)
    states = bbc.groupby('CLUSTER')[clone_cols].first()
    rows = []
    for cluster, row in states.iterrows():
        alleles = [tuple(int(v.split('|')[k]) for v in row) for k in (0, 1)]
        ampdel = any(min(a) < base < max(a) for a in alleles)
        rows.append((cluster, ampdel, ' '.join(row.values)))
    df = pd.DataFrame(rows, columns=['cluster', 'ampdel', 'profile']).set_index('cluster')
    df['weight'] = weight
    profiles[config] = df
    used = df[df.ampdel]
    print(
        f'{config:9s} clusters using amp+del: {len(used):3d}  '
        f'genome fraction: {used.weight.sum():.4f}  ids: {list(used.index)[:20]}'
    )

for config in configs[1:]:
    changed = profiles[config].profile != profiles['default'].profile
    print(
        f'{config:9s} clusters whose clone profile '
        f'differs from default (clone labels may permute): {int(changed.sum())}'
    )
print('cluster 25:', {c: profiles[c].loc[25, 'profile'] for c in configs})
