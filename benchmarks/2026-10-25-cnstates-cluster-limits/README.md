# Synthetic benchmark: per-cluster limits on tumor copy-number states

This benchmark measures what `-d/--cnstates D` costs and what it changes, alone and with `--cnstates-exempt` or
`--cluster-cnstates`:
- runtime and memory
- model size
- which number of clones compute-cn selects
- accuracy against a known truth

The implementation and its rationale are in `ai-logs/2026-10-25-CNSTATES_CLUSTER_LIMITS.md`.

## Files

| File | Purpose |
|---|---|
| `scenarios.json` | Scenario parameters (`defaults`, one override block per scenario, `smoke` shrinks any scenario) |
| `simulate.py` | Writes `bulk.bbc`/`bulk.seg`, the truth (`truth/`), and the runs of a scenario; `--self-check` verifies the simulator |
| `bench_cnstates.sh` | SLURM array driver: one task per scenario; simulates, builds, runs compute-cn per configuration under `benchmarks/common/measure.py` |
| `summarize_cnstates.py` | Markdown summary of every run in `WORK` (also written to `WORK/summary.md` by each task) |

## Scenarios

The parameters are in `scenarios.json`. The noise levels are cluster- and bin-level standard deviations.

| Index | Scenario | Clusters | Samples | True tumor clones | Truth | Question |
|---|---|---|---|---|---|---|
| 0 | `lowhet` | 120 | 12 | 3 | 90% clonal, 7% ≤2 states, 3% hot (3 states) | Does the limit remove spurious clones in the motivating case? Formulation and `D` sweep, modes 2 and 0. |
| 1 | `lowhet_noisy` | 120 | 12 | 3 | as `lowhet`, 2× noise | The same, with more pressure to overfit |
| 2 | `moderate` | 150 | 10 | 4 | 30% ≤2 states, 10% hot (≤4) | How it behaves when real heterogeneity is wider |
| 3 | `nohot` | 120 | 10 | 3 | all ≤2 states | Overhead when there is nothing to exempt |
| 4 | `wgd` | 120 | 10 | 3 | as `lowhet`, tetraploid (cmax 12, 4 bits) | Tetraploid path |
| 5 | `scale` | 300 | 20 | 5 | 15% ≤3 states, 5% hot | Memory and time stress, and old vs new formulation; modes 2 and 0 |

## Configurations

`D` is the scenario's `background_states` (2, or 3 for `scale`). Hot clusters are simulated with more than `D` states.

| Config | compute-cn arguments | What it shows |
|---|---|---|
| `free` | none | Baseline |
| `free_old` | none, run with the code at `9083d8be` | No regression without `-d` (same objectives; similar time) |
| `d` | `-d D` | A global limit forces the hot clusters to compromise |
| `d_exempt` | `-d D --cnstates-exempt <hot>` | The intended use |
| `d_exempt_partial` | `-d D --cnstates-exempt <half of hot>` | A misspecified list |
| `map_oracle` | `--cluster-cnstates <true count per cluster>` | Upper bound on the benefit |
| `d_exempt_pairwise` | as `d_exempt`, with the old pairwise formulation (native only) | Old vs new formulation, everything else equal |
| `d1_exempt`, `d3_exempt` | `-d 1` / `-d 3` with the true exemptions | Sensitivity to `D` |

All runs use the following settings (in `scenarios.json` and the driver):
- `--diploid -eD 6`, or `--tetraploid -eT 12` for `wgd`
- `-ts 0`, so no cluster is dropped for size
- `-td 0.03`
- `-n 2,6`, or `2,8` for `scale`
- 100 seeds, with `-j` set to the task's CPUs

## Running

Run these from this directory. Nothing here submits jobs by itself.

```
PREFLIGHT_ONLY=1 bash bench_cnstates.sh          # minutes on a login node: simulate, build both solvers, check arguments
SMOKE=1 sbatch --array=0-5 bench_cnstates.sh     # tiny inputs, 4 seeds, WORK=...-smoke: end-to-end check
sbatch --array=0-5 bench_cnstates.sh             # full benchmark (8 CPUs, 32 GB, up to 2 days per task)
sbatch --array=5 --mem=64G bench_cnstates.sh     # resubmit 'scale' alone with more memory if needed
```

Useful overrides, passed with `--export=ALL,...`:
- `CONFIGS="d d_exempt"`
- `MODES=2`
- `SEEDS=400`
- `SOLVER=gurobi` for the Pyomo path; it skips `d_exempt_pairwise`
- `WORK=/path`

Each task appends nothing to the repo. Its outputs are:
- `WORK/<scenario>/input/`: inputs, truth, `runs.tsv`
- `WORK/<scenario>/mode<M>/<config>/`: compute-cn results, `log.err`, `time.json`, `coordinate_descence.tsv`
- `WORK/summary.md`

## Verification and follow-up

1. `PREFLIGHT_ONLY=1` must end with `pre-flight OK`.
2. After the smoke run, `WORK-smoke/summary.md` must have one row per configuration with `exit 0`. In the hot scenarios,
   the `d` row's `hot states` should be at most `D` (1 in smoke), and `d_exempt` should show more.
3. In the full summary, check these:
   - **Resource cost:** `mean C step`, `mean C rows / discrete` and `peak RSS` of `d` against `free`, and of
     `d_exempt_pairwise` against `d_exempt`. The pairwise formulation should need roughly 10× the C-step rows.
   - **Model selection:** in `lowhet`/`lowhet_noisy`, `chosen`/`tumor clones` of `d_exempt` and `map_oracle` should be
     closer to the true count than `free` if spurious clones are the issue.
   - **Accuracy:** `frac MAE`, `exact CN`, and `over`. Over-counted genome fraction is the signature of spurious
     subclonal states.
   - **Misspecification:** `d_exempt_partial` against `d_exempt` shows the cost of missing exemptions: `under` at the
     non-exempt hot clusters.
   - **Regression:** `free_old` against `free` should give the same objectives per clone count. Time may differ slightly.
4. Re-summarize without rerunning: `$PY summarize_cnstates.py $WORK`.
