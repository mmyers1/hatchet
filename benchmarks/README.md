# Benchmarks and verification scripts

Each subdirectory belongs to one task and is named after that task's log in `ai-logs/`, which records the plan, the
results and the commands. The scripts are meant to be run by hand or submitted to SLURM by the user. None of them
submits jobs by itself.

| Directory | Task log | Contents |
|---|---|---|
| `common/` | | `measure.py`: runs a command and records wall/CPU time and peak RSS (process tree, Python processes, largest process) as JSON |
| `2026-09-24-preprocessing-efficiency/` | `ai-logs/2026-09-24-preprocessing-efficiency.md` | `env.sh` (shared settings), `compare.sh`, `run_tests.sh`, `hg38_phase_check.sh` and `wgs_compare.sh`: old-vs-new preprocessing output and performance comparisons |
| `2026-09-29-compute-cn-efficiency/` | `ai-logs/2026-09-29-compute-cn-efficiency.md` | `bench_compute_cn.sh` and `compare_compute_cn.py`: compute-cn on real data, old vs new solver |
| `2026-10-22-ampdel-exempt-clusters/` | `ai-logs/2026-10-22-AMPDEL_EXEMPT_CLUSTERS.md` | `ampdel_usage.py`: per-cluster amp+del usage in compute-cn outputs |
| `2026-10-25-cnstates-cluster-limits/` | `ai-logs/2026-10-25-CNSTATES_CLUSTER_LIMITS.md` | Synthetic benchmark of `--cnstates`, `--cnstates-exempt` and `--cluster-cnstates`; see its `README.md` |

Job logs and other run outputs stay next to their task log in `ai-logs/`.

The scripts use absolute paths, so that a copy made by `sbatch` still finds its files. Edit `REPO`, `PY` and the data
paths at the top of a script if the repo moves.
