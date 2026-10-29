# Per-Cluster Limits on Tumor Copy-Number States (`--cnstates-exempt`, `--cluster-cnstates`)

## Status

- Overall: implemented and verified here wherever no Gurobi licence is needed. The native solve tests and the
  synthetic benchmark are for you to run (see "Your next steps").
- Last updated: 2026-10-25
- Branch `osteo`, nothing committed.
- Approved plan: `/home/myersm2/.claude/plans/constraining-number-of-copy-number-encapsulated-curry.md`

## Goal

The problem: for patients with many samples but little heterogeneity, compute-cn adds clones that differ only at a few
clusters. The fix is to limit most clusters to a few distinct tumor copy-number states (`-d/--cnstates D`) while
designated clusters may use more:
- `--cnstates-exempt CLUSTERS` gives the listed clusters no limit, up to n−1 tumor states.
- `--cluster-cnstates ID:D,...` gives each listed cluster its own limit.

Only one of the two may be used.

Decisions confirmed with you:
- Implement both forms, and make it an error to use them together.
- Fix the existing `--cnstates` bugs and reformulate the constraint as slot states.
- Add heavier synthetic benchmarks for you to run.
- Put all benchmark code in `benchmarks/` at the repo root, one subdirectory per task, including the code moved from
  ai-logs.

## Method

### Relation to `--cnstates`

`-d/--cnstates D` already limited the distinct (A, B) states across the tumor clones of every cluster. The normal clone
is not counted, so `D = 1` makes every cluster clonal. The new options make that limit per-cluster. The "deprecated"
label came from upstream HATCHet (the first Sphinx docs commit `aa1b5306`). It is now removed, because the bugs below
were the likely reason the option was abandoned.

### Existing bugs (file:line refer to commit `9083d8be`)

1. **Out-of-bounds hot start.** `firstHotStart` (`src/ilp-min.cpp:403-460`) and `first_hot_start`
   (`ilp_subset.py:566-599`) indexed length-`d` target arrays with `i % min(n,k)`. Whenever `d < min(n,k)`, which is
   the normal case, that was an out-of-bounds read in C++ and an `IndexError` in Python. Python's
   `np.random.choice(..., d, replace=False)` also raised when `d > k`.
2. **Python over-constraint.** The symmetry-breaking row `_sum_l <= _sum_l1` (`ilp_subset.py:412`, and its copy in
   `ILPSubsetSplit`) was added inside the clone loop. The resulting prefix-sum rows force tumor clones 1 and 2 into the
   same state in every cluster, which effectively removes one clone.
3. **Cost.** The pairwise bit-equality formulation added about `4·M·d·(n−1)(n−2)/2` rows per cluster. In CARCH it also
   built all `2·M·m·n` copy-number bit binaries. Together that made each coordinate-descent C-step about 11× larger.

### New design

- **One code path.** C++ uses a per-segment array `cnstates[s]` (`-1` means no limit), built by `getCnstates` in
  `solve.cpp`. Python uses the global `d` plus an override map `cluster_cnstates` (cluster ID → limit, or None for no
  limit), built by `resolve_cluster_cnstates` in `solve/__init__.py`. Both replace the scalar `d` everywhere: model,
  coordinate descent, hot start and validators.
- **When a limit applies.** `stateLimit(s)` / `state_limit(cluster_id)` returns 0 for clonal clusters and for a limit
  of at least n−1, which cannot bind. The small-n runs of a clone sweep therefore get no extra variables.
- **Slot formulation.** For a limited cluster with limit D:
  - continuous slot states `sA`/`sB` in `[0, upperCN]`, one per slot
  - an assignment binary `z[t][l]` for each tumor clone t and each slot `l < min(t, D)`
  - `Σ_l z = 1`
  - `|c − s_l| ≤ U·(1 − z)` for A and B
  - restricted-growth symmetry breaking, `z[t][l] ≤ Σ_{t'<t} z[t'][l−1]`, which numbers slots by first use and is exact
  - copy-number bits only in FULL mode
- **Hot start.**
  - `firstHotStart` gives a limited cluster at most D states: clone i starts from the ((i−1) % D)-th of D random
    samples. Unlimited clusters keep the old behaviour exactly, and draw no random numbers.
  - `hotStart` / `hot_start` set the `z` and slot starts by first use. They skip a cluster whose start has more than D
    states, which a fixed-matrix overlay can cause.
- **Benchmark-only native switch.** `--cnstates-formulation pairwise|slot` (default `slot`, global
  `g_cnstatesFormulation`) keeps the old pairwise rows and bits behind the fixed hot start. It exists for benchmarking
  and the equivalence test, is not exposed by `hatchet compute-cn`, and can be deleted later.
- **Model sizes.** `coordinate_descence.tsv` rows gain three trailing columns: rows, integer and binary variables, and
  nonzeros. The FULL ILP size is logged at verbosity 2 ("ILP model size: ...").

## Results

### Tests run here

The environment was `hatchet-test` with `PYTHONPATH=src`, `module load gurobi`, `HATCHET_TEST_MIP_SOLVER=cbc` and
`HATCHET_TEST_CPP_SOLVER=build-native/solve`.

- `tests/test_cnstates_limits.py` (new, 46 tests): 40 pass. The 6 native solve tests fail only because the Gurobi
  licence server is unreachable from the sandbox.
- **Watched failing first (TDD):**
  - All new Python tests failed before the change, with `TypeError: unexpected keyword argument 'cluster_cnstates'`.
  - The bug-2 regression scenario (d=2, tumor clones 1 and 2 must differ) was also run on the old code exported from
    `HEAD`: objective 0.9, with c2 forced to A=1 in both clones. The new code reaches objective 0.
- **C++ hot start, checked with a scratch driver under `-D_GLIBCXX_ASSERTIONS`** (no ASan runtime on this host):
  - The old `firstHotStart` with d=1, n=4, k=9 aborts on an out-of-bounds vector access (exit 134).
  - The new one gives 1, 3 and 2 tumor states for limits 1, none and 2.
  - The driver is at
    `/tmp/claude-164077095/-data1-shahs3-users-myersm2-repos-hatchet/5a20d02f-48af-40fa-a94c-a4fe1ae6f155/scratchpad/hotstart/check.cpp`.
- **Wider run** of `test_cnstates_limits`, `test_ampdel_exempt`, `test_fixed_matrix`, `test_sample_clonal`,
  `test_solver` and `test_config`: 225 passed, 21 skipped, 42 failed. All 42 failures were checked by script, and every
  one is a Gurobi licence or token-server error. None is an assertion about the changed behaviour.
- **Not run:** `test_steps*.py`, `test_phase.py` and `test_counting_equivalence.py`. They exercise preprocessing, not
  touched here, and are heavy.
- **Native build:** `cmake --build build-native -j 4` succeeds, and `solve -h` lists the new options.

### Benchmark tooling

- `simulate.py --self-check` passes for all six scenarios at full and smoke size. Without noise, the SEG values
  reproduce the truth exactly.
- A simulated `lowhet` input goes through compute-cn's parsing, filtering and neutral-cluster detection: all 120
  clusters are kept, cluster 1 is neutral, and `-d 2 --cnstates-exempt …` is forwarded to the native solver. It stops
  only at the licence.
- `summarize_cnstates.py` was checked on fake runs:
  - Ground truth with permuted clone labels gives MAE 0, exact CN 1.000 and no over- or under-count.
  - Two clones collapsed into one gives under-count 0.088 (the hot cluster's weight), hot states 1/2, and exact CN
    0.956.
- `PREFLIGHT_ONLY=1 bash bench_cnstates.sh` with WORK in the scratchpad ends in `pre-flight OK`. It built both
  solvers and validated 54 argument sets.
- `PREFLIGHT_ONLY=1` of the moved `bench_compute_cn.sh`, with the repo's `tests/data/fw` input, exits 0.
- All moved scripts pass `bash -n`.

### Model-size impact

The worked example has 162 clusters, 9 samples, n=8 (7 tumor clones), cmax 6 (3 bits) and D=3. The counts are rows and
discrete variables added to each solve.

| | FULL ILP (modes 0/1) | CARCH step of CD (mode 2, default) |
|---|---|---|
| Baseline, no `-d` | ~210k rows, ~8.2k binaries | ~15.5k rows, 2.6k int + 324 bin |
| Old `-d 3` | +159k rows (+75%), +3.4k bin | +161k rows (~11×), +11.2k bin (bits + z) |
| New `-d 3`, no exemptions | +14.6k rows (+7%), +2.8k bin, 1k continuous slots | +14.6k rows (~2×), same variables |
| Each exempt cluster (new) | −90 rows, −17 bin | same |

Runtime depends on the data, so the benchmark measures it. The new `coordinate_descence.tsv` columns let it compare
measured C-step sizes with these estimates.

## Progress

- [x] F. `benchmarks/` directory, with the existing ai-logs benchmark code moved into it and paths fixed.
  `.gitignore` re-includes `benchmarks/**/*.sh`.
- [x] A–D. Native solver: options, per-segment limits, slot formulation, pairwise switch, hot start, validators,
  model-size output
- [x] A–D. Python/Pyomo: per-cluster limits, slot formulation (shared by both classes), hot start, validators,
  `solve()` label resolution, CLI and ini
- [x] Tests: the Python tests, the C++ validator-driver tests, and the native rejection tests pass. The native solve
  tests are written but need your licensed run.
- [x] E. Docs and ini
- [x] G. Synthetic benchmark suite: simulator, scenarios, driver, summarizer, README; preflight passes
- [x] Self-review of the whole diff (no subagent review; see the rulings)
- [ ] Your runs: native solve tests, smoke benchmark, full benchmark

## Your next steps

1. **Native and Gurobi tests**, on a node that reaches the licence server. About 2 minutes.

   ```
   cd /data1/shahs3/users/myersm2/repos/hatchet
   module load gurobi
   cmake --build build-native -j 4
   HATCHET_TEST_CPP_SOLVER=$PWD/build-native/solve PYTHONPATH=src \
     conda run -n hatchet-test python -m pytest tests/test_cnstates_limits.py tests/test_ampdel_exempt.py \
     tests/test_fixed_matrix.py tests/test_sample_clonal.py tests/test_solver.py -q
   ```

   Expected: everything passes. These six tests were not run here:
   - `test_native_solver_enforces_the_global_limit[*]`
   - `test_native_solver_frees_exempt_and_bounds_per_cluster_limits`
   - `test_native_pairwise_and_slot_formulations_reach_the_same_optimum`
   - `test_native_coordinate_descent_reports_model_sizes`

   If any fails, send me the `-x -vv` output.
2. **Smoke benchmark** (minutes per task):

   ```
   cd benchmarks/2026-10-25-cnstates-cluster-limits
   PREFLIGHT_ONLY=1 bash bench_cnstates.sh
   SMOKE=1 sbatch --array=0-5 bench_cnstates.sh
   ```

   Check `/data1/shahs3/users/myersm2/hatchet-cnstates-bench-smoke/summary.md`:
   - every row exits 0
   - in the `d` rows, hot states are at most 1
3. **Full benchmark:**

   ```
   sbatch --array=0-5 bench_cnstates.sh
   ```

   Rerun `scale` (index 5) with `--mem=64G` if it runs out of memory. What to read in `summary.md` is listed in that
   directory's README, under "Verification and follow-up".
4. **On real data**, for example the osteo patient that motivated this:

   ```
   hatchet compute-cn ... -d 2 --cnstates-exempt <IDs>
   ```

   or `--cluster-cnstates ID:D,...`. Then compare the chosen clone count and the per-cluster states against your
   current run.
5. **Reinstall** hatchet in the env your pipeline uses (`pip install .`), so that `hatchet compute-cn` and its packaged
   `solve` binary pick up the change.

## Rulings / decisions made during implementation

- Work happens in the current checkout on branch `osteo`, not in a new worktree. The installed `hatchet-test` package
  is non-editable, so tests already run with `PYTHONPATH=src` against this checkout, and earlier features were done the
  same way. Nothing is committed. Cost if wrong: none.
- Python keeps the global `d` plus the `cluster_cnstates` override map instead of a materialized per-segment array. The
  Pyomo classes and the validators are keyed by cluster ID, and existing callers pass `d`, so the behaviour matches C++
  with a smaller change. Cost if wrong: none functionally.
- Under a limit, the hot start takes A and B for each slot from one random sample. The old code shuffled A and B
  independently. Pairs from one sample are coherent allele states. Cost if wrong: different, equally valid start values.
- Slot states are continuous, not integer. A used slot equals an integer copy number through the big-M rows, and an
  unused slot then needs no branching. The plan said "integer". Cost if wrong: none for feasibility; possibly a slightly
  different branching path.
- The pairwise switch is a process-wide global set in `solve.cpp`, not a constructor argument. That keeps the
  benchmark-only code out of four constructor signatures, so it is easy to delete. Cost if wrong: a little more coupling
  until it is removed.
- The pairwise formulation sets no `z` start, as before. It reproduces the old formulation faithfully, apart from the
  fixed hot start. Cost if wrong: the pairwise benchmark rows look slightly worse than an "optimized old" version would.
- The solving tests take their Pyomo solver from `HATCHET_TEST_MIP_SOLVER` (default `gurobi`), because the sandbox
  cannot reach the Gurobi token server. I ran them with CBC. Cost if wrong: none.
- Test fixture: samples hold 10% normal cells, not 0%. With a zero normal proportion, the normal clone's copy-number
  variables do not appear in the LP, and Pyomo returns `None` for them. That quirk predates this work. Cost if wrong:
  none.
- Simulator: there is no rejection sampling for cluster separation. compute-cn's defaults `mergerdr = mergebaf = 0`
  merge only identical clusters, and noise makes every cluster unique. The benchmark passes `-ts 0`, so no cluster is
  dropped for size, and `-td 0.03`. The `lowhet` preflight check confirmed all 120 clusters were kept. Cost if wrong: a
  scenario run with non-default merge settings could merge clusters.
- Benchmark ploidy is fixed per scenario: `--diploid`, or `--tetraploid` for `wgd`. This isolates clone-count selection
  from ploidy selection. Cost if wrong: the benchmark does not test ploidy choice. Edit `compute_cn_args` in
  `scenarios.json` to include it.
- Benchmark run directories, logs and summaries are kept outside the repo (default
  `/data1/shahs3/users/myersm2/hatchet-cnstates-bench`, as for the 2026-09-29 benchmark). Only code is under
  `benchmarks/`. Cost if wrong: none.
- Not moved to `benchmarks/`:
  - job logs, which stay in ai-logs next to their task log
  - `ai-logs/2026-09-24-preprocessing-efficiency/temp/`, a gitignored scratch copy of the tests
  - `ai-logs/2026-10-19-READ_FLAG_FILTERS/regenerate_fixtures.sh`, which is fixture tooling, not a benchmark

  Cost if wrong: one `git mv` each.
- There is no fresh-context subagent review, because subagents are not spawned unless you ask for one. The final review
  was a self-review: the C++ diff, the restricted-growth index ranges, the hot-start paths in all modes, and the CLI
  error paths. Cost if wrong: an author's blind spot. Ask for a `/code-review` if you want a second pass.

## Summary

- **Best way to implement it:** make the existing `--cnstates` limit per-cluster. `--cnstates-exempt` frees the listed
  clusters, and `--cluster-cnstates ID:D` gives them their own limit. Both work in both solvers and in every mode.
- **Bugs fixed:** two `--cnstates` bugs made the option unusable. One was an out-of-bounds hot start in both solvers;
  the other was a Python over-constraint that merged tumor clones 1 and 2.
- **Reformulation:** the constraint is now slot states with exact symmetry breaking. That is about 9–10× fewer rows
  than before, and coordinate-descent steps no longer need bit variables.
- **Runtime and memory:** with the new formulation, `-d` adds about 7% to the exact ILP and roughly doubles a
  coordinate-descent step's row count. The old formulation added about 75% and about 11×. Exempt clusters cost
  nothing.
- **Benchmarks:** all benchmark code now lives in `benchmarks/`. The new synthetic suite measures time, memory, model
  size, clone-count selection and accuracy against a known truth, and is ready to submit.

## Updates made (exhaustive, with reasons)

### Native solver

1. `src/utils.h` / `src/utils.cpp`:
   - Added `CNSTATES_FORMULATION_t`, the global `g_cnstatesFormulation` (default slot), and `segmentCnstates()` (the
     per-segment limit or -1).
   - Added `#include <numeric>` for `std::iota`.
   - Reason: per-segment limits and the benchmark switch.
2. `src/argparse.{h,cpp}`:
   - New options: `--cnstates-exempt`, `--cluster-cnstates` (ID:INT pairs, positive, no duplicates) and the hidden
     `--cnstates-formulation pairwise|slot`.
   - New errors: both forms used together, and `--cnstates-exempt` without `-d`.
   - Updated the help and `toString`. The `-d` help no longer claims a default of 3.
   - Reason: the user-facing options.
3. `src/solve.cpp`:
   - Added `getCnstates()`, which builds the per-segment array from `-d` and the new options and rejects unknown
     cluster IDs.
   - Sets `g_cnstatesFormulation`, and passes the array everywhere `argparse.d()` was used.
   - Logs the ILP model size at verbosity 2.
   - Reason: the single code path and model-size reporting.
4. `src/ilp-min.{h,cpp}`:
   - `d` became `cnstates`. Added `stateLimit()`, `anyStateLimit()` and `needsBits()`.
   - New slot variables and constraints (`buildOptionalVariables`, `buildOptionalConstraints`, `buildSlotConstraints`).
   - Moved the old rows into `buildPairwiseConstraints`, kept for the switch.
   - Added slot starts (`startSlots`) and the rewritten bug-free `firstHotStart`.
   - Added the model-size getters `numConstrs`, `numDiscreteVars` and `numNonzeros`.
   - Reason: bug 1, the reformulation and the size reporting.
5. `src/worker.{h,cpp}` and `src/coordinate_descent.{h,cpp}`:
   - Take the per-segment array.
   - `coordinate_descence.tsv` gets rows, discrete variables and nonzeros per step.
   - The `snprintf` size argument was 2014 for a 1024-byte buffer; it is now 1024.
   - Reason: plumbing, benchmark metrics, and an overflow-safe buffer.
6. `src/fixed_matrix.{h,cpp}` and `src/sample_clonal.{h,cpp}`: the validators take the per-segment array and check
   each cluster against its own limit. Reason: section D.

### Python

7. `src/hatchet/utils/solve/ilp_subset.py`:
   - Added `cluster_cnstates` to `ILPSubset` and `ILPSubsetSplit`, including `__copy__`.
   - Added `state_limit()`, plus a shared `build_cnstates_constraints()` and `start_cnstates()`.
   - Removed the duplicated buggy z blocks from both classes.
   - Bits are built only in FULL mode.
   - Rewrote `first_hot_start`.
   - Reason: bugs 1–2, the reformulation, and removing the duplication.
8. `src/hatchet/utils/solve/__init__.py`:
   - Added `resolve_cluster_cnstates()`, which matches labels to seg IDs as strings. Its errors: unknown ID, both
     forms, exempt without `d`, non-positive limit.
   - `solve()` takes `cnstates_exempt` and `cluster_cnstates`, and passes the resolved map to the validators and to
     every ILP and coordinate-descent construction.
   - Reason: label resolution, done as for `ampdel_exempt`.
9. `src/hatchet/utils/solve/cd.py`: `CoordinateDescent` and `CoordinateDescentSplit` accept and forward
   `cluster_cnstates`. Reason: plumbing.
10. `src/hatchet/utils/fixed_matrix.py` and `src/hatchet/utils/sample_clonal.py`: the validators take
    `cluster_cnstates`, and each cluster uses its own limit. Reason: section D.
11. `src/hatchet/bin/HATCHet.py`:
    - New `--cnstates-exempt` and `--cluster-cnstates` arguments, with parsing and validation.
    - A SEG-ID check, forwarding to the native solver in `makeBaseCMD`, and pass-through to the Pyomo path in
      `execute_python`.
    - Reason: the CLI.
12. `src/hatchet/hatchet.ini`: added `cnstatesexempt =` and `clustercnstates =`. Reason: ini defaults and the
    `HATCHET_COMPUTE_CN_*` env overrides.

### Docs and tests

13. `docs/source/doc_compute_cn.md`:
    - Removed "deprecated" and documented the tumor-clone semantics.
    - Added both new options to the table.
    - Added a section on when limits apply, the INI keys, solver coverage, the resource impact, and model selection.
    - Reworded the two later `--cnstates` mentions.
    - Reason: section E.
14. `tests/test_cnstates_limits.py` (new): model structure, hot start, CBC/Gurobi solves, validators, `solve()`, CLI,
    `makeBaseCMD`, docs/ini, C++ validator drivers, native rejections and native solves. Reason: verification.
15. `tests/fixed_matrix_cpp_test.cpp` and `tests/sample_clonal_cpp_test.cpp`: adapted to the per-segment API, with new
    `cnstates` and `cnstates-exempt` scenarios. Reason: the signature change and coverage.

### Benchmarks

16. Created `benchmarks/` and moved the earlier code into it:
    - `git mv` of `measure.py` (to `benchmarks/common/`), `compare_compute_cn.py` and `ampdel_usage.py`.
    - `mv` of the untracked `env.sh`, `compare.sh`, `run_tests.sh`, `hg38_phase_check.sh`, `wgs_compare.sh` and
      `bench_compute_cn.sh`.
    - Updated the hardcoded `HERE`, `MEASURE` and `source .../env.sh` paths. `env.sh` gained `MEASURE`.
    - Reason: your request to organize benchmarks by task.
17. `benchmarks/2026-09-29-compute-cn-efficiency/compare_compute_cn.py` reads only the first 6 columns of
    `coordinate_descence.tsv`. Reason: the file now has 9 columns.
18. `.gitignore`: added `!benchmarks/**/*.sh`. Reason: the global `*.sh` rule would leave the benchmark scripts
    untracked, which is why the earlier ones never were tracked.
19. `ai-logs/2026-09-24-preprocessing-efficiency.md`, `ai-logs/2026-09-29-compute-cn-efficiency.md` and
    `ai-logs/2026-10-22-AMPDEL_EXEMPT_CLUSTERS.md`: updated the script paths and added a "moved on 2026-10-25" note.
    Reason: keep the logs' commands runnable.
20. `benchmarks/README.md` (new): an index of the subdirectories and their logs.
21. `benchmarks/2026-10-25-cnstates-cluster-limits/` (new): `scenarios.json`, `simulate.py`, `bench_cnstates.sh`,
    `summarize_cnstates.py` and `README.md`. Reason: section G and your request for heavier synthetic benchmarks.
22. This log.

## Tasks prevented by sandbox-limited file access

- Entries such as `.bashrc`, `.gitconfig`, `.idea` and `.mcp.json` appear at the repo root inside Bash commands. They
  are `/dev/null` character devices that the sandbox bind-mounts over deny-listed paths, not real files. Workaround:
  none needed. They are not added, read, or deleted.
- Gurobi solves (Pyomo `gurobi` and the native solver) fail with "Failed to connect to token server
  'islogin01.mskcc.org' (port 41954)". The sandbox network blocks the licence server, including a retry with that host
  declared. Workarounds:
  - The Pyomo solving tests ran with CBC.
  - The bug-2 regression was confirmed on the old code, exported with `git archive` into `$TMPDIR`, also with CBC.
  - The C++ hot-start fix was checked with a licence-free driver.
  - The native solve tests and the benchmarks are on your run list.
- The default benchmark WORK directories (`/data1/shahs3/users/myersm2/hatchet-*-bench*`) are not writable from the
  sandbox. Workaround: both preflights ran with WORK in the session scratchpad. `bench_compute_cn.sh` used the repo's
  `tests/data/fw` input, because the default patient BBC is outside the readable paths.
- No ASan runtime (`libasan.so`) is installed. Workaround: bounds checking with `-D_GLIBCXX_ASSERTIONS`.
- flake8 and blue (the repo's pre-commit linters) are not installed in any env. Workaround: a manual 120-column check
  and `py_compile`.

## Side effect to note

The existing native-solver tests (`test_sample_clonal.py`, `test_solver.py`) run the solver from the repo root. The
native solver truncates `./coordinate_descence.tsv` at the start of every coordinate-descent run. The untracked
`coordinate_descence.tsv` in the repo root is now empty, last written 2026-10-25 14:50. If it held output from an
earlier run that you need, it is lost, because it was untracked and not backed up. `tests/out/solver/` was also
rewritten by `test_solver.py`, as on every run of that test.
