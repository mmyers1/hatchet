# compute_cn Solver Efficiency

> 2026-10-25: the benchmark and verification scripts of this log moved to `benchmarks/2026-09-29-compute-cn-efficiency/` (shared `measure.py`: `benchmarks/common/measure.py`); job logs stay in `ai-logs/`. See `ai-logs/2026-10-25-CNSTATES_CLUSTER_LIMITS.md`.

## Status

- Overall: Implemented and verified on test fixtures; real-data benchmark script ready, run pending (user)
- Last updated: 2026-09-30

## Problem

`hatchet compute-cn` became much slower after the fixed-matrix (`aaf95e0`), sample-clonal (`d825609`) and
max-homdel (`1d6ab0b`) commits, even when those features are unused.

## Method

Diffed `99db3b4..1d6ab0b` in both backends (C++ `src/ilp-min.cpp`, `src/worker.cpp`, `src/solve.cpp`; Python
`src/hatchet/utils/solve/`).

Findings:

1. New feature code mostly exits early when unused (`fixMatrixEntries`, `fixSampleClonal`, Python `apply_*`).
2. **Root cause:** the symmetry-breaking "repair" turned a no-op row (`sumI <= sumI`) into an active constraint in both
   FULL and CARCH modes. In CARCH (the coordinate-descent C-step) U is fixed, so clones are not interchangeable: the row
   is invalid there and couples every cluster into one dense row, destroying the per-cluster separability of the
   C-step. The C-step runs seeds x iterations times, and `--mode 2` (CD only) is the default.
3. **Latent bug exposed:** hot starts applied the inverse of the sorting permutation. With >= 3 tumor clones the start
   can violate the now-active symmetry row, so Gurobi rejects the MIP start.

Decisions (user): drop symmetry breaking from both models and canonicalize clone order after solving (skipped when the
fixed matrix addresses tumor clones); also add C++ C-step warm starts, complete CD->ILP MIP starts, variable bounds in
place of equality rows, and the fixed-matrix symmetry cleanup.

## Changes

- **No clone-order constraint** in FULL or CARCH models (C++ `IlpSubset::init`, Python `create_model` for standard and
  binwise models). `buildSymmetryBreaking`, `build_symmetry_breaking`, and `checkRank` removed. `symmCoeff` is kept for
  the `cnstates` state-slot ordering, which is a different (valid) symmetry.
- **Identity hot starts** (C++ `hotStart`, Python `hot_start`): columns are no longer permuted. Starts are now complete
  for integer variables: `cA/cB`, bit expansions, amp/del indicators, and, when U is supplied in FULL mode, `u` and the
  presence indicators `x`. Python `first_hot_start` applies the fixed-matrix overlay once instead of `m*n` times.
- **CD -> ILP handoff** passes the CD proportions (`solve.cpp` BOTH mode; Python `'both'` for standard and binwise).
  The C++ worker already warm-started each C-step from the previous iteration's copy numbers, so no change was needed.
- **Canonical clone order after solving:** C++ `canonicalizeCloneOrder` (`src/utils.cpp`) and Python
  `canonicalize_clone_order` (`solve/utils.py`) stable-sort tumor clones by `sum_s (s+1)(A+B)`, permuting A/B columns
  and U rows together. Skipped when `FixedMatrix::addressesTumorClones()` /
  `ResolvedFixedMatrix.addresses_tumor_clones`. Python applies it in `_finalize_solution` (was `_verified_solution`).
- **Removed** the "fully fixed clone burdens violate canonical clone order" validation (C++ and Python); such files are
  now valid.
- **Bounds instead of equality rows:** C++ `computeFixedBounds` and Python `compute_fixed_bounds` derive variable bounds
  from the normal clone (1/1), `--clonal` states, known purities, fixed-matrix entries, sample-clonal CARCH states, and
  UARCH exclusions. Conflicting values fall back to explicit equality rows (`addFixedValueConflicts`,
  `model.fixed_value_conflicts`) so infeasibility behavior is unchanged. Omitted rows: `cA+cB <= upper` when bounds
  imply it (including the former duplicate in CARCH with `cnstates`), and homdel rows when bounds imply them. FULL-mode
  sample-clonal big-M rows are unchanged.
- **Docs:** `docs/source/doc_compute_cn.md` now describes post-solve ordering instead of the removed constraint.

## Tests

- `tests/test_fixed_matrix.py`: fixed-entry checks now inspect variable bounds; new tests for no clone-order rows
  (Python model inspection plus a C++ source guard), identity/complete hot starts, clonal clusters as bounds with omitted
  rows, conflict fallback rows, canonical ordering (Python and the C++ harness `--canonicalize` mode), and
  `_finalize_solution` skipping addressed clones. The burden-order rejection case became an acceptance test. The
  "passes local validation but infeasible" test now uses an all-zero clone with a fixed positive proportion, because
  the old scenario was infeasible only because of the removed row.
- `tests/test_sample_clonal.py`: CARCH/UARCH enforcement checks now inspect bounds.
- `tests/test_solver.py`: the n=3 fixed-matrix solve no longer asserts canonical order (the fixture addresses clone1);
  new licensed test `test_solve_reports_tumor_clones_in_canonical_order` (cpp and gurobi, ILP and CD+ILP, no fixed
  matrix).

## Results

Verification (`hatchet-test` env, `module load gurobi`, native solver rebuilt with CMake and installed to
`src/hatchet/solve`):

- `tests/test_fixed_matrix.py tests/test_sample_clonal.py`: 185 passed.
- `tests/test_solver.py`: 24 passed, 1 failed. The failure is `test_solve_command_cbc`: CBC segfaults (return code -11).
  It fails identically on unchanged HEAD, so it is pre-existing and unrelated.
- All touched C++ files pass `g++ -std=c++11 -Wall -fsyntax-only` with no new warnings.

Small-fixture comparison (`tests/data/fw/bulk`: 10 clusters, 3 samples; `-n 4 -p 20 -j 4 -r 6700 -c 5:1:1`), pre-change
binary vs new binary:

| Measure | Pre-change | New |
|---|---|---|
| CD-only wall time (`-M 2`) | 13.5 s | 12.9 s |
| CD+ILP wall time (`-M 0`) | 28.2 s | 24.3 s |
| Total C-step Gurobi time (`-M 2`, 20 seeds) | 2.55 s (74 steps) | 0.98 s (67 steps) |
| Mean C-step time | 34.5 ms | 14.6 ms |
| CD iterations per seed | 3.70 | 3.35 |
| CD objective | 8.970899 | 8.421963 |
| Full ILP runtime after CD | 14.8 s | 11.4 s |
| Final ILP objective (gap 0) | 8.354420 | 8.354420 |

Wall times on this fixture are dominated by a fixed 2 s sleep per seed in `CoordinateDescent::runWorker` (see follow-up).
Both binaries' outputs have tumor clones in nondecreasing burden order.

## Follow-up opportunity (not implemented)

`CoordinateDescent::runWorker` (`src/coordinate_descent.cpp:336`) sleeps 2 s before every seed; this dates from the
original 2018 code. With the packaged defaults (`seeds = 400`, `jobs = 1`, `clones = 2,8`, both ploidy branches), that
is about 800 s of idle time per clone count and branch, roughly 3 hours per default run. It may have been added to
stagger Gurobi environment/license startup, so it should be tested before removal.

## User benchmark (pending)

The first benchmark command in this log (2026-09-29) did not run. `compute-cn` requires the `-x` running directory to
exist ("Running directory not found!"), the `<BBC_PREFIX>` placeholder was never filled in, and `hatchet` in the
`hatchet-test` env is a stale non-editable install (site-packages 2.1.0), so it would not have run the workspace code
anyway. It is replaced by `benchmarks/2026-09-29-compute-cn-efficiency/bench_compute_cn.sh`, which follows the pattern of
`benchmarks/2026-09-24-preprocessing-efficiency/wgs_compare.sh`:

- The dataset is the `BBC` variable, defaulting to
  `/data1/shahs3/users/myersm2/repos/osteo/results/hatchet2/IID_H134708/bbc/bulk.bbc` (165 clusters, 9 samples). The
  script derives the `-i` prefix and checks that `bulk.seg` exists next to it.
- The defaults mirror `IID_H134708.ini`: `clones=3,5`, `-eD 9 -eT 15`. `SEEDS=100` (packaged default 400) and
  `JOBS=$SLURM_CPUS_PER_TASK`. `SOLVER=cpp` by default; `SOLVER=gurobi` benchmarks the Pyomo backend that the ini uses.
  `MODES=2` (the default CD-only mode); `MODES="2 0"` adds CD + ILP, together with a `TIMELIMIT`.
- The old code is `git archive 4baaf666` into `<WORK>/code_old`; the new code is the working tree. Both native solvers
  are built out of source into `<WORK>/build_old` and `<WORK>/build_new`. Each run uses `PYTHONPATH` on its own code
  and an explicit solver path.
- Each run gets a fresh `<WORK>/mode<M>/<version>` that is both `-x` and the working directory. This keeps a stray
  `hatchet.ini` from being picked up and keeps each version's `coordinate_descence.tsv` separate.
  `HATCHET_COMPUTE_CN_SOLVER` forces the solver.
- The preflight validates the exact arguments with `compute-cn`'s own `parsing_arguments` without solving.
- `compare_compute_cn.py` writes `<WORK>/summary.md`. It reports wall/CPU/peak memory (via the existing
  `benchmarks/common/measure.py`), C/U step counts and times, the objective per ploidy and clone
  count, the chosen model, and whether the chosen solution's clones are in canonical order.

Commands (from `benchmarks/2026-09-29-compute-cn-efficiency/`):

```bash
PREFLIGHT_ONLY=1 bash bench_compute_cn.sh      # ~20 s: builds both solvers, validates arguments, no solving
sbatch bench_compute_cn.sh                      # the benchmark; job log: hatchet-cn-bench-<jobid>.log
# variants
sbatch --export=ALL,SOLVER=gurobi bench_compute_cn.sh
sbatch --export=ALL,MODES="2 0",TIMELIMIT=600 bench_compute_cn.sh
sbatch --export=ALL,BBC=/path/to/other/bbc/bulk.bbc bench_compute_cn.sh
```

Verification and follow-up:

- The job log ends with the summary tables (also in `/data1/shahs3/users/myersm2/hatchet-cn-bench/IID_H134708/summary.md`).
  Both runs should have exit 0.
- `new` should have a lower wall time and a lower mean C-step time than `old`.
- Objectives per ploidy and clone count should be equal or lower for `new`. CD is a heuristic, so an occasional small
  increase for one clone count is possible, but the typical case should improve.
- The "clone order of best" column should say `sorted` for `new`.
- If the chosen model differs between versions, compare the objective tables: model selection depends on the
  objectives, so a better objective can change the choice.
- To re-summarize without rerunning: `$PY compare_compute_cn.py <WORK> <SOLVER>`, where `$PY` is the `hatchet-test`
  Python.

Checks run (2026-09-30):

- `PREFLIGHT_ONLY=1` with the osteo default: both solvers built and both argument sets were accepted (19 s).
- Full smoke run on `tests/data/fw/bulk.bbc` (`SEEDS=4 CLONES=3,4 EXTRA_ARGS="-eD 6 -eT 12 -g 0.35 -l 0.6"`): both
  versions exited 0 and the summary was complete. Both chose diploid with 3 clones, and both outputs are in canonical
  order. The mean C-step time was 43.7 ms old vs 16.8 ms new. Objectives were identical except diploid n4 (4.689353 old
  vs 4.678702 new).
- Wall times from that smoke run (265 s vs 242 s) are not meaningful. The host (`isca082`, 4 cores, load average ~27)
  was oversubscribed: a direct solver call spent 44.6 s in system CPU and 0.17 s in user CPU. Run the benchmark with
  `sbatch`.

## Update log (every change and its reasoning)

### 2026-09-29: solver changes

C++ solver

1. `src/ilp-min.cpp`, `IlpSubset::init`: stopped calling `buildSymmetryBreaking()`.
   *Why:* this was the root cause. In CARCH (fixed U) clones are not interchangeable, so the row is invalid and couples
   all clusters into one dense row. In FULL mode it replaced the no-op the code had always effectively run with.
2. `src/ilp-min.cpp`: deleted `buildSymmetryBreaking()` and `checkRank()`; kept `symmCoeff()`.
   *Why:* they are unused once the constraint is gone. `symmCoeff()` still orders the `cnstates` state slots, which is a
   valid, different symmetry.
3. `src/ilp-min.cpp`, `IlpSubset::hotStart`: removed the rank/`sort_indexes` column remap and used the columns as given.
   *Why:* the remap applied the inverse permutation, which could produce infeasible starts. It is also meaningless
   without an order constraint, and wrong in CARCH, where U fixes clone identity.
4. `src/ilp-min.cpp` and `src/ilp-min.h`, `hotStart`: added the optional `hU` argument. The start now also sets bit
   variables, amp/del indicators, and in FULL mode `u` and `x`.
   *Why:* a start that fixes every integer variable is completed by an LP instead of a sub-MIP, so Gurobi accepts it
   faster.
5. `src/ilp-min.cpp`, `src/ilp-min.h`: added `computeFixedBounds()`, the `_lbA/_ubA/_lbB/_ubB/_lbU/_ubU` members, and
   the `FixedValue` struct. `buildVariables()` creates `cA`, `cB` and `u` with those bounds.
   *Why:* fixed cells (normal clone, `--clonal`, purities, fixed-matrix entries, sample-clonal states and exclusions)
   become variable bounds instead of equality rows, which presolve removes for free.
6. `src/ilp-min.cpp`: added `addFixedValueConflicts()`.
   *Why:* when two sources fix one cell to different values, the old equality rows made the model infeasible. Re-adding
   the conflicting equality keeps that behavior exactly and never throws, which would kill worker threads.
7. `src/ilp-min.cpp`: deleted `fixGivenCN()` and `fixMatrixEntries()`. `fixSampleClonal()` now keeps only the
   FULL-mode big-M rows.
   *Why:* those rows are now bounds (item 5).
8. `src/ilp-min.cpp`, `buildConstraints`: removed the `cA[s][0]==1` and `cB[s][0]==1` rows. The `cA+cB <= upper` row is
   now added once, in FULL and CARCH, and only when the bounds do not already imply it. The homdel row `cA+cB >= 1` is
   skipped when the lower bounds imply it.
   *Why:* these rows are trivially satisfied by fixed cells; CARCH with `cnstates` used to add the upper row twice.
9. `src/ilp-min.cpp`, `src/ilp-min.h`: added `upperCN(s)`.
   *Why:* the same upper-bound expression had been repeated in four places.
10. `src/solve.cpp`: in CD+ILP mode, `hotStart` now receives `cd.getU()`.
    *Why:* this gives the full ILP a complete start from the CD result.
11. `src/solve.cpp`: after any solve mode, calls `canonicalizeCloneOrder` unless `fixedMatrix.addressesTumorClones()`.
    *Why:* the user chose to drop the order constraint and sort afterwards; clones addressed by the fixed matrix keep
    their labels.
12. `src/utils.h`, `src/utils.cpp`: added `canonicalizeCloneOrder`, a stable sort of tumor clones by
    `sum_s (s+1)(A+B)` that permutes the A/B columns and U rows together.
    *Why:* this is the post-solve replacement for the constraint.
13. `src/fixed_matrix.h`: added `FixedMatrix::addressesTumorClones()`.
    *Why:* it decides when sorting would break user-assigned clone labels.
14. `src/fixed_matrix.cpp`: removed the "fully fixed clone burdens violate canonical clone order" validation.
    *Why:* without an order constraint, such files are valid.
15. `src/worker.cpp`, `src/coordinate_descent.cpp`: no change.
    *Why:* the worker already warm-starts each C-step from the previous iteration's copy numbers.

Python solver

16. `src/hatchet/utils/solve/ilp_subset.py`, `create_model` in both `ILPSubset` and `ILPSubsetSplit`: removed the
    `build_symmetry_breaking` call.
    *Why:* same as item 1.
17. Same file: added `compute_fixed_bounds()` and `copy_number_upper_bound()`. `cA`, `cB` and `u` are created with the
    resulting bounds, and the `u[0] == 1 - purity` rows, `cA[m][0] == 1` rows and `fix_given_cn` rows are gone.
    *Why:* this mirrors items 5, 7 and 8.
18. Same file: added `apply_fixed_value_conflicts()` (the `model.fixed_value_conflicts` list), replacing
    `apply_fixed_matrix_constraints()`.
    *Why:* this mirrors item 6.
19. Same file: added `build_copy_number_limits()`, a single upper row that is skipped when implied (it replaces the rows
    in the bit block and the CARCH block). `build_homdel_constraints()` skips rows the bounds imply.
    *Why:* this mirrors item 8.
20. Same file: `apply_sample_clonal_constraints()` now only builds the FULL-mode big-M rows.
    *Why:* the CARCH states and UARCH exclusions are bounds.
21. Same file: `hot_start(f_a, f_b, u=None)` uses identity mapping, and `get_rank` was removed. It sets bits, amp/del
    indicators, and in FULL mode with `u`, the `u` and `x` values. `create_model` keeps `self._bitcA/_bitcB/_adA/_adB/_x`
    for this.
    *Why:* this mirrors items 3 and 4.
22. Same file: `first_hot_start` applies the fixed-matrix overlay once, after the loops, instead of `m*n` times.
    *Why:* it was redundant repeated work with the same result.
23. Same file: removed `build_symmetry_breaking`, `apply_fixed_matrix_constraints`, `fix_given_cn` and the unused
    `purities` local.
    *Why:* they were dead after the changes above.
24. `src/hatchet/utils/solve/utils.py`: added `canonicalize_clone_order`.
    *Why:* this mirrors item 12.
25. `src/hatchet/utils/solve/__init__.py`: `_verified_solution` became `_finalize_solution`, which verifies the user
    constraints and then sorts clones unless the fixed matrix addresses tumor clones. The CD+ILP paths pass the CD `u`
    to `hot_start`.
    *Why:* this mirrors items 10 and 11. The CD result is handed to the ILP unsorted, so C and U stay consistent.
26. `src/hatchet/utils/fixed_matrix.py`: added `ResolvedFixedMatrix.addresses_tumor_clones` and removed the
    fully-fixed burden-order validation.
    *Why:* this mirrors items 13 and 14.

Documentation

27. `docs/source/doc_compute_cn.md`: replaced the paragraph describing the enforced clone order and its validation.
    It now describes post-solve sorting, which is skipped when a tumor clone is addressed.
    *Why:* the documented behavior no longer existed.

Tests

28. `tests/test_fixed_matrix.py`:
    - Fixed-entry assertions now read variable bounds (`_fixed_entry_values`).
      *Why:* the entries are bounds now.
    - The two symmetry tests were replaced by a Python model check (no constraint touches copy numbers of two tumor
      clones in two clusters) and a C++ source guard.
      *Why:* they asserted the removed row.
    - The burden-order rejection case became an acceptance test.
      *Why:* item 14.
    - Hot-start tests now cover identity mapping and the `u`, `x` and bit starts.
      *Why:* items 3 and 4.
    - New tests cover clonal clusters as bounds with omitted rows, and the conflict fallback rows.
      *Why:* item 5 changes behavior that needs coverage.
    - New canonical-order tests (Python, C++ harness, `_finalize_solution`).
      *Why:* items 11 and 12 are new code paths.
    - The "passes local validation but infeasible" test now uses an all-zero clone with a fixed positive proportion.
      *Why:* the old scenario was infeasible only because of the removed row.
29. `tests/test_sample_clonal.py`: the CARCH and UARCH enforcement tests read bounds.
    *Why:* item 20.
30. `tests/test_solver.py`:
    - The n=3 fixed-matrix solve no longer asserts canonical order and was renamed.
      *Why:* its fixture addresses clone1, so labels are kept.
    - New licensed test `test_solve_reports_tumor_clones_in_canonical_order`.
      *Why:* it covers item 11 end to end for both backends.
    - `run_fixed_matrix_solve` and `assert_solve_output` accept no fixed matrix.
      *Why:* the new test needs them to.
31. `tests/fixed_matrix_cpp_test.cpp`: added a `--canonicalize` mode and an `addressesTumorClones()` check.
    *Why:* these are native unit checks for items 12 and 13.

Build artifact

32. `src/hatchet/solve` (gitignored): rebuilt from the working tree and installed. The pre-change binary was backed up
    in the session scratchpad.
    *Why:* the licensed tests use the packaged binary.

### 2026-09-30: benchmark fix

33. Added `benchmarks/2026-09-29-compute-cn-efficiency/bench_compute_cn.sh`.
    *Why:* the benchmark command in this log failed (missing running directory, placeholder input, stale `hatchet`
    install). The script takes the dataset as `BBC` (default: the osteo `IID_H134708` BBC), uses that dataset's
    `compute_cn` settings, builds both solvers, validates arguments before solving, and runs both versions on one node.
34. Added `benchmarks/2026-09-29-compute-cn-efficiency/compare_compute_cn.py`.
    *Why:* it turns the runs into one comparable summary (time, memory, step times, objectives, chosen model, clone order).
35. This log: rewrote "User benchmark", updated "Status", and added this section.
    *Why:* the old instructions were wrong, and the user's instructions require an exhaustive change list.

## Summary

The slowdown came from the symmetry-breaking fix activating a constraint in the coordinate-descent C-step, where it is
invalid and couples all clusters, and from an inverted hot-start permutation that the active constraint exposed. Both
backends now solve without clone-order constraints, warm-start with complete starts, enforce fixed values through
variable bounds, and report interchangeable clones in canonical order after solving. On the test fixture, C-steps are
about 2.4x faster, CD converges in fewer iterations to a better objective, and the full ILP reaches the same optimum
faster.
