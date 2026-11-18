# Sample-Clonal Coordinate Descent Never Uses More Than k Tumor Clones

## Status

- Overall: Implemented; native (cpp) tests and the osteo rerun pending (user to run)
- Last updated: 2026-10-09

## Goal

With `sampleclonal`, solutions that allow more clones do not use them even when they could improve the fit. In the
osteo run `IID_H135990` (3 samples, `clones=4,8`, CD-only `M: 2`, Pyomo `solver=gurobi`), the objective is
35.80704747074… for every clone count from 4 to 8. Every solution has the same three tumor proportions
(T1: 0.7048 + 0.1806, T03: 0.7868, T04: 0.5556), and the extra clone columns have `u = 0` in every sample
(`results/results.diploid.n{4..8}.seg.ucn.tsv`).

## Method

### Root cause

With sample-clonal input, coordinate-descent seeds come from `ILPSubset.build_sample_clonal_u` (Python) and
`CoordinateDescent::buildSampleClonalU` (C++). Both give each sample's entire residual purity to exactly one tumor
clone (`U[assignment[s]][s] = residual[s]`). As a result:

1. Every seed has at most k present tumor clones (k = number of samples), regardless of n.
2. CARCH ties copy numbers to the data only for present clones (`_fixed_u > tol` in the fA/fB sums). Absent clones'
   A/B are left arbitrary and never get optimized.
3. In UARCH, `compute_fixed_bounds` forces `U = 0` for any clone whose arbitrary A/B break a sample's states. Even when
   such a clone is allowed, its unoptimized profile does not help the fit, so it never enters.

CD is therefore limited to the best solution with at most k tumor clones, and n = 4…8 converge to the same answer.
Without sample-clonal input, `build_random_u` spreads each sample over up to `size_bubbles` random clones, so all
clones are exercised. The 2026-09-03 sample-clonal implementation log gives no reason for the one-clone design.

The run used the Python/Pyomo path: the ini sets `solver=gurobi`, and `solver: None` in the log is the unused native
solver path argument. The native backend has the same seed design.

### Fix

Keep the existing backtracking assignment, which guarantees one compatible clone per sample. After it, spread each
sample's residual over a random number of compatible tumor clones:

- The target part count is `max(randint(1, n-1), randint(1, n-1))`, which mirrors `build_random_u`. It is capped at
  `floor((residual + tol) / mu)` so every share can be at least `mu`.
- Extra clones use the same candidate predicate as the assignment step, so a clone is never shared by incompatible
  samples.
- The shares are split by placing uniform cut points over the residual above the `mu` floor.

## Progress

- [x] Root cause identified from code and run output.
- [x] Plan approved.
- [x] Failing Python seed test written and observed failing.
- [x] Python fix implemented; seed test passes.
- [x] C++ fix implemented; compiles and builds.
- [x] End-to-end regression test (both backends) added. The Gurobi/Pyomo variant was run and fails before the fix.
- [x] Cheap verification run (Python backends).
- [ ] User runs the native (cpp) tests. The sandbox cannot reach the Gurobi token server.
- [ ] User reinstalls and reruns osteo IID_H135990.

## Results

All runs use the `hatchet-test` environment with `PYTHONPATH=src`.

- **RED:** `pytest tests/test_sample_clonal.py -k seed`, before the fix.
  - The new `test_sample_clonal_seeds_spread_purity_over_additional_clones` failed with `assert 1 >= 2`: across
    50 seeds, no sample ever had more than one tumor clone.
- **RED, end-to-end:** the unfixed Python seed code was checked out from HEAD into a scratchpad copy of the package.
  - The new `test_compute_cn_coordinate_descent_mixes_tumor_clones_with_sample_clonal_states[gurobi]` failed with
    `u_clone1 = 0.8, u_clone2 = 0.0`. The second clone is never used even though only a two-clone mixture fits
    cluster 3, which reproduces the osteo symptom.
- **GREEN:** after the fix:
  - the seed tests pass (4 passed, including the two existing sample-clonal seed tests and the fixed-matrix seed
    test);
  - the end-to-end mixture test passes on the Gurobi/Pyomo backend, with calibrated gamma 3.0 and purity 0.8 as
    designed;
  - the existing `conflicting_sample_clonal_states[gurobi-{two,three}-states]` tests still pass.
- **Suites:** `pytest tests/test_sample_clonal.py tests/test_fixed_matrix.py tests/test_solver.py -k "not cpp and not
  native"` gave **170 passed, 1 skipped, 48 deselected** (36 s).
  - The skip is `test_solve_binary`, whose `GRB_LICENSE_FILE` check fails.
  - The deselected tests need the native binary, which needs the token-server license (see Sandbox-Limited Access).
- **C++:**
  - `g++ -std=c++11 -Wall -Wextra -fsyntax-only` on `src/coordinate_descent.cpp` reports no new warnings; the
    sign-compare warnings are all on untouched lines.
  - A full CMake Release build into the scratchpad succeeded.
- **Real-data seed coverage**, solver-free: seeds were built from the osteo `bbc/bulk.seg` with the run's
  `sampleclonal` spec, n = 8, mu = 0.1, 400 seeds. The derived purities 0.8853 / 0.7868 / 0.5556 match the run log.

  | | Tumor clones per seed (min / median / max) | Max clones in one sample | Seeds using all 7 tumor clones |
  |---|---|---|---|
  | HEAD | 1 / 3 / 3 | 1 | 0 / 400 |
  | Fixed | 5 / 7 / 7 | 7 | 346 / 400 |

- `git diff --check` is clean, and no added line exceeds 120 characters.

### Commands for the user

Native backend: rebuild and install first, because the test fixture prefers the packaged `src/hatchet/solve`.

```bash
cd /data1/shahs3/users/myersm2/repos/hatchet
module load gurobi
GUROBI_ROOT="$GUROBI_HOME"; case "$(basename "$GUROBI_ROOT")" in linux64|mac64) GUROBI_ROOT="$(dirname "$GUROBI_ROOT")";; esac
cmake -S . -B build-native -DCMAKE_BUILD_TYPE=Release -DGUROBI_HOME="$GUROBI_ROOT"
cmake --build build-native --config Release
install -m 755 build-native/solve src/hatchet/solve
# expect 4 passed: mixture [gurobi, cpp] and conflicting [cpp-two-states, cpp-three-states]
conda run -n hatchet-test env PYTHONPATH="$PWD/src" pytest -q tests/test_sample_clonal.py \
  -k "mixes_tumor_clones or (conflicting_sample_clonal_states and cpp)"
conda run -n hatchet-test env PYTHONPATH="$PWD/src" pytest -q \
  tests/test_sample_clonal.py tests/test_fixed_matrix.py tests/test_solver.py
```

Osteo rerun: the `hatchet` env holds a non-editable install, so reinstall first.

```bash
module load gurobi && conda run -n hatchet pip install /data1/shahs3/users/myersm2/repos/hatchet
sbatch /data1/shahs3/users/myersm2/repos/osteo/results/hatchet2/2sample/IID_H135990/run_IID_H135990_heavy.sh
# after it finishes:
grep "Best objective found" /data1/shahs3/users/myersm2/repos/osteo/results/hatchet2/2sample/IID_H135990/hatchet.log
```

Expected on the rerun:

- The objectives for n = 5..8 are no longer identical, and are typically lower than n = 4. The n = 4 objective may also
  improve, because a sample can now be split across clones in more seeds.
- More than 3 tumor clones can have `u > 0` in `results/results.diploid.n*.seg.ucn.tsv`.
- Model selection may pick a different n.
- CD per-seed runtime rises somewhat, because more present clones means more free integer variables in each C-step.

## Summary

With `sampleclonal`, every coordinate-descent seed put each sample's whole tumor purity on one clone. Copy numbers are
optimized only for present clones, and the U-step excludes clones whose arbitrary copy numbers break a sample's
states. CD therefore could never use more than k tumor clones (k = number of samples), and n = 4..8 converged to the
same solution.

Seeds now spread each sample over a random number of compatible clones, as non-sample-clonal seeds already do:

- a share of at least `mu` per clone;
- a clone is never shared by incompatible samples;
- fixed U entries are kept.

Both backends were changed, and a seed unit test and a two-backend end-to-end test guard the behavior.

## Updates and Reasoning

1. `tests/test_sample_clonal.py::test_sample_clonal_seeds_spread_purity_over_additional_clones` (new).
   - Setup: 50 seeds from `_conflicting_sample_clonal_ilp(n=6)`.
   - Per seed, asserts: fixed normal proportion; unit column sums; nonzero tumor shares ≥ mu - tol; no clone shared by
     the incompatible samples.
   - Across seeds, asserts: some sample uses ≥ 2 clones; all 5 tumor clones are used.
   - Reason: a solver-free reproduction of the root cause. It failed before the fix.
2. `tests/test_sample_clonal.py::_write_mixture_solver_inputs` (new) and
   `test_compute_cn_coordinate_descent_mixes_tumor_clones_with_sample_clonal_states[gurobi|cpp]` (new).
   - Data: a one-sample, three-cluster input calibrated to gamma 3 / purity 0.8. Cluster 3 is a 0.6 / 0.2 mixture of
     (1,1) and (3,1), so its f_A = 1.4 has no single-clone fit.
   - The test runs CD-only (`--mode 2`, `-n 3`) through the public CLI and requires both tumor clones to be present
     and different at cluster 3.
   - Reason: an end-to-end regression on both backends. The Gurobi variant was confirmed to fail before the fix.
3. `src/hatchet/utils/solve/ilp_subset.py::ILPSubset.build_sample_clonal_u`: after the existing backtracking
   assignment (unchanged), each sample, in shuffled order, adds up to `n_parts - 1` extra clones.
   - `n_parts = max(randint(1, n), randint(1, n))`, capped at `floor((residual + tol) / mu)` so every share can be at
     least mu.
   - Extra clones are drawn from the same candidate set as the assignment step (not fixed for that sample, compatible
     with every sample already on the clone).
   - Reason: the root-cause fix. It mirrors the clone-count distribution of `build_random_u`, and feasibility is
     preserved because every sample already has its base clone.
4. `src/hatchet/utils/solve/ilp_subset.py::ILPSubset._split_sample_purity` (new).
   - Each part gets `mu`. The remaining purity is cut at sorted uniform points, and the last share absorbs rounding so
     the column sums exactly.
   - Reason: random shares that `fix_u` and the minimum-proportion rule accept. It is inherited by `ILPSubsetSplit`,
     so binwise mode is covered.
5. `src/coordinate_descent.cpp::CoordinateDescent::buildSampleClonalU`:
   - The per-clone candidate test of `assignSample` moved into a `cloneAvailable` lambda, so the assignment and the new
     extension step share one rule.
   - Added the same extension step as (3), using `g_rng`.
   - Reason: native parity with (3). Native seeds were also limited to one clone per sample.
6. `src/coordinate_descent.cpp::CoordinateDescent::splitSamplePurity` (new) and its declaration in
   `src/coordinate_descent.h`. Reason: native parity with (4).
7. This log. Reason: the CLAUDE.md logging requirement.

No production files other than (3) to (6) changed. Scratch artifacts are outside the repository, in the session
scratchpad:

- a HEAD copy of the package, used for the RED end-to-end run;
- the native build;
- `seed_coverage.py`.

## Sandbox-Limited Access

1. At first, Bash could not list `/data1/shahs3/users/myersm2/repos/osteo/results/hatchet2/2sample/IID_H135990/`
   ("No such file or directory"), because the read sandbox denies `/data1` outside the repository.
   - Workaround: the directory was then added as a working directory.
   - Files read: the run's `hatchet.log`, `IID_H135990.ini` and `run_IID_H135990_heavy.sh`; the `u_*` columns of
     `results/results.diploid.n{4,5,6,8}.seg.ucn.tsv`; and `bbc/bulk.seg` (cluster RDR/BAF only, for the
     solver-free seed coverage check).
2. Native (cpp) solve tests could not run.
   - Cause: the native binary links the site Gurobi 9.5.2 library, whose license
     (`/admin/software/gurobi/gurobi.lic`) needs the token server `islogin01.mskcc.org:41954`, and the sandbox blocks
     that connection ("Failed to connect to token server").
   - Also, `module load gurobi` does not set `GUROBI_HOME` or `GRB_LICENSE_FILE` inside the sandbox.
   - Workaround: the binary was built against `/admin/software/gurobi/gurobi952` in the scratchpad, to confirm it
     compiles and links. The cpp tests are left for the user (commands above).
   - The Pyomo/Gurobi tests ran without `GRB_LICENSE_FILE`, which reaches no token server.

## Out of Scope (noted, not changed)

- The Python seed candidate predicate lacks the C++ `cloneCompatible` check against fixed-matrix A/B entries. This gap
  is pre-existing; seeds that hit it are already skipped as infeasible.
