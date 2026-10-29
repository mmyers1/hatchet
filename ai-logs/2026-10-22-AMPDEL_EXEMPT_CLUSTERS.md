# Per-Cluster Exemptions From the Amp-Del Constraint (`--ampdel-exempt`)

> 2026-10-25: the benchmark and verification scripts of this log moved to `benchmarks/2026-10-22-ampdel-exempt-clusters/` (shared `measure.py`: `benchmarks/common/measure.py`); job logs stay in `ai-logs/`. See `ai-logs/2026-10-25-CNSTATES_CLUSTER_LIMITS.md`.

## Status

- Overall: Implemented and tested in both solvers. The user ran the resource benchmark on real data with the
  reinstalled native solver, and the results are below. The native pytest solve tests (which need the licence) have not
  been reported yet.
- Last updated: 2026-10-22

## Goal

`compute-cn -f/--noampdel` turns the amp-del assumption off for every cluster or for none. The assumption requires
each allele of each non-clonal cluster to be all ≥ base or all ≤ base across tumor clones. This change adds
`--ampdel-exempt CLUSTERS` (ini key `ampdelexempt`): the assumption stays on everywhere except the listed clusters,
which may then hold both amplified and deleted states.

Decisions confirmed with the user:
- A new comma-separated flag; `-f/--noampdel` stays a plain boolean.
- `--ampdel-exempt` together with `--noampdel` is an error.

## Method

A cluster is subject to amp-del only if all three hold: `ampdel` is on, the cluster is not a clonal (`-c`) cluster,
and it is not exempt. Exempting a clonal cluster is allowed and does nothing. Labels are matched to the seg `#ID` as
strings, and an unknown label is an error.

The exemption applies at every site where the constraint is used, in both solvers:
- Python/Pyomo (`ILPSubset`, `ILPSubsetSplit`):
  - `adA`/`adB` variable creation
  - amp-del constraints
  - `first_hot_start`
  - `validate_fixed_matrix`
  - `validate_sample_clonal_constraints`
- Native C++/Gurobi (`IlpSubset`):
  - `fixC` assertions
  - `hotStart` start values for `adA`/`adB`
  - `buildVariables`
  - `buildConstraints`
  - `firstHotStart`
  - `FixedMatrix::validate`
  - `SampleClonal::validate`

Plumbing follows commit `1d6ab0b7` (max-homdel-size).

## Resource impact (analysis)

Each non-clonal cluster under amp-del has 2 binaries (`adA`, `adB`) and `4·(n−1)` linear constraints. Exempting a
cluster removes exactly those and adds nothing; exempting every cluster is the same as `--noampdel`.

Worked example: IID_H134708 has 162 clusters and 9 samples. With n=4 (normal + 3 tumor clones), diploid, cmax=6
(M=3 bits):

| Model (per solve) | Discrete vars | of which `ad*` | Constraints (approx.) | of which amp-del |
|---|---|---|---|---|
| FULL ILP (modes 0/1) | ~4.2k binary + 1.3k int | 322 (≈7.6% of binaries) | ~117k | 1,932 (≈1.6%) |
| CARCH step of CD (mode 2, U fixed) | 1.3k int + 322 binary | 322 (≈20%) | ~11k | 1,932 (≈17%) |
| Per exempted cluster | −2 | | −12 | |

Search space per exempted cluster, counting one allele's values across 3 tumor clones (ignoring A+B ≤ cmax):
`(base+1)^3 + (cmax−base+1)^3 − 1` with amp-del versus `(cmax+1)^3` without it.
- Diploid: 223 → 343 (×1.5).
- Tetraploid with cmax=12: 1357 → 2197 (×1.6).

Expected consequences:
- **Memory:** essentially unchanged; the model shrinks by ≤ 2% even with every cluster exempt.
- **CD (mode 2, the default):** with U fixed, CARCH splits into independent per-cluster blocks, so an exemption only
  enlarges that one cluster's small subproblem. Runtime should stay between the default and `--noampdel`. UARCH has
  no amp-del terms and is unaffected.
- **ILP (modes 0/1):** clusters are coupled through U, so the feasible set grows across exempt clusters. The big-M
  amp-del LP relaxation is weak, though, and there are fewer binaries to branch on. The net effect can only be
  measured empirically.
- **Objective:** an exemption relaxes the problem, so the optimal objective is ≤ the default (CD is a heuristic, so
  usually there too). That can shift the chosen clone count.

## Progress

- [x] Python CLI / driver (`HATCHet.py`, both ini files)
- [x] Python solver (`solve/__init__.py`, `ilp_subset.py`, `cd.py`, both validators)
- [x] Native solver (argparse, solve.cpp, ilp-min, worker, coordinate_descent, both validators)
- [x] Tests: `tests/test_ampdel_exempt.py` plus one C++-validator test each in `test_fixed_matrix.py` and
  `test_sample_clonal.py`
- [x] Docs (`docs/source/doc_compute_cn.md`)
- [x] Verification (Python, Pyomo/Gurobi and native validation paths)
- [ ] Native end-to-end pytest solve tests (user run; needs the Gurobi licence). Native solves with exemptions did
  run successfully in the benchmark.
- [x] Benchmark (user run, 2026-10-22)

## Results / Verification

Environment: `hatchet-test` with `PYTHONPATH=src`, because the env has a non-editable install. The Pyomo/Gurobi tests
use gurobipy's bundled size-limited licence.

- `tests/test_ampdel_exempt.py`: 17 passed and 3 skipped. The skipped tests are native solves that need
  `GRB_LICENSE_FILE`.
  - Against the old installed package (no `PYTHONPATH=src`), 11 of the original 13 tests failed, as expected. The
    other 2 assert that the flag is absent.
  - The Pyomo end-to-end test ran with Gurobi. With `ampdel=True`, fixed entries `A[cluster1] = (3, 0)` make the
    model infeasible; with `cluster1` exempt it solves and returns `(3, 0)`.
  - The native validation tests ran against the rebuilt `build-native/solve`. They cover: amp+del fixed entries
    rejected without an exemption, an unknown exempt ID rejected, and `--ampdel-exempt` without `-f` rejected.
- C++ validator drivers (compiled with g++ by pytest): the `ampdel-exempt` scenario passes for the fixed-matrix
  driver and the sample-clonal driver. The sample-clonal `ampdel` scenario still rejects the same input.
- Regression set (`test_ampdel_exempt.py`, `test_fixed_matrix.py`, `test_sample_clonal.py`, `test_config.py`,
  `test_solver.py`): 217 passed, 5 skipped, 20 failed.
  - All 20 failures are native (`cpp`) solver runs that abort with `GRBException` because no Gurobi licence is
    available in the sandbox.
  - The same failure sets occur on a pristine `git archive HEAD` copy using the same prebuilt
    `src/hatchet/solve` binary: 10/10 for `test_sample_clonal.py` and 10/10 for `test_solver.py`. They are not
    caused by this change.
- The native build (`cmake --build build-native -j 4`) succeeded, and `solve -h` lists `--ampdel-exempt CLUSTERS`.
- Formatting: `black -S -l 120 --diff` (a stand-in for the repo's `blue` hook) reports no new differences in any
  touched Python file.

### Resource benchmark (user run, 2026-10-22)

Setup: IID_H134708 (162 clusters, 9 samples), `--diploid -n 4 -p 20 -j 4 -r 1 -M 2` (coordinate descent only),
native solver from the reinstalled `hatchet-test` package, and one run per configuration.
- `exempt1` exempts cluster 25.
- `exempt16` exempts clusters 1–16, which are the 16 largest because the IDs are already ordered by bin count.

| Config | Wall time | vs default | Peak RSS (MB) | Objective | Clusters using amp+del | Genome fraction using amp+del |
|---|---|---|---|---|---|---|
| default | 4:32.6 (272.6 s) | — | 451.8 | 114.607383 | 0 | 0 |
| exempt1 (25) | 4:21.9 (261.9 s) | −3.9% | 453.0 | 114.607383 | 0 | 0 |
| exempt16 (1–16) | 3:48.7 (228.7 s) | −16.1% | 452.6 | 114.590693 | 7 (IDs 5, 8, 9, 10, 12, 13, 14) | 0.134 |
| noampdel | 3:34.5 (214.5 s) | −21.3% | 453.6 | 104.194705 | 120 | 0.499 |

Amp+del usage was counted with `benchmarks/2026-10-22-ampdel-exempt-clusters/ampdel_usage.py`: an allele uses amp+del if
its values across tumor clones fall both below and above base 1. Command:
`conda run -n hatchet-test python -I benchmarks/2026-10-22-ampdel-exempt-clusters/ampdel_usage.py <bench dir> 1`.

Interpretation:
- **Memory:** unchanged, as predicted. Peak RSS spans 451.8–453.6 MB (0.4%).
- **Runtime:**
  - As predicted, every exemption run fell between default and `--noampdel`.
  - Runtime dropped steadily as more clusters were exempted, so in coordinate-descent mode removing the constraint
    made the solver faster, not slower.
  - The −3.9% for one cluster is probably within run-to-run noise; there was one replicate and the node was shared.
  - −16% for removing about 10% of the amp-del binaries (16 of 161 constrained clusters) is more than proportional.
    Fewer coordinate-descent iterations to converge is a plausible cause, but the logs at this verbosity cannot show
    it, and replicates would be needed to separate it from noise.
- **Objective:** the predicted order `noampdel ≤ exempt16 ≤ exempt1 ≤ default` held.
- **Correctness on real data:** in `exempt16`, every cluster that uses amp+del is in the exempt set, and no
  non-exempt cluster violates the constraint.
- **Cluster 25 (exempt1):**
  - The solution is identical to the default (same objective, same clone profiles). With every other cluster
    constrained, the fit did not need amp+del for cluster 25 at n=4.
  - In the unconstrained run, cluster 25 does use amp+del (A allele 0, 0, 2 across clones), but there the whole clone
    structure is different.
  - If cluster 25 is expected to be amp+del, try more seeds or ILP refinement (`-M 0`) with the exemption.
- **Model selection:** `--noampdel` uses amp+del in 120 of 162 clusters (half the genome) for a 9% lower objective.
  Exempting everything therefore gives the factorization much more freedom; a targeted exemption list keeps the
  solution close to the constrained one (exempt16 changed the objective by 0.015%).

### Follow-up for the user

1. Native end-to-end tests, which need the licence:
   ```bash
   module load gurobi
   cmake --build build-native -j 4
   HATCHET_TEST_CPP_SOLVER=$PWD/build-native/solve PYTHONPATH=src \
     conda run -n hatchet-test python -m pytest tests/test_ampdel_exempt.py tests/test_sample_clonal.py tests/test_solver.py -q
   ```
   Expected: everything passes, including the 3 `test_native_solver_allows_amplification_and_deletion_in_exempt_cluster`
   modes and the 20 native tests that fail here only for lack of a licence.
2. Reinstall hatchet (`pip install .` in the env your pipeline uses) so that `hatchet compute-cn` and the packaged
   `hatchet/solve` binary pick up `--ampdel-exempt`.
3. Benchmark done (see above). Optional: run 2–3 replicates with different `-r` seeds to measure timing noise, and
   repeat `default`/`exempt16` in ILP mode (`-M 1 -s 900`) if ILP modes are used.

## Summary

`compute-cn --ampdel-exempt 3,17` (or `ampdelexempt = 3,17` in `[compute_cn]`) keeps the amp-del assumption for
every cluster except the listed ones, in both the Pyomo and native solvers. For listed clusters:
- the `adA`/`adB` binaries and their `4·(n−1)` constraints are not built
- hot starts do not clamp them toward one direction
- fixed-matrix and sample-clonal validation no longer reject amp+del requirements for them

Invalid use fails early:
- unknown IDs, in the driver against the SEG file and in `solve()`/`solve.cpp` against the cluster index
- an empty list
- combining the option with `--noampdel` (or, natively, using it without `-f`)

## Exhaustive list of updates and reasoning

1. `src/hatchet/bin/HATCHet.py`:
   - New `--ampdel-exempt` argument (default from config). Parsed into a de-duplicated label list that keeps order.
     `str()` is applied first because the ini reader turns a single numeric value into an int.
   - Errors for an empty list and for combining with `--noampdel` (the user's choice: error, not warning).
   - `main()` checks the labels against the SEG keys before any solving.
   - `makeBaseCMD` forwards `--ampdel-exempt` to the native solver, but only inside the `-f` branch.
   - `execute_python` forwards the labels to `solve()`.
2. `src/hatchet/hatchet.ini` and `script/hatchet.ini`: new `ampdelexempt =` key, so that ini and environment
   (`HATCHET_COMPUTE_CN_AMPDELEXEMPT`) configuration work like other options.
3. `src/hatchet/utils/solve/ilp_subset.py`:
   - `ampdel_exempt` is stored as a frozenset on `ILPSubset` and `ILPSubsetSplit` and carried through both
     `__copy__` methods, because coordinate descent builds its CARCH/UARCH models by copying.
   - New `ampdel_applies(cluster_id)` helper, so that variable creation, constraints and `first_hot_start` share a
     single rule.
   - The constraint loop iterates over the created `adA` keys, so the variable and constraint sets cannot drift
     apart.
   - Removed the now-unused local `ampdel` variables.
4. `src/hatchet/utils/solve/cd.py`: forwards `ampdel_exempt` into the ILP that coordinate descent wraps.
5. `src/hatchet/utils/solve/__init__.py`:
   - `solve(..., ampdel_exempt=None)` resolves labels to the SEG `#ID` values by string match. The `#ID` column is
     numeric, so `'3'` must become `3`.
   - Errors for unknown labels, and for exemptions with `ampdel=False` (consistent with the CLI and native solver).
   - Passes the resolved set to both validators and all 8 model constructions.
6. `src/hatchet/utils/fixed_matrix.py`, `src/hatchet/utils/sample_clonal.py`: the validators skip exempt clusters in
   their amp+del conflict check, so user-fixed amp+del entries are accepted only where the model allows them.
7. `src/argparse.h/.cpp`:
   - Native `--ampdel-exempt` option stored as a `StringSet`, with help, usage and `toString` lines.
   - Errors for an empty list, a repeated option, and use without `-f`. Pattern copied from `--max-homdel-size`.
8. `src/solve.cpp`: `getAmpdelExempt()`, modeled on `getCNMap()`, maps labels to cluster indices and throws for
   unknown IDs. The resulting `IntSet` is passed to the validators, `CoordinateDescent`, both `IlpSubset`
   constructions and `firstHotStart`.
9. `src/ilp-min.h/.cpp`:
   - New `const IntSet& _ampdelExempt` member, stored by reference like `_cn`.
   - New `ampdelApplies(s)` helper, used for `adA/adB` variable creation, constraints and hot-start values.
   - Exempt clusters are skipped in the `fixC` assertions and in the static `firstHotStart`.
   - Without these skips the hot start would set the start attribute on a non-existent `GRBVar` for an exempt
     cluster.
10. `src/worker.h/.cpp`, `src/coordinate_descent.h/.cpp`: pass the exempt set through to the per-seed `IlpSubset`s.
11. `src/fixed_matrix.h/.cpp`, `src/sample_clonal.h/.cpp`: the native validators skip exempt clusters, mirroring
    item 6.
12. `tests/fixed_matrix_cpp_test.cpp`, `tests/sample_clonal_cpp_test.cpp`:
    - Updated for the new validator signatures.
    - Added `ampdel-exempt` scenarios. The sample-clonal driver also gets an `ampdel` scenario with a second sample,
      because an amp+del conflict needs two states of the same cluster.
13. `tests/test_ampdel_exempt.py` (new): tests for
    - model structure (−2 variables and −4·(n−1) constraints per exempt cluster, for both model classes)
    - clonal-plus-exempt being a no-op
    - hot-start clamping
    - both validators
    - `solve()` label resolution and errors
    - CLI and config parsing
    - SEG membership check in `main()`
    - `makeBaseCMD` forwarding
    - a Pyomo/Gurobi infeasible→feasible solve
    - native rejection paths
    - native end-to-end solves in modes 0, 1 and 2
    - the docs and ini keys

    These are in a new file, not in `test_fixed_matrix.py` as the plan said, to keep the feature's tests together.
    The C++-driver tests went into the existing modules so they reuse those modules' compile fixtures.
14. `tests/test_fixed_matrix.py`, `tests/test_sample_clonal.py`: one native validator test each for the new
    scenario.
15. `docs/source/doc_compute_cn.md`: a table row and a paragraph covering semantics, the ini key, both solvers and
    the resource or model-selection impact.
16. `benchmarks/2026-10-22-ampdel-exempt-clusters/ampdel_usage.py`: reproducible script behind the benchmark's amp+del
    usage columns.

## Tasks prevented by sandbox-limited file access

1. During the repo-wide `grep -rn ampdel .` search, `grep` hit "Permission denied" on `./.claude/loop.md`. The file
   is Claude Code configuration, not hatchet code, so no workaround was needed.
2. `module load gurobi` had no effect because Lmod's module path `/admin/software/lmod/modulefiles` is outside the
   sandbox's readable paths. Only `/admin/software/gurobi` is allowed. As a result `GRB_LICENSE_FILE` and
   `GUROBI_HOME` stay unset, and every native solve aborts with `GRBException`.
   - The native build still worked, because `build-native/CMakeCache.txt` already holds the Gurobi library paths
     under `/admin/software/gurobi`.
   - Workaround: native code paths that run before Gurobi starts (argument parsing, label resolution, validators)
     were tested directly. The pre-existing native test failures were shown to be environmental by comparing
     against pristine HEAD. Native end-to-end solves and the benchmark are left to the user (see Follow-up).
