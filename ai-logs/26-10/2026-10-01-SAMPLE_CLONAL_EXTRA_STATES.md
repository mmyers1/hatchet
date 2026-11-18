# Sample-Clonal: More Than Two States per Sample

## Status

- Overall: Implemented; licensed solve verification pending (user to run)
- Last updated: 2026-10-01

## Goal

Allow `--sample-clonal` / `--sample-clonal-file` (`sampleclonal` / `sampleclonalfile`) to list more than two
sample-clonal states for each sample, adding as little code as possible.

Agreed semantics:

- The first two states listed for each sample (TSV row order or inline order) are the calibration pair. They
  must use distinct clusters and distinct total copy numbers, and they alone enter the HATCHet Theorem 1 system
  that derives the sample's scaling factor and purity.
- Any further states for that sample use distinct clusters but have no total-copy-number restriction. They are
  **not** used to derive scaling or purity. They are fixed on every tumor clone present in that sample, exactly like
  the first two.
- Additional states are validated by the same checks as the calibration pair (copy-number bound, amp/del,
  `--cnstates`, `--max-homdel-size`, `--fixed-matrix` compatibility) and by the same post-solve verification.

## Method

Discovery showed that every sample-clonal enforcement path already iterates over all entries of a sample in both
backends:

| Concern | Python | C++ |
|---------|--------|-----|
| CARCH (fixed U) | `ILPSubset.compute_fixed_bounds` fixes A/B with the same `fix()` bounds used by `--fixed-matrix` | `IlpSubset::computeFixedBounds` uses `fixCN` |
| UARCH (fixed A/B) | incompatible clones get `U = 0` | same, via `fixProportion` |
| FULL | big-M constraints in `apply_sample_clonal_constraints` | `IlpSubset::fixSampleClonal` |
| CD seeds | `build_sample_clonal_u` | `CoordinateDescent::buildSampleClonalU` |
| Validation | `validate_sample_clonal_constraints`, `validate_fixed_matrix(sample_clonal=...)` | `SampleClonal::validate`, fixed-matrix validation |
| Verification | `verify_sample_clonal_solution` | `SampleClonal::verify` |

Entry order is preserved from parsing through resolution, the Python-to-native handoff file
(`write_sample_clonal_file`) and `entries_for_sample` / `entriesForSample`, so "first two" is well defined in both
backends.

Only five sites hard-coded two entries: the Python TSV/inline validation, the Python inline parser, the Python
calibration, the C++ reader, and the C++ calibration. The change relaxes these sites to "at least two" and makes the
calibration use the first two. All enforcement code is reused unchanged.

## Progress

- [x] Discovery and plan approved by the user.
- [x] Python parser/validation/calibration changes.
- [x] C++ reader/calibration changes.
- [x] CLI help text.
- [x] Documentation.
- [x] Tests (parser messages, calibration, verification, C++ harness, end-to-end solve case, docs contract).
- [x] Cheap verification run.
- [ ] Licensed solve verification (to be run by the user; commands below).

## Results

- 2026-10-01: `conda run -n hatchet-test env PYTHONPATH="$PWD/src" pytest -q tests/test_sample_clonal.py
  -k "parser or calibration or verify or resolution or cpp_sample_clonal"`: 29 passed, 69 deselected (2.98 s).
  This covers the new calibration and verification tests and both C++ harness parametrizations (the harness was
  compiled and executed, not skipped).
- 2026-10-01: Documentation contract tests (`-k documentation` over `test_sample_clonal.py` and
  `test_fixed_matrix.py`): 3 passed.
- 2026-10-01: `compute-cn --help` renders the updated `--sample-clonal` / `--sample-clonal-file` help text.
- 2026-10-01: The new end-to-end cases are collected as
  `test_compute_cn_separates_samples_with_conflicting_sample_clonal_states[{gurobi,cpp}-{two,three}-states]`.
  They have not been run here because they perform licensed Gurobi solves and build the native solver.
- `git diff --check` is clean for the touched files. No touched line exceeds 120 characters; the only long line in
  `src/sample_clonal.cpp` (L349) predates this change.
- 2026-10-01 (user run): `pytest -q tests/test_sample_clonal.py tests/test_fixed_matrix.py tests/test_solver.py`
  gave 213 passed, 2 failed in 87.8 s.
  - `test_compute_cn_separates_samples_with_conflicting_sample_clonal_states[cpp-three-states]`: **stale native
    binary, not a code defect.** The `sample_clonal_cpp_solver` fixture uses the packaged `src/hatchet/solve` whenever
    its `-h` output advertises `-S SAMPLE-CLONAL` and `--max-homdel-size`, and builds from source only otherwise.
    The packaged binary dates from 2026-09-29 14:59, before `src/sample_clonal.cpp` was edited (2026-10-01 17:19), and
    `strings src/hatchet/solve` still contains `' must have exactly two sample-clonal entries`. The old binary
    therefore rejected the 3-state file; `cpp-two-states` passed on the same binary.
  - `test_solve_command_cbc`: **pre-existing and environmental.** `cbc -quit` segfaults by itself in `hatchet-test`.
    `2026-09-29-compute-cn-efficiency.md` records the same failure on unchanged HEAD. The test doesn't use
    sample-clonal input, and none of the files touched here are on its path.
  - The user chose not to change the fixture; the native solver must instead be rebuilt and installed before the
    licensed tests run.

### Commands for the user

Rebuild and install the native solver **first**, because the test fixture prefers the packaged binary:

```bash
cd /data1/shahs3/users/myersm2/repos/hatchet
module load gurobi
GUROBI_ROOT="$GUROBI_HOME"
case "$(basename "$GUROBI_ROOT")" in linux64|mac64) GUROBI_ROOT="$(dirname "$GUROBI_ROOT")";; esac
cmake -S . -B build-native -DCMAKE_BUILD_TYPE=Release -DGUROBI_HOME="$GUROBI_ROOT"
cmake --build build-native --config Release
install -m 755 build-native/solve src/hatchet/solve
strings src/hatchet/solve | grep -c 'at least two sample-clonal entries'   # expect 1
```

These flags match the fixture's own fallback build. `build-native/` is an untracked scratch directory. (An earlier
version of this log gave `cmake --build build --target solve`, which is wrong: `build/` is a setuptools directory with
no CMake cache.)

Then run:

```bash
# expect 4 passed
conda run -n hatchet-test env PYTHONPATH="$PWD/src" pytest -q tests/test_sample_clonal.py \
  -k conflicting_sample_clonal_states

# expect 0 failures (214 passed); the CBC test is excluded because CBC segfaults in this environment
conda run -n hatchet-test env PYTHONPATH="$PWD/src" pytest -q \
  tests/test_sample_clonal.py tests/test_fixed_matrix.py tests/test_solver.py \
  --deselect tests/test_solver.py::test_solve_command_cbc
```

The `three-states` cases must report `2|0` for cluster 3 in the clone present in `TumorA`, with `u_normal = 0.5` in
both samples. If `cpp-three-states` still fails after the `strings` check prints `1`, rerun it with
`-x -vv --basetemp=/tmp/hatchet-sc` and inspect `best.bbc.ucn` and the solver stderr under that directory.

## Summary

`--sample-clonal`, `--sample-clonal-file`, `sampleclonal` and `sampleclonalfile` now accept two or more states per
sample. The first two states per sample, in input order, still form the calibration pair (distinct clusters and
distinct totals) and alone determine the scaling factor and purity. Additional states are fixed on every tumor clone
present in that sample through the existing sample-clonal enforcement in both backends, which already iterated over
all entries per sample. In coordinate-descent CARCH steps this enforcement uses the same `fix()` / `fixCN` bound
mechanism as `--fixed-matrix`. No solver, seed, wrapper-flow, INI or output-format code changed.

## Updates and Reasoning

1. `src/hatchet/utils/sample_clonal.py::_validate_entries`: changed `len(sample_entries) != 2` to `< 2`, with the
   message `must have at least two entries`. Reason: allow additional states while still requiring a calibration
   pair.
2. Same function: the distinct-total check now applies to `sample_entries[:2]` and reports the line of the second
   entry. Reason: only the calibration pair needs distinct totals (the Theorem 1 system is singular otherwise).
   Additional states may share totals with the pair. The duplicate (sample, cluster) check is unchanged, so all
   clusters within a sample remain distinct.
3. `src/hatchet/utils/sample_clonal.py::parse_sample_clonal_inline`: changed `len(states) != 2` to `< 2`, with the
   same message. Reason: same as (1), for the inline grammar.
4. `src/hatchet/utils/sample_clonal.py::derive_sample_clonal_calibration`: unpacks
   `entries_for_sample(...)[:2]` and adds a comment. Reason: the explicit requirement that states beyond the first two
   are not used to calibrate.
5. `src/sample_clonal.cpp::SampleClonal::read`: changed `entries.size() != 2` to `< 2` and the message to
   `at least two`, and added a comment above the unchanged `entries[0]` / `entries[1]` distinct-total check. Reason:
   native parity with (1) and (2).
6. `src/sample_clonal.cpp::SampleClonal::calibrate`: changed `entries.size() != 2` to `< 2`; the function already
   used only `entries[0]` and `entries[1]`. Reason: native parity with (4).
7. `src/hatchet/bin/HATCHet.py`: `--sample-clonal` and `--sample-clonal-file` help text now says "two or more" and
   that only the first two calibrate. Reason: the user asked for explicit acknowledgement of these semantics.
8. `docs/source/doc_compute_cn.md` (sample-clonal section): changed "two" to "two or more", added a paragraph
   stating that only the first two states are references and that additional states are fixed but not used for
   scaling or purity, extended the inline grammar, added a third state to the inline/INI/TSV examples, changed the TSV
   rule from "exactly twice" to "at least twice" with first-two-in-file-order semantics, and added constraint
   bullets. Reason: user-facing contract and explicit acknowledgement.
9. `tests/test_sample_clonal.py`: changed three expected messages from `exactly two` to `at least two`. Reason: the
   messages changed in (1), (3) and (5).
10. `tests/test_sample_clonal.py::test_calibration_ignores_states_after_first_two` (new): a 3-state spec, inline and
    TSV, whose third state repeats the first total and whose cluster has an extreme RDR still parses, preserves
    order, and yields gamma 3 and purity 0.5. Reason: covers (2), (4) and order preservation.
11. `tests/test_sample_clonal.py::test_verify_sample_clonal_solution_checks_additional_states` (new): post-solve
    verification rejects a present clone that violates only the third state. Reason: solver-free proof that
    additional states are enforced by the existing generic verification.
12. `tests/test_sample_clonal.py::test_cpp_sample_clonal_parser_and_calibration_match_python`: parametrized with a
    3-row file. `tests/sample_clonal_cpp_test.cpp` gained cluster `c3`, a 3-row RDR matrix (c3 RDR 5.0, which must be
    ignored) and an `entries().size() < 2` sanity check. Reason: covers (5) and (6) natively.
13. `tests/test_sample_clonal.py::_run_conflicting_sample_clonal_solve` gained an `extra_tumor_a_state` argument, and
    `test_compute_cn_separates_samples_with_conflicting_sample_clonal_states` is parametrized with `,3:2:0` for
    TumorA, asserting `2|0` for cluster 3 in TumorA's clone. Reason: the end-to-end proof through the public CLI on
    both backends. The observed BAF is 0.5, so `2|0` appears only when the state is enforced, and amp/del remains
    feasible (A: 2/1, B: 0/1 relative to base 1).
14. `tests/test_sample_clonal.py::test_compute_cn_documentation_covers_sample_clonal_contract`: added an assertion for
    `not used to derive scaling or purity`. Reason: keeps the documented acknowledgement from regressing.
15. This log: recorded the user's licensed run and the diagnosis of both failures, and replaced the incorrect
    rebuild command with a configure/build/install sequence that must run before the licensed tests. Reason: the
    `cpp-three-states` failure came from the fixture selecting the stale packaged binary, and the earlier command
    targeted a non-CMake directory. No source or test files changed, per the user's choice.
