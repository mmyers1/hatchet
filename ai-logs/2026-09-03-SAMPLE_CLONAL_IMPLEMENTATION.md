# Sample-Clonal Override Implementation

## Status

- Overall: Complete
- Current phase: Complete
- Current task: None
- Last updated: 2026-09-03

## Agreed Behavior

- Support an inline command-line/INI specification and a versioned TSV file.
- Require exactly two sample-clonal reference states for every input sample.
- Require the user to select either `--diploid` or `--tetraploid`.
- Reject use with the existing global `--clonal` option.
- Allow `--purities` only when the values match the derived purities within tolerance.
- Use the sample-clonal states to derive scaling factors and purities.
- Enforce each supplied A/B state on every tumor clone present in that sample.
- Allow one cluster to have different states in different samples.
- Do not reinterpret sample-clonal states as global A/B constraints.

## Completed Discovery

- [x] Read the local HATCHet manuscript PDF and supplement.
- [x] Read the local HATCHet2 manuscript PDF.
- [x] Verify the per-sample scaling and purity equations from HATCHet Supplementary Theorem 1.
- [x] Trace automatic and manual global-clonal handling in the Python wrapper.
- [x] Trace scaling, purity, factorization base, and global A/B constraints in the Python solver.
- [x] Trace scaling, coordinate descent, ILP modes, and global A/B constraints in the C++ solver.
- [x] Identify current Python/C++ purity-constraint and sample-ordering differences.
- [x] Decide the public option precedence and ploidy requirements.
- [x] Decide that downstream sample clonality must be enforced conditionally on clone presence.

## Phase 1: Parsing and Calibration

- [x] Add immutable sample-clonal entry, resolved-entry, and calibration types.
- [x] Add versioned TSV parser.
- [x] Add inline parser.
- [x] Require exactly two distinct clusters and distinct total copy numbers per sample.
- [x] Resolve sample and cluster identifiers exactly against the SEG input.
- [x] Require complete coverage of all SEG samples.
- [x] Implement order-independent Theorem 1 scaling and purity derivation.
- [x] Validate singular systems, finite values, positive scaling, and purity bounds.
- [x] Add focused parser, resolution, and algebra tests.

## Phase 2: Wrapper and Configuration

- [x] Add `--sample-clonal` inline option.
- [x] Add `--sample-clonal-file` TSV option.
- [x] Add `sampleclonal` and `sampleclonalfile` INI settings.
- [x] Reject simultaneous inline and file inputs.
- [x] Reject sample-clonal input with `--clonal`.
- [x] Require exactly one explicit ploidy flag.
- [x] Compare optional `--purities` with named derived purities.
- [x] Bypass neutral/global-clonal discovery and parsing.
- [x] Preserve all legacy branches when no sample-clonal input is supplied.

## Phase 3: Python Solver

- [x] Pass an explicit factorization base independently of global clonal copy numbers.
- [x] Keep the global copy-number map empty for sample-clonal runs.
- [x] Apply derived gamma values to RDR by sample name.
- [x] Fix normal proportions to `1 - purity` in FULL and UARCH modes.
- [x] Add presence variables when sample-clonal enforcement is active, including `minprop = 0`.
- [x] Add conditional A/B constraints in FULL mode.
- [x] Add direct A/B constraints for present clones in CARCH mode.
- [x] Make incompatible clones unavailable in UARCH mode.
- [x] Apply equivalent changes to the binwise models.
- [x] Verify sample-clonal conditions after solving.

## Phase 4: Coordinate Descent

- [x] Generate purity-aware U seeds.
- [x] Generate seeds that respect sample-state compatibility.
- [x] Combine sample-clonal seed requirements with fixed U entries.
- [x] Skip infeasible seeds without aborting remaining work.
- [x] Report a clear error when no feasible seed exists.

## Phase 5: Fixed-Matrix Interactions

- [x] Reject fixed normal proportions that conflict with derived purity.
- [x] Reject fixed-positive tumor clones with incompatible reference A/B states.
- [x] Treat incompatible fixed A/B states as making a clone unavailable when U is not fixed positive.
- [x] Validate copy-number bounds and large-cluster nonzero constraints.
- [x] Validate amp/del and `cnstates` lower-bound interactions.
- [x] Add focused Python fixed-matrix interaction tests.

## Phase 6: C++ Backend

- [x] Add native versioned-TSV parser and named ID resolver.
- [x] Add native Theorem 1 calibration with Python parity.
- [x] Add internal sample-clonal file and explicit-base arguments.
- [x] Bypass automatic/global calibration in native sample-clonal mode.
- [x] Add named derived purity constraints.
- [x] Add FULL, CARCH, and UARCH sample-clonal enforcement.
- [x] Add purity- and compatibility-aware native seeds.
- [x] Handle infeasible native coordinate-descent seeds.
- [x] Extend native fixed-matrix validation.
- [x] Add final native solution verification.
- [x] Add the new native source to CMake.

## Phase 7: Documentation and Compatibility

- [x] Document both input formats and examples.
- [x] Document the expert-use warning and mathematical conditions.
- [x] Document option precedence and ploidy requirements.
- [x] Document interactions with purity, fixed matrices, copy bounds, amp/del, and clone ranges.
- [x] Add legacy automatic/global-clonal regression tests.
- [x] Confirm output formats and filenames remain unchanged.

## Public Parameter Solve Coverage

- [x] Solve through `hatchet compute-cn --sample-clonal` with the Python backend.
- [x] Solve through `hatchet compute-cn --sample-clonal-file` with the Python backend.
- [x] Solve through `hatchet compute-cn --sample-clonal` with the C++ backend.
- [x] Solve through `hatchet compute-cn --sample-clonal-file` with the C++ backend.
- [x] Exercise CD-plus-ILP, ILP-only, and CD-only for both public input forms and both backends.
- [x] Assert `best.bbc.ucn` and `best.seg.ucn` contain the enforced states and derived 50% purity.
- [x] Build the native test solver from the Gurobi module when an updated binary is not supplied.

## Documentation Coverage

- [x] Test sample-clonal ploidy requirements and both-ploidy rejection.
- [x] Test inline/file mutual exclusion with an existing file.
- [x] Test global-clonal incompatibility.
- [x] Test purity agreement at, inside, and outside the inclusive `0.05` tolerance.
- [x] Test that derived purity remains authoritative in public output.
- [x] Test present/absent clone enforcement in FULL, CARCH, and UARCH models.
- [x] Test a real multi-sample solve with conflicting states for one cluster on Python and C++.
- [x] Test clone-count ranges and clear insufficient-clone reporting.
- [x] Test successful sample-clonal plus fixed-matrix solving on Python and C++.
- [x] Test fixed-normal purity, fixed-positive state, and remaining-purity conflicts.
- [x] Test active diploid/tetraploid bounds, amp/del, `cnstates`, and large-cluster nonzero rules.
- [x] Test fixed-matrix UTF-8 BOM, invalid UTF-8, blank lines, row order, sparse semantics, and duplicates.
- [x] Test fixed-matrix CLI/INI/environment precedence and clone-range rejection.
- [x] Test fixed-matrix `minprop` boundaries and native validation parity.
- [x] Test a positive-residual constrained solve and partial global infeasibility.
- [x] Test one fixed matrix is routed to both automatic model-selection branches and branch-specific failures are explicit.
- [x] Retain static offline-editor packaging checks; interactive editor behavior tests were explicitly excluded by user request.

## Validation Log

- 2026-09-03: `conda run -n hatchet-test env PYTHONPATH="$PWD/src" pytest -q tests/test_sample_clonal.py` - 20 passed.
- 2026-09-03: Same focused command after CLI/config integration - 26 passed.
- 2026-09-03: Same focused command after Python solver calibration routing - 28 passed.
- 2026-09-03: `conda run -n hatchet-test env PYTHONPATH="$PWD/src" pytest -q tests/test_fixed_matrix.py tests/test_sample_clonal.py` - 84 passed after explicit-base constructor propagation.
- 2026-09-03: Same combined command after FULL/CARCH/UARCH and binwise enforcement - 94 passed.
- 2026-09-03: Focused sample-clonal suite after conflict-aware seed generation - 38 passed.
- 2026-09-03: Focused sample-clonal suite after early/fixed-matrix/post-solve validation - 45 passed.
- 2026-09-03: Focused Python plus standalone C++ parser suite - 50 passed, including normalized native handoff cleanup.
- 2026-09-03: Combined sample-clonal and fixed-matrix suites - 107 passed.
- 2026-09-03: Real Gurobi optimization on an exact enforced fixture - ILP, CD, and CD-plus-ILP all passed with objective 0.
- 2026-09-03: `module load gurobi` provides Gurobi 9.5.2; the full native optimizer builds successfully and solves the exact fixture in ILP, CD, and CD-plus-ILP modes.
- 2026-09-03: Six public Python `compute-cn` solves passed: inline/file sample-clonal input across all three solve modes.
- 2026-09-03: Six public C++ `compute-cn` solves passed: inline/file sample-clonal input across all three solve modes.
- 2026-09-03: Complete `tests/test_sample_clonal.py` run with the module-built native solver - 66 passed in 19.87s.
- 2026-09-03: Native integration fixture independently configured and built the solver from `module load gurobi`, then completed a public `--sample-clonal` ILP solve - 1 passed in 15.04s.
- 2026-09-03: Expanded sample-clonal and fixed-matrix feature suites - 162 passed.
- 2026-09-03: Legacy fixed-matrix solver scenarios against the selected current native binary - 16 passed.
- 2026-09-03: Final combined `test_sample_clonal.py`, `test_fixed_matrix.py`, and `test_solver.py` run - 183 passed in 88.82s.
- The initial base-environment run could not collect because `pandas` was unavailable; all implementation tests pass in the required `hatchet-test` environment.
