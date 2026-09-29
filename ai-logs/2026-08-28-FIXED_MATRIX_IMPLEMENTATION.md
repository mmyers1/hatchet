# Fixed Matrix Entries Implementation

## Goal

Add an optional `compute_cn.fixedmatrix` / `--fixed-matrix` interface that lets users fix arbitrary entries of the HATCHet factorization matrices:

- `A[cluster, clone]`: integer A-allele copy number.
- `B[cluster, clone]`: integer B-allele copy number.
- `U[clone, sample]`: clone proportion in a bulk sample.

Unlisted cells remain solver variables. Fixed cells are hard equality constraints in every solving mode and in both the C++ and Pyomo backends.

The work also repairs the existing symmetry-breaking defect in both backends. The intended tumor-clone ordering is

$$
\sum_s (s + 1)(A_{s,i} + B_{s,i})
\leq
\sum_s (s + 1)(A_{s,i+1} + B_{s,i+1}).
$$

## Compatibility Requirements

- Runs without a fixed-matrix file retain the existing interface and output schema.
- The default compiled C++ solver and all supported Pyomo solvers honor the same file.
- ILP-only, coordinate-descent-only, and coordinate-descent-plus-ILP modes enforce the same cells.
- Standard and binwise Python models use the same semantics.
- The symmetry repair can reorder equivalent `clone1`, `clone2`, ... columns for solutions with at least two tumor clones. Objective values and represented mixtures should remain equivalent.
- Existing `clonal`, `purities`, `minprop`, `cnstates`, copy-number bounds, and amp/del constraints remain authoritative; contradictory fixed entries fail before solving when the conflict can be determined locally.

## Agreed Interface

### Configuration

Add an empty optional setting to every distributed `[compute_cn]` section:

```ini
fixedmatrix =
```

Expose the same setting on the command line:

```text
hatchet compute-cn --fixed-matrix PATH ...
```

The command-line value overrides the INI value. The existing `HATCHET_COMPUTE_CN_FIXEDMATRIX` environment override continues to work through the config wrapper.

A fixed-matrix run requires one exact clone count, such as `clones = 3`. Clone ranges are rejected because clone-specific cells need stable identities. The count includes the normal clone.

If neither ploidy is forced, the same constraints are applied to both diploid and tetraploid branches. A branch-specific conflict is reported rather than silently dropping constraints or skipping the branch.

### File Format

Use a versioned, sparse, long-form TSV:

```tsv
# HATCHET_FIXED_ENTRIES 1
# clones 3
matrix	row	column	value
A	1	clone1	2
B	1	clone1	1
U	normal	Tumor1	0.35
U	clone1	Tumor1	0.40
```

Rules:

- `A` and `B` rows use an exact SEG cluster ID in `row` and `normal|clone1|...` in `column`.
- `U` rows use the clone label in `row` and an exact SEG sample ID in `column`.
- Only listed cells are fixed.
- Comments and arbitrary data-row order are allowed.
- Duplicate coordinates are errors, including duplicates with the same value.
- `A` and `B` values are nonnegative integers.
- `U` values are finite numbers in `[0, 1]`.
- `normal` is clone 0; fixed normal A/B values must be 1.

## Implementation Checklist

### 1. Symmetry Repair

- [x] Correct `IlpSubset::buildSymmetryBreaking()` in `src/ilp-min.cpp` to use clone `i + 1` on the right-hand side.
- [x] Correct `ILPSubset.build_symmetry_breaking()` in `src/hatchet/utils/solve/ilp_subset.py` identically.
- [x] Add a solver-independent Pyomo model test that inspects/evaluates the adjacent-clone constraint.
- [x] Keep the current coefficient and loop bounds unchanged.

### 2. Python File Contract

- [x] Add `src/hatchet/utils/fixed_matrix.py` with an immutable sparse representation.
- [x] Parse the version marker, clone-count metadata, header, and rows with source line numbers.
- [x] Reject malformed rows, unsupported matrices, duplicate coordinates, non-finite values, and invalid clone labels.
- [x] Resolve cluster and sample names only after the SEG context is available.
- [x] Validate normal clone, bounds, fixed A+B totals, U sums/residual capacity, `minprop`, purities, clonal states, amp/del direction, and locally decidable `cnstates` conflicts.
- [x] Do not require exact reconstruction of fractional copy numbers because the HATCHet objective permits residual error.

### 3. Python Solver Integration

- [x] Add `fixedmatrix` to packaged and user-facing INI files and `--fixed-matrix` to `compute-cn`.
- [x] Require an exact clone count and parse the file after SEG/sample discovery.
- [x] Pass resolved entries through `hatchet.utils.solve.solve`.
- [x] Preserve entries in `ILPSubset`, `ILPSubsetSplit`, their copy methods, coordinate-descent workers, and multiprocessing payloads.
- [x] Add hard equalities in FULL mode.
- [x] In CARCH, overlay fixed U cells on the current full U matrix and constrain sparse A/B cells.
- [x] In UARCH, overlay fixed A/B cells on the current full copy-number matrices and constrain sparse U cells.
- [x] Generate U seeds around fixed cells while maintaining sum-to-one and `minprop` constraints.
- [x] Apply fixed A/B cells to first starts without moving user-addressed columns during rank normalization.
- [x] Verify requested cells in returned solutions before writing output files.

### 4. C++ Solver Integration

- [x] Add `src/fixed_matrix.h` and `src/fixed_matrix.cpp` with the same TSV contract and defensive validation.
- [x] Add the new source to `CMakeLists.txt`.
- [x] Add an internal/direct `-F PATH` option to `ArgParse` and forward it from the Python wrapper.
- [x] Resolve textual labels through `InputInstance` after reading the SEG file.
- [x] Pass entries through `solve.cpp`, `CoordinateDescent`, `Worker`, and `IlpSubset`.
- [x] Enforce sparse equalities and constant overlays in FULL, CARCH, and UARCH modes.
- [x] Condition C++ random U seeds and A/B starts on fixed cells.
- [x] Verify fixed cells after solving.

### 5. Offline Editor

- [x] Add a `hatchet fixed-matrix-editor` command with a small `importlib.resources` / `webbrowser` launcher.
- [x] Include a no-open option for remote or headless systems.
- [x] Store a self-contained HTML asset under `src/hatchet/data/`, which is already included by the package-data wildcard.
- [x] Support local SEG import, exact clone count, A/B/U tabs, keyboard editing, paste, filtering, import of an existing fixed TSV, live validation, clear actions, and canonical TSV download.
- [x] Keep normal A/B cells visibly locked to 1 and all other blank cells free.
- [x] Use no server, external JavaScript, network request, or new runtime dependency.

### 6. Documentation and Tests

- [x] Document the format, matrix dimensions, workflow, validation, and option interactions in `docs/source/doc_compute_cn.md`.
- [x] Update all distributed INI examples.
- [x] Add parser/config/editor tests in `tests/test_fixed_matrix.py`.
- [x] Add a small constraints fixture under `tests/data/fw/`.
- [x] Add model-level FULL/CARCH/UARCH, copy-persistence, seed, hot-start, binwise, and symmetry tests that do not require a solver license.
- [x] Add availability-gated C++ and Pyomo end-to-end tests for all solve modes.
- [x] Apply fixed-matrix fixture and parser coverage to both fixed-width and variable-width binning scenarios.
- [x] Retain no-fixed-file regression comparisons.

## End-to-End Solve Test Plan

### Goals

The licensed solve tests verify the behavior that model inspection and parser tests cannot prove:

- each backend honors hard A, B, and U entries through actual optimization;
- coordinate-descent subproblems preserve constraints between C and U steps;
- combined mode preserves constraints when the CD result becomes an ILP hot start;
- fixed-width and variable-width binning inputs resolve cluster/sample identifiers correctly;
- output BBC files retain the existing schema and contain valid mixtures; and
- solutions with at least two tumor clones follow the repaired canonical clone order.

### Test Matrix

The primary matrix contains 12 cases in `test_fixed_matrix_solve_binning_scenario`:

| Dimension | Values |
|-----------|--------|
| Backend | compiled `cpp`, Pyomo `gurobi` |
| Binning | fixed-width `tests/data/fw`, variable-width `tests/data/vw/bbc` |
| Solve mode | `0` CD plus ILP, `1` ILP only, `2` CD only |
| Clone count | 2 total clones (`normal`, `clone1`) |

These cases use `tests/data/{fw,vw}/fixed-matrix.tsv`. Each fixture contains at least one A entry, one B entry, and fixed U entries for a real sample from the matching SEG file.

A second four-case matrix in `test_fixed_matrix_solve_preserves_canonical_clone_order` covers:

| Dimension | Values |
|-----------|--------|
| Backend | compiled `cpp`, Pyomo `gurobi` |
| Binning | fixed-width, variable-width |
| Solve mode | `1` ILP only |
| Clone count | 3 total clones (`normal`, `clone1`, `clone2`) |

These cases use `fixed-matrix-n3.tsv` fixtures and assert nondecreasing weighted copy-number burden across the two tumor clones. Solver-independent model tests continue to cover the symmetry constraint in CARCH and FULL model construction, so the licensed symmetry matrix is intentionally limited to one solve mode.

### Runtime Profile

The integration helper uses deterministic, bounded test settings:

- one coordinate-descent seed;
- one worker, with one thread passed to the compiled solver;
- random seed 6700;
- at most two coordinate-descent iterations;
- a 60-second limit for each ILP invocation;
- forced diploid mode to avoid unrelated WGD model selection; and
- exact clone counts declared by each fixed-matrix fixture.

The production objective and normal `compute-cn` output path remain in use. Tests do not mock solver/model behavior.

### Assertions

Every successful primary-matrix solve checks:

1. `best.bbc.ucn` is produced.
2. Original BBC columns are unchanged and followed by the expected interleaved CN/U columns.
3. Every `cn_normal` value is `1|1`.
4. Clone proportions sum to 1 in every sample within tolerance 0.001.
5. Every A/B/U coordinate listed in the scenario fixture has the requested output value.

The three-clone matrix additionally reconstructs each tumor clone's burden in SEG cluster order and verifies that burdens are nondecreasing.

### Availability and Failure Policy

`solver_available()` can report an installed executable even when its Gurobi license is expired. A session-scoped fixture therefore performs a real one-variable Gurobi optimization before any integration case runs. If that probe fails, all licensed cases skip with one explicit reason.

After the runtime probe:

- C++ cases also require the packaged `solve` executable and the existing C++ availability check.
- Pyomo cases also require `solver_available('gurobi')`.
- Once a case starts a real HATCHet solve, any parser, solver, output, or constraint failure fails the test; it is not converted to a skip.

### Compatibility Baselines

The pre-existing no-fixed-file tests remain unchanged. They compare `best.bbc.ucn` against the fixed-width reference for C++, gurobipy, external Gurobi, and CBC when available. This preserves exact-output regression coverage independently of the new constraint tests.

### Execution

Collect the drafted matrix without running a solver:

```bash
PYTHONPATH=src $PYTHON -m pytest \
	tests/test_solver.py::test_fixed_matrix_solve_binning_scenario \
	tests/test_solver.py::test_fixed_matrix_solve_preserves_canonical_clone_order \
	--collect-only -q
```

Run all 16 licensed cases:

```bash
PYTHONPATH=src $PYTHON -m pytest \
	tests/test_solver.py::test_fixed_matrix_solve_binning_scenario \
	tests/test_solver.py::test_fixed_matrix_solve_preserves_canonical_clone_order \
	-xvs
```

Useful slices:

```bash
# C++ only
PYTHONPATH=src $PYTHON -m pytest tests/test_solver.py -k 'fixed_matrix_solve and cpp' -xvs

# Pyomo Gurobi only
PYTHONPATH=src $PYTHON -m pytest tests/test_solver.py -k 'fixed_matrix_solve and gurobi' -xvs

# Coordinate descent modes only
PYTHONPATH=src $PYTHON -m pytest tests/test_solver.py -k 'fixed_matrix_solve and (cd or cd_plus_ilp)' -xvs
```

## Validation Gates

Use the requested environment for all Python checks:

```bash
PYTHON=/home/myersm2/miniforge3/envs/hatchet/bin/python3.9
```

Focused gates:

```bash
PYTHONPATH=src $PYTHON -m pytest tests/test_fixed_matrix.py -xvs
PYTHONPATH=src $PYTHON -m pytest tests/test_solver.py -xvs
PYTHONPATH=src $PYTHON -m pytest tests/test_config.py tests/test_steps.py::test_compute_cn -xvs
```

`PYTHONPATH=src` is required until the workspace package is rebuilt or installed into the environment; otherwise this environment can import an older installed HATCHet copy.

Additional gates:

- Rebuild the packaged C++ solver against the configured Gurobi installation.
- Run repository formatting and lint hooks on touched files.
- Build the Sphinx documentation.
- Exercise editor import, edit, paste, clear, download, and re-import on desktop and narrow layouts.
- Feed an editor-generated file through both the INI and CLI interfaces.
- Compare fixed cells, U tolerance, canonical clone burden order, and output schemas across C++ and Pyomo results.

## Out of Scope

- Changing BBC or SEG result formats.
- Supporting JSON or CSV constraint variants in the first implementation.
- Separate sections for each candidate clone count or ploidy branch.
- Silently selecting a different solver backend.
- Reworking unrelated solver formulations or `cnstates` behavior.
- Changing the symmetry coefficient or ordering beyond the confirmed `i` versus `i + 1` defect.

## Progress Log

### 2026-08-28

- [x] Traced `.ini` configuration through `HATCHet.py` into both compiled and Python solver dispatch.
- [x] Confirmed matrix dimensions and identifiers: A/B are cluster-by-clone, U is clone-by-sample, and clone 0 is normal.
- [x] Confirmed `fix_c()` and `fix_u()` currently replace whole coordinate-descent subproblem matrices and do not implement sparse persistent user constraints.
- [x] Traced model rebuilds through FULL, CARCH, and UARCH modes in both backends.
- [x] Confirmed the default backend is `cpp`, so backend parity is required for a usable default interface.
- [x] Confirmed binwise solving uses `ILPSubsetSplit` and must inherit the Python constraint behavior.
- [x] Selected a dependency-free, packaged, offline HTML editor architecture.
- [x] Selected exact clone-count semantics and application of one file to both ploidy branches.
- [x] Confirmed the symmetry-breaking bug is duplicated in C++ and Python: both sides currently reference clone `i`, while C++ `checkRank()` demonstrates the intended `i + 1` comparison.
- [x] Established the versioned sparse TSV contract and validation boundaries.
- [x] Created this repository-local implementation tracker.
- [x] Corrected the adjacent-clone symmetry expressions in the C++ and Python implementations.
- [x] Added solver-independent regression coverage for the Python model coefficients and a C++ source guard.
- [x] Verified the symmetry slice with `PYTHONPATH=src /home/myersm2/miniforge3/envs/hatchet/bin/python3.9 -m pytest tests/test_fixed_matrix.py -xvs` (`2 passed`).
- [x] Identified that the requested environment imports a stale installed package without `PYTHONPATH=src`; recorded workspace routing in the validation commands.
- [x] Added frozen sparse and resolved fixed-matrix representations with exact textual SEG identifier lookup.
- [x] Added the versioned TSV parser with source-line errors for malformed metadata, headers, values, clone labels, and duplicate coordinates.
- [x] Added context-aware validation for CN bounds, clonal assignments, A+B totals, U feasibility, `minprop`, purities, amp/del direction, large `(0,0)` states, `cnstates`, and fully fixed symmetry order.
- [x] Verified the expanded parser/validation/symmetry suite with the requested Python 3.9 environment (`37 passed`).
- [x] Confirmed the new Python files contain no lines over the repository's 120-character limit and have no VS Code diagnostics.
- [ ] The requested environment does not currently provide the `flake8` module; run the repository lint hook when its development dependencies are available.
- [x] Added `fixedmatrix` to all distributed INI files and exposed `--fixed-matrix` in `compute-cn` help.
- [x] Added file-existence and exact clone-count checks while preserving the empty-option path.
- [x] Added post-SEG loading that checks declared clone metadata and exact cluster/sample identifiers.
- [x] Verified parser, validation, symmetry, config, and front-door behavior (`41 passed` focused tests plus `test_config.py`).
- [x] Propagated resolved entries through all regular and binwise Python ILP/CD construction and copy paths.
- [x] Added named hard equalities in FULL/CARCH/UARCH and overlaid sparse cells onto whole fixed subproblem matrices.
- [x] Conditioned coordinate-descent U seeds and A/B starts while preserving user-addressed clone columns.
- [x] Added defensive whole-matrix shape, U-sum, and `minprop` checks.
- [x] Added post-solve verification for direct ILP, CD, combined CD+ILP, and binwise results.
- [x] Verified Python enforcement, copying, binwise behavior, starts, result checks, parsing, and symmetry (`46 passed`).
- [x] Compiled all touched Python modules successfully and confirmed no touched Python line exceeds 120 characters.
- [x] Added a temporary fail-closed dispatch check so the default C++ backend cannot silently ignore fixed entries while native support is pending (`47 passed` focused tests).
- [ ] The existing Gurobi end-to-end regression is blocked because the installed license expired on 2023-10-25, even though Pyomo reports the solver as available.
- [x] Added the native C++ TSV parser, resolved sparse representation, local validation, overlays, seed conditioning, and result verification.
- [x] Added direct native `-F` parsing and quoted Python-to-binary forwarding, then removed the temporary fail-closed dispatch guard.
- [x] Propagated immutable entries through `solve.cpp`, `CoordinateDescent`, `Worker`, and every `IlpSubset` mode.
- [x] Added native hard equalities for FULL/CARCH/UARCH, conditioned A/B starts and U seeds, and disabled rank remapping for user-addressed clone columns.
- [x] Added a solver-independent C++ harness compiled by pytest to exercise parsing, validation, overlays, U completion, verification, and malformed clone labels.
- [x] Added `tests/data/fw/fixed-matrix.tsv` for native and future licensed integration tests.
- [x] Verified the complete focused suite (`49 passed`), all native translation units with C++11 syntax checks, and the linked CMake `solve` target.
- [x] Verified the rebuilt binary advertises `-F` and reports a source-line clone-count mismatch before solver initialization.
- [x] Re-ran Python compilation, `test_config.py`, native diagnostics, and `git diff --check` successfully.
- [ ] Run constrained C++ and Pyomo optimization tests when a valid Gurobi license is available; the current license expired on 2023-10-25.
- [x] Replaced the brittle native command assertion with semantic `-F` token lookup, independent of trailing command options.
- [x] Added separate `tests/data/fw/fixed-matrix.tsv` and `tests/data/vw/fixed-matrix.tsv` fixtures with scenario-specific cluster IDs and A/B/U entries.
- [x] Parameterized Python resolution/validation and the standalone C++ parser harness over both fixed-width and variable-width SEG scenarios (`53 passed` focused tests).
- [x] Added and collected four licensed end-to-end cases covering `{cpp, gurobi} x {fixed-width, variable-width}` and asserting every fixture entry in `best.bbc.ucn`.
- [x] Fixed quiet C++ execution: `execute()` now includes stderr captured by `communicate()` in final-objective and failure parsing when outer verbosity is below 2.
- [x] Added a mocked quiet-execution regression proving a native `Final objective` line returns the parsed objective instead of a false missing-objective error.
- [x] Re-verified semantic `-F` forwarding and Python/native fixture handling for both fixed-width and variable-width scenarios (`54 passed` focused tests).
- [ ] Execute the four licensed dual-width integration cases in an environment with a valid Gurobi license.
- [x] Added and registered `hatchet fixed-matrix-editor` with browser launch and headless `--no-open` behavior.
- [x] Added the packaged, self-contained editor with local SEG/TSV import, A/B/U grids, filtering, keyboard navigation, multi-cell paste, live validation, clear actions, and canonical TSV download.
- [x] Added launcher/resource/feature tests and verified the full fixed-matrix suite (`57 passed`).
- [x] Verified the embedded JavaScript parses successfully and `hatchet.data` package discovery includes the 33 KB editor asset.
- [x] Rendered and inspected nonblank, overlap-free desktop (1440x900) and mobile (390x844) screenshots with headless Firefox.
- [ ] Automated loaded-SEG browser interaction remains pending because Playwright, Selenium, and geckodriver are unavailable in this environment.
- [x] Documented A/B/U dimensions and semantics, the versioned sparse TSV contract, editor workflow, CLI/INI/environment usage, exact clone-count policy, backend support, and interactions with existing constraints.
- [x] Added a documentation contract test covering the marker, editor command, CLI/INI keys, exact clone count, environment override, and dual-ploidy behavior (`58 passed` focused tests).
- [x] Built Sphinx HTML successfully and confirmed the generated compute-cn page contains the new guide. The build retains 58 pre-existing repository-wide theme, toctree, and link warnings.
- [ ] Complete licensed all-mode integration tests and automated loaded-SEG browser interaction when the required license and browser tooling are available.
