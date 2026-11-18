# Maximum Homozygous-Deletion Cluster Size Implementation

Expose the maximum cluster-weight fraction that may receive a homozygous deletion as `--max-homdel-size`. The default remains `0.005`, and clusters with weight fractions greater than or equal to the threshold must have nonzero total copy number in every tumor clone.

## Decisions

- CLI: `--max-homdel-size`
- INI: `maxhomdelsize`
- Environment: `HATCHET_COMPUTE_CN_MAXHOMDELSIZE`
- Python API: `max_homdel_size`
- Native API: `maxHomdelSize`
- Valid range: `[0, 1]`
- Cutoff: inclusive (`cluster_fraction >= max_homdel_size`)
- Coverage: standard and binwise Pyomo ILP/CD, native C++ ILP/CD, fixed-matrix validation, and sample-clonal validation

## Progress

- [x] **Phase 1: Public Python plumbing**
  - [x] Add packaged, script, and documentation INI defaults.
  - [x] Add CLI parsing and range validation.
  - [x] Forward the option to Python and native execution.
- [x] **Phase 2: Python enforcement and validation**
  - [x] Propagate the API through solve and coordinate-descent paths.
  - [x] Share model enforcement across standard and binwise models.
  - [x] Parameterize fixed-matrix and sample-clonal validation.
  - [x] Add focused Python tests.
- [x] **Phase 3: Native C++ parity**
  - [x] Parse and validate the native option.
  - [x] Propagate it through solve, coordinate descent, workers, and ILP models.
  - [x] Parameterize native validators and add model enforcement.
  - [x] Add focused native helper tests.
- [x] **Phase 4: Documentation and test coverage**
  - [x] Document units, defaults, cutoff semantics, and configuration names.
  - [x] Add focused Python and native regression tests.
- [x] **User verification**
  - [x] Build the native solver after `module load gurobi`.
  - [x] Run focused and full test suites in the `hatchet-test` Conda environment.

## Verification Log

Commands are run after `module load gurobi` using the `hatchet-test` Conda environment.

- Public parser/config precedence: passed.
- Standard and binwise FULL/CARCH model constraints: passed.
- Fixed-matrix and sample-clonal custom cutoff validation: passed.
- Focused Python modules: 159 passed; 8 native cases await the Phase 3 rebuild.
- Rebuilt Python/native parity suite: 170 passed.
- Final API-boundary test additions are pending user verification; no tests were run after the user requested manual verification.
