# CBC Test Segfault and Sample-Clonal Pyomo W1002 Warnings

## Status

- Overall: Implemented; full-suite and osteo verification pending (user to run)
- Last updated: 2026-10-01

## Goal

1. Fix `tests/test_solver.py::test_solve_command_cbc`, which fails with
   `Solver (cbc) returned non-zero return code (-11)` (`test_log.txt`).
2. Remove the Pyomo warnings that flood `compute_cn` logs whenever `sampleclonal` is used, even with only two states
   per sample. The osteo IID_H135990 `hatchet.log` (gurobi, CD-only `M: 2`, `clones=4,8`, 3 samples) contains 6,870
   warnings such as `WARNING (W1002): Setting Var 'cA_4_2' to a numeric value 1 outside the bounds (3, 3).`

## Method

### CBC segfault (environment)

- `cbc -quit` alone segfaults in `hatchet-test` (return code 139), so HATCHet code is not involved.
- `LD_LIBRARY_PATH` is empty, and `ldd` resolves every library inside the env.
- `LD_DEBUG=libs` shows every library initializer completing, so the crash happens after control passes to `main`.
- `catchsegv` backtrace: `free()` on an invalid pointer inside `CbcMain1` (`libCbcSolver.so.3`), which indicates an
  ABI mismatch between COIN-OR libraries.
- `conda-meta/history`:
  - 2026-08-28: `mamba install -c conda-forge coincbc` installed a matched stack: cbc 2.10.13, utils 2.11.13,
    osi 0.108.12, clp 1.17.11, cgl 0.60.10.
  - 2026-09-25: `mamba install bcftools=1.7 samtools=1.7` downgraded it to cbc 2.10.12 (build `h8b142ea_1`, which
    requires libstdcxx>=13, so a recent build) on top of 2022-era utils 2.11.6, osi 0.108.7, clp 1.17.6 and cgl 0.60.6.
    The downgrade was allowed because cbc pins only `coin-or-utils >=2.11,<2.12` (and similar for its other COIN-OR
    dependencies).
  - bcftools and samtools were removed later, but the downgraded COIN-OR stack stayed. This is why the 2026-09-29 log
    saw the same failure on unchanged HEAD.
- A `mamba install --dry-run` of the matched versions resolved cleanly. It upgrades the five COIN-OR packages and
  `coincbc`, plus zlib, sqlite and tk (the other pieces the same transaction had downgraded), and installs libsqlite.

### W1002 warnings (code)

- Coordinate descent: `Worker.run` in `src/hatchet/utils/solve/cd.py` calls `fix_u(seed)`, then `create_model()`, then
  `hot_start(_cA, _cB)`. On the first iteration `_cA/_cB` are `cd.hcA/hcB`, built once by `first_hot_start()` and
  shared by all seeds. That start includes the normal clone, global clonal states and fixed-matrix cells, but not
  sample-clonal states, which depend on which clones each seed's U makes present.
- In CARCH mode, `compute_fixed_bounds` narrows present clones' `cA/cB` bounds to the sample-clonal state.
  `hot_start` then writes out-of-bound values. Pyomo still sets them but logs W1002, and Gurobi receives an
  infeasible MIP start.
- Reproduced without a solver using the `_sample_clonal_ilp()` test helper: 4 `pyomo.core` records with
  `id='W1002'`.
- Later CD iterations are clean: UARCH excludes clones whose `_cA/_cB` violate a sample's states, so present clones
  already match.
- The native backend has the same start gap, but Gurobi's C++ `Start` attribute discards it silently, so that path
  is out of scope.

## Progress

- [x] Diagnosis of both problems; plan approved.
- [x] Failing regression test written.
- [x] `hot_start` clamp implemented; test passes.
- [x] `hatchet-test` COIN-OR stack repaired; `cbc -quit` and a Pyomo CBC solve succeed.
- [ ] User runs the CBC regression, the full suites and the osteo rerun.

## Results

- New regression test before the fix: failed and captured `pyomo.core` W1002 records, e.g.
  `Setting Var 'cB_2_2' to a numeric value 0 outside the bounds (1, 1)`.
- After the fix: `pytest -q tests/test_sample_clonal.py tests/test_fixed_matrix.py -k "hot_start or carch or uarch or
  seeds or binwise or split"` gave 8 passed.
- After the fix: solver-free tests across both feature suites
  (`-k "not solve and not solves and not compute_cn and not gurobi and not native and not cpp_solver"`) gave 120 passed.
- `git diff --check` is clean, and no touched line exceeds 120 characters.
- Env repair: before the change `hatchet-test` was at revision 6. `mamba install` linked coin-or-cbc 2.10.13,
  coin-or-utils 2.11.13, coin-or-osi 0.108.12, coin-or-clp 1.17.11, coin-or-cgl 0.60.10, coincbc 2.10.13, zlib/libzlib
  1.3.2, sqlite/libsqlite 3.53.4 and tk 8.6.13. Afterwards, `cbc -quit` returns 0 (it previously segfaulted with 139),
  and a 1-variable Pyomo integer model solved with `SolverFactory('cbc')` is `optimal` (x = 3). To revert:
  `conda install -n hatchet-test --revision 6`.

### Commands for the user

```bash
cd /data1/shahs3/users/myersm2/repos/hatchet
# CBC regression. The earlier test_log.txt run imported the stale site-packages copy; PYTHONPATH=src avoids that.
conda run -n hatchet-test env PYTHONPATH="$PWD/src" pytest -q tests/test_solver.py::test_solve_command_cbc
# full suites (after `module load gurobi` and installing the rebuilt native solver to src/hatchet/solve)
conda run -n hatchet-test env PYTHONPATH="$PWD/src" pytest -q \
  tests/test_sample_clonal.py tests/test_fixed_matrix.py tests/test_solver.py
```

If the CBC test now runs but fails in `assert_frame_equal` against `tests/data/fw/best.bbc.ucn`, that is a
solution-level difference, not the segfault.

The osteo run uses the `hatchet` env, which holds a non-editable install of this repo, so reinstall before resubmitting:

```bash
module load gurobi && conda run -n hatchet pip install /data1/shahs3/users/myersm2/repos/hatchet
sbatch /data1/shahs3/users/myersm2/repos/osteo/results/hatchet2/2sample/IID_H135990/run_IID_H135990_heavy.sh
grep -c "W1002" /data1/shahs3/users/myersm2/repos/osteo/results/hatchet2/2sample/IID_H135990/hatchet.log   # expect 0
```

## Summary

- CBC: the segfault came from the `hatchet-test` env, not HATCHet. A 2026-09-25 bcftools/samtools install
  downgraded the COIN-OR libraries under a newer CBC build, an ABI mismatch that crashes `CbcMain1`. Restoring the
  matched 2026-08-28 stack fixes the binary. No code changed.
- W1002: in coordinate descent, every seed's first C-step was warm-started from one shared start that does not know
  the sample-clonal states of the clones that seed's U makes present. `hot_start` now clamps copy-number starts to the
  model's variable bounds, which already encode those states. Without sample-clonal input this is a no-op, because all
  other fixed bounds were already satisfied by the starts. With sample-clonal input, CD warm starts are now valid
  (feasible for Gurobi), so results from sample-clonal runs may differ slightly from earlier runs.

## Updates and Reasoning

1. `tests/test_sample_clonal.py`: added `import logging` and `test_carch_hot_start_respects_present_clone_states`.
   The test builds a CARCH model with clone 1 present, warm-starts it with values that violate that clone's
   sample-clonal states, and asserts no W1002 records, the clamped present-clone states, and unchanged absent-clone
   values. Reason: reproduces the reported warning without a solver, and was confirmed to fail before the fix.
2. `src/hatchet/utils/solve/ilp_subset.py::ILPSubset.hot_start` (also used by `ILPSubsetSplit`): clamps each `cA`/`cB`
   start to `self._bounds`, sets the bit variables from the clamped values, and computes the `_adA`/`_adB` starts from
   the clamped tumor values. Reason: removes the W1002 warnings and gives Gurobi a start consistent with the bounds.
   A comment explains why the shared first hot start can violate those bounds.
3. Env `hatchet-test`: `mamba install -c conda-forge "coin-or-cbc>=2.10.13" "coin-or-utils>=2.11.13"
   "coin-or-osi>=0.108.12" "coin-or-clp>=1.17.11" "coin-or-cgl>=0.60.10"`. Reason: restores the matched COIN-OR stack
   that the 2026-09-25 bcftools/samtools transaction downgraded. Revert with `--revision 6`.
