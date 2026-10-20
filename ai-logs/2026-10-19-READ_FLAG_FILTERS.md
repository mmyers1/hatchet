# Exclude Unmapped/Secondary/QC-fail/Duplicate/Supplementary Reads From All Read Counting

## Status

- Overall: Implemented. All 6 reported integration failures have been explained (2 are expected effects of the
  change, 4 are unrelated). The `np.nan` fix is done. Fixture regeneration (`regenerate_fixtures.sh`) is waiting on
  the user's run.
- Last updated: 2026-10-19

## Goal

Make every read-counting path use the same flag filter, `-F 3844` (UNMAP 4 + SECONDARY 256 + QCFAIL 512 + DUP 1024 +
SUPPLEMENTARY 2048), or `--ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY` for bcftools mpileup.

Before this change each tool used its own default:

| Call site | Previous filter | Counted SUPPLEMENTARY? | Counted DUP/SECONDARY/QCFAIL? |
|---|---|---|---|
| `count_reads.count_chromosome` (`samtools view -q`) | none | yes | yes |
| `count_reads.run_mosdepth` | mosdepth default 1796 | yes | no |
| `count_alleles` / `genotype_snps` (`bcftools mpileup`) | UNMAP,SECONDARY,QCFAIL,DUP | yes | no |
| `TotalCounting` (`samtools view -c -q`, fixed-width path) | none | yes | yes |
| `BAMBinning` (pysam `count` + `read_callback`) | MAPQ only (a callable replaces pysam's default filter) | yes | yes |

## Method

- Added two shared constants to `src/hatchet/utils/Supporting.py`: `EXCLUDE_FLAGS = 3844` and
  `MPILEUP_EXCLUDE_FLAGS = 'UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY'`. Every call site uses them, so the five paths
  cannot drift apart.
- The scope (genotype_snps and the fixed-width path in addition to count_reads/mosdepth/count_alleles) was confirmed
  with the user.

## Results / Verification

- `tests/test_counting_equivalence.py`: 319 passed. Command:
  `PATH=/home/myersm2/miniforge3/envs/hatchet/bin:$PATH PYTHONPATH=src conda run -n hatchet-test pytest tests/test_counting_equivalence.py -q`
  - `hatchet-test` has a **non-editable** install of hatchet in site-packages, so `PYTHONPATH=src` is needed for the
    tests to exercise the working tree.
  - `hatchet-test` has no samtools, so the BAM test was skipped until the `hatchet` env's bin was added to PATH.
  - The extended `test_count_chromosome_writes_read_starts` fails against the old installed code (it returns
    `10,100,150,200`) and passes with the change (`10,200`).
- Synthetic BAM smoke test (5 reads at the same position with flags 0, 0x800, 0x400, 0x100, 0x200; tools from the
  `hatchet` env: samtools 1.21, mosdepth, bcftools):
  - `samtools view -c`: 5 by default → 1 with `-F 3844`
  - mosdepth depth: 2 with `-F 1796` (the default) → 1 with `-F 3844`
  - `bcftools mpileup` DP: 2 by default → 1 with `--ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY`

## Follow-up for the user

- Integration tests (not run; they may be slow):
  `PYTHONPATH=src conda run -n hatchet-test pytest tests/test_steps.py tests/test_steps_vw.py -q`
  If the test BAMs contain supplementary or duplicate reads, the expected counts will shift.
- Reinstall hatchet into `hatchet-test` (or `pip install -e .`) so that `hatchet` commands pick up the change.
- Existing preprocessing outputs are reused when present (`*.starts.gz`, `*.per-base.bed.gz.csi`, chromosome arrays).
  Delete them before rerunning count-reads or count-alleles on samples that were already processed, or the old
  unfiltered counts will remain.
- Sanity check on a real BAM: `samtools view -c -F 3844 -q <readquality> <bam> chr1` should equal the number of lines
  in `<sample>.chr1.starts.gz`.

## Integration test failures (user run, `source iris_setup.sh`, from the repo root)

The user reported 6 failures. Each was traced to a cause:

| Test | Cause | From this change? |
|---|---|---|
| `test_steps.py::test_count_reads` (fw) | `bulk.1bed` has fewer reads in 2 bins: Tumor1 14576 vs 14582, Tumor2 14930 vs 14935. Normal is unchanged. The `DataFrame.columns` message appears because `.1bed` has no header, so the changed first row is read as column names. | **Yes, expected.** In the test region, `samtools-1.7 view -c -q 11 -f 2048` gives 6 supplementary reads in `bulk_03clone1_06clone0_01normal.sorted.bam` and 0 in `normal.bam`. No reads in either have flags 4/256/512/1024. |
| `test_steps_vw.py::test_count_reads` | Same BAMs: mosdepth and samtools now drop supplementary reads. | **Yes, expected** (inferred; the outputs were deleted before I could diff them, see below). |
| `test_steps.py::test_genotype_snps` | `HATCHET_PATHS_REFERENCE=reference/hg19.fa.gz` is plain gzip (header FLG=00), not bgzip, and has no `.fai`. bcftools fails with `[E::fai_build3] Cannot index files compressed with gzip, please use bgzip` (`tests/out/snps/Normal_chr22:..._bcftools.log`). | No. bcftools-1.7 fails the same way without `--ff`. |
| `test_steps.py::test_count_alleles_normal_snps` | Same reference (mpileup `-f`). Even with a valid reference, this test and `test_genotype_snps` had fixture mismatches with bcftools 1.19; see `ai-logs/2026-09-24-preprocessing-efficiency/run_tests.sh`. | No (inferred: same mpileup `-f` reference). |
| `test_steps_vw.py::test_combine_counts`, `..._phasing` | `combine_counts.py:112` uses `np.NAN`, which NumPy 2.0 removed; `hatchet-test` has numpy 2.0.2. These tests read static fixtures in `tests/data/vw/rdr`, so read counting does not affect them. | No. |

Also checked: bcftools-1.7 accepts `--ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY` (synthetic BAM DP 2 → 1), and
samtools-1.7 accepts `-F 3844`.

**Incident:** while getting the combine_counts traceback, I ran `pytest tests/test_steps_vw.py::test_combine_counts` in
the repo. Its module fixture runs `rmtree` on `tests/out/vw`, which deleted the user's vw outputs from their run
(regenerable test outputs). Later runs used a copy of `tests/` in the session scratchpad.

## Fixture regeneration (user chose: regenerate the whole chain)

Dependency chain of the fixtures (input → command → expected output, which is the next step's input):

- fw: `fw/bin/*.1bed` → combine-counts-fw → `fw/bb/bulk.bb` → cluster-bins-gmm → `fw/bbc/bulk.{seg,bbc}` → compute-cn →
  `fw/results/best.bbc.ucn` (also read by plot_bins/plot_cn)
- fw → vw: `fw/bb/bulk.bb` → cluster-bins → `vw/bbc/{bulk,bulk_onecl}.{seg,bbc}` → `test_solver.py` (with the
  hand-written `vw/fixed-matrix*.tsv`, which name clusters by ID)
- vw: `vw/rdr/*` → combine-counts → `vw/bb/bulk_{no,yes}phase.bb`

Not regenerated, because no test produces them: `fw/baf`, `vw/baf` (allele counts; the original count-alleles
parameters are unknown), `snps/*`, the top-level `fw/bulk.*` and `fw/best.bbc.ucn`, and `*/fixed-matrix*.tsv`.

`ai-logs/2026-10-19-READ_FLAG_FILTERS/regenerate_fixtures.sh` runs the stages `count combine cluster solve verify`
in order (the user runs it: mosdepth needs ~4 GB and compute_cn needs Gurobi). Checked so far:

- `bash -n` passes.
- The guards stop the script when HATCHET_PATHS_REFERENCE is unset or not faidx-indexed.
- The fw combine invocation reproduces the current `fw/bb/bulk.bb` byte for byte from the current inputs.

The fixtures have not been changed yet.

Steps for the user:

```bash
cd /data1/shahs3/users/myersm2/repos/hatchet
# 1. an indexed reference (bcftools/bedtools cannot read plain gzip); ~3 GB
zcat reference/hg19.fa.gz > reference/hg19.fa
/home/myersm2/miniforge3/envs/hatchet/bin/samtools faidx reference/hg19.fa
#    then in iris_setup.sh: export HATCHET_PATHS_REFERENCE="/data1/shahs3/users/myersm2/repos/hatchet/reference/hg19.fa"
# 2. on a compute node
source iris_setup.sh
ai-logs/2026-10-19-READ_FLAG_FILTERS/regenerate_fixtures.sh 2>&1 | tee ai-logs/2026-10-19-READ_FLAG_FILTERS/regenerate.log
```

Expected and to check: the verify stage passes. If the `variable-width` fixed-matrix tests in `test_solver.py` fail,
the regenerated `vw/bbc` clusters no longer match the cluster IDs or values in `vw/fixed-matrix*.tsv`, and those
files must be edited by hand. `test_genotype_snps` and `test_count_alleles_normal_snps` should pass with an indexed
reference and bcftools-1.7: `normal.bam` has no reads with flags 4/256/512/1024/2048 in the test region.

## Exhaustive list of changes

1. `src/hatchet/utils/Supporting.py`: added `EXCLUDE_FLAGS` and `MPILEUP_EXCLUDE_FLAGS`, so all call sites share one
   definition.
2. `src/hatchet/utils/count_reads.py`
   - import `EXCLUDE_FLAGS`
   - `count_chromosome`: `samtools view -F 3844`, so read starts drop unmapped/secondary/QC-fail/dup/supplementary
     reads (previously none were filtered).
   - `run_mosdepth`: `-F 3844`, so per-base coverage also excludes supplementary reads (default was 1796).
   - `total_counts_from_arrays` docstring: now gives the `samtools view` equivalent with `-F 3844`.
3. `src/hatchet/utils/count_alleles.py`: mpileup `--ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY`, so allele counts
   exclude supplementary reads.
4. `src/hatchet/utils/genotype_snps.py`: the same `--ff`, so SNP calling and allele counting see the same reads.
5. `src/hatchet/utils/TotalCounting.py`: `samtools view -c -F 3844`, so fixed-width total counts match.
6. `src/hatchet/utils/BAMBinning.py`: `read_callback` also requires `not (r.flag & EXCLUDE_FLAGS)`. A callable
   `read_callback` replaces pysam's default filter, so dups/secondary/QC-fail reads were previously counted here.
8. `src/hatchet/utils/combine_counts.py:112`: `np.NAN` → `np.nan` (approved by the user). NumPy 2.0 removed `np.NAN`,
   which broke both vw combine_counts tests. Verified: those tests now get past line 112 and stop at line 147
   (`bedtools nuc`), which fails on the plain-gzip reference.
9. `ai-logs/2026-10-19-READ_FLAG_FILTERS/regenerate_fixtures.sh` (new): staged regeneration of every fixture that
   depends on read counting, for the user to run (it is too expensive to run here).
7. `tests/test_counting_equivalence.py`: `test_count_chromosome_writes_read_starts` now writes a supplementary and a
   duplicate high-MAPQ read on chr1, and the expected output is unchanged (`10\n200\n`). This covers the new filter
   while keeping the BAM tiny (6 reads).

## Tasks prevented by sandbox / environment limits

- None blocked by the sandbox.
- `flake8` is not installed in `hatchet-test`. Workaround: checked the edited files for lines over 120 characters
  (the `.flake8` limit) with awk, and rewrapped one docstring.
- The first synthetic-BAM attempt failed because `conda run` does not forward a stdin heredoc. Workaround: wrote the
  generator script to the session scratchpad and ran it as a file.
