# Plan: faster, leaner count_reads / count_alleles / genotype_snps / phase_snps with identical output

Status: **done and verified locally; job 4 (WGS) fixed and ready to submit; job 5 (rerun on a compute node) pending** (branch `osteo`, not committed). Progress and results are logged at the bottom of this file.

## Record

Before any code changes, create **`ai-logs/2026-09-24-preprocessing-efficiency.md`** in the hatchet repo. It holds this plan, a status line and a dated `## Progress log`. The `ai-logs/` directory is new in this repo and follows the style of `read-overlap-hmmcopy/ai-logs`.

- Append an entry after each numbered part and each verification step: what changed, test counts, measured numbers (wall time, max RSS), and any `Ruling:` lines for judgment calls made along the way.
- Keep the status line current.
- Finish with a **`## Final summary`** section with these parts:
  - changes per file (line counts before → after);
  - old-vs-new equivalence results;
  - a time/memory table per step, old vs new, on the `data/` BAMs and the synthetic benchmark;
  - intentional deviations, if any;
  - outstanding items.
- Commit only if asked.

## Context

You asked to refactor the four preprocessing submodules for time and memory while keeping their output semantically identical. Exploration found these hotspots:

1. **genotype_snps and count_alleles read the whole BAM in every per-chromosome job.** They run `bcftools mpileup -T <targets>` without `-r`. `bcftools mpileup --help` (1.19) says `-T` "streams rather than index-jumps", so every (sample, chromosome) job decompresses the entire BAM.
2. **count_reads loops in pure Python:**
   - It runs a Python loop over every read start, all loaded into a Python `list[int]`, about 2 GB per sample for 30× chr1.
   - It runs another Python loop over every mosdepth per-base record, written to a temporary per-chromosome tabix dump first.
   - Each chromosome worker re-reads the whole BAF file.
   - A second full `samtools view -c` pass over every BAM exists only to fill `total.tsv`.
3. **count_alleles selects het SNPs in the parent with one scalar `beta.ppf` call per SNP**, about 85 µs each. In hatchet-test, the array call runs about 55× faster and is **bitwise identical**, checked on 64k count pairs × 5 gamma values. The parent also keeps nested tuples and dicts for every sample and chromosome in memory.
4. **phase_snps spends almost all its time in shapeit and picard.** Only small savings are available: a redundant copy and an extra compress/decompress round trip.

Decisions (user):
- Derive `total.tsv` from the count arrays instead of the extra samtools pass. Same numbers; on a resume it reflects whatever `-q` built the arrays.
- Verify on the BAMs in `data/`: `normal.bam` plus three `bulk_*` tumors, with `chr` contig names.
- Develop in the `hatchet-test` conda env; `module load gurobi` before building or installing.

## Semantics to preserve (all parts)

- Output files keep the same bytes and line order: VCF records, `.1bed` files, `*.total.gz`, `*.thresholds.gz`, `samples.txt` and `total.tsv`. Only VCF header lines that record bcftools command strings may differ.
- The public API is unchanged:
  - `count_alleles.counting()` keeps its return format; `tests/test_steps.py` compares it with a pickle.
  - `selectHetSNPs`, `isHet`, `checkShift`, `count_reads.read_snps`, `form_counts_array` and `genotype_snps.call` keep their signatures and return values.
- Every assert and error reachable from valid input stays the same:
  - The read-start binning assert depends on read order.
  - Coverage quirks for sample 0: `arr[n,0] = end_reads`, and samples i>0 use the extended thresholds.
  - `checkShift` raises `ZeroDivisionError` when A+B==0 on a site that passes `isHet`.
- count_alleles key and ordering rules stay the same:
  - Results are keyed by (sample, CHROM of the first output line), and later tasks override earlier ones.
  - A duplicate position keeps its first occurrence in the normal output, which is sorted. In the tumor output it keeps the first occurrence's position and the last occurrence's value.
  - The tumor sample `set` is iterated in the same order.

## 1. `src/hatchet/utils/count_reads.py`

- **`form_counts_array`, read starts (even columns):**
  - Read with `pd.read_csv(header=None, dtype=int64, chunksize≈2M)`, which infers gzip.
  - Assign each read to a row with `np.searchsorted(T, r, side='right') - 1` and count rows with `np.bincount`. Reads ≥ `T[-1]` go to row n-2.
  - The assert holds exactly as before when, up to the first such read, row indices are ≥0 and non-decreasing. That state carries across chunks.
  - Handle an empty file, and drop the unused second `open`.
  - Pure helper for this: `_bin_read_starts(chunks, thresholds)`.
- **`form_counts_array`, coverage (odd columns):**
  - Pipe `tabix <perpos> <chrom>` straight into chunked `pd.read_csv`, replacing the temporary `<perpos>.<chrom>` file.
  - Each threshold q=T-1 goes to the first record with end>q, found by `searchsorted` on the running max of `end` carried across chunks. This matches the sequential loop even when ends are not monotonic.
  - Keep the `start <= q` assert, the sample-0 checks (`idx == len(T)`, `chr_end > T[-1]`), the threshold extension and the `arr[idx,0] = end_reads` quirk.
  - Pure helper: `_coverage_at_thresholds(chunks, thresholds)`.
  - Make `arr` int64. `savetxt(fmt='%d')` writes the same text.
- **Read the BAF file once in `main`:**
  - `_load_baf_table` reads it once, with CHR and SAMPLE as categories.
  - `_snps_for_chromosome(table, ch, names)` converts CHR and SAMPLE back to object, then runs read_snps's existing validation and pivot.
  - Positions are built only for chromosomes that are not X/Y and whose outputs are missing. The table is freed before the Pool starts, and `run_chromosome` takes `positions` instead of `baf_file`.
  - `read_snps()` stays as a thin wrapper.
- **total.tsv:** replace `tc.tcount` with per-sample sums of `arr[:-1, 2i]` over each `<ch>.total.gz`, read with pandas. These equal `samtools view -c -q Q bam ch`, because the starts files come from `samtools view -q Q bam ch | cut -f4`. `TotalCounting.py` stays for `count_reads_fw`.
- **Minor fixes:**
  - `tabix -f` indexing runs in parallel across samples (ThreadPool, ≤ j).
  - `count_chromosome` closes its gzip output handle and checks gzip's return code.
  - `get_chr_end` waits for its `zcat | tail` processes.

## 2. `src/hatchet/utils/count_alleles.py`

- **`AlleleCounter`:**
  - Add `-r <chrom>` whenever the contig is in that BAM's header, read once per BAM with `pysam.AlignmentFile(...).references`. pysam is already a dependency.
  - When a tumor BAM lacks the contig, keep the old streaming command, which gives empty output. For a region name `chr:a-b`, check the part before the last `:`.
  - Parse the `bcftools query` stdout line by line instead of `communicate()` plus a split of the whole string. Close the parent's copy of `mpileup.stdout` and wait on both processes. The log file and errors stay the same.
  - Shared `_iter_rows(...)`; `countAlleles` returns `list(_iter_rows(...))`, so `counting()` is unchanged.
- **`_het_mask(A, B, gamma, maxshift)`**, vectorized:
  - It calls `beta.ppf` on arrays, then applies the shift check.
  - It raises `ZeroDivisionError('float division by zero')` where the original short-circuit `and` would.
  - `selectHetSNPs` is rebuilt on it with the same return value. The scalar `isHet` and `checkShift` stay.
- **Normal pass in the workers (`HetSelector`):**
  - Each row reduces to (pos, a0, c0, a1, c1), where a1/c1 is the first alt allele with the highest count.
  - `_het_mask` runs per task. The task returns `(sample, first_CHROM, HetSites | None)`, where HetSites is numpy arrays sorted by position, with first occurrences kept via `np.unique(return_index=True)`.
  - `main` applies the same override-then-drop-empty order, and writes the TMP target files and `outputNormal` with the same format strings as before.
- **Tumor pass in the workers (`TumorCounter`):**
  - Each worker inherits the normal's per-chromosome pos→(REF, ALT) as an attribute, via fork.
  - Rows map through `dict(r[3]).get(REF/ALT, 0)`. The result is compact arrays (pos, refcount, altcount) in first-occurrence order, with last-occurrence values.
  - `main` applies the same (sample, CHROM) override and `close()` message, and writes in `for sample in samples: for chro in chromosomes` order.
  - The raw `rcounts`/`counts` structures are never built.

## 3. `src/hatchet/utils/genotype_snps.py`

- With a snplist, add `-r {chromosome}` next to `-T tgtfile`, so each job index-jumps instead of streaming the whole BAM.
- Compress the targets file with `gzip -1` instead of `-9`. The content is identical, and the file is deleted after use.
- Close the parent's copies of the pipes and wait on every process.
- Replace the hand-rolled `Caller`/`JoinableQueue` with the shared `hatchet.utils.multiprocessing.Worker`, as count_alleles uses. This also fixes a hang: an exception in a task never calls `task_done`, so `tasks.join()` blocks forever. `call()` still returns the sorted list of VCF paths.
- Count the final SNPs in parallel (ThreadPool, ≤ j) instead of one file at a time.

## 4. `src/hatchet/utils/phase_snps.py`

The gains here are small but free.
- hg19 without `chr` names: drop the `stage_vcfs` copy and pass the input VCF straight to `biallelic`.
- hg19 with `-N`: pipe `bcftools annotate --rename-chrs … -Ou` into `bcftools view --max-alleles 2 --exclude-types indels -Oz -o {chr}_filtered.vcf.gz`. This is a Popen pair with return-code checks, as in `liftover`, and drops one compressed intermediate. `test_phase` covers it.
- `print_log` closes the files it opens.
- The hg38 liftover path is untouched, because no hg38 test data is available to verify it.

## 5. New tests: `tests/test_counting_equivalence.py`

These need no BAMs. Verbatim copies of the original loops and functions serve as references, compared against the new code on randomized inputs:
- `_bin_read_starts`: duplicate thresholds, reads ≥ the last threshold, an empty file, chunk sizes 1 and 3, and unsorted input, where both versions raise AssertionError.
- `_coverage_at_thresholds`: sample 0 and samples >0, plus non-monotonic ends.
- `selectHetSNPs` / `_het_mask`: tied alt counts, duplicate positions, A=B=0 (ZeroDivisionError in both), and gamma 0 and 1.
- The tumor mapping, against the original `form`/dict logic.

## Verification

Environment:
- Tools on PATH: `PATH=/home/myersm2/miniforge3/envs/hatchet/bin:$PATH`, which provides samtools 1.21, bcftools 1.19, tabix, bgzip, mosdepth, shapeit, picard and java. Always run `/home/myersm2/miniforge3/envs/hatchet-test/bin/python` explicitly.
- Reference: `hg19.fa.gz`, which is plain gzip, unpacked into the scratchpad, then `samtools faidx`.
- Environment variables: `HATCHET_TESTS_BAM_DIRECTORY=<repo>/data`, `HATCHET_PATHS_REFERENCE=<scratch>/hg19.fa`, and `PYTHONHASHSEED=0` so the tumor sample `set` is iterated in a reproducible order.

1. **Baseline, before any edit.** The installed hatchet-test package matches HEAD for these modules (`diff -rq` clean). Run a golden script from the repo root under `/usr/bin/time -v`, saving outputs plus time/RSS to `scratch/baseline/`:
   - `genotype_snps` on chr22, both without a snplist and with one. The snplist is the baseline chr22 VCF, tabix-indexed, which exercises `-r`+`-T`.
   - `count_alleles` with the normal and three tumors, with `-L` set to the genotype VCF.
   - `count_reads` with `-b tests/data/vw/baf/bulk.1bed -V hg19 --chromosomes chr22`.
   - `phase_snps` with `-V hg19 -N` on `tests/data/vw/snps/chr22.vcf.gz`, using `1000GP_Phase3/` in the repo root.
2. **Implement parts 1-5**, logging each in the progress log.
3. **Install:** `bash -lc 'module load gurobi && /home/myersm2/miniforge3/envs/hatchet-test/bin/pip install --no-deps .'`
4. **Golden diff:** rerun into `scratch/new/` and compare with the baseline. All must be identical:
   - VCF records, excluding `#` lines;
   - `normal.1bed` and `bulk.1bed`, byte for byte;
   - `zcat` of `chr22.total.gz` and `chr22.thresholds.gz`, plus `samples.txt` and `total.tsv`, byte for byte;
   - `phased.vcf.gz` records.

   Record wall time and max RSS per step, old vs new.
5. **pytest in hatchet-test:** `tests/test_counting_equivalence.py`, `tests/test_steps_vw.py::test_count_reads`, `tests/test_steps.py::test_genotype_snps`, `tests/test_steps.py::test_count_alleles_normal_snps` and `tests/test_phase.py`. Tests that also fail on the baseline, such as the genotype md5 pinned to bcftools 1.7 in CI, are reported as pre-existing.
6. **Synthetic benchmark** of the Python-bound kernels:
   - 20M sorted starts and 40M per-base records, old `form_counts_array` loops vs new;
   - 3M-SNP het selection, old vs new.

   Report time and peak memory.
7. **Write the `## Final summary`** in the ai-log.

## Progress log

- 2026-09-24: plan approved; log created.
- **Verification 1 done: baseline.** Golden script (`scratchpad/golden.sh`, `measure.py`) run with the unmodified hatchet-test install, `-j 4` on this 4-core node, `PYTHONHASHSEED=0`.
  - Inputs: the `data/` BAMs (normal + 3 tumors, mostly chr22); `genotype`/`alleles` on `chr21 chr22`; `reads` on `chr22 chrX` (so the sex-chromosome path runs) with `tests/data/vw/baf/bulk.1bed`; `phase` with `-V hg19 -N` on `tests/data/vw/snps/chr22.vcf.gz`. The `-R`/`-L` SNP lists for `genotype_tgt` and `alleles` are frozen copies of the baseline genotype output (`scratchpad/fixed_snps/`), so old and new code read identical inputs.
  - Reference: `/data1/shahs3/users/myersm2/reference/hg19_UCSC/hg19.fa`. The repo-root `hg19.fa.gz` is **truncated** (gzip "unexpected end of file"; the sequence stops inside chr19, so chr21/chr22 are missing), and the repo-root `hg19_no_chr.fa` holds a single partial contig `>1`. Both come from `download-panel` and would break the hg38 liftover path of phase-snps; see Outstanding.
  - Wall time and peak memory (peak RSS summed over the process tree / largest single process):

    | step | wall | CPU | peak tree RSS | max single RSS |
    |---|---|---|---|---|
    | genotype-snps | 202 s | 195 s | 0.61 GB | 0.23 GB |
    | genotype-snps with `-R` targets | 30 s | 38 s | 0.68 GB | 0.22 GB |
    | count-alleles | 89 s | 148 s | 2.04 GB | 0.99 GB |
    | count-reads | 177 s | 382 s | 14.13 GB | 4.03 GB |
    | phase-snps | 476 s | 470 s | 2.10 GB | 1.74 GB |
- **Part 1 done: count_reads.**
  - `form_counts_array` now calls `bin_read_starts` (searchsorted + bincount per chunk, carrying the "previous bin" state so the original order assertion fires in exactly the same cases) and `coverage_at_thresholds` (searchsorted on the running max of record ends, which equals the original sequential assignment even for non-monotonic ends). Files are read in chunks of 1e6 lines (`read_starts_chunks`; `read_perpos_chunks` streams `tabix` output directly, so the temporary `<perpos>.<chr>` dump is gone). `arr` is int64.
  - The BAF file is read once in `main` (`load_baf_table`, categorical CHR/SAMPLE) and validated per chromosome (`snps_for_chromosome`, the old `read_snps` body); workers receive the positions. `read_snps` is a thin wrapper.
  - `total.tsv` comes from `total_counts_from_arrays` (per-sample sums of the read-start columns, excluding the chromosome-end row). The `TotalCounting` pass over every BAM is gone from count-reads (`count_reads_fw` still uses it).
  - `tabix -f` indexing runs in a ThreadPool; `count_chromosome` closes its output handle and checks gzip's return code; `get_chr_end` reaps its processes.
  - Ruling: a malformed BAF file (e.g. a chromosome missing from it) now raises in `main` before any chromosome's arrays are written, instead of in its worker after the other chromosomes' arrays were written. Same error and message; only the partial outputs left behind differ. Cost if wrong: none for valid input.
  - `tests/test_counting_equivalence.py` (count_reads part): 300 random cases for `bin_read_starts` (including 119 where the original asserts, and disorder after the tail that the original does not check), 300 for `coverage_at_thresholds` (contiguous mosdepth-like records, 74 overlapping cases with non-monotonic ends, arbitrary intervals that assert, and thresholds past the last record), each at chunk sizes 1, 3 and 1000; and 20 end-to-end `form_counts_array` cases against a verbatim copy of the old function on real gzip/bgzip+tabix files, including an empty starts file and other chromosomes around the target. **620 passed.**
- **Part 2 done: count_alleles.**
  - `AlleleCounter.iterAlleles` streams `bcftools query` line by line (no `communicate()` of the whole chromosome) and adds `-r <chromosome>` when `Supporting.indexed_region` finds the BAM indexed with that contig (pysam header + `has_index()`, cached per process; ~2 ms). Otherwise the command is unchanged, so an unindexed BAM or a tumor missing the contig behaves as before. `countAlleles` returns `list(iterAlleles(...))`, so `counting()` is unchanged.
  - `het_mask` vectorizes `isHet and checkShift` (raising the same `ZeroDivisionError` for a heterozygous site with no reads); `het_candidates` + `select_het_sites` reproduce `selectHetSNPs` (first alt allele among ties; first heterozygous record of a repeated position), and `selectHetSNPs` is now built on them.
  - `main` uses `counting_het_sites` (`HetSiteCounter`: selection inside the workers, only `HetSites` arrays come back) and `counting_tumor_alleles` (`TumorAlleleCounter`: maps tumor counts to the normal's REF/ALT inside the workers and returns `TumorCounts` arrays; position of the first record, counts of the last). `collect_by_first_site` keeps `counting`'s keying by (sample, chromosome of the first site), later tasks replacing earlier ones.
  - Ruling: `-r` is added only when the BAM is indexed and has the contig, because the original `-T`-only command works on unindexed BAMs (count-alleles never checks for an index). Cost if wrong: none.
  - Ruling: a tumor site not among the normal's heterozygous SNPs (impossible with `-T`) now fails inside the worker (reported as "WORK i FAILED" with the KeyError traceback) rather than as a bare KeyError in `main`.
- **Part 3 done: genotype_snps.**
  - `-r <chromosome>` added next to `-T <targets>` (same `indexed_region` guard), targets file compressed with `gzip -1`, parent copies of the pipes closed, all processes waited.
  - `Caller` is now a `hatchet.utils.multiprocessing.Worker` (as in count-alleles), so a failing task raises instead of hanging `tasks.join()` forever; `call()` still returns the sorted VCF paths. The final SNP count runs in a ThreadPool.
- **Part 4 done: phase_snps.** hg19 without `chr`: the input VCF goes straight to the biallelic filter (no copy). hg19 with `-N`: `bcftools annotate --rename-chrs ... -Ou | bcftools view --max-alleles 2 --exclude-types indels -Oz` in one pipe with both exit codes checked (`biallelic(..., rename=...)`), so `{chr}_toFilter.vcf.gz` is no longer written. `stage_vcfs` removed; `print_log` closes its files. hg38 liftover untouched.
- **Part 5 done: tests.** `tests/test_counting_equivalence.py` adds, against verbatim copies of the original functions and `main` logic: 200 `selectHetSNPs` cases (gamma 0/0.05/0.2/1, maxshift 0-0.5; 40 raising ZeroDivisionError in both; ~1.7k repeated candidate positions, ~2.2k tied alt counts), 200 `het_mask` cases (30k count pairs; 99 raising), and 100 end-to-end normal + tumor passes (the worker classes driven with stubbed `iterAlleles`, compared on the normal output lines and TMP lists, the selected keys, and the tumor output lines in order). **1,120 passed** in total (620 count_reads + 500 count_alleles), run against `src/` with `PYTHONPATH=src`.
- **Verification 3 done: install.** `module load gurobi && pip install --no-deps .` into hatchet-test (GUROBI_HOME `/admin/software/gurobi/gurobi952/linux64`); the `solve` binary was rebuilt and the installed modules match `src/`.
- **Verification 4 done: golden diff, all identical** (`scratchpad/compare.sh baseline new`):
  - genotype-snps: chr21 (25) and chr22 (276,686) records identical; with `-R` targets, chr21 (25) and chr22 (276,685) identical. Headers differ only in the recorded commands (`-r chr22` now appears) and dates.
  - count-alleles: `normal.1bed` (262,126 lines) and `bulk.1bed` (785,668 lines) byte-identical.
  - count-reads: `chr22.total.gz`, `chr22.thresholds.gz` (251,467 lines), `chrX.total.gz`, `chrX.thresholds.gz` (29,522 lines, the sex-chromosome path), `samples.txt` and `total.tsv` byte-identical after decompression. So `total.tsv` from the arrays equals the old `samtools view -c` pass.
  - phase-snps: `phased.vcf.gz` records (276,157) and `phased.log` identical.
  - The same files are left in every output directory.
  - Timing on a shared node (load average ~68 from other users' jobs; this session's cgroup has 4 cores), so ±10% is noise:

    | step | wall old → new | CPU old → new | peak tree RSS old → new | max single RSS old → new |
    |---|---|---|---|---|
    | genotype-snps | 202 s → 194 s | 195 s → 202 s | 0.60 → 0.62 GB | 0.22 → 0.22 GB |
    | genotype-snps `-R` | 30 s → 30 s | 38 s → 32 s | 0.67 → 0.60 GB | 0.22 → 0.22 GB |
    | count-alleles | 89 s → 53 s (1.7×) | 148 s → 96 s | 1.99 → 1.23 GB | 0.97 → 0.25 GB |
    | count-reads | 177 s → 112 s (1.6×) | 382 s → 313 s | 13.80 → 13.80 GB | 3.94 → 3.94 GB |
    | phase-snps | 476 s → 502 s | 470 s → 492 s | 2.05 → 2.05 GB | 1.70 → 1.70 GB |
  - Reading the table: these BAMs are almost entirely chr22, so `-r` saves little here: only the chr21 task stopped streaming the whole BAM (genotype `-R` CPU 38 → 32 s). On WGS every chromosome's task streams the whole BAM today. count-reads' peak memory is mosdepth's (4 concurrent processes, ~3.5 GB each), not the Python code's; it is broken down per program below. phase-snps is one shapeit run (~7.5 min) whose time varies with node load; the changed part is isolated below.
- **count-reads, old vs new back to back, twice** (old code run from a `git archive HEAD` copy via `PYTHONPATH`; `measure.py` now also records the Python processes' summed RSS and each program's largest process):

  | run | wall | CPU | Python processes, summed peak RSS | largest Python process | mosdepth |
  |---|---|---|---|---|---|
  | old 1 / old 2 | 184 s / 173 s | 380 s / 388 s | 1.63 / 1.62 GB | 1.25 GB | 4.03 GB |
  | new 1 / new 2 | 109 s / 112 s | 308 s / 320 s | 0.95 / 0.96 GB | 0.36 GB | 4.03 GB |

  All four runs' arrays and `total.tsv` are identical to the baseline. The largest Python process shrinks 3.5×. The run's overall peak (14 GB) is mosdepth, 4 samples × ~3.5 GB, unchanged.
- **`-r` isolated** (one count-alleles mpileup + query on `normal.bam`, targets = the frozen SNP list): chr21 **5.5 s → 0.3 s**; chr22 25.1 s → 25.6 s. Output identical in both cases (7 and 276,598 lines). The old per-task time is set by decompressing the whole BAM whatever the chromosome, so on WGS each of ~24 tasks per BAM pays for the whole BAM; with `-r` each pays only for its chromosome.
- **phase-snps preprocessing isolated** (hg19 `-N` path on `tests/data/vw/snps/chr22.vcf.gz`, 3 repeats): old annotate `-Oz` then view: 2.86-2.95 s; new annotate `-Ou | view`: 1.20-1.66 s. Filtered records (276,614) and the `#CHROM` line identical. The 26 s difference in the full phase-snps run is shapeit's variance on the loaded node.
- **Verification 6 done: synthetic kernels** (`scratchpad/bench.py`, each case in its own process under `measure.py`, old = the verbatim copies in the test file):

  | kernel | old | new | outputs |
  |---|---|---|---|
  | `form_counts_array`, 1 sample, 200 Mb chromosome: 20M read starts, 36.3M per-base records, 200k thresholds | 45.7 s, 0.91 GB | **12.5 s, 0.28 GB** (3.7×, 3.3× less memory) | arrays and thresholds equal |
  | het selection (`selectHetSNPs`), 3M candidate sites, 1.43M selected | 270 s | **10.2 s** (26×) | dicts equal |

  - The count-reads kernel scales per sample per chromosome: at 30× WGS (chr1 ≈ 50M reads, >100M per-base records) the old code holds a ~2 GB Python list per sample in each worker, and the new code stays at the chunk size (1e6 lines).
  - The het-selection memory figures (1.77 → 1.90 GB) are dominated by the benchmark's own 3M-row input list, which `selectHetSNPs` takes as its argument. In `count-alleles` itself that list is never built: `HetSiteCounter` reduces rows to columns as bcftools streams them, inside the workers.
- **Verification 5 done: existing tests**, run from a scratch copy of `tests/` with a scratch working directory (the step fixtures `rmtree` `tests/out`, and `download_panel` regenerates chain/rename files in its cwd, so the repo was left untouched), for the new code (installed) and the old code (`PYTHONPATH` = HEAD copy):

  | test | old code | new code |
  |---|---|---|
  | `test_steps_vw.py::test_count_reads` | passed | passed |
  | `test_phase.py::test_script` | passed | passed |
  | `test_steps.py::test_genotype_snps` | failed, md5 `d8706bd7…` | failed, same md5 `d8706bd7…` |
  | `test_steps.py::test_count_alleles_normal_snps` | failed | failed |
  | `test_counting_equivalence.py` (1,120 cases) | n/a | passed |

  The two failures are pre-existing: their fixtures (an md5 and `normal_snps.pik`) were generated with bcftools 1.7 in CI, and this environment has 1.19. Checked directly: `counting()` on the test's inputs returns the same 15,575 rows with old and new code (`old == new: True`); both differ from the pickle only in allele depths (e.g. `('A', 21)` in the pickle vs `('A', 23)` from bcftools 1.19).

## Final summary

### Changes (line counts HEAD → now)

- **`src/hatchet/utils/count_reads.py` (593 → 693).** Vectorized, chunked `form_counts_array` (`bin_read_starts`, `coverage_at_thresholds`, `read_starts_chunks`, `read_perpos_chunks`); tabix output streamed instead of dumped to a temporary file. BAF file read once (`load_baf_table`, `snps_for_chromosome`; `read_snps` kept as a wrapper). `total.tsv` from the arrays (`total_counts_from_arrays`) instead of a `samtools view -c` pass over every BAM. Parallel tabix indexing; pipe/handle fixes in `count_chromosome` and `get_chr_end`.
- **`src/hatchet/utils/count_alleles.py` (286 → 437).** `-r <chromosome>` alongside `-T` when the BAM is indexed; bcftools output streamed; het selection vectorized (`het_mask`, `het_candidates`, `select_het_sites`) and done in the workers (`HetSiteCounter`); tumor counts mapped in the workers (`TumorAlleleCounter`), returning compact arrays. `counting()`, `selectHetSNPs`, `isHet`, `checkShift` keep their behavior.
- **`src/hatchet/utils/genotype_snps.py` (307 → 193).** `-r <chromosome>` alongside `-T`; `gzip -1` for the temporary targets file; shared `Worker` instead of hand-rolled queues (fixes a hang on task failure); parallel SNP count.
- **`src/hatchet/utils/phase_snps.py` (383 → 388).** No copy of the input VCF (hg19); annotate piped into the biallelic filter (hg19 `-N`); `print_log` closes its files.
- **`src/hatchet/utils/Supporting.py` (221 → 240).** `indexed_region()` helper (pysam header + index check, cached per process).
- **`tests/test_counting_equivalence.py` (new, 471 lines).** 1,120 randomized old-vs-new equivalence cases against verbatim copies of the original code; needs no BAMs (the 20 file-level cases need tabix/bgzip on PATH and skip otherwise).

### Equivalence

On the `data/` BAMs, every output of all four commands is identical to the baseline: genotype-snps VCF records (with and without `-R`), count-alleles `normal.1bed`/`bulk.1bed` byte for byte, count-reads arrays for chr22 and chrX plus `samples.txt` and `total.tsv`, and phase-snps `phased.vcf.gz` records and `phased.log`. The old and new test outcomes are identical, and `counting()` returns identical rows.

### Time and memory (old → new)

| | old | new |
|---|---|---|
| count-alleles, data/ BAMs, chr21+chr22 | 89 s, largest process 0.97 GB | **53 s, 0.25 GB** |
| count-reads, data/ BAMs, chr22+chrX (2 repeats) | 173-184 s; Python peak 1.62 GB, largest 1.25 GB | **109-112 s; 0.95 GB, largest 0.36 GB** |
| one mpileup task for a chromosome that isn't most of the BAM (chr21) | 5.5 s | **0.3 s** |
| phase-snps preprocessing per chromosome (hg19 `-N`) | 2.9 s | **1.2-1.7 s** |
| `form_counts_array`, 1 sample, 20M starts / 36M per-base records | 45.7 s, 0.91 GB | **12.5 s, 0.28 GB** |
| het selection, 3M candidate sites | 270 s | **10.2 s** |
| genotype-snps (no `-R`), phase-snps (shapeit) | 202 s, 476 s | unchanged within node noise (194 s, 502 s) |

The biggest win is not visible on these test BAMs, which are almost all chr22: with only `-T`, **each (sample, chromosome) task of genotype-snps `-R` and count-alleles decompressed the entire BAM**. On WGS that is ~24 whole-BAM passes per BAM and per step; with `-r` each task reads only its chromosome. Likewise count-reads no longer makes its extra full pass over every BAM for `total.tsv`.

### Intentional deviations (none change outputs)

1. A malformed BAF file now raises in `main` before any chromosome's arrays are written (same error), instead of after the other chromosomes' arrays were written.
2. A tumor site missing from the normal's heterozygous SNPs (impossible with `-T`) fails inside the worker ("WORK i FAILED" + KeyError traceback) rather than as a bare KeyError in `main`.
3. genotype-snps: a failing task now raises instead of hanging forever; its progress bar is the shared `Worker`'s.
4. On a resume where the arrays exist but `total.tsv` does not, `total.tsv` reflects the `-q` the arrays were built with (user decision).
5. VCF headers record `-r <chromosome>` in the mpileup command line; phase-snps (hg19) no longer writes `{chr}_toFilter.vcf.gz`.

### Outstanding

- **Broken panel reference files in the repo root** (`hg19.fa.gz` truncated, `hg19_no_chr.fa` partial): superseded on 2026-09-25, since this machine uses `/data1/shahs3/users/myersm2/reference/phasing-panel`, which is complete. The untracked repo-root copies are still broken and could be deleted.
- **Test fixtures pinned to bcftools 1.7:** `test_genotype_snps` (md5) and `test_count_alleles_normal_snps` (pickle) fail with bcftools 1.19 with old and new code alike. Regenerate them if 1.19 is the target.
- **count-reads peak memory is mosdepth** (~3.5-4 GB per concurrent sample for hg19 contig lengths); with the Python side reduced, `-j` for count-reads is bounded by mosdepth.
- **Pre-existing, left alone:** `phase_snps.cleanup()` only removes intermediates for chromosomes named 1-22 (X's remain). (`count_chromosome`'s skip check, listed here before, was fixed on 2026-09-25.)
- **hg38 liftover path of phase-snps:** optimized on 2026-09-25 (`filter -Ou | norm`), checked at the bcftools level; end-to-end check is job 3 of the Continuation section.
- **Not yet measured on a real WGS dataset**; all numbers above are from the `data/` BAMs and synthetic inputs. See job 4 of the Continuation section.
- Not committed. The verification scripts now live in `ai-logs/2026-09-24-preprocessing-efficiency/`; the 2026-09-24 run outputs were in the (temporary) session scratchpad.

## Continuation (2026-09-25)

User instructions for this round: use `refpaneldir=/data1/shahs3/users/myersm2/reference/phasing-panel` on this machine; do not run heavy jobs here, but propose jobs and verification steps for the user to run.

### Method

- The phasing-panel directory is complete, unlike the repo-root files: full-size `hg19.fa.gz` (949 MB), `hg19_no_chr.fa` (3.2 GB) with its `.dict`, all four chain files, both rename files, and the chrX panel files. With it the hg38 path of phase-snps (liftover hg38 → hg19 panel → hg38) can be exercised, so that path gets its `-Ou` change now, plus a job that verifies it end to end.
- Only light work runs here: unit tests, a seconds-long bcftools check, script syntax checks, the reinstall. The heavier checks are scripts in `ai-logs/2026-09-24-preprocessing-efficiency/` for the user to run.

### Progress log

- **phase-snps hg38 liftover:** `bcftools filter` now writes uncompressed BCF (`--output-type u`) into `bcftools norm`, as the hg19 path already pipes annotate into view. Light check on an indexed copy of the chr22 test VCF standing in for picard's output: `filter -Oz | norm` and `filter -Ou | norm` give identical records (276,614) and both write bgzip; only the recorded `##bcftools_*Command` lines differ. The full path (picard included) is job 3 below.
- **Tests trimmed to minimal size** (project CLAUDE.md): `tests/test_counting_equivalence.py` now runs **316 cases in 4.7 s** (was 1,120 in ~20 s). Every edge class is still covered: 32 read-start cases where the original asserts plus 20 with unchecked disorder past the last threshold; 20 each of contiguous, non-monotonic, asserting and partial coverage cases; 16 selection and 14 het-mask cases raising ZeroDivisionError; the empty starts file in `test_form_counts_array`.
- **hatchet-test reinstalled** (`module load gurobi && pip install --no-deps .`); installed modules match `src/`.
- **Verification scripts** in `ai-logs/2026-09-24-preprocessing-efficiency/` (syntax-checked; `compare.sh` checked against this session's scratch runs, where it reproduces the all-identical result):
  - `env.sh`: paths (hatchet-test python, tools from the `hatchet` env, phasing-panel, hg19 and GRCh38 references), `PYTHONHASHSEED=0`, `code_src old|new` (old code = `git archive 1d6ab0b src`, run via `PYTHONPATH`; new = `src/`), `measure`.
  - `measure.py`: wall/CPU time, peak RSS of the process tree, of the Python processes, and of the largest process per program.
  - `compare.sh <old> <new>`: every output file (VCF records + `#CHROM` line; gzipped files decompressed; others byte for byte; logs, indexes and shapeit's `*_alignments.log` skipped), then the time/memory table. Exit 1 on any difference.
  - `run_tests.sh <workdir> [old|new]`, `hg38_phase_check.sh <workdir>`, `wgs_compare.sh` (sbatch-ready; settings at the top).
- **Finding (pre-existing, not fixed):** `download_panel.mk_rename_file` writes `rename_chrs1.txt`/`rename_chrs2.txt` for chr1-22 only, and so do the phasing-panel copies. With hg19 BAMs using `chr` names (`-N`), chrX reaches shapeit as `chrX` instead of `X` (the panel's name), and is not renamed back. The "support phasing X" commit (99db3b4) changed only `phase_snps.py`. Fix: add the X line to `mk_rename_file` and to the two files. The hg38 path is unaffected (its chain files drop `chr` for every contig).

### Jobs for the user to run

All paths are absolute; the work directories below are suggestions (anywhere with a few GB free).

1. **Unit tests** (light, any node; 5 s):
   ```bash
   cd /data1/shahs3/users/myersm2/repos/hatchet
   PATH=/home/myersm2/miniforge3/envs/hatchet/bin:$PATH ~/miniforge3/envs/hatchet-test/bin/python -m pytest -q tests/test_counting_equivalence.py
   ```
   Expect `316 passed`.
2. **Existing step tests with the phasing panel**, new and old code (~15 min each, 4 cores):
   ```bash
   V=/data1/shahs3/users/myersm2/repos/hatchet/ai-logs/2026-09-24-preprocessing-efficiency
   sbatch -p cpushort -c 4 --mem 16G -t 1:30:00 -o tests-new-%j.log --wrap "bash $V/run_tests.sh /data1/shahs3/users/myersm2/hatchet-preproc-verify/tests new"
   sbatch -p cpushort -c 4 --mem 16G -t 1:30:00 -o tests-old-%j.log --wrap "bash $V/run_tests.sh /data1/shahs3/users/myersm2/hatchet-preproc-verify/tests old"
   ```
   Expect, for both: `test_steps_vw.py::test_count_reads` and `test_phase.py::test_script` pass; `test_steps.py::test_genotype_snps` fails with md5 `d8706bd7…` and `test_count_alleles_normal_snps` fails (pre-existing: bcftools 1.7 fixtures); new also `316 passed` equivalence cases.
3. **hg38 phase-snps, old vs new** (~30 min, 1-2 cores, picard `-Xmx8g`):
   ```bash
   sbatch -p cpushort -c 2 --mem 16G -t 2:00:00 -o hg38-phase-%j.log --wrap "bash $V/hg38_phase_check.sh /data1/shahs3/users/myersm2/hatchet-preproc-verify/hg38"
   ```
   It builds a `.dict` for the GRCh38 FASTA inside the work directory, lifts the chr22 test VCF to hg38 with picard, then runs `phase-snps -V hg38 -N -D <phasing-panel>` with both code versions. Expect `IDENTICAL  phase/phased.vcf.gz` and `IDENTICAL  phase/phased.log`, no `DIFFERENT`/`MISSING`.
4. **Real WGS, old vs new** (hours; 8 cores, 64 GB): edit the settings block of `$V/wgs_compare.sh` (BAMs, sample names, reference, `CHROMS`, work directory), then `sbatch $V/wgs_compare.sh`.
   - With a SNP list it first extracts the chosen chromosomes from it into the work directory (genotype-snps needs a `.tbi`; the lists in `reference/snp-lists` only have `.csi`).
   - Expect every file `IDENTICAL`. genotype-snps and count-alleles should drop from roughly one whole-BAM pass per chromosome task to one pass over each chromosome; count-reads' Python memory should drop to about the chunk size per worker (mosdepth remains the peak); phase-snps about unchanged.

**What to check in each log:** `grep -E 'DIFFERENT|MISSING' <log>` prints nothing, the pytest summary matches the expectations above, and the time table at the end of the compare output. Paste the tables back to have them added to this log.

## Round 3 (2026-09-25): count_chromosome skip check; the user's test run

### Progress log

- **`count_chromosome` skip check fixed** (user request). It checked `<sample>.<chr>.starts` but writes `<sample>.<chr>.starts.gz`, so a rerun never skipped finished files. It now checks `.starts.gz`.
  - Ruling: with a working skip, a partial file from an interrupted run would be taken as complete. So the pipeline now writes to `<sample>.<chr>.starts.gz.tmp`, renames it with `os.replace` once samtools, cut and gzip all exit 0, and deletes it when any of them fails. Cost if wrong: none; the finished file is the same.
  - Caveat: a partial `.starts.gz` left by an interrupted run of the *old* code (which wrote the final name directly) would now be skipped. Delete leftover `*.starts.gz` files after an interrupted old run.
  - Three tests added to `tests/test_counting_equivalence.py`: an existing `.starts.gz` is skipped without running samtools (`false` stands in for samtools), a failing samtools leaves no output or temporary file, and a real run on a 4-read BAM built with pysam writes exactly the MAPQ-passing chr1 starts. Both behavioral tests fail against the pre-refactor code. **319 passed** in 4.9 s; hatchet-test reinstalled with `module load gurobi`.
- **User's `run_tests.sh new` run: 4 failures.**
  - `test_steps.py::test_genotype_snps` (md5 `d8706bd7…`) and `test_count_alleles_normal_snps`: the expected, pre-existing bcftools 1.7 fixtures.
  - `test_steps_vw.py::test_count_reads`: `subprocess.CalledProcessError` from `mosdepth -t 1 -Q 11 …`. mosdepth runs before any code this refactor changed, with an unchanged command; this test passed with old and new code on 2026-09-24.
  - `test_phase.py::test_script`: `WORK 0 FAILED`. The surviving shapeit log shows the input was right (276,157 SNPs, the same as on 2026-09-24, so the new annotate | view pipe worked) and that shapeit's phasing stopped mid "Burn-in iteration [1/1]" without an error message, about 4 minutes into a ~7.5 minute run. The shapeit command is unchanged. The failing process was killed or crashed.
  - The error output that would explain the two was deleted: in the old `run_tests.sh`, `test_steps.py` ran last, and its fixture deletes `tests_<version>/out`, including the `vw/` outputs where mosdepth's and shapeit's error logs were. That was a flaw in the script.
  - Checked and ruled out, here: the environment (`module load gurobi` only sets `LD_LIBRARY_PATH` to gurobi's lib directory; mosdepth runs normally with it); OOM kills of the user's processes on isca020 (the kernel log has none; this session's cgroup counts `oom_kill 0`). The run did not go through SLURM (no job in `sacct`), so it ran in one of the interactive sessions (isca020, iscb006, isca023, isca081, iscb004); in any other of those, a memory limit or other load cannot be checked from here. mosdepth (~4 GB) and shapeit (~2 GB) are the two largest processes of these tests, which fits memory pressure, but that is not established.
  - Also noted: the repo-root panel files (`1000GP_Phase3/`, `hg19*`, chain and rename files) are gone since 13:48, before the test run and while this session was idle; the tests now use the phasing-panel directory, so they are not needed.
- **`run_tests.sh` rewritten** so a failure can be diagnosed: `test_steps.py` → `test_steps_vw.py` → `test_phase.py` order (no fixture deletes another test's outputs); the full pytest output goes to `<workdir>/pytest_<version>.log` together with the host, the session cgroup's memory limit, peak usage and `oom_kill` count before and after; an optional third argument is a pytest `-k` selection. Smoke-tested on isca020 with the three `count_chromosome` tests.

### Jobs for the user to run

5. **Rerun the two unexplained tests**, new then old, in a fresh work directory (~12 min each; give the session ≥ 8 GB):
   ```bash
   V=/data1/shahs3/users/myersm2/repos/hatchet/ai-logs/2026-09-24-preprocessing-efficiency
   W=/data1/shahs3/users/myersm2/hatchet-preproc-verify/tests2
   bash $V/run_tests.sh $W new 'test_count_reads or test_script'
   bash $V/run_tests.sh $W old 'test_count_reads or test_script'
   ```
   Or through SLURM: `sbatch -p cpushort -c 4 --mem 16G -t 1:00:00 -o tests2-%j.log --wrap "bash $V/run_tests.sh $W new 'test_count_reads or test_script'"`.
   - If both pass: the first failure was environmental.
   - If a test fails: `$W/pytest_new.log` has the mosdepth or phase-snps error; the phase-snps logs are in `$W/tests_new/out/vw/phase/phase/*.log`, and the last lines of `$W/pytest_new.log` give the session's memory peak and `oom_kill` count. If `oom_kill` went up, rerun with more memory.
- **User's rerun (`tests2`, samtools and bcftools 1.7): 321 passed, 2 failed; cause found, environmental.**
  - With bcftools 1.7, `test_genotype_snps` and `test_count_alleles_normal_snps` now **pass**: their fixtures were the only problem, as expected.
  - `test_count_reads` and `test_script` (phase) still failed. The new log shows why: the run was on **`islogin01`, a login node**. The session's memory peaked at exactly **2.0 GB**, and its cgroup's `oom_kill` counter went from 99 to **104** during the run: all four mosdepth processes died with `SIGKILL` (mosdepth needs ~4 GB), and shapeit was killed in its burn-in (~2 GB), leaving the same truncated log as in the first run. Nothing in hatchet failed; the old code fails the same way there. The 2026-09-24 runs that passed were on isca020 (64 GB SLURM session).
  - `run_tests.sh` now refuses to start when the smallest memory limit on its cgroup path is under 8 GB (message names the host; `ALLOW_LOW_MEMORY=1` overrides). Checked on isca020 (64 GB: runs) and on a simulated login-node cgroup capped at 2 GB on the user slice (stops); islogin01's real cgroup layout could not be inspected from here.
  - Next step for the user: rerun job 5 on a compute node, e.g. `sbatch -p cpushort -c 4 --mem 16G -t 1:00:00 -o tests2-%j.log --wrap "bash $V/run_tests.sh $W new 'test_count_reads or test_script'"`, and the same with `old`.

## Round 4 (2026-09-28): WGS comparison job (job 4)

### Diagnosis

- Both submissions failed in the first step of the old-code pass: job 14909123 (Sep 25, 2 m 51 s, of which ~3 min were the SNP-list extraction) and job 15255781 (Sep 28, 11 s). `wgs/old/log_genotype.err`: genotype-snps' argument parsing (`extractChromosomes`, `ArgParsing.py`) stopped with "The given reference cannot be used because the chromosome names are inconsistent!".
- The BAMs are GRCh37d5 (hs37d5; header `UR:…/GRCh37d5/gr37.fasta`): 86 contigs `1…22, X, Y, MT`, GL contigs, `NC_007605`, `hs37d5`. The settings had `REF=$HG19`, the UCSC hg19 FASTA, whose contigs are `chr1…`. This is not a hatchet problem; old and new code check the same.
- Two more problems the job would have hit next:
  - `SNPLIST` pointed to `snp-lists/hg19_chr/00-All.vcf.gz`, with `chr` names.
  - The script reused any existing `$WORK/snplist.vcf.gz`. The one there was extracted on Sep 25 with `chr20 chr21 chr22` (`tabix -l`), so with `CHROMS="20 21 22"` genotype-snps would have had no targets. Deleted.
- Replacements that exist here:
  - `reference/GRCh37-lite.fa` (`.fai`, `.dict`): same names; the MD5s of 20, 21, 22 and X in its `.dict` equal the M5 tags in the T04 BAM header, i.e. identical sequences. No hs37d5 FASTA was found under the shared reference directories.
  - `snp-lists/hg19_nochr/00-All.vcf.gz`: contigs `1…22 X Y MT`, with a `.tbi`.
- Size, for the estimates: normal 61 GB, tumors 145 GB and 144 GB; chromosomes 20-22 carry 4.3% of the normal's 1.1 billion mapped reads.

### Changes to `wgs_compare.sh`

- Settings: `REF` → GRCh37-lite, `SNPLIST` → `hg19_nochr`; the user's BAMs, names, `CHROMS="20 21 22"`, `CHRNOT=` and `componc_cpu` partition are kept.
- **Pre-flight** (0.3 s), before anything else; `PREFLIGHT_ONLY=1` runs only this. It prints the host, the samtools/bcftools hatchet will use (honoring `HATCHET_PATHS_*`) and the settings, then checks, reporting every problem:
  - BAMs present and indexed; `NAMES` count;
  - each chromosome of `CHROMS` in every BAM header, in `REF.fai` and in the SNP list;
  - every standard chromosome of the normal BAM in `REF.fai` (the check that failed);
  - M5 of each BAM header vs `REF`'s `.dict`, where both exist;
  - `CHRNOT` consistent with the names; the panel present.
- **SNP-list extract** named `snplist.<md5 of SNPLIST and CHROMS>.vcf.gz`, written to a temporary name and renamed with its `.tbi` last, then checked to hold exactly the `CHROMS` contigs. It is made only when the genotype step runs.
- **Resumable steps**: `VERSIONS` (default `old new`) and `STEPS` (default `genotype alleles reads phase`). Each step checks that its inputs exist, replaces only its own outputs, and on failure prints the last 30 lines of its stderr and stops the job. `compare.sh` runs only when both versions finished all four steps.

### Verification (light, here)

- `bash -n`; pre-flight on the real settings: **OK**.
- Pre-flight on a copy with the failing settings (`REF=$HG19`, `hg19_chr` list, `CHRNOT=-N`): **FAILED**, listing chromosomes 20-22 missing from the reference and the SNP list, all 24 standard chromosomes of the normal missing from the reference, and the `CHRNOT` mismatch. That is the error of job 15255781, in seconds.
- Pre-flight with a fake `.dict` whose chromosome 20 MD5 was altered: FAILED on T04's M5 (the other two BAM headers carry no M5).
- On the small `data/` BAMs (a scratch copy of the script: chr22, a 276k-site SNP list):
  - alleles without genotype stops with "needs …/snps/chr22.vcf.gz: run the genotype step";
  - genotype runs (31 s), creating `snplist.f695ad1d.vcf.gz`, and a second run reuses it;
  - a chromosome missing from the SNP list fails the pre-flight;
  - a failing step (`J=0`) prints genotype-snps' error and exits 1.

### Job for the user

```bash
cd /data1/shahs3/users/myersm2/repos/hatchet/ai-logs/2026-09-24-preprocessing-efficiency
PREFLIGHT_ONLY=1 bash wgs_compare.sh      # seconds; the last line must be "pre-flight OK"
sbatch wgs_compare.sh                     # componc_cpu, 8 cores, 64 GB, 2 days
```

- **Expected:** a few minutes of SNP-list extraction (20-22 of dbSNP), then the old pass (~3-5 h: each old genotype and count-alleles task streams its whole BAM, 61 / 145 / 144 GB), then the new pass (~1-1.5 h, mostly mosdepth, the same in both), then `compare.sh`. Peak memory ~15 GB (mosdepth ~4 GB per sample, 3 at a time).
- **Verify:** `grep -E 'FAILED|DIFFERENT|MISSING' hatchet-preproc-compare-<jobid>.log` prints nothing; the log ends with the time/memory table.
- **If a step fails:** the log shows `FAILED: <version> <step>` and that step's stderr tail. After fixing, resume without redoing finished steps, e.g. `sbatch --export=ALL,VERSIONS=new,STEPS="reads phase" wgs_compare.sh`.
- **Follow-up:** paste the table (and any failure) back to be recorded here.
