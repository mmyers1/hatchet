# phase_snps_standalone.py

A standalone script for phasing germline SNPs using a reference panel (SHAPEIT + 1000 Genomes Phase 3). This is a self-contained version of the `hatchet phase-snps` step with no dependencies on the hatchet Python package.

## Dependencies

### Python

- **Python 3.7+** (uses f-strings, `os.makedirs(..., exist_ok=True)`, multiprocessing)
- No third-party Python packages required — only the standard library is used.

### External Tools

The following command-line tools must be installed and accessible (either on `$PATH` or with paths provided via arguments):

| Tool | Purpose | Install |
|------|---------|---------|
| **bcftools** | VCF filtering, annotation, indexing, concatenation | `conda install -c bioconda bcftools` |
| **bgzip** | Block-gzip compression of VCF files | `conda install -c bioconda htslib` |
| **shapeit** (v2) | Haplotype phasing using reference panel | `conda install -c dranew shapeit` |
| **picard** | Liftover of VCF coordinates (only needed for hg38) | `conda install -c bioconda picard` |

### Reference Data

- **1000 Genomes Phase 3 reference panel** — must be downloaded before running. The panel directory should contain a `1000GP_Phase3/` subdirectory with `.hap.gz`, `.legend.gz`, `.sample`, and genetic map files.
- **Reference genome FASTA** (matching BAM alignments) — e.g., `hg19.fa` or `hg38.fa`.
- If using **hg38**: liftover chain files (`hg38ToHg19.*.chain`, `hg19ToHg38.*.chain`) and the hg19 reference (`hg19_no_chr.fa`) must be present in the reference panel directory.
- If using **hg19 with chr notation**: chromosome renaming files (`rename_chrs1.txt`, `rename_chrs2.txt`) must be present in the reference panel directory.

## Usage

```bash
python phase_snps_standalone.py \
    -D /path/to/refpanel_dir \
    -g /path/to/reference.fa \
    -V hg19 \
    -o /path/to/output_dir \
    -L chr1.vcf.gz chr2.vcf.gz chr3.vcf.gz ... \
    [-j 4] \
    [-N] \
    [-si /path/to/shapeit/dir] \
    [-pc /path/to/picard/dir] \
    [-bt /path/to/bcftools/dir] \
    [-bg /path/to/bgzip/dir] \
    [--picard-java-flags "-Xmx8g"]
```

### Arguments

| Flag | Long | Required | Description |
|------|------|----------|-------------|
| `-D` | `--refpaneldir` | Yes | Path to directory containing 1000GP reference panel |
| `-g` | `--refgenome` | Yes | Path to reference genome FASTA |
| `-V` | `--refversion` | Yes | Reference genome version (`hg19` or `hg38`) |
| `-o` | `--outdir` | Yes | Output directory for phased VCF files |
| `-L` | `--snps` | Yes | Per-chromosome VCF files of SNPs to phase (space-separated) |
| `-j` | `--processes` | No | Number of parallel processes (default: 2) |
| `-N` | `--chrnotation` | No | Flag: chromosomes use "chr" prefix notation |
| `-si` | `--shapeit` | No | Directory containing `shapeit` executable |
| `-pc` | `--picard` | No | Directory containing `picard` or `picard.jar` |
| `-bt` | `--bcftools` | No | Directory containing `bcftools` executable |
| `-bg` | `--bgzip` | No | Directory containing `bgzip` executable |
| | `--picard-java-flags` | No | Java flags for picard (default: `-Xmx8g`) |

### Input SNP Files

Input VCF files (`-L`) should be named by chromosome (e.g., `1.vcf.gz`, `chr1.vcf.gz`, `22.vcf.gz`). The script derives chromosome identifiers from the filenames by stripping the `.vcf.gz` extension and any `chr` prefix.

## Output

- `<outdir>/phased.vcf.gz` — concatenated phased VCF across all chromosomes
- `<outdir>/phased.log` — per-chromosome summary of phased vs. original SNP counts

## Pipeline Overview

For each chromosome, the script:

1. **Preprocesses** the input VCF:
   - hg38: lifts over coordinates to hg19 (where the reference panel is defined)
   - hg19 with chr: renames chromosomes to remove "chr" prefix
   - hg19 without chr: copies as-is

2. **Filters** to biallelic SNPs only (removes indels and multi-allelic sites)

3. **Phases** with SHAPEIT using the 1000GP Phase 3 reference panel

4. **Postprocesses** the phased VCF:
   - hg38: lifts back to hg38 coordinates
   - hg19 with chr: adds "chr" prefix back
   - hg19 without chr: uses as-is

5. **Concatenates** all per-chromosome phased VCFs into a single output file

Chromosomes are processed in parallel (controlled by `-j`).
