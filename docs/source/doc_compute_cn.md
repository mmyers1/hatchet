# compute-cn

This step infers allele-specific copy numbers for each tumor clone as well as the mixture proportions of the tumor clones in each sample. It consists of the following components:

1. Compute allele-specific fractional copy numbers.
2. Solve a constrained distance-based simultaneous factorization to compute allele and clone-specific copy numbers and clone proportions.
3. Deploy a model-selection criterion to select the number of clones by explicitly considering the trade-off between subclonal copy-number aberrations and whole-genome duplication.

The step offers some parameters to control each of these features.

## Input

This step requires in input two tab-separated files that must have the same prefix `${PRE}` and the only difference is the extension:

1. The name of the first file must be `${PRE}.bbc` and is a tab-separated file describing the RDR and BAF of clustered genomic bins with the following fields

| Field | Description |
|-------|-------------|
| `CHR` | Name of a chromosome |
| `START` | Starting genomic position of a genomic bin in `CHR` |
| `END` | Ending genomic position of a genomic bin in `CHR` |
| `SAMPLE` | Name of a tumor sample |
| `RD` | RDR of the bin in `SAMPLE` |
| `#SNPS` | Number of SNPs present in the bin in `SAMPLE` |
| `COV` | Average coverage in the bin in `SAMPLE` |
| `ALPHA` | Alpha parameter related to the binomial model of BAF for the bin in `SAMPLE`, typically total number of reads from A allele |
| `BETA` | Beta parameter related to the binomial model of BAF for the bin in `SAMPLE`, typically total number of reads from B allele |
| `BAF` | BAF of the bin in `SAMPLE` |
| `CLUSTER` | The name of a cluster |

2. The name of the second file must be `${PRE}.seg` and is a tab-separated file describing the RDR and BAF of clusters of genomic regions with the following fields

| Field | Description |
|-------|-------------|
| `ID` | The name of a cluster |
| `SAMPLE` | The name of a sample |
| `#BINS` | The number of bins included in `ID` |
| `RD` | The RDR of the cluster `ID` in `SAMPLE` |
| `#SNPS` | The total number of SNPs in the cluster `ID` |
| `COV` | The average coverage in the cluster `ID` |
| `ALPHA` | The alpha parameter of the binomial model for the BAF of the cluster `ID` |
| `BETA` | The beta parameter of the binomial model for the BAF of the cluster `ID` |
| `BAF` | The BAF of the cluster `ID` in `SAMPLE` |

## Output

This step produces 2 files, a BBC result file and a SEG result file, for every considered number of clones and for both the assumption of the presence and absence of a WGD.
For every considered number of clones `${N}` and assuming `${P}` is either equal to `diploid` or `tetraploid`, the steps generates a pair of files as follows:

1. The first BBC file has the following name `results.${P}.n${N}.bbc.ucn.tsv` and is a copy of the BBC input file with `2 (${N} + 1)` additional fields according to the following:

| Field | Description |
|-------|-------------|
| ... | ... |
| `cn_normal` | The copy number state of the normal diploid clone equal to <code>1&#124;1</code> |
| `u_normal` | The normal admixture of the normal diploid cells in the corresponding sample |
| `cn_clone${n}` | The copy number state of the `${n}` tumor clone in the format <code>A&#124;B</code> where `A` and `B` are the two allele-specific copy numbers of the corresponding genomic bin |
| `u_clone${n}` | The clone proportion of the `${n}` tumor clone in the corresponding sample |

2. The second SEG file has the following name `results.${P}.n${N}.seg.ucn.tsv` and is a file of genomic segments with the following fields:

| Field | Description |
|-------|-------------|
| `CHR` | The name of a chromosome |
| `START` | The genomic position that starts the corresponding genomic segment |
| `END` | The genomic position that ends the corresponding genomic segment |
| `SAMPLE` | The name of a sample |
| `cn_normal` | The copy number state of the normal diploid clone equal to <code>1&#124;1</code> |
| `u_normal` | The normal admixture of the normal diploid cells in the corresponding sample |
| `cn_clone${n}` | The copy number state of the `${n}` tumor clone in the format <code>A&#124;B</code> where `A` and `B` are the two allele-specific copy numbers of the corresponding genomic bin |
| `u_clone${n}` | The clone proportion of the `${n}` tumor clone in the corresponding sample |

In addition, the two files for the best solution assuming both the absence or presence of a WGD are copied in the files `chosen.diploid.bbc.ucn`, `chosen.diploid.seg.ucn` and `chosen.tetraploid.bbc.ucn`, `chosen.tetraploid.seg.ucn`, respectively.
Last, hatchet copies the best solution according to the prediction of WGD in the two files `best.bbc.ucn` and `best.seg.ucn`.

## Fractional copy numbers

The first feature of hatchet is the explicit estimation of allele-specific fractional copy numbers.
hatchet performs this estimation based on rigorous criterion stating that allele-specific fractional copy numbers can be obtained from RDR and BAF of clusters of genomic regions.
More specifically, the criterion states that:
- when there is no WGD, the identification of the diploid cluster with copy-number state (1, 1) is sufficient
- when there is a WGD, the identification of the tetraploid cluster with copy-number state (2, 2) and of a second clonal cluster with a copy-number state (A, B) in all tumor clones are sufficient.

As such, hatchet identifies a heuristic to identifies the required clusters and copy-number states.
This heuristic can be controlled by the following parameters:

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-c`, `--clonal` | The required clusters and corresponding copy number states | User can directly specifies the required clusters and the corresponding copy number states to compute the allele-specific fractional copy number states. These must be speficied in the format `IDX-1:A-1:B-1, ..., IDX-M:A-M:B-M` where `IDX-S` is the name of cluster `S` and `(A-S, B-S)` is the corresponding copy-number state. Moreover, user can specify addittional clusters and copy numbers beyond the required ones. The copy-numbers for these clusters will be fixed during the computation. This can be an usefule feature for especially noisy datasets. | None |
| `-ts`, `--minsize` | Threshold for size of clusters | The minimum size of the clusters to consider for the heuristic that identifies clonal clusters. The non-selected clusters will not be considered as potential tumor-clonal clusters. The threshold must be expressed as a fraction of the entire genome. | 0.02, e.g. `2%` of genome |
| `-tc`, `--minchrs` | Threshold for number of chromosomes | The minimum number of chromosomes covered by the clusters to consider for the heuristic that identifies clonal clusters. The non-selected clusters will not be considered as potential tumor-clonal clusters. | 1 |
| `-td`, `--maxneutralshift` | Maximum BAF shift allowed for diploid cluster | The maximum expected shift from 0.5 for BAF for a diploid or tetraploid cluster (i.e. with copy-number states (1, 1) or (2, 2)). This threshold is used for two goals to identify the diploid or tetraploid cluster. | 0.1 |
| `-tR`, `--toleranceRDR` | Maximum RDR tolerance | The maximum RDR tolerance used by the heuristic when estimating the position of all clonal clusters | 0.04 |
| `-tB`, `--toleranceBAF` | Maximum BAF tolerance | The maximum BAF tolerance used by the heuristic when estimating the position of all clonal clusters | 0.03 |
| `--merge` | Activate merging of clusters | When activated, the heuristic will merge together clusters that appear to have the same values of RDR and BAF, according to the values below. This procedure can help the heuristic by refining the clustering and merging clusters that are likely to have the same copy-number states and unlikely to be clonal. | False, not used |
| `-mR`, `--mergeRDR` | RDR merging threhsold | The maximum difference in RDR considered by the merging procedure | 0 |
| `-mB`, `--mergeBAF` | BAF merging threhsold | The maximum difference in BAF considered by the merging procedure | 0 |

## Simultaneous factorization

The second feature corresponds to the simultaneous factorization of allele-specific copt numbers and clone proportions.
hatchet solves a constrained and distance-based variant of the factorization that can be controlled by the following parameters.

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-n`, `--clones` | Interval for number of clones | User can specify the minimum and maximum number of clones to consider | `"2,8"` |
| `-eD`, `--diploidcmax` | Maximum copy number with no WGD | The value of maximum copy number that is considered when assuming no WGD. When `0`-value is specified the maximum copy number is directly inferred from the data by rounding the maximum fractional copy number | 8 |
| `-eT`, `--tetraploidcmax` | Maximum copy number with a WGD | The value of maximum copy number that is considered when assuming there is a WGD. When `0`-value is specified the maximum copy number is directly inferred from the data by rounding the maximum fractional copy number | 8 |
| `-u`, `--minprop` | Minimum clone proportion | In every sample, each clone either is non present (clone proportion equal to 0.0) or has a clone proportion higher than this threshold | 0.03 |
| `-f`, `--noampdel` | Activate clone evolutionary contraints | User can decide whether to enable or not constrained about the evolution of tumor clones. These constrained force each allele to be either amplified or deleted across all tumor clones | Activated |
| `-d`, `--cnstates` | Maximum number of distinct copy-number states per cluster | When enabled, the maximum number of distinct copy-number states per cluster is fixed. This option is deprecated | Not used |
| `--fixed-matrix` | Fixed factorization entries | Path to a sparse TSV file containing entries of the copy-number matrices `A` and `B` or mixture matrix `U` that must remain fixed during solving | None |

HATCHet implements two methods to solve the constrained and distance-based simultaneous factorization: (1) a integer-linear programming (ILP) and (2) a coordinate-descent method (CD).
These methods can be combined in 3 different modes:
- (0) CD + ILP: the solution found by CD is used to start the ILP. As such, ILP attempts to improve the solution found by CD.
- (1) ILP only
- (2) CD only

The mode can be specified through the flag `--mode` whose default is mode 2.
In addition, the solving methods can be controlled by the following parameters.

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-j`, `--jobs` | Number of parallel jobs | Every job will run on parallel a different seed of the CD algorithm and all the jobs can be used for solving the ILP | Number of available processes on machine |
| `-p`, `--seeds` | Number of CD seeds | Number of random restarts executed by the CD algorithm. In general, the higher the number of restarts, the better the result | 400 |
| `-r`, `--randomseed` | Seed for random generator | Every execution started from the same seed have a deterministic behaviour and can be exactly replicated | None, random seed |
| `-s`, `--timelimit` | Time limit | The time limit, expressed in seconds, is imposed to every step of the CD algorithm or to the whole ILP | None |
| `-m`, `--timelimit` | Memory limit | The memory limit, expressed in megabytes, is imposed to every step of the CD algorithm or to the whole ILP. The execution will not be interrupted when reaching the threshold but disk is used | None |
| `--maxiterations` | Maximum number of iteration per seed | This number is imposed as the maximum number of iterations executed for every restart of the CD algorithm | 40 |
| `--diploid` | Assume no WGD | When enabled, HATCHet assumes the absence of a WGD | Not used |
| `--tetraloid` | Assume a WGD | When enabled, HATCHet assumes the occurrence of a WGD | Not used |

## Model selection

This steps have two main parameters to control the model-selection criterion:

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-l`, `--limitinc` | Sensitivity level | The sensitivity level is used to control the confidence in evaluating the presence of tumor clones characterized by small CNAs. By decreasing the value of the sensitivity, HATCHet is more sensible to the presence of small CNAs and small clusters or with small shifts in RDR/BAF are more likely considered as the signal of an additional tumor clone. The possible values of this parameter are between 1.0 and 0.0 and specifically corresponds to an upper bound for the left relative improvement of the objective function. | None, reasonable values to use can be 0.6, 0.5, 0.4, 0.3, 0.2, ... according to the values of the objective function. |
| `-g`, `--ghostprop` | Confidence in the presence of a single tumor clone | This value expresses the confidence of HATCHet when evaluating the presence of a single tumor clone. The higher the value the more likely the presence of a single clone is considered | 0.2 |

## [EXPERT FEATURE] Fixing entries of A, B, and U

Experts may provide hard constraints for selected entries of the simultaneous factorization while leaving all other entries free. This can be useful when copy-number states are known for particular clones or when external evidence supports mixture proportions for a subset of samples. **WARNING: This feature is not recommended for users who lack considerable experience refining HATCHet2 solutions, as forcing specific values can easily produce poorly supported results.**

For `m` segmental clusters, `n` total clones, and `k` bulk samples, HATCHet uses:

- `A`, an `m` by `n` matrix of integer A-allele copy numbers;
- `B`, an `m` by `n` matrix of integer B-allele copy numbers; and
- `U`, an `n` by `k` matrix of clone proportions.

The first clone is the normal clone. It is addressed as `normal`; tumor clones are addressed as `clone1`, `clone2`, and so on. The normal copy numbers are always `A = 1` and `B = 1`. Every sample column of `U` sums to 1.

The observed allele-specific fractional copy numbers are approximated by the products `A U` and `B U`. Fixed entries are hard equalities in this factorization. They do not require an exact fit to the observed fractional copy numbers because the HATCHet objective permits residual error.

### Creating a file in the offline editor

HATCHet includes a self-contained graphical editor that runs locally in a web browser and makes no network requests:

```bash
hatchet fixed-matrix-editor
```

On a remote or headless system, print the packaged HTML path without trying to open a browser:

```bash
hatchet fixed-matrix-editor --no-open
```

In the editor:

1. Open the `${PRE}.seg` file used by `compute-cn`.
2. Set the exact total number of clones, including the normal clone.
3. Enter fixed cells in the A allele, B allele, or U mixture tabs. Blank cells remain free.
4. Review the validation status and download the TSV file.

The editor supports filtering, keyboard entry, rectangular spreadsheet paste, clearing one matrix or all matrices, and importing an existing fixed-entry TSV. Normal A/B cells are displayed as locked values of 1 and do not need to be written to the downloaded file.

The editor's maximum-copy-number and minimum-proportion fields provide immediate checks while editing. The corresponding `compute-cn` options remain authoritative when the file is used for solving.

### File format

The auxiliary file is UTF-8, tab-separated, sparse, and versioned. A complete example for three total clones is:

```text
# HATCHET_FIXED_ENTRIES 1
# clones 3
matrix	row	column	value
A	1	clone1	2
B	1	clone1	1
A	4	clone2	3
B	4	clone2	1
U	normal	Tumor1	0.35
U	clone1	Tumor1	0.40
```

The first line must be the format marker shown above. The clone-count metadata must appear before the header and must equal the exact `compute-cn` clone count. `n` includes the normal clone, so `# clones 3` defines `normal`, `clone1`, and `clone2`.

Each data row has four fields:

| Field | A or B entry | U entry |
|-------|--------------|---------|
| `matrix` | `A` or `B` | `U` |
| `row` | Exact cluster ID from the `#ID` column of the SEG file | `normal` or `cloneN` |
| `column` | `normal` or `cloneN` | Exact sample name from the `SAMPLE` column of the SEG file |
| `value` | Nonnegative integer copy number | Finite number in `[0, 1]` |

Only listed entries are fixed. Data rows may appear in any order, blank lines and comment lines are allowed, and identifiers are matched exactly. Duplicate matrix coordinates are rejected rather than overwritten. A file must contain at least one data row.

### Running compute-cn

Specify the auxiliary file on the command line and use one exact total clone count:

```bash
hatchet compute-cn \
	-i /path/to/bulk \
	-n 3 \
	--fixed-matrix /path/to/patient.fixed-matrix.tsv
```

The equivalent `.ini` interface is:

```ini
[compute_cn]
clones = 3
fixedmatrix = /path/to/patient.fixed-matrix.tsv
```

The command-line option overrides the `.ini` value. The setting can also be overridden through `HATCHET_COMPUTE_CN_FIXEDMATRIX` environment variable.

Clone ranges such as `clones = 2,6` are rejected when a fixed-matrix file is present because `clone1`, `clone2`, and subsequent labels must have stable meanings. Output filenames and BBC/SEG columns are unchanged.

### Interactions with other constraints

HATCHet validates the auxiliary file before optimization and reports the source line for errors when possible. Important interactions include:

- `--clonal`: a clonal copy-number assignment fixes that cluster in every tumor clone. A conflicting A or B entry is rejected.
- `--purities`: in the Pyomo solver path, fixing purity sets `U[normal, sample] = 1 - purity`. A conflicting normal proportion is rejected.
- `--minprop`: a fixed tumor proportion must be 0 or at least the configured minimum, within numerical tolerance.
- `--diploidcmax` and `--tetraploidcmax`: fixed copy numbers and fully specified A+B totals must satisfy the active branch's bounds.
- Amp/del constraints: when enabled, fixed values cannot force the same allele to be both amplified and deleted across tumor clones for one cluster.
- `--cnstates`: fully specified states count toward the configured limit, although this option remains deprecated.
- Large clusters: a fully fixed tumor state of `(0, 0)` is not allowed where the existing large-cluster constraint requires at least one copy.
- Mixture sums: fixed U entries for a sample may sum to at most 1. If the entire column is fixed, it must sum to 1. Any remaining mass must be assignable under `--minprop`.

When neither `--diploid` nor `--tetraploid` is specified, the same file is applied to both model-selection branches. A constraint that is incompatible with either branch causes an explicit error. Force a ploidy mode when the expert constraints are intended for only that branch.

Tumor clone labels follow HATCHet's canonical weighted-copy-number order. For adjacent tumor clones, HATCHet enforces:

```text
sum_s (s + 1) (A[s, i] + B[s, i]) <= sum_s (s + 1) (A[s, i + 1] + B[s, i + 1])
```

Fully specified entries that violate this ordering are rejected when the conflict can be determined before solving. Partial entries may still produce a globally infeasible model; in that case the solver reports that no constrained solution is feasible.


## Additional parameters
