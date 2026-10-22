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
| `--max-homdel-size` | Maximum cluster size with a homozygous deletion | Maximum fraction of total covered cluster weight that may receive tumor copy number `(0, 0)`. Clusters at or above the threshold must retain at least one copy in every tumor clone. | 0.005 |
| `-f`, `--noampdel` | Activate clone evolutionary contraints | User can decide whether to enable or not constrained about the evolution of tumor clones. These constrained force each allele to be either amplified or deleted across all tumor clones | Activated |
| `--ampdel-exempt` | Clusters exempt from the amp-del constraint | Comma-separated cluster IDs (as in the SEG `#ID` column) that may be both amplified and deleted across tumor clones while the amp-del constraint stays active for every other cluster. Cannot be combined with `--noampdel` | None |
| `-d`, `--cnstates` | Maximum number of distinct copy-number states per cluster | When enabled, the maximum number of distinct copy-number states per cluster is fixed. This option is deprecated | Not used |
| `--fixed-matrix` | Fixed factorization entries | Path to a sparse TSV file containing entries of the copy-number matrices `A` and `B` or mixture matrix `U` that must remain fixed during solving | None |

HATCHet implements two methods to solve the constrained and distance-based simultaneous factorization: (1) a integer-linear programming (ILP) and (2) a coordinate-descent method (CD).
These methods can be combined in 3 different modes:
- (0) CD + ILP: the solution found by CD is used to start the ILP. As such, ILP attempts to improve the solution found by CD.
- (1) ILP only
- (2) CD only

The mode can be specified through the flag `--mode` whose default is mode 2.
In addition, the solving methods can be controlled by the following parameters.

`--max-homdel-size` is expressed as a fraction in `[0, 1]`, not as a percentage, and its cutoff is inclusive. For example, the default `0.005` prevents a cluster covering exactly 0.5% or more of the total cluster weight from receiving `(0, 0)` in a tumor clone. This constraint is independent of the amp/del option and is applied by both the Pyomo and native solvers, including coordinate-descent copy-number steps and the Pyomo bin-wise objective.

The equivalent `[compute_cn]` INI setting is `maxhomdelsize = 0.005`. It can also be overridden with the `HATCHET_COMPUTE_CN_MAXHOMDELSIZE` environment variable.

`--ampdel-exempt` relaxes the amp-del constraint for the listed clusters only. Each allele of an exempt cluster may then be above its baseline in some tumor clones and below it in others (the baseline is 1 copy per allele without a WGD and 2 with a WGD). Listing a clonal cluster has no effect, since its copy numbers are already fixed. An ID that is not in the SEG file is an error. The option can also be set as `ampdelexempt` in the `[compute_cn]` section of `hatchet.ini`, and it applies to both the Pyomo and native solvers, including fixed-matrix and sample-clonal validation. Each exempt cluster removes 2 binary variables and `4 × (number of tumor clones)` linear constraints from the factorization, so memory use is essentially unchanged. Coordinate-descent runtime is expected to fall between the runs with and without `--noampdel`. Because exemptions relax the problem, the objective at a given number of clones can only stay the same or improve (up to the heuristic nature of coordinate descent), and this can change which number of clones is selected.

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

# Expert options

**WARNING: These options are not recommended for users who lack considerable experience refining HATCHet2 solutions, as they restrict the factorization and may produce implausible results if used improperly.**

## [EXPERT OPTION] Sample-clonal override

An expert may bypass both automatic clonal-cluster identification and the global `--clonal` override by supplying two or more **sample-clonal** states for every sample. A sample-clonal cluster has one copy-number state in every tumor clone that is present in that sample; the state and cluster may differ between samples.

Only the **first two** states listed for each sample are reference states used for solving the sample's scaling factor and purity (below). Any additional states for that sample are fixed in every tumor clone present in that sample, exactly like the reference states, but are not used to derive scaling or purity and are not checked against the sample's RDR or BAF.

This option is intended for experts who have visually identified reliable clonal clusters. Incorrect reference states impose hard constraints on the factorization and can produce infeasible or misleading results.

For each sample `p`, the two reference states provide cluster RDRs `r1` and `r2` and distinct total copy numbers `w1` and `w2`. HATCHet derives the sample scaling factor and purity from the HATCHet Theorem 1 system:

```text
gamma[p] * r1 = 2 * (1 - purity[p]) + w1 * purity[p]
gamma[p] * r2 = 2 * (1 - purity[p]) + w2 * purity[p]
```

The two total copy numbers must differ and the system must be nonsingular. The resulting scaling factor must be positive and the purity must be in `[0, 1]`. Scaling uses total copy number, while factorization enforcement uses the supplied allele-specific `cn_a` and `cn_b` values.

### Inline form

Use `--sample-clonal` with one semicolon-separated specification per sample:

```bash
hatchet compute-cn \
	-i /path/to/bulk \
	--diploid \
	--sample-clonal 'Tumor1=3:1:1,5:2:1,7:1:0;Tumor2=4:1:1,8:2:0'
```

Here `Tumor1` uses clusters `3` and `5` as reference states and additionally fixes cluster `7` to `(1, 0)` in every tumor clone present in `Tumor1`.

The grammar is:

```text
SAMPLE=CLUSTER:CN_A:CN_B,CLUSTER:CN_A:CN_B[,CLUSTER:CN_A:CN_B...][;SAMPLE=...]
```

The first two `CLUSTER:CN_A:CN_B` states for a sample are its reference states.

The identifiers must exactly match the `SAMPLE` and `#ID` columns of the SEG input. The characters `=`, `;`, `,`, and `:` are reserved delimiters in this form. Use the TSV form when an identifier contains one of these characters.

The equivalent INI setting is:

```ini
[compute_cn]
sampleclonal = Tumor1=3:1:1,5:2:1,7:1:0;Tumor2=4:1:1,8:2:0
diploid = True
```

It can also be supplied through `HATCHET_COMPUTE_CN_SAMPLECLONAL`.

### Versioned TSV form

Use `--sample-clonal-file` for larger or programmatically generated specifications:

```bash
hatchet compute-cn \
	-i /path/to/bulk \
	--tetraploid \
	--sample-clonal-file /path/to/patient.sample-clonal.tsv
```

The UTF-8, tab-separated file is versioned and has four columns:

```text
# HATCHET_SAMPLE_CLONAL 1
sample	cluster	cn_a	cn_b
Tumor1	3	2	2
Tumor1	5	2	1
Tumor2	4	2	2
Tumor2	8	2	0
Tumor1	7	2	0
```

Blank lines and comment lines are allowed after the format marker. Every SEG sample must appear at least twice, the cluster IDs for one sample must be distinct, and all allele copy numbers must be nonnegative integers. Rows for different samples may be interleaved in any order, but within each sample the first two rows in file order are the reference states; later rows for that sample (such as `Tumor1 7` above) are additional states.

The equivalent INI and environment settings are:

```ini
[compute_cn]
sampleclonalfile = /path/to/patient.sample-clonal.tsv
tetraploid = True
```

```text
HATCHET_COMPUTE_CN_SAMPLECLONALFILE=/path/to/patient.sample-clonal.tsv
```

### Constraints and option interactions

- Exactly one of `--diploid` or `--tetraploid` is required. The sample-clonal override does not run automatic WGD model selection.
- `--sample-clonal` and `--sample-clonal-file` are mutually exclusive.
- Sample-clonal input cannot be combined with global `--clonal` assignments.
- `--purities` may be supplied as an independent check. Every value must agree with the purity derived from the sample-clonal references within absolute tolerance `0.05`; the derived values are used by the solver.
- The normal proportion is fixed to `1 - purity[p]` in every sample.
- For each supplied state, reference or additional, every tumor clone with positive proportion in that sample is constrained to the supplied A/B state. Clones absent from the sample are not constrained by that sample.
- Only the two reference states of a sample must have distinct total copy numbers. Additional states have no such requirement and do not affect the derived scaling factor or purity.
- The same cluster may have different states in different samples. In that case, a shared tumor clone cannot be present in both samples unless it satisfies both states.
- Clone-count ranges remain supported because the input does not name tumor clone columns. The solver reports an error when too few tumor clones are available to satisfy incompatible sample states.
- A `--fixed-matrix` file may be used simultaneously. Fixed-positive tumor clones must satisfy the corresponding sample-clonal states, fixed normal proportions must match derived purity, and remaining tumor purity must be assignable to a compatible clone.
- Required states, including additional states, must satisfy the active `--diploidcmax` or `--tetraploidcmax`, amp/del constraints, the deprecated `--cnstates` limit, and the existing nonzero constraint for large clusters.


## [EXPERT OPTION] Fixing entries of A, B, and U

Experts may provide hard constraints for selected entries of the simultaneous factorization while leaving all other entries free. This can be useful when copy-number states are known for particular clones or when external evidence supports mixture proportions for a subset of samples.

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

Each data row has four fields. There are two types of rows, one for A/B entries and the other for U entries.

| `matrix` | `row` | `column` | `value` |
|----------|-------|----------|---------|
|`A` or `B` | Exact cluster ID from the `#ID` column of the SEG file | `normal` or `cloneN` | Nonnegative integer copy number |
| `U` | `normal` or `cloneN` | Exact sample name from the `SAMPLE` column of the SEG file | Finite number in `[0, 1]` |

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
- `--max-homdel-size`: a tumor state of `(0, 0)` is not allowed for a cluster whose fraction of total cluster weight is at or above the configured cutoff.
- Mixture sums: fixed U entries for a sample may sum to at most 1. If the entire column is fixed, it must sum to 1. Any remaining mass must be assignable under `--minprop`.

When neither `--diploid` nor `--tetraploid` is specified, the same file is applied to both model-selection branches. A constraint that is incompatible with either branch causes an explicit error. Force a ploidy mode when the expert constraints are intended for only that branch.

Tumor clones are interchangeable in the model, so HATCHet does not constrain their order while solving. Without a fixed-matrix file, or when the file addresses only the normal clone, the reported tumor clones are sorted after solving by nondecreasing weighted copy number, where `s` is the cluster's position in the SEG file:

```text
sum_s (s + 1) (A[s, i] + B[s, i])
```

When the file addresses a tumor clone (an A or B entry in a `cloneN` column, or a U entry in a `cloneN` row), clone labels keep the meaning given by the file and are not reordered. Local validation cannot detect every conflict: partial entries may still produce a globally infeasible model, and in that case the solver reports that no constrained solution is feasible.
