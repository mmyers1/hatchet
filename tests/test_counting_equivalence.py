"""
Equivalence tests for the vectorized counting code in count_reads and count_alleles.

Each new function is compared against a verbatim copy of the loop it replaced (the `_old_*` functions below), on
randomized inputs that include the edge cases of the original code. These tests do not need any BAM files.
"""
import gzip
import os
import shutil
import subprocess

import numpy as np
import pytest
from scipy.stats import beta

from hatchet.utils import count_reads, count_alleles


# ---------------------------------------------------------------------------------------------------------------------
# count_reads: reference implementations (loop bodies copied from the original form_counts_array)
# ---------------------------------------------------------------------------------------------------------------------


def _old_bin_read_starts(read_starts, thresholds):
    arr = np.zeros(len(thresholds) + 1)
    end = thresholds[0]
    j = 0
    for idx in range(1, len(thresholds)):
        start = end
        end = thresholds[idx]

        my_count = 0
        while j < len(read_starts) and read_starts[j] < end:
            assert read_starts[j] >= start
            my_count += 1
            j += 1

        arr[idx - 1] = my_count

    arr[idx - 1] += len(read_starts) - j
    return arr[:-2]


def _old_coverage(lines, thresholds):
    coverage = np.zeros(len(thresholds))
    idx = 0
    last_record = None
    for line in lines:
        tokens = line.split()
        if len(tokens) == 4:
            start = int(tokens[1])
            end = int(tokens[2])
            nreads = int(tokens[3])

            while idx < len(thresholds) and thresholds[idx] - 1 < end:
                assert thresholds[idx] - 1 >= start
                coverage[idx] = nreads
                idx += 1
            last_record = line
    return coverage, idx, last_record


def _old_form_counts_array(starts_files, perpos_files, thresholds, chromosome, tabix, chunksize=1e5):
    arr = np.zeros((thresholds.shape[0] + 1, len(starts_files) * 2))   # add one for the end of the chromosome

    for i in range(len(starts_files)):
        # populate starts in even entries
        fname = starts_files[i]
        if fname.endswith('.gz'):
            f = gzip.open(fname)
        else:
            f = open(fname)

        read_starts = [int(a) for a in gzip.open(fname, 'r')]
        end = thresholds[0]
        j = 0
        for idx in range(1, len(thresholds)):
            start = end
            end = thresholds[idx]

            my_count = 0
            while j < len(read_starts) and read_starts[j] < end:
                assert read_starts[j] >= start
                my_count += 1
                j += 1

            arr[idx - 1, i * 2] = my_count

        arr[idx - 1, i * 2] += len(read_starts) - j
        f.close()

    for i in range(len(perpos_files)):
        # populate threshold coverage in odd entries
        fname = perpos_files[i]

        chr_sample_file = os.path.join(fname[:-3] + '.' + chromosome)

        if not os.path.exists(chr_sample_file):
            with open(chr_sample_file, 'w') as f:
                subprocess.run([tabix, fname, chromosome], stdout=f)

        with open(chr_sample_file, 'r') as records:
            idx = 0
            last_record = None
            for line in records:
                tokens = line.split()
                if len(tokens) == 4:
                    start = int(tokens[1])
                    end = int(tokens[2])
                    nreads = int(tokens[3])

                    while idx < len(thresholds) and thresholds[idx] - 1 < end:
                        assert thresholds[idx] - 1 >= start
                        arr[idx, (2 * i) + 1] = nreads
                        idx += 1
                    last_record = line

            if i == 0:
                # identify the (effective) chromosome end as the last well-formed record
                assert idx == len(thresholds)
                _, _, chr_end, end_reads = last_record.split()
                chr_end = int(chr_end)
                end_reads = int(end_reads)

                assert chr_end > thresholds[-1]
                assert len(thresholds) == len(arr) - 1

                # add the chromosome end to thresholds
                thresholds = np.concatenate([thresholds, [chr_end]])

                # count the number of reads covering the chromosome end
                arr[idx, 0] = end_reads

        if os.path.exists(chr_sample_file):
            os.remove(chr_sample_file)

    return arr, thresholds


# ---------------------------------------------------------------------------------------------------------------------
# count_reads: random inputs
# ---------------------------------------------------------------------------------------------------------------------


def _random_thresholds(rng, n, length):
    # non-decreasing, starting at 1, with some duplicates (as happen around the centromere)
    t = np.sort(rng.integers(2, length, size=n - 1))
    if n > 3:
        dup = rng.integers(1, n - 1, size=max(1, n // 20))
        t[dup] = t[dup - 1]
    return np.concatenate([[1], t]).astype(np.int64)


def _chunks(a, size):
    return [a[i : i + size] for i in range(0, len(a), size)]


def _outcome(fn, *args):
    try:
        return fn(*args)
    except AssertionError:
        return 'AssertionError'


def _random_perpos(rng, length, kind='contiguous'):
    """
    Per-base coverage records: 'contiguous' like mosdepth (covering [0, length)), 'overlapping' (sorted starts but
    non-monotonic ends, still covering every position), or 'arbitrary' intervals.
    """
    if kind == 'contiguous':
        breaks = np.unique(np.concatenate([[0, length], rng.integers(1, length, size=rng.integers(1, 200))]))
        starts, ends = breaks[:-1], breaks[1:]
    elif kind == 'overlapping':
        starts = np.unique(np.concatenate([[0], rng.integers(1, length, size=rng.integers(1, 200))]))
        gaps = np.diff(np.concatenate([starts, [length]]))
        overhang = rng.integers(0, 3, size=len(starts)) * rng.integers(0, length // 10 + 1, size=len(starts))
        ends = starts + gaps + overhang
    else:
        starts = rng.integers(0, length, size=rng.integers(1, 100))
        ends = starts + rng.integers(1, length // 4, size=len(starts))
    depths = rng.integers(0, 50, size=len(starts))
    return starts.astype(np.int64), ends.astype(np.int64), depths.astype(np.int64)


@pytest.mark.parametrize('seed', range(80))
def test_bin_read_starts(seed):
    rng = np.random.default_rng(seed)
    length = int(rng.integers(50, 5000))
    thresholds = _random_thresholds(rng, int(rng.integers(2, 60)), length)
    n_reads = int(rng.integers(0, 400))
    starts = np.sort(rng.integers(1, length + 200, size=n_reads))
    case = seed % 4
    if case == 1 and n_reads > 2:
        # unsorted reads: the original asserts unless the disorder is confined to the tail past thresholds[-1]
        i, j = rng.integers(0, n_reads, size=2)
        starts[i], starts[j] = starts[j], starts[i]
    elif case == 2 and n_reads > 0:
        # a read before the first threshold
        starts[0] = 0
    elif case == 3 and n_reads > 5:
        # disorder after the first read past the last threshold is not checked by the original
        starts[-3:] = starts[-3:][::-1]

    expected = _outcome(_old_bin_read_starts, starts.tolist(), thresholds)
    for size in (1, 3, 1000):
        got = _outcome(count_reads.bin_read_starts, _chunks(starts, size), thresholds)
        if isinstance(expected, str):
            assert got == expected
        else:
            assert np.array_equal(got, expected)


@pytest.mark.parametrize('seed', range(80))
def test_coverage_at_thresholds(seed):
    rng = np.random.default_rng(seed)
    length = int(rng.integers(50, 5000))
    thresholds = _random_thresholds(rng, int(rng.integers(2, 60)), length)
    if seed % 4 == 3:
        # a threshold past the last record (as for the chromosome end of samples after the first)
        thresholds = np.concatenate([thresholds, [length + 10]])
    kind = ['contiguous', 'overlapping', 'arbitrary', 'contiguous'][seed % 4]
    starts, ends, depths = _random_perpos(rng, length, kind=kind)
    lines = [f'chr1\t{s}\t{e}\t{d}\n' for s, e, d in zip(starts, ends, depths)]

    expected = _outcome(_old_coverage, lines, thresholds)
    for size in (1, 3, 1000):
        records = list(zip(_chunks(starts, size), _chunks(ends, size), _chunks(depths, size)))
        got = _outcome(count_reads.coverage_at_thresholds, records, thresholds)
        if isinstance(expected, str):
            assert got == expected
        else:
            cov, idx, last = got
            assert np.array_equal(cov, expected[0])
            assert idx == expected[1]
            _, _, chr_end, end_reads = expected[2].split()
            assert last == (int(chr_end), int(end_reads))


@pytest.mark.skipif(
    shutil.which('tabix') is None or shutil.which('bgzip') is None,
    reason='tabix and bgzip are needed to index per-base coverage files',
)
@pytest.mark.parametrize('seed', range(6))
def test_form_counts_array(seed, tmp_path):
    rng = np.random.default_rng(seed)
    length = int(rng.integers(1000, 20000))
    thresholds = _random_thresholds(rng, int(rng.integers(4, 100)), length - 1)
    starts_files, perpos_files = [], []
    for sample in range(3):
        n_reads = 0 if (seed % 5 == 0 and sample == 2) else int(rng.integers(0, 3000))
        starts = np.sort(rng.integers(1, length, size=n_reads))
        fname = str(tmp_path / f's{sample}.chr1.starts.gz')
        with gzip.open(fname, 'wt') as f:
            f.write(''.join(f'{s}\n' for s in starts))
        starts_files.append(fname)

        bed = tmp_path / f's{sample}.per-base.bed'
        with open(bed, 'w') as f:
            # another chromosome before and after, to check that only chr1 is used
            f.write(f'chr0\t0\t{length}\t7\n')
            s, e, d = _random_perpos(rng, length)
            f.write(''.join(f'chr1\t{a}\t{b}\t{c}\n' for a, b, c in zip(s, e, d)))
            f.write(f'chr2\t0\t{length}\t9\n')
        subprocess.run(['bgzip', '-f', str(bed)], check=True)
        subprocess.run(['tabix', '-f', '-p', 'bed', str(bed) + '.gz'], check=True)
        perpos_files.append(str(bed) + '.gz')

    arr_old, thr_old = _old_form_counts_array(starts_files, perpos_files, thresholds, 'chr1', 'tabix')
    for chunksize in (7, 1e6):
        arr_new, thr_new = count_reads.form_counts_array(
            starts_files, perpos_files, thresholds, 'chr1', 'tabix', chunksize=chunksize
        )
        assert np.array_equal(arr_new, arr_old)
        assert np.array_equal(thr_new, thr_old)
        assert thr_new.dtype == thr_old.dtype
    # no temporary per-chromosome files are left behind
    assert sorted(os.listdir(tmp_path)) == sorted(
        [os.path.basename(f) for f in starts_files]
        + [os.path.basename(f) for f in perpos_files]
        + [os.path.basename(f) + '.tbi' for f in perpos_files]
    )


# ---------------------------------------------------------------------------------------------------------------------
# count_reads.count_chromosome: skipping existing outputs, and not leaving partial ones
# ---------------------------------------------------------------------------------------------------------------------


def test_count_chromosome_skips_existing_output(tmp_path):
    starts = tmp_path / 'S.chr1.starts.gz'
    with gzip.open(starts, 'wt') as f:
        f.write('5\n')
    # 'false' as samtools would fail if it were run
    count_reads.count_chromosome('chr1', str(tmp_path), 'false', 'unused.bam', 'S', 0)
    with gzip.open(starts, 'rt') as f:
        assert f.read() == '5\n'


def test_count_chromosome_failure_leaves_no_output(tmp_path):
    with pytest.raises(ValueError, match='samtools'):
        count_reads.count_chromosome('chr1', str(tmp_path), 'false', 'unused.bam', 'S', 0)
    assert os.listdir(tmp_path) == []


@pytest.mark.skipif(shutil.which('samtools') is None, reason='samtools is needed to read the BAM')
def test_count_chromosome_writes_read_starts(tmp_path):
    import pysam

    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': 'chr1', 'LN': 1000}, {'SN': 'chr2', 'LN': 1000}]}
    bam = str(tmp_path / 'tiny.bam')
    reads = [('chr1', 10, 30), ('chr1', 10, 5), ('chr1', 200, 60), ('chr2', 50, 60)]   # (contig, start, MAPQ)
    with pysam.AlignmentFile(bam, 'wb', header=header) as f:
        for i, (contig, start, mapq) in enumerate(reads):
            read = pysam.AlignedSegment(f.header)
            read.query_name = f'r{i}'
            read.reference_name = contig
            read.reference_start = start - 1
            read.mapping_quality = mapq
            read.cigarstring = '20M'
            read.query_sequence = 'A' * 20
            f.write(read)
    pysam.index(bam)

    count_reads.count_chromosome('chr1', str(tmp_path), 'samtools', bam, 'S', 20)
    with gzip.open(tmp_path / 'S.chr1.starts.gz', 'rt') as f:
        assert f.read() == '10\n200\n'   # MAPQ >= 20 on chr1 only
    assert not os.path.exists(tmp_path / 'S.chr1.starts.gz.tmp')


# ---------------------------------------------------------------------------------------------------------------------
# count_alleles: reference implementations (copied from the original count_alleles)
# ---------------------------------------------------------------------------------------------------------------------


def _old_selectHetSNPs(counts, gamma, maxshift):
    hetSNPs = {
        c: [(r[0], r[1], r[2], r[3][0], max(r[3][1:], key=(lambda x: x[1]))) for r in counts[c] if len(r[3]) >= 2]
        for c in counts
    }
    check = lambda r: _old_isHet(r[3][1], r[4][1], gamma) and _old_checkShift(r[3][1], r[4][1], maxshift)
    hetSNPs = {c: list(filter(check, hetSNPs[c])) for c in hetSNPs if len(hetSNPs[c]) > 0}
    return {c: {int(r[2]): (r[3], r[4]) for r in reversed(hetSNPs[c])} for c in hetSNPs if len(hetSNPs[c]) > 0}


def _old_isHet(countA, countB, gamma):
    p_lower = gamma / 2.0
    p_upper = 1.0 - p_lower
    [c_lower, c_upper] = beta.ppf([p_lower, p_upper], countA + 1, countB + 1)
    return c_lower <= 0.5 <= c_upper


def _old_checkShift(countA, countB, maxshift):
    return (0.5 - (float(min(countA, countB)) / float(countA + countB))) <= maxshift


def _old_counting_results(task_rows):
    # what `counting` does with the workers' results
    return {(v[0][0], v[0][1]): v for v in task_rows if v}


def _old_tumor_counts(task_rows, hetSNPs, normal):
    # the tumor part of the original main()
    rcounts = _old_counting_results(task_rows)
    rcounts = {c: dict(map(lambda r: (int(r[2]), dict(r[3])), rcounts[c])) for c in rcounts}
    het = lambda chro: hetSNPs[normal, chro]
    form = lambda REF, ALT, T: (
        (REF, T[REF] if REF in T else 0),
        (ALT, T[ALT] if ALT in T else 0),
    )
    return {c: {o: form(het(c[1])[o][0][0], het(c[1])[o][1][0], rcounts[c][o]) for o in rcounts[c]} for c in rcounts}


# ---------------------------------------------------------------------------------------------------------------------
# count_alleles: random inputs
# ---------------------------------------------------------------------------------------------------------------------


def _random_rows(rng, sample, chrom, n, positions=None):
    """Rows as reported by AlleleCounter: (sample, chromosome, position string, [(allele, count), ...])."""
    rows = []
    for _ in range(n):
        if positions is not None:
            pos = int(rng.choice(positions))
        else:
            # few distinct positions so that some repeat
            pos = int(rng.integers(1, max(2, n)))
        n_alleles = int(rng.integers(0, 5))
        alleles = list(rng.choice(list('ACGT'), size=n_alleles, replace=False))
        scale = [0, 3, 30, 300][int(rng.integers(0, 4))]
        # ties between the other alleles, and sites without reads
        counts = [int(rng.integers(0, scale + 1)) for _ in alleles]
        if n_alleles >= 3 and rng.random() < 0.3:
            counts[2] = counts[1]
        rows.append((sample, chrom, str(pos), list(zip(alleles, counts))))
    return rows


def _no_empty_sites(rows):
    # the original raises ZeroDivisionError for a heterozygous site without reads; keep most cases away from it
    return [r for r in rows if len(r[3]) < 2 or r[3][0][1] + max(c for _, c in r[3][1:]) > 0]


class _StubRows:
    """Replaces AlleleCounter.iterAlleles with given rows, keyed by (sample, chromosome)."""

    def __init__(self, rows_by_task):
        self.rows_by_task = rows_by_task

    def __call__(self, bamfile, samplename, chromosome):
        yield from self.rows_by_task[samplename, chromosome]


_COUNTER_ARGS = dict(
    bcftools='bcftools', reference='ref.fa', q=0, Q=11, mincov=0, dp=1000, E=False, verbose=False, outdir='.'
)


@pytest.mark.parametrize('seed', range(80))
def test_select_het_snps(seed):
    rng = np.random.default_rng(seed)
    gamma = [0.05, 0.0, 1.0, 0.2][seed % 4]
    maxshift = [0.5, 0.1, 0.0, 0.3][(seed // 4) % 4]
    counts = {}
    for chrom in ['chr1', 'chr2', 'chr3']:
        rows = _random_rows(rng, 'Normal', chrom, int(rng.integers(0, 80)))
        if seed % 5 != 0:
            rows = _no_empty_sites(rows)
        counts['Normal', chrom] = rows

    expected = _outcome_exc(_old_selectHetSNPs, counts, gamma, maxshift)
    got = _outcome_exc(count_alleles.selectHetSNPs, counts, gamma, maxshift)
    assert got == expected


def _outcome_exc(fn, *args, **kwargs):
    try:
        return fn(*args, **kwargs)
    except ZeroDivisionError as e:
        return ('ZeroDivisionError', str(e))


@pytest.mark.parametrize('seed', range(30))
def test_het_mask(seed):
    rng = np.random.default_rng(seed)
    gamma = float(rng.choice([0.0, 1.0, 0.05, rng.random()]))
    maxshift = float(rng.choice([0.0, 0.5, rng.random() / 2]))
    a = rng.integers(0, [2, 20, 2000][seed % 3], size=300)
    b = rng.integers(0, [2, 20, 2000][seed % 3], size=300)
    expected = []
    for x, y in zip(a.tolist(), b.tolist()):
        expected.append(_outcome_exc(lambda: _old_isHet(x, y, gamma) and _old_checkShift(x, y, maxshift)))
    if any(isinstance(e, tuple) for e in expected):
        assert _outcome_exc(count_alleles.het_mask, a, b, gamma, maxshift) == (
            'ZeroDivisionError',
            'float division by zero',
        )
    else:
        assert count_alleles.het_mask(a, b, gamma, maxshift).tolist() == [bool(e) for e in expected]


@pytest.mark.parametrize('seed', range(40))
def test_normal_and_tumor_passes(seed):
    """
    The normal pass (HetSiteCounter) and the tumor pass (TumorAlleleCounter), with the collection done in main(),
    against the original main(): selected SNPs, the normal output, and the tumor output.
    """
    rng = np.random.default_rng(seed)
    gamma, maxshift = 0.05, 0.5
    chromosomes = ['chr1', 'chr2', 'chr3']
    normal = 'Normal'

    # normal pass: tasks per chromosome; some report nothing, some report no heterozygous SNPs
    normal_rows = {}
    for chrom in chromosomes:
        n = int(rng.choice([0, 1, 5, 60]))
        normal_rows[normal, chrom] = _no_empty_sites(_random_rows(rng, normal, chrom, n))
    tasks = [(normal, chrom) for chrom in chromosomes]

    old_het = _old_selectHetSNPs(_old_counting_results([normal_rows[t] for t in tasks]), gamma, maxshift)

    counter = count_alleles.HetSiteCounter(gamma=gamma, maxshift=maxshift, snplist=None, **_COUNTER_ARGS)
    counter.iterAlleles = _StubRows(normal_rows)
    results = count_alleles.collect_by_first_site([counter.work('bam', *t) for t in tasks])
    new_het = {c: het for c, het in results.items() if het is not None}

    assert set(new_het) == set(old_het)
    for c in old_het:
        # normal output lines and TMP target lists: sorted positions
        expected = [(snp, v[0][0], v[0][1], v[1][0], v[1][1]) for snp, v in sorted(old_het[c].items())]
        got = list(zip(*(column.tolist() for column in new_het[c])))
        assert got == expected

    # tumor pass: the sites reported are the normal's heterozygous SNPs (possibly repeated), any alleles
    tumor_rows = {}
    for sample in ['Tumor1', 'Tumor2']:
        for chrom in chromosomes:
            if (normal, chrom) in old_het and rng.random() < 0.8:
                positions = list(old_het[normal, chrom])
                tumor_rows[sample, chrom] = _random_rows(rng, sample, chrom, int(rng.integers(1, 40)), positions)
            else:
                tumor_rows[sample, chrom] = []
    tasks = list(tumor_rows)

    old_counts = _old_tumor_counts([tumor_rows[t] for t in tasks], old_het, normal)

    counter = count_alleles.TumorAlleleCounter(
        het_sites={chro: het for (_, chro), het in new_het.items()}, snplist=None, **_COUNTER_ARGS
    )
    counter.iterAlleles = _StubRows(tumor_rows)
    new_counts = count_alleles.collect_by_first_site([counter.work('bam', *t) for t in tasks])

    assert set(new_counts) == set(old_counts)
    for c in old_counts:
        # tumor output lines, in order
        expected = [(snp, v[0][1], v[1][1]) for snp, v in old_counts[c].items()]
        got = list(zip(*(column.tolist() for column in new_counts[c])))
        assert got == expected
