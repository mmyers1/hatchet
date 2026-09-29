import sys
import os
from multiprocessing import Pool
from multiprocessing.pool import ThreadPool
import subprocess as sp
import numpy as np
import pandas as pd
import subprocess
import traceback
from importlib.resources import path
import hatchet.data

from hatchet.utils.ArgParsing import parse_count_reads_args
from hatchet.utils.Supporting import log, logArgs, error


def main(args=None):
    log(msg='# Parsing and checking input arguments\n', level='STEP')
    args = parse_count_reads_args(args)
    logArgs(args, 80)

    bams = args['bams']
    names = args['names']

    tbams, tnames = zip(*sorted(zip(*(bams[1:], names[1:])), key=lambda x: x[1]))
    bams = [bams[0]] + list(tbams)
    names = [names[0]] + list(tnames)

    chromosomes = args['chromosomes']
    samtools = args['samtools']
    processes = args['j']
    outdir = args['outdir']
    mosdepth = args['mosdepth']
    tabix = args['tabix']
    readquality = args['readquality']

    if len(check_array_files(outdir, chromosomes)) == 0:
        log(
            msg='# Found all array files, skipping to total read counting. \n',
            level='STEP',
        )
    else:

        if len(check_counts_files(outdir, chromosomes, names)) == 0:
            log(
                msg='# Found all count files, skipping to forming arrays\n',
                level='STEP',
            )

        else:
            params = zip(
                np.repeat(chromosomes, len(bams)),
                [outdir] * len(bams) * len(chromosomes),
                [samtools] * len(bams) * len(chromosomes),
                bams * len(chromosomes),
                names * len(chromosomes),
                [readquality] * len(bams) * len(chromosomes),
            )

            n_workers_samtools = min(processes, len(bams) * len(chromosomes))
            with Pool(n_workers_samtools) as p:   # divide by 2 because each worker starts 2 processes
                p.map(count_chromosome_wrapper, params)

            n_workers_mosdepth = min(processes, len(bams))

            # compute number of decompression threads to use for each call to mosdepth
            if processes > len(bams):
                base = min(int(processes / n_workers_mosdepth), 4)
                threads_per_task = [base] * len(bams)

                if base < 4:
                    remainder = processes % n_workers_mosdepth
                    i = 0
                    while remainder > 0:
                        threads_per_task[i] += 1
                        i += 1
                        remainder -= 1
            else:
                threads_per_task = [1] * len(bams)

            # Note: These function calls are the only section that uses mosdepth
            mosdepth_params = [
                (
                    outdir,
                    names[i],
                    bams[i],
                    threads_per_task[i],
                    mosdepth,
                    readquality,
                )
                for i in range(len(bams))
            ]
            with Pool(n_workers_mosdepth) as p:
                p.map(mosdepth_wrapper, mosdepth_params)

            if len(check_counts_files(outdir, chromosomes, names)) > 0:
                raise ValueError(error('Missing some counts files!'))

        ### Aggregate count files into count arrays for adaptive binning ###
        # (formerly formArray)
        use_chr = args['use_chr']

        n_workers = min(len(chromosomes), processes)

        # Read in centromere locations table
        with path(hatchet.data, f'{args["refversion"]}.centromeres.txt') as centromeres:
            centromeres = pd.read_table(
                centromeres,
                header=None,
                names=['CHR', 'START', 'END', 'NAME', 'gieStain'],
            )

        chr2centro = {}
        for ch in centromeres.CHR.unique():
            my_df = centromeres[centromeres.CHR == ch]
            assert (my_df.gieStain == 'acen').all()
            # Each centromere should consist of 2 adjacent segments
            assert len(my_df == 2)
            assert my_df.START.max() == my_df.END.min()
            if use_chr:
                if ch.startswith('chr'):
                    chr2centro[ch] = my_df.START.min(), my_df.END.max()
                else:
                    chr2centro['chr' + ch] = my_df.START.min(), my_df.END.max()
            else:
                if ch.startswith('chr'):
                    chr2centro[ch[3:]] = my_df.START.min(), my_df.END.max()
                else:
                    chr2centro[ch] = my_df.START.min(), my_df.END.max()

        for ch in chromosomes:
            if ch not in chr2centro:
                raise ValueError(
                    error(f'Chromosome {ch} not found in centromeres file. Inspect file provided as -C argument.')
                )

        # Use Tabix to index per-position coverage bed files for each sample
        perpos_files = [os.path.join(outdir, name + '.per-base.bed.gz') for name in names]
        with ThreadPool(min(processes, len(perpos_files))) as p:
            p.map(lambda perpos_file: subprocess.run([tabix, '-f', perpos_file]), perpos_files)

        # Read the SNP positions once (rather than once per chromosome worker), and only for the chromosomes that
        # still need arrays; sex chromosomes do not use SNPs
        todo = [ch for ch in chromosomes if not is_sex_chromosome(ch) and not chromosome_arrays_exist(outdir, ch)]
        positions = {}
        if len(todo) > 0:
            log(msg=f'Reading SNPs file {args["baf_file"]}\n', level='INFO')
            all_snps = load_baf_table(args['baf_file'])
            for ch in todo:
                positions[ch], _, _ = snps_for_chromosome(all_snps, ch, names)
            del all_snps

        # form parameters for each worker
        params = [
            (
                outdir,
                names,
                ch,
                chr2centro[ch][0],
                chr2centro[ch][1],
                positions.get(ch),
                tabix,
            )
            for ch in chromosomes
        ]

        # dispatch workers
        with Pool(n_workers) as p:
            p.map(run_chromosome_wrapper, params)

        np.savetxt(os.path.join(outdir, 'samples.txt'), names, fmt='%s')

        if len(check_array_files(outdir, chromosomes)) > 0:
            raise ValueError(error('Missing some output arrays!'))
        else:
            log(
                msg='# Array forming completed successfully, removing intermediate count files. \n',
                level='STEP',
            )
            [os.remove(f) for f in expected_counts_files(outdir, chromosomes, names)]

    totals_file = os.path.join(outdir, 'total.tsv')
    if os.path.exists(totals_file):
        log(msg='# Found total reads file, exiting. \n', level='STEP')
        return
    else:
        log(
            msg='# Counting total number of reads for normal and tumor samples\n',
            level='STEP',
        )
        total = total_counts_from_arrays(outdir, chromosomes, names)

        log(
            msg='# Writing the total read counts for all samples in {}\n'.format(totals_file),
            level='STEP',
        )
        with open(totals_file, 'w') as f:
            for name in names:
                f.write('{}\t{}\n'.format(name, total[name]))


def mosdepth_wrapper(params):
    run_mosdepth(*params)


def run_mosdepth(outdir, sample_name, bam, threads, mosdepth, readquality):
    try:
        last_file = os.path.join(outdir, sample_name + '.per-base.bed.gz.csi')
        if os.path.exists(last_file):
            log(
                f'Skipping mosdepth on sample {sample_name} (output file {last_file} exists)\n',
                level='STEP',
            )
            return

        log(
            f'Starting mosdepth on sample {sample_name} with {threads} threads\n',
            level='STEP',
        )
        sys.stderr.flush()

        md = sp.run(
            [
                mosdepth,
                '-t',
                str(threads),
                '-Q',
                str(readquality),
                os.path.join(outdir, sample_name),
                bam,
            ]
        )
        md.check_returncode()
        log(f'Done mosdepth on sample {sample_name}\n', level='STEP')

    except Exception as e:
        log('Exception in countPos: {}\n'.format(e), level='ERROR')
        raise e


def count_chromosome(ch, outdir, samtools, bam, sample_name, readquality, compression_level=6):
    try:
        outfile = os.path.join(outdir, f'{sample_name}.{ch}.starts.gz')
        if os.path.exists(outfile):
            log(
                f'Skipping sample {sample_name} chromosome {ch} (output file exists)\n',
                level='STEP',
            )
            return

        log(f'Sample {sample_name} -- Starting chromosome {ch}\n', level='STEP')
        sys.stderr.flush()

        # Get start positions. They are written to a temporary file that is renamed once complete, so that an
        # interrupted run leaves no partial output file for a later run to skip.
        tmpfile = outfile + '.tmp'
        with open(tmpfile, 'w') as out:
            st = sp.Popen((samtools, 'view', '-q', str(readquality), bam, ch), stdout=sp.PIPE)
            cut = sp.Popen(('cut', '-f', '4'), stdin=st.stdout, stdout=sp.PIPE)
            st.stdout.close()   # so that samtools gets SIGPIPE if cut exits early
            gz = sp.Popen(('gzip', '-{}'.format(compression_level)), stdin=cut.stdout, stdout=out)
            cut.stdout.close()
            st.wait()
            cut.wait()
            gz.wait()
        for name, p in (('samtools', st), ('cut', cut), ('gzip', gz)):
            if p.returncode != 0:
                os.remove(tmpfile)
                raise ValueError('{} subprocess returned nonzero value: {}'.format(name, p.returncode))
        os.replace(tmpfile, outfile)

        log(f'Sample {sample_name} -- Done chromosome {ch}\n', level='STEP')

    except Exception as e:
        log('Exception in countPos: {}\n'.format(e), level='ERROR')
        raise e


def count_chromosome_wrapper(param):
    count_chromosome(*param)


def read_snps(baf_file, ch, all_names):
    """
    Read and validate SNP data for this patient (TSV table output from HATCHet deBAF.py).
    """
    return snps_for_chromosome(load_baf_table(baf_file), ch, all_names)


def load_baf_table(baf_file):
    """
    Read the HATCHet BAF table for all chromosomes (the categorical CHR and SAMPLE columns keep it compact).
    """
    return pd.read_table(
        baf_file,
        names=['CHR', 'POS', 'SAMPLE', 'REF', 'ALT'],
        dtype={
            'CHR': 'category',
            'POS': np.uint32,
            'SAMPLE': 'category',
            'ALT': np.uint32,
            'REF': np.uint32,
        },
    )


def snps_for_chromosome(all_snps, ch, all_names):
    """
    Validate the SNP data of one chromosome from a table read by load_baf_table.
    Returns the positions present in all tumor samples, their total counts, and the table of those SNPs.
    """
    all_names = all_names[1:]   # remove normal sample -- not looking for SNP counts from normal

    # Keep only SNPs on this chromosome
    snps = all_snps[all_snps.CHR == ch].astype({'CHR': object, 'SAMPLE': object}).sort_values(by=['POS', 'SAMPLE'])
    snps = snps.reset_index(drop=True)

    if len(snps) == 0:
        raise ValueError(
            error(f'Chromosome {ch} not found in SNPs file (chromosomes in file: {np.asarray(all_snps.CHR.unique())})')
        )

    n_samples = len(all_names)
    if n_samples != len(snps.SAMPLE.unique()):
        raise ValueError(
            error(f'Expected {n_samples} samples, found {len(snps.SAMPLE.unique())} samples in SNPs file.')
        )

    if set(all_names) != set(snps.SAMPLE.unique()):
        raise ValueError(
            error(
                f'Expected sample names did not match sample names in SNPs file.\n\
                Expected: {sorted(all_names)}\n  Found:{sorted(snps.SAMPLE.unique())}'
            )
        )

    # Add total counts column
    snpsv = snps.copy()
    snpsv['TOTAL'] = snpsv.ALT + snpsv.REF

    # Create counts array and find SNPs that are not present in all samples
    snp_counts = snpsv.pivot(index='POS', columns='SAMPLE', values='TOTAL')
    missing_pos = snp_counts.isna().any(axis=1)

    # Remove SNPs that are absent in any sample
    snp_counts = snp_counts.dropna(axis=0)
    snpsv = snpsv[~snpsv.POS.isin(missing_pos[missing_pos].index)]

    # Pivot table for dataframe should match counts array and have no missing entries
    check_pivot = snpsv.pivot(index='POS', columns='SAMPLE', values='TOTAL')
    assert np.array_equal(check_pivot, snp_counts), 'SNP file reading failed'
    assert not np.any(check_pivot.isna()), 'SNP file reading failed'
    assert np.array_equal(all_names, list(snp_counts.columns)), (
        all_names,
        list(snp_counts.columns),
    )   # make sure that sample order is the same

    return np.array(snp_counts.index), np.array(snp_counts), snpsv


def form_counts_array(starts_files, perpos_files, thresholds, chromosome, tabix, chunksize=1e6):
    """
    NOTE: Assumes that starts_files[i] corresponds to the same sample as perpos_files[i]
    Parameters:
        starts_files: list of <sample>.<chromosome>.starts.gz files each containing a list of start positions
        perpos_files: list of <sample>.per-base.bed.gz files containing per-position coverage from mosdepth
        thresholds: list of potential bin start positions (thresholds between SNPs)
        chromosome: chromosome to extract read counts for
        chunksize: number of lines read at a time from each file

    Returns: <n> x <2d> np.ndarray
        entry [i, 2j] contains the number of reads starting in (starts[i], starts[i + 1]) in sample j
        entry [i, 2j + 1] contains the number of reads covering position starts[i] in sample j
    """
    # add one for the end of the chromosome
    arr = np.zeros((thresholds.shape[0] + 1, len(starts_files) * 2), dtype=np.int64)
    chunksize = int(chunksize)

    for i in range(len(starts_files)):
        # populate starts in even entries
        arr[:-2, i * 2] = bin_read_starts(read_starts_chunks(starts_files[i], chunksize), thresholds)

    for i in range(len(perpos_files)):
        # populate threshold coverage in odd entries
        coverage, idx, last_record = coverage_at_thresholds(
            read_perpos_chunks(perpos_files[i], chromosome, tabix, chunksize), thresholds
        )
        arr[: len(thresholds), (2 * i) + 1] = coverage

        if i == 0:
            # identify the (effective) chromosome end as the last well-formed record
            assert idx == len(thresholds)
            chr_end, end_reads = last_record

            assert chr_end > thresholds[-1]
            assert len(thresholds) == len(arr) - 1

            # add the chromosome end to thresholds
            thresholds = np.concatenate([thresholds, [chr_end]])

            # count the number of reads covering the chromosome end
            arr[idx, 0] = end_reads

    return arr, thresholds


def bin_read_starts(chunks, thresholds):
    """
    Count read starts (given as chunks of positions, in file order) between consecutive thresholds.
    Returns an array c of length len(thresholds) - 1, where c[k] is the number of reads starting in
    [thresholds[k], thresholds[k + 1]); the reads from the first one starting at or after thresholds[-1] onward are
    all added to c[-1]. Before that read, reads may not start before the bin of the previous read (asserted).
    """
    n_bins = len(thresholds) - 1
    counts = np.zeros(n_bins, dtype=np.int64)
    last_bin = 0   # bin of the previous read
    n_tail = 0   # number of reads from the first read at or after thresholds[-1] onward
    for starts in chunks:
        if n_tail > 0:
            n_tail += len(starts)
            continue

        bins = np.searchsorted(thresholds, starts, side='right') - 1
        in_tail = bins >= n_bins
        first_tail = np.argmax(in_tail) if in_tail.any() else len(bins)

        body = bins[:first_tail]
        if len(body) > 0:
            assert body[0] >= last_bin and np.all(body[1:] >= body[:-1])
            counts += np.bincount(body, minlength=n_bins)
            last_bin = body[-1]
        n_tail += len(bins) - first_tail

    counts[-1] += n_tail
    return counts


def coverage_at_thresholds(records, thresholds):
    """
    Find the depth at each threshold (0-based position thresholds[k] - 1) from per-base coverage records, given as
    chunks of (start, end, depth) arrays in file order.
    Each threshold takes the depth of the first record whose end is past its position; that record must start at or
    before the position (asserted).
    Returns the depths (0 for thresholds past the last record), the number of thresholds that were assigned a record,
    and the (end, depth) of the last record (None if there are no records).
    """
    positions = thresholds - 1
    coverage = np.zeros(len(thresholds), dtype=np.int64)
    idx = 0
    last_record = None
    for starts, ends, depths in records:
        if len(ends) == 0:
            continue
        last_record = ends[-1], depths[-1]

        # The first record whose end is past a position is also the first whose running maximum end is past it.
        # Positions not assigned in earlier chunks are at or past the ends of all earlier records.
        max_ends = np.maximum.accumulate(ends)
        n_here = np.searchsorted(positions[idx:], max_ends[-1], side='left')
        if n_here == 0:
            continue

        my_positions = positions[idx : idx + n_here]
        my_records = np.searchsorted(max_ends, my_positions, side='right')
        assert np.all(starts[my_records] <= my_positions)
        coverage[idx : idx + n_here] = depths[my_records]
        idx += n_here

    return coverage, idx, last_record


def read_starts_chunks(fname, chunksize):
    """
    Yield the read start positions in a gzipped <sample>.<chromosome>.starts.gz file as arrays of up to chunksize.
    """
    try:
        reader = pd.read_csv(fname, header=None, dtype=np.int64, compression='gzip', chunksize=chunksize)
    except pd.errors.EmptyDataError:
        # no reads on this chromosome
        return
    with reader:
        for chunk in reader:
            yield chunk[0].to_numpy()


def read_perpos_chunks(perpos_file, chromosome, tabix, chunksize):
    """
    Yield (start, end, depth) arrays of up to chunksize per-base coverage records on this chromosome, streamed from
    tabix.
    """
    query = subprocess.Popen([tabix, perpos_file, chromosome], stdout=subprocess.PIPE)
    try:
        try:
            reader = pd.read_csv(
                query.stdout,
                sep='\t',
                header=None,
                usecols=[1, 2, 3],
                dtype={1: np.int64, 2: np.int64, 3: np.int64},
                chunksize=chunksize,
            )
        except pd.errors.EmptyDataError:
            return
        with reader:
            for chunk in reader:
                yield chunk[1].to_numpy(), chunk[2].to_numpy(), chunk[3].to_numpy()
    finally:
        query.stdout.close()
        query.wait()


def get_chr_end(stem, all_names, chromosome):
    starts_files = []
    for name in all_names:
        starts_files.append(os.path.join(stem, name + '.' + chromosome + '.starts.gz'))

    last_start = 0
    for sfname in starts_files:
        zcat = subprocess.Popen(('zcat', sfname), stdout=subprocess.PIPE)
        tail = subprocess.Popen(('tail', '-1'), stdin=zcat.stdout, stdout=subprocess.PIPE)
        zcat.stdout.close()
        my_last = int(tail.communicate()[0].decode('utf-8').strip())
        zcat.wait()

        if my_last > last_start:
            last_start = my_last

    return last_start


def is_sex_chromosome(ch):
    return ch.endswith('X') or ch.endswith('Y')


def chromosome_arrays_exist(outdir, ch):
    return os.path.exists(os.path.join(outdir, f'{ch}.total.gz')) and os.path.exists(
        os.path.join(outdir, f'{ch}.thresholds.gz')
    )


def total_counts_from_arrays(outdir, chromosomes, names):
    """
    Total number of reads of each sample over the given chromosomes, from the count arrays. Their read-start columns
    count every read reported by `samtools view -q <readquality> <bam> <chromosome>` (the last row holds the coverage
    at the chromosome end instead), so they sum to `samtools view -c -q <readquality> <bam> <chromosome>`.
    """
    with open(os.path.join(outdir, 'samples.txt')) as f:
        array_names = [line.rstrip('\n') for line in f]

    total = dict.fromkeys(array_names, 0)
    for ch in chromosomes:
        arr = pd.read_csv(os.path.join(outdir, f'{ch}.total.gz'), sep=' ', header=None, dtype=np.int64).to_numpy()
        for name, n_reads in zip(array_names, arr[:-1, 0::2].sum(axis=0).tolist()):
            total[name] += n_reads

    try:
        return {name: total[name] for name in names}
    except KeyError:
        raise KeyError(error('Either a chromosome or a sample has not been considered in the total counting!'))


def run_chromosome(
    outdir,
    all_names,
    chromosome,
    centromere_start,
    centromere_end,
    positions,
    tabix,
):
    """
    Construct arrays that contain all counts needed to perform adaptive binning for a single chromosome
    (across all samples).
    positions: SNP positions on this chromosome (from snps_for_chromosome); not used for sex chromosomes
    """

    try:
        totals_out = os.path.join(outdir, f'{chromosome}.total.gz')
        thresholds_out = os.path.join(outdir, f'{chromosome}.thresholds.gz')

        if chromosome_arrays_exist(outdir, chromosome):
            log(
                msg=f'Output files already exist, skipping chromosome {chromosome}\n',
                level='INFO',
            )
            return

        log(msg=f'Loading chromosome {chromosome}\n', level='INFO')
        # Per-position coverage bed files for each sample
        perpos_files = [os.path.join(outdir, name + '.per-base.bed.gz') for name in all_names]

        # Identify the start-positions files for this chromosome
        starts_files = []
        for name in all_names:
            starts_files.append(os.path.join(outdir, name + '.' + chromosome + '.starts.gz'))

        if is_sex_chromosome(chromosome):
            log(
                msg='Running on sex chromosome -- ignoring SNPs and min SNP reads\n',
                level='INFO',
            )

            # TODO: do this procedure only for XY
            last_start = get_chr_end(outdir, all_names, chromosome)
            positions = np.arange(5000, last_start, 5000)

        thresholds = np.trunc(np.vstack([positions[:-1], positions[1:]]).mean(axis=0)).astype(np.uint32)
        last_idx_p = np.argwhere(thresholds > centromere_start)[0][0]
        first_idx_q = np.argwhere(thresholds > centromere_end)[0][0]
        all_thresholds = np.concatenate(
            [
                [1],
                thresholds[:last_idx_p],
                [centromere_start],
                [centromere_end],
                thresholds[first_idx_q:],
            ]
        )

        log(msg=f'Loading counts for chromosome {chromosome}\n', level='INFO')
        # Load read count arrays from file (also adds end of chromosome as a threshold)
        total_counts, complete_thresholds = form_counts_array(
            starts_files, perpos_files, all_thresholds, chromosome, tabix=tabix
        )

        np.savetxt(totals_out, total_counts, fmt='%d')
        np.savetxt(thresholds_out, complete_thresholds, fmt='%d')

        log(msg=f'Done chromosome {chromosome}\n', level='INFO')
    except Exception as e:
        print(f'Error in chromosome {chromosome}:')
        print(e)
        traceback.print_exc()
        raise e


def run_chromosome_wrapper(param):
    run_chromosome(*param)


def expected_counts_files(dcounts, chrs, all_names):
    expected = []

    for chr in chrs:
        for name in all_names:
            fname = os.path.join(dcounts, '.'.join([name, chr, 'starts.gz']))
            expected.append(fname)

    # mosdepth output
    for name in all_names:
        fname = os.path.join(dcounts, name + '.mosdepth.global.dist.txt')
        expected.append(fname)

        fname = os.path.join(dcounts, name + '.mosdepth.summary.txt')
        expected.append(fname)

        fname = os.path.join(dcounts, name + '.per-base.bed.gz')
        expected.append(fname)

        fname = os.path.join(dcounts, name + '.per-base.bed.gz.csi')
        expected.append(fname)

    return expected


def check_counts_files(dcounts, chrs, all_names):
    return [a for a in expected_counts_files(dcounts, chrs, all_names) if not os.path.exists(a)]


def expected_arrays(darray, chrs):
    expected = []
    # formArray (abin/<chr>.<total/thresholds> files))

    fname = os.path.join(darray, 'samples.txt')
    expected.append(fname)

    for ch in chrs:
        fname = os.path.join(darray, f'{ch}.total.gz')
        expected.append(fname)

        fname = os.path.join(darray, f'{ch}.thresholds.gz')
        expected.append(fname)

    return expected


def check_array_files(darray, chrs):
    return [a for a in expected_arrays(darray, chrs) if not os.path.exists(a)]


if __name__ == '__main__':
    main()
