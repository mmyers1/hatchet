import sys
import os
import os.path
import shlex
import subprocess as pr
from array import array
from collections import namedtuple
import numpy as np
from scipy.stats import beta
import tempfile

import hatchet.utils.ArgParsing as ap
from hatchet.utils.Supporting import log, logArgs, error, close, indexed_region, MPILEUP_EXCLUDE_FLAGS
from hatchet.utils.multiprocessing import Worker

ALLELES = {'A', 'C', 'G', 'T'}

# Heterozygous SNPs of one sample and chromosome, as arrays sorted by position
HetSites = namedtuple('HetSites', ['pos', 'ref', 'refcount', 'alt', 'altcount'])
# Read counts of a tumor sample for the two alleles of heterozygous SNPs, as arrays in the order they were reported
TumorCounts = namedtuple('TumorCounts', ['pos', 'refcount', 'altcount'])


def main(args=None):
    log(
        msg=(
            '# Parsing the input arguments, checking the consistency of given files, and extracting required '
            'information\n'
        ),
        level='STEP',
    )
    args = ap.parse_count_alleles_arguments(args)
    logArgs(args, 80)

    counter_args = dict(
        bcftools=args['bcftools'],
        reference=args['reference'],
        q=args['q'],
        Q=args['Q'],
        mincov=args['mincov'],
        dp=args['maxcov'],
        E=args['E'],
        verbose=args['verbose'],
        outdir=args['outputSnps'],
    )
    normal = args['normal'][1]

    log(
        msg='# Counting SNPs alleles from the matched-normal sample\n',
        level='STEP',
    )
    log(msg='# Selecting heterozygous SNPs\n', level='STEP')
    hetSNPs = counting_het_sites(
        samples=[args['normal']],
        chromosomes=args['chromosomes'],
        num_workers=args['j'],
        snplist=args['snps'],
        gamma=args['gamma'],
        maxshift=args['maxshift'],
        **counter_args,
    )
    if not hetSNPs:
        close('No heterozygous SNPs found in the selected regions of the normal!\n')

    log(
        msg='# Writing the list of selected SNPs, covered and heterozygous in the normal sample\n',
        level='STEP',
    )
    # put lists of SNPs in uniquely named temporary dir to prevent concurrent hatchet runs from overwriting
    # each others temp files
    with tempfile.TemporaryDirectory(dir=args['outputSnps']) as tmpdirname:
        hetsnpsfiles = {}
        for chro in args['chromosomes']:
            hetsnpsfiles[chro] = os.path.join(tmpdirname, 'TMP_{}.tsv'.format(chro))
            with open(hetsnpsfiles[chro], 'w') as f:
                if (normal, chro) in hetSNPs:
                    f.writelines('{}\t{}\n'.format(chro, snp) for snp in hetSNPs[normal, chro].pos.tolist())

        log(
            msg='# Writing the allele counts of the normal sample for selected SNPs\n',
            level='STEP',
        )
        handle = open(args['outputNormal'], 'w') if args['outputNormal'] is not None else sys.stdout
        for chro in args['chromosomes']:
            if (normal, chro) in hetSNPs:
                het = hetSNPs[normal, chro]
                handle.writelines(
                    '{}\t{}\t{}\t{}\t{}\n'.format(chro, snp, normal, refcount, altcount)
                    for snp, refcount, altcount in zip(het.pos.tolist(), het.refcount.tolist(), het.altcount.tolist())
                )
        if handle is not sys.stdout:
            handle.close()

        log(msg='# Counting SNPs alleles from tumour samples\n', level='STEP')
        counts = counting_tumor_alleles(
            samples=args['samples'],
            chromosomes=args['chromosomes'],
            num_workers=args['j'],
            snplist=hetsnpsfiles,
            het_sites={chro: het for (_, chro), het in hetSNPs.items()},
            **counter_args,
        )
        if not counts:
            close('The selected SNPs are not covered in the tumors!\n')

    log(
        msg='# Writing the allele counts of tumor samples for selected SNPs\n',
        level='STEP',
    )
    handle = open(args['outputTumors'], 'w') if args['outputTumors'] is not None else sys.stdout
    for sample in args['samples']:
        for chro in args['chromosomes']:
            if (sample[1], chro) in counts:
                count = counts[sample[1], chro]
                handle.writelines(
                    '{}\t{}\t{}\t{}\t{}\n'.format(chro, snp, sample[1], refcount, altcount)
                    for snp, refcount, altcount in zip(
                        count.pos.tolist(), count.refcount.tolist(), count.altcount.tolist()
                    )
                )
    if handle is not sys.stdout:
        handle.close()


def selectHetSNPs(counts, gamma, maxshift):
    selected = {}
    for c in counts:
        het = select_het_sites(*het_candidates(counts[c]), gamma=gamma, maxshift=maxshift)
        if het is not None:
            selected[c] = {
                snp: ((ref, refcount), (alt, altcount))
                for snp, ref, refcount, alt, altcount in zip(*(column.tolist() for column in het))
            }
    return selected


def het_candidates(rows):
    """
    Columns (position, first allele, its count, most covered other allele, its count) of the rows reported by
    AlleleCounter that have at least two alleles. The first of the equally covered other alleles is taken.
    """
    pos, ref, refcount, alt, altcount = array('q'), [], array('q'), [], array('q')
    for r in rows:
        if len(r[3]) >= 2:
            a = max(r[3][1:], key=(lambda x: x[1]))
            pos.append(int(r[2]))
            ref.append(r[3][0][0])
            refcount.append(r[3][0][1])
            alt.append(a[0])
            altcount.append(a[1])
    return pos, ref, refcount, alt, altcount


def select_het_sites(pos, ref, refcount, alt, altcount, gamma, maxshift):
    """
    HetSites of the candidates (as returned by het_candidates) that are heterozygous (see het_mask), or None if there
    are none. A position reported more than once keeps its first heterozygous record.
    """
    if len(pos) == 0:
        return None
    refcount = np.asarray(refcount, dtype=np.int64)
    altcount = np.asarray(altcount, dtype=np.int64)
    het = np.flatnonzero(het_mask(refcount, altcount, gamma, maxshift))
    if len(het) == 0:
        return None
    pos, first = np.unique(np.asarray(pos, dtype=np.int64)[het], return_index=True)
    keep = het[first]
    return HetSites(pos, np.asarray(ref)[keep], refcount[keep], np.asarray(alt)[keep], altcount[keep])


def het_mask(countA, countB, gamma, maxshift):
    """
    Vectorized `isHet(countA, countB, gamma) and checkShift(countA, countB, maxshift)` over arrays of counts.
    Like checkShift, raises ZeroDivisionError for a heterozygous site with no reads.
    """
    countA = np.asarray(countA, dtype=np.int64)
    countB = np.asarray(countB, dtype=np.int64)
    p_lower = gamma / 2.0
    p_upper = 1.0 - p_lower
    c_lower = beta.ppf(p_lower, countA + 1, countB + 1)
    c_upper = beta.ppf(p_upper, countA + 1, countB + 1)
    het = (c_lower <= 0.5) & (0.5 <= c_upper)

    total = countA + countB
    if np.any(het & (total == 0)):
        raise ZeroDivisionError('float division by zero')
    with np.errstate(divide='ignore', invalid='ignore'):
        shift = 0.5 - np.minimum(countA, countB) / total
    return het & (shift <= maxshift)


def isHet(countA, countB, gamma):
    p_lower = gamma / 2.0
    p_upper = 1.0 - p_lower
    [c_lower, c_upper] = beta.ppf([p_lower, p_upper], countA + 1, countB + 1)
    return c_lower <= 0.5 <= c_upper


# def isHet(countA, countB, gamma):
#     lb, ub = proportion_confint(min(countA, countB), countA+countB, alpha=gamma, method='beta')
#     return lb <= 0.5 <= ub

# def isHet(countA, countB, gamma):
#     lb, ub = proportion_confint(min(countA, countB), countA+countB, alpha=gamma, method='jeffreys')
#     return lb <= 0.5 <= ub


def checkShift(countA, countB, maxshift):
    return (0.5 - (float(min(countA, countB)) / float(countA + countB))) <= maxshift


def counting(
    bcftools,
    reference,
    samples,
    chromosomes,
    num_workers,
    snplist,
    q,
    Q,
    mincov,
    dp,
    E,
    verbose,
    outdir,
):

    work = []
    for bam in samples:
        for chro in chromosomes:
            work.append((bam[0], bam[1], chro))

    worker = AlleleCounter(
        bcftools,
        reference,
        q,
        Q,
        mincov,
        dp,
        E,
        snplist,
        verbose,
        outdir,
    )

    results = worker.run(work=work, n_instances=num_workers)
    results = {(v[0][0], v[0][1]): v for v in results if v}
    return results


def counting_het_sites(samples, chromosomes, num_workers, gamma, maxshift, **counter_args):
    """
    Equivalent to `selectHetSNPs(counting(...), gamma, maxshift)`, with HetSites values: the selection is done by the
    workers, so that only the heterozygous SNPs are passed back.
    """
    work = [(bam[0], bam[1], chro) for bam in samples for chro in chromosomes]
    worker = HetSiteCounter(gamma=gamma, maxshift=maxshift, **counter_args)
    results = collect_by_first_site(worker.run(work=work, n_instances=num_workers))
    return {c: het for c, het in results.items() if het is not None}


def counting_tumor_alleles(samples, chromosomes, num_workers, het_sites, **counter_args):
    """
    Read counts of the tumor samples for the two alleles of the normal's heterozygous SNPs (het_sites, by chromosome),
    as TumorCounts keyed like `counting`.
    """
    work = [(bam[0], bam[1], chro) for bam in samples for chro in chromosomes]
    worker = TumorAlleleCounter(het_sites=het_sites, **counter_args)
    return collect_by_first_site(worker.run(work=work, n_instances=num_workers))


def collect_by_first_site(results):
    """
    Key the results of tasks that reported sites by (sample, chromosome of the first site); as in `counting`, a later
    task replaces an earlier one with the same key.
    """
    collected = {}
    for result in results:
        if result is not None:
            key, value = result
            collected[key] = value
    return collected


class AlleleCounter(Worker):
    def __init__(
        self,
        bcftools,
        reference,
        q,
        Q,
        mincov,
        dp,
        E,
        snplist,
        verbose,
        outdir,
    ):
        self.bcftools = bcftools
        self.reference = reference
        self.q = q
        self.Q = Q
        self.mincov = mincov
        self.dp = dp
        self.E = E
        self.snplist = snplist
        self.verbose = verbose
        self.outdir = outdir

    def work(self, bamfile, samplename, chromosome):
        return self.countAlleles(bamfile=bamfile, samplename=samplename, chromosome=chromosome)

    def countAlleles(self, bamfile, samplename, chromosome):
        return list(self.iterAlleles(bamfile=bamfile, samplename=samplename, chromosome=chromosome))

    def iterAlleles(self, bamfile, samplename, chromosome):
        """
        Yield (samplename, chromosome, position, [(allele, count), ..]) for each site reported by bcftools, keeping
        only the A/C/G/T alleles.
        """
        cmd_mpileup = '{} mpileup {} -Ou -f {} --skip-indels -a INFO/AD --ff {} -q {} -Q {} -d {} -T {}'.format(
            self.bcftools,
            bamfile,
            self.reference,
            MPILEUP_EXCLUDE_FLAGS,
            self.q,
            self.Q,
            self.dp,
            self.snplist[chromosome],
        )
        # With only -T, mpileup streams the whole BAM; restricting it to the chromosome makes it jump there instead
        region = indexed_region(bamfile, chromosome)
        if region is not None:
            cmd_mpileup += ' -r {}'.format(region)
        cmd_query = "{} query -f '%CHROM\\t%POS\\t%REF,%ALT\\t%AD\\n' -i 'SUM(AD)<={} & SUM(AD)>={}'".format(
            self.bcftools, self.dp, self.mincov
        )
        if self.E:
            cmd_mpileup += ' -E'
        errname = os.path.join(self.outdir, '{}_{}_bcftools.log'.format(samplename, chromosome))

        with open(errname, 'w') as err:
            mpileup = pr.Popen(shlex.split(cmd_mpileup), stdout=pr.PIPE, stderr=err)
            query = pr.Popen(
                shlex.split(cmd_query),
                stdin=mpileup.stdout,
                stdout=pr.PIPE,
                stderr=err,
                universal_newlines=True,
            )
            mpileup.stdout.close()   # so that mpileup gets SIGPIPE if query exits early
            try:
                for line in query.stdout:
                    line = line.rstrip('\n')
                    if line != '':
                        p = line.strip().split()
                        yield (
                            samplename,
                            p[0],
                            p[1],
                            [(a, int(c)) for a, c in zip(p[2].split(','), p[3].split(',')) if a in ALLELES],
                        )
            finally:
                query.stdout.close()
                codes = [mpileup.wait(), query.wait()]
        if any(c != 0 for c in codes):
            raise ValueError(
                error(
                    'Allele counting failed on {} of {}, please check errors in {}!'.format(
                        chromosome, samplename, errname
                    )
                )
            )
        else:
            os.remove(errname)


class HetSiteCounter(AlleleCounter):
    """
    Returns ((sample, chromosome of the first site), HetSites or None) for the heterozygous SNPs of each task, or None
    for a task without sites.
    """

    def __init__(self, gamma, maxshift, **counter_args):
        super().__init__(**counter_args)
        self.gamma = gamma
        self.maxshift = maxshift

    def work(self, bamfile, samplename, chromosome):
        rows = self.iterAlleles(bamfile=bamfile, samplename=samplename, chromosome=chromosome)
        first = next(rows, None)
        if first is None:
            return None
        candidates = het_candidates(_chain(first, rows))
        return (first[0], first[1]), select_het_sites(*candidates, gamma=self.gamma, maxshift=self.maxshift)


class TumorAlleleCounter(AlleleCounter):
    """
    Returns ((sample, chromosome of the first site), TumorCounts) with the read counts of the two alleles of the
    normal's heterozygous SNPs (het_sites: chromosome -> HetSites) for each task, or None for a task without sites.
    A position reported more than once keeps its place but takes the counts of its last record.
    """

    def __init__(self, het_sites, **counter_args):
        super().__init__(**counter_args)
        self.het_sites = het_sites

    def work(self, bamfile, samplename, chromosome):
        rows = self.iterAlleles(bamfile=bamfile, samplename=samplename, chromosome=chromosome)
        first = next(rows, None)
        if first is None:
            return None
        het = self.het_sites[first[1]]
        alleles = dict(zip(het.pos.tolist(), zip(het.ref.tolist(), het.alt.tolist())))

        counts = {}
        for r in _chain(first, rows):
            snp = int(r[2])
            ref, alt = alleles[snp]
            observed = dict(r[3])
            counts[snp] = (observed.get(ref, 0), observed.get(alt, 0))

        n = len(counts)
        return (first[0], first[1]), TumorCounts(
            np.fromiter(counts.keys(), dtype=np.int64, count=n),
            np.fromiter((c[0] for c in counts.values()), dtype=np.int64, count=n),
            np.fromiter((c[1] for c in counts.values()), dtype=np.int64, count=n),
        )


def _chain(first, rest):
    yield first
    yield from rest


if __name__ == '__main__':
    main()
