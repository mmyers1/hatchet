import os
import os.path
import shlex
import subprocess as pr
from multiprocessing.pool import ThreadPool

from hatchet.utils.ArgParsing import parse_genotype_snps_arguments
from hatchet.utils.Supporting import log, logArgs, error, indexed_region, MPILEUP_EXCLUDE_FLAGS
from hatchet.utils.multiprocessing import Worker


def main(args=None):
    log(
        msg=(
            '# Parsing the input arguments, checking the consistency of given files, and extracting required ',
            'information\n',
        ),
        level='STEP',
    )
    args = parse_genotype_snps_arguments(args)
    logArgs(args, 80)

    log(msg='# Inferring SNPs from the normal sample\n', level='STEP')
    snps = call(
        bcftools=args['bcftools'],
        reference=args['reference'],
        samples=[args['normal']],
        chromosomes=args['chromosomes'],
        num_workers=args['j'],
        q=args['q'],
        Q=args['Q'],
        mincov=args['mincov'],
        dp=args['maxcov'],
        E=args['E'],
        snplist=args['snps'],
        outdir=args['outputsnps'],
        verbose=args['verbose'],
    )

    log(msg='# Counting number of identified SNPs\n', level='STEP')

    def count(f):
        cmd_bcf = "{} query -f '%CHROM,%POS\n' {}".format(args['bcftools'], f)
        cmd_wcl = 'wc -l'
        bcf = pr.Popen(
            shlex.split(cmd_bcf),
            stdout=pr.PIPE,
            stderr=pr.PIPE,
            universal_newlines=True,
        )
        number = pr.Popen(
            shlex.split(cmd_wcl),
            stdin=bcf.stdout,
            stdout=pr.PIPE,
            stderr=pr.PIPE,
            universal_newlines=True,
        ).communicate()[0]
        number = ''.join([_l for _l in number if _l.isdigit()])
        return int(number) if len(number) > 0 else 0

    with ThreadPool(max(1, min(args['j'], len(snps)))) as p:
        number_snps = sum(p.map(count, snps))

    if number_snps == 0:
        raise ValueError(error('No SNPs found in the normal!\n'))
    else:
        log(
            msg='{} SNPs have been identified in total\n'.format(number_snps),
            level='INFO',
        )

    log(msg='# SNP Calling is concluded\n', level='STEP')
    log(
        msg='## Called SNPs have been written per chromosome in:\n{}\n'.format('\n'.join(snps)),
        level='INFO',
    )


def call(
    bcftools,
    reference,
    samples,
    chromosomes,
    num_workers,
    q,
    Q,
    mincov,
    dp,
    E,
    outdir,
    snplist=None,
    verbose=False,
):
    work = [(bam[0], bam[1], chro) for bam in samples for chro in chromosomes]
    caller = Caller(bcftools, reference, q, Q, mincov, dp, E, outdir, snplist, verbose)
    return sorted(caller.run(work=work, n_instances=num_workers))


class Caller(Worker):
    def __init__(
        self,
        bcftools,
        reference,
        q,
        Q,
        mincov,
        dp,
        E,
        outdir,
        snplist,
        verbose,
    ):
        self.bcftools = bcftools
        self.reference = reference
        self.q = q
        self.Q = Q
        self.mincov = mincov
        self.dp = dp
        self.E = E
        self.outdir = outdir
        self.snplist = snplist
        self.verbose = verbose

    def work(self, bamfile, samplename, chromosome):
        return self.callSNPs(bamfile=bamfile, samplename=samplename, chromosome=chromosome)

    def callSNPs(self, bamfile, samplename, chromosome):
        errname = os.path.join(self.outdir, '{}_{}_bcftools.log'.format(samplename, chromosome))

        outfile = os.path.join(self.outdir, '{}.vcf.gz'.format(chromosome))

        if self.snplist is not None:
            cmd_tgt = "{} query -f '%CHROM\t%POS\n' -r {} {}".format(self.bcftools, chromosome, self.snplist)
            # temporary file, removed after use: fast compression
            cmd_gzip = 'gzip -1 -'
            tgtfile = os.path.join(self.outdir, 'target_{}.pos.gz'.format(chromosome))
            with open(tgtfile, 'w') as tout, open(errname, 'w') as err:
                tgt = pr.Popen(shlex.split(cmd_tgt), stdout=pr.PIPE, stderr=err)
                gzip = pr.Popen(shlex.split(cmd_gzip), stdin=tgt.stdout, stdout=tout, stderr=err)
                tgt.stdout.close()
                codes = [tgt.wait(), gzip.wait()]
            if any(c != 0 for c in codes):
                raise ValueError(
                    error(
                        'SNP Calling failed on {} of {}, please check errors in {}!'.format(
                            chromosome, samplename, errname
                        )
                    )
                )
            else:
                os.remove(errname)

        cmd_mpileup = '{} mpileup {} -Ou -f {} --skip-indels -a INFO/AD,AD,DP --ff {} -q {} -Q {} -d {}'.format(
            self.bcftools, bamfile, self.reference, MPILEUP_EXCLUDE_FLAGS, self.q, self.Q, self.dp
        )

        if self.snplist is not None:
            assert os.path.isfile(tgtfile)
            cmd_mpileup += ' -T {}'.format(tgtfile)
            # With only -T, mpileup streams the whole BAM; restricting it to the chromosome makes it jump there instead
            region = indexed_region(bamfile, chromosome)
            if region is not None:
                cmd_mpileup += ' -r {}'.format(region)
        else:
            cmd_mpileup += ' -r {}'.format(chromosome)
        if self.E:
            cmd_mpileup += ' -E'
        cmd_call = '{} call -mv -Ou'.format(self.bcftools)
        cmd_filter = "{} view -i 'FMT/DP>={}' -Oz -o {}".format(self.bcftools, self.mincov, outfile)

        with open(errname, 'w') as err:
            mpileup = pr.Popen(shlex.split(cmd_mpileup), stdout=pr.PIPE, stderr=err)
            call = pr.Popen(shlex.split(cmd_call), stdin=mpileup.stdout, stdout=pr.PIPE, stderr=err)
            mpileup.stdout.close()
            filter = pr.Popen(shlex.split(cmd_filter), stdin=call.stdout, stderr=err)
            call.stdout.close()
            pcss = [mpileup, call, filter]
            codes = [p.wait() for p in pcss]
        if any(c != 0 for c in codes):
            raise ValueError(
                error('SNP Calling failed on {} of {}, please check errors in {}!').format(
                    chromosome, samplename, errname
                )
            )
        else:
            os.remove(errname)
            if self.snplist is not None:
                os.remove(tgtfile)
        return outfile


if __name__ == '__main__':
    main()
