#!/bin/bash
# hg38_phase_check.sh <workdir>
# phase-snps on hg38 input (-V hg38 -N), old code then new code, and the comparison of their outputs.
# This covers the liftover path (hg38 -> hg19 panel -> hg38), which the data/ BAMs (hg19) do not reach.
# About 25 minutes on 1 core (shapeit ~8 min per run, plus picard); picard runs with -Xmx8g, so ask for ~12 GB.
#
# The hg38 input is the chr22 test VCF (tests/data/vw/snps/chr22.vcf.gz, hg19) lifted over to hg38 with picard.
set -euo pipefail
source /data1/shahs3/users/myersm2/repos/hatchet/benchmarks/2026-09-24-preprocessing-efficiency/env.sh
W=$1
mkdir -p $W/ref $W/input

# hg38 reference with the sequence dictionary picard needs, made here rather than next to the shared FASTA
ln -sf $HG38 $W/ref/hg38.fa
ln -sf $HG38.fai $W/ref/hg38.fa.fai
[ -s $W/ref/hg38.dict ] || samtools dict $W/ref/hg38.fa -o $W/ref/hg38.dict

# hg38 SNPs; phase-snps takes the chromosome from the file name, so it has to be chr22.vcf.gz
if [ ! -s $W/input/chr22.vcf.gz ]; then
    picard -Xmx8g LiftoverVcf -I $REPO/tests/data/vw/snps/chr22.vcf.gz -O $W/input/lifted.vcf.gz \
        -CHAIN $REFPANELDIR/hg19ToHg38.over.chain.gz -R $W/ref/hg38.fa -REJECT $W/input/rejected.vcf.gz \
        --WARN_ON_MISSING_CONTIG true > $W/input/liftover.log 2>&1
    bcftools view -t chr22 -Oz -o $W/input/chr22.vcf.gz $W/input/lifted.vcf.gz
    bcftools index -t $W/input/chr22.vcf.gz
fi
echo "hg38 input: $(bcftools view -H $W/input/chr22.vcf.gz | wc -l) records"

for version in old new; do
    out=$W/$version
    rm -rf $out
    mkdir -p $out/cwd
    (
        cd $out/cwd   # phase-snps deletes shapeit* logs in its working directory
        PYTHONPATH=$(code_src $version $W) measure phase_hg38 $out $PY -m hatchet phase-snps -D $REFPANELDIR \
            -g $W/ref/hg38.fa -V hg38 -N -o $out/phase -L $W/input/chr22.vcf.gz -j 1
    )
done

bash $HERE/compare.sh $W/old $W/new
