#!/bin/bash
#SBATCH --job-name=hatchet-preproc-compare
#SBATCH --partition=componc_cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=hatchet-preproc-compare-%j.log
# wgs_compare.sh: genotype-snps, count-alleles, count-reads and phase-snps on real WGS BAMs with the old code, then the
# new code, in one job (same node), followed by the comparison of every output and a time/memory table.
#
# Edit the settings below, check them (seconds, fine on a login node), then submit from this directory:
#   PREFLIGHT_ONLY=1 bash wgs_compare.sh
#   sbatch wgs_compare.sh
# To resume after a failure without redoing the steps that finished, select versions and steps, e.g.:
#   sbatch --export=ALL,VERSIONS=new,STEPS="reads phase" wgs_compare.sh
# Rough cost for chromosomes 20-22 of a 61 GB normal and two ~145 GB tumors: old code ~3-5 h (every count-alleles task
# streams a whole BAM), new code ~1-1.5 h (mostly mosdepth, the same in both); peak memory ~15 GB (mosdepth ~4 GB per
# sample, 3 at a time).
set -euo pipefail
source /data1/shahs3/users/myersm2/repos/hatchet/benchmarks/2026-09-24-preprocessing-efficiency/env.sh

# ---- settings
WORK=/data1/shahs3/users/myersm2/hatchet-preproc-verify/wgs   # a few GB per sample
NORMAL=/data1/papaemme/isabl/data/analyses/32/73/163273/I-H-135387-N1-1-D1-1.bam
TUMORS="/data1/papaemme/isabl/data/analyses/32/72/163272/I-H-135387-T1-1-D1-1.bam /data1/papaemme/isabl/data/analyses/85/91/408591/IID_H135387_T04_01_WG01.bam"
NAMES="Normal T01 T04"   # normal first, then the tumors in the order above
# The reference the BAMs were aligned to, with the same chromosome names. These BAMs are GRCh37d5 (hs37d5: 1..22, X, Y,
# MT and decoys); GRCh37-lite has the same names, and the same sequences for 20, 21, 22 and X (its MD5s equal the M5
# tags of the BAM headers). The pre-flight checks names, and MD5s where the BAM header has them.
REF=/data1/shahs3/users/myersm2/reference/GRCh37-lite.fa
REFVERSION=hg19
CHRNOT=                  # -N if chromosome names have the chr prefix, empty otherwise
# known SNPs (indexed VCF, chromosome names as in the BAMs); empty to call SNPs without a list
SNPLIST=/data1/shahs3/users/myersm2/reference/snp-lists/hg19_nochr/00-All.vcf.gz
# A few chromosomes are enough to see the per-task cost (with only -T, every old task streams the whole BAM).
# Keep chrX out with CHRNOT=-N: the panel's rename_chrs files have no chrX line (a separate, pre-existing issue).
CHROMS="20 21 22"
J=${SLURM_CPUS_PER_TASK:-8}
VERSIONS=${VERSIONS:-"old new"}
STEPS=${STEPS:-"genotype alleles reads phase"}
# ----

header_contigs() {  # contig names in a BAM header
    samtools view -H "$1" | awk -F'\t' '$1 == "@SQ" { for (i = 2; i <= NF; i++) if ($i ~ /^SN:/) print substr($i, 4) }'
}
header_md5() {  # header_md5 <BAM or .dict> <contig>: the M5 tag of the contig, if any
    grep -P "^@SQ\tSN:$2\t" "$1" 2> /dev/null | grep -oP 'M5:\K[0-9a-f]+' || true
}

# Check that the BAMs, the reference, the SNP list and the settings agree; list every problem found
preflight() {
    local problems=() bam c n contigs ref_contigs normal_contigs bam_m5 ref_m5 missing=""
    local dict=${REF%.*}.dict
    # the samtools and bcftools hatchet will use (HATCHET_PATHS_* override the PATH)
    local samtools=${HATCHET_PATHS_SAMTOOLS:+$HATCHET_PATHS_SAMTOOLS/}samtools
    local bcftools=${HATCHET_PATHS_BCFTOOLS:+$HATCHET_PATHS_BCFTOOLS/}bcftools
    echo "host: $(hostname); $($samtools --version | head -1); $($bcftools --version | head -1)"
    echo "normal: $NORMAL"
    echo "tumors: $TUMORS"
    echo "names: $NAMES; reference: $REF ($REFVERSION); CHRNOT='$CHRNOT'; chromosomes: $CHROMS"
    echo "SNP list: ${SNPLIST:-none}; panel: $REFPANELDIR; work: $WORK; versions: $VERSIONS; steps: $STEPS; J=$J"

    if [ -s $REF.fai ]; then ref_contigs=$(cut -f1 $REF.fai); else problems+=("no index $REF.fai"); ref_contigs=""; fi
    [ $(echo $NAMES | wc -w) -eq $(echo $NORMAL $TUMORS | wc -w) ] ||
        problems+=("NAMES must name the normal and each tumor")
    for bam in $NORMAL $TUMORS; do
        if [ ! -s $bam ]; then
            problems+=("BAM not found: $bam")
            continue
        fi
        [ -s $bam.bai ] || [ -s ${bam%.bam}.bai ] || problems+=("no index (.bai) for $bam")
        contigs=$(header_contigs $bam)
        for c in $CHROMS; do
            if ! grep -qx -- "$c" <<< "$contigs"; then
                problems+=("chromosome $c is not in the header of $bam")
                continue
            fi
            # same sequence as the reference, where both carry an MD5
            bam_m5=$(samtools view -H $bam | header_md5 /dev/stdin $c)
            ref_m5=$(header_md5 $dict $c)
            if [ -n "$bam_m5" ] && [ -n "$ref_m5" ] && [ "$bam_m5" != "$ref_m5" ]; then
                problems+=("chromosome $c of $bam (M5 $bam_m5) differs from the reference's (M5 $ref_m5)")
            fi
        done
    done
    for c in $CHROMS; do
        grep -qx -- "$c" <<< "$ref_contigs" || problems+=("chromosome $c is not in the reference $REF")
    done
    # genotype-snps and count-alleles require every standard chromosome of the normal BAM to be in the reference
    if [ -s $NORMAL ]; then
        normal_contigs=$(header_contigs $NORMAL)
        for c in $(seq 1 22) X Y; do
            for n in $c chr$c; do
                if grep -qx -- "$n" <<< "$normal_contigs" && ! grep -qx -- "$n" <<< "$ref_contigs"; then
                    missing="$missing $n"
                fi
            done
        done
        [ -z "$missing" ] || problems+=("chromosomes of the normal BAM missing from the reference $REF:$missing")
    fi
    if [ -n "$SNPLIST" ]; then
        if [ ! -s $SNPLIST ] || { [ ! -s $SNPLIST.tbi ] && [ ! -s $SNPLIST.csi ]; }; then
            problems+=("SNP list or its index not found: $SNPLIST")
        else
            contigs=$(tabix -l $SNPLIST)
            for c in $CHROMS; do
                grep -qx -- "$c" <<< "$contigs" || problems+=("chromosome $c is not in the SNP list $SNPLIST")
            done
        fi
    fi
    case "$(echo $CHROMS | cut -d' ' -f1)" in
        chr*) [ "$CHRNOT" == -N ] || problems+=("chromosome names have the chr prefix: set CHRNOT=-N") ;;
        *) [ -z "$CHRNOT" ] || problems+=("chromosome names have no chr prefix: leave CHRNOT empty") ;;
    esac
    [ -s $REFPANELDIR/1000GP_Phase3/1000GP_Phase3.sample ] || problems+=("no reference panel in $REFPANELDIR")

    if [ ${#problems[@]} -gt 0 ]; then
        echo "pre-flight FAILED:"
        printf '  - %s\n' "${problems[@]}"
        return 1
    fi
    echo "pre-flight OK"
}

# The chosen chromosomes of the SNP list, with the .tbi genotype-snps requires. The file name depends on the SNP list
# and the chromosomes, so changing either makes a new extract; the .tbi is written last and marks it complete.
snplist_extract() {
    local key tmp got
    key=$(echo "$SNPLIST $CHROMS" | md5sum | cut -c1-8)
    SNPX=$WORK/snplist.$key.vcf.gz
    if [ ! -s $SNPX.tbi ]; then
        echo "== extracting chromosomes $CHROMS from $SNPLIST"
        tmp=$WORK/snplist.$key.tmp.vcf.gz
        bcftools view -r "$(echo $CHROMS | tr ' ' ',')" -Oz -o $tmp $SNPLIST
        tabix -f -p vcf $tmp
        mv $tmp $SNPX
        mv $tmp.tbi $SNPX.tbi
    fi
    got=$(tabix -l $SNPX | sort | tr '\n' ' ')
    if [ "$got" != "$(printf '%s\n' $CHROMS | sort | tr '\n' ' ')" ]; then
        echo "The SNP list extract $SNPX has chromosomes '$got', not '$CHROMS'"
        exit 1
    fi
}

# run_step <old|new> <step>: run one step, replacing only that step's outputs; stop the job if it fails
run_step() {
    local version=$1 step=$2 out=$WORK/$1 dir status=0 vcfs v
    vcfs=$(for c in $CHROMS; do echo $out/snps/$c.vcf.gz; done)
    case $step in
        genotype) dir=snps ;;
        alleles) dir=baf ;;
        reads) dir=rdr ;;
        phase) dir=phase ;;
        *) echo "unknown step: $step" && exit 1 ;;
    esac
    # inputs made by earlier steps
    case $step in
        alleles | phase) for v in $vcfs; do
            [ -s $v ] || { echo "$version $step needs $v: run the genotype step" && exit 1; }
        done ;;
        reads) [ -s $out/baf/bulk.1bed ] || { echo "$version reads needs $out/baf/bulk.1bed: run the alleles step" &&
            exit 1; } ;;
    esac
    rm -rf $out/$dir $out/time_$step.json $out/log_$step.out $out/log_$step.err
    mkdir -p $out/$dir $out/cwd

    case $step in
        genotype)
            measure genotype $out $PY -m hatchet genotype-snps -N $NORMAL -r $REF $SNPARG -o $out/snps -j $J \
                --chromosomes $CHROMS || status=$?
            ;;
        alleles)
            measure alleles $out $PY -m hatchet count-alleles -N $NORMAL -T $TUMORS -S $NAMES -r $REF -j $J -L $vcfs \
                -O $out/baf/normal.1bed -o $out/baf/bulk.1bed -l $out/baf --chromosomes $CHROMS || status=$?
            ;;
        reads)
            measure reads $out $PY -m hatchet count-reads -N $NORMAL -T $TUMORS -S $NAMES -j $J \
                -b $out/baf/bulk.1bed -O $out/rdr -V $REFVERSION --chromosomes $CHROMS || status=$?
            ;;
        phase)
            # phase-snps deletes shapeit* logs in its working directory
            (cd $out/cwd && measure phase $out $PY -m hatchet phase-snps -D $REFPANELDIR -g $REF -V $REFVERSION \
                $CHRNOT -o $out/phase -L $vcfs -j $J) || status=$?
            ;;
    esac
    if [ $status -ne 0 ]; then
        echo "FAILED: $version $step (exit $status). Last lines of $out/log_$step.err:"
        tail -30 $out/log_$step.err | sed 's/\x1b\[[0-9;]*m//g'
        exit 1
    fi
}

finished() {  # finished <old|new>: all four steps ran successfully
    local s
    for s in genotype alleles reads phase; do
        grep -qs '"returncode": 0' $WORK/$1/time_$s.json || return 1
    done
}

preflight
[ -z "${PREFLIGHT_ONLY:-}" ] || exit 0

mkdir -p $WORK
SNPARG=""
if [ -n "$SNPLIST" ] && [[ " $STEPS " == *" genotype "* ]]; then
    snplist_extract
    SNPARG="-R $SNPX"
fi
for version in $VERSIONS; do
    export PYTHONPATH=$(code_src $version $WORK)
    for step in $STEPS; do
        run_step $version $step
    done
done
unset PYTHONPATH

if finished old && finished new; then
    bash $HERE/compare.sh $WORK/old $WORK/new
else
    echo "Not comparing: the old and new runs have not both finished all four steps."
fi
