# Shared settings for the verification scripts in this directory; source it. Edit paths here.
REPO=/data1/shahs3/users/myersm2/repos/hatchet
HERE=$REPO/benchmarks/2026-09-24-preprocessing-efficiency
MEASURE=$REPO/benchmarks/common/measure.py   # wall/CPU time and peak memory of a command
PY=/home/myersm2/miniforge3/envs/hatchet-test/bin/python   # hatchet-test env: python, pandas, pysam, psutil, pytest
TOOLS=/home/myersm2/miniforge3/envs/hatchet/bin            # samtools, bcftools, tabix, bgzip, mosdepth, shapeit, picard, java
REFPANELDIR=/data1/shahs3/users/myersm2/reference/phasing-panel
HG19=/data1/shahs3/users/myersm2/reference/hg19_UCSC/hg19.fa
HG38=/data1/shahs3/users/myersm2/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa
OLD_COMMIT=1d6ab0b   # last commit before the refactor

export PATH=$TOOLS:$PATH
# count-alleles writes the tumor samples in the iteration order of a set of strings; fix the hash seed so that old and
# new runs write them in the same order
export PYTHONHASHSEED=0

# code_src old|new <workdir>: the source directory to put on PYTHONPATH (the old code is exported from git into
# <workdir>/old_code the first time)
code_src() {
    if [ "$1" == old ]; then
        [ -d "$2/old_code/src" ] || { mkdir -p "$2/old_code" && git -C $REPO archive $OLD_COMMIT src | tar -x -C "$2/old_code"; }
        echo "$2/old_code/src"
    else
        echo "$REPO/src"
    fi
}

# measure <name> <outdir> <command...>: run a command, saving its stdout/stderr and its wall/CPU time and peak memory
# (<outdir>/time_<name>.json)
measure() {
    local name=$1 out=$2
    shift 2
    echo "== $name"
    $PY $MEASURE $out/time_$name.json $out/log_$name.out $out/log_$name.err -- "$@"
}
