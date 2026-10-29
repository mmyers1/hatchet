#!/bin/bash
# run_tests.sh <workdir> [old|new] [pytest -k expression]   (default: new, all the tests below)
# Runs the tests that cover count-reads, genotype-snps, count-alleles and phase-snps on the data/ BAMs.
# About 15 minutes and ~5 GB at peak (mosdepth ~4 GB, shapeit ~2 GB); run it on a compute node.
# The full output goes to <workdir>/pytest_<old|new>.log, with the host and the memory limit and peak usage of the
# session's cgroup; the tests' outputs stay in <workdir>/tests_<old|new>/out.
#
# The tests run from a copy of tests/ in <workdir>. Their fixtures delete their output folders (test_steps.py: out/,
# test_steps_vw.py: out/vw, test_phase.py: out/vw/phase), so they run in that order and every output survives.
# The reference panel is the phasing-panel directory, seen through a directory of symlinks, because download-panel
# (which test_phase.py calls first) rewrites the chain and rename files in its refpaneldir.
# Expected with bcftools 1.19: test_genotype_snps and test_count_alleles_normal_snps fail for old and new code alike
# (their md5 / pickle fixtures were made with bcftools 1.7); everything else passes.
# Example, rerunning two tests only: run_tests.sh <workdir> new 'test_count_reads or test_script'
set -euo pipefail
source /data1/shahs3/users/myersm2/repos/hatchet/benchmarks/2026-09-24-preprocessing-efficiency/env.sh
W=$1
VERSION=${2:-new}
SELECT=${3:-}
mkdir -p $W
SRC=$(code_src $VERSION $W)
T=$W/tests_$VERSION
LOG=$W/pytest_$VERSION.log
rm -rf $T $W/panel_$VERSION
mkdir -p $T $W/panel_$VERSION
rsync -a --exclude out --exclude __pycache__ $REPO/tests/ $T/
for f in 1000GP_Phase3 hg19_no_chr.fa hg19_no_chr.dict hg19ToHg38.over.chain.gz hg38ToHg19.over.chain.gz; do
    ln -s $REFPANELDIR/$f $W/panel_$VERSION/$f
done

export HATCHET_DOWNLOAD_PANEL_REFPANELDIR=$W/panel_$VERSION
export HATCHET_TESTS_BAM_DIRECTORY=$REPO/data
export HATCHET_PATHS_REFERENCE=$HG19
TESTS="$T/test_steps.py::test_genotype_snps $T/test_steps.py::test_count_alleles_normal_snps"
TESTS="$TESTS $T/test_steps_vw.py::test_count_reads $T/test_phase.py"
if [ $VERSION == new ]; then
    TESTS="$T/test_counting_equivalence.py $TESTS"
fi

# memory cgroup of this session (SLURM job), to tell whether a killed process hit its limit
CGROOT=/sys/fs/cgroup/memory
CG=$CGROOT$(grep ':memory:' /proc/self/cgroup | cut -d: -f3 | sed 's|/step_.*||')
cgroup_memory() {
    for f in memory.limit_in_bytes memory.max_usage_in_bytes; do
        [ -r $CG/$f ] && awk -v f=$f '{printf "%s: %.1f GB\n", f, $1 / 2^30}' $CG/$f
    done
    [ -r $CG/memory.oom_control ] && grep oom_kill $CG/memory.oom_control || true
}
# the smallest memory limit from this session's cgroup up to the root, in GB (login nodes cap each user's slice)
memory_limit_gb() {
    local d=$CG min=$((1 << 62)) limit
    while [ ${#d} -gt ${#CGROOT} ]; do
        if [ -r $d/memory.limit_in_bytes ]; then
            limit=$(cat $d/memory.limit_in_bytes)
            [ $limit -lt $min ] && min=$limit
        fi
        d=$(dirname $d)
    done
    echo $((min >> 30))
}
{
    echo "host: $(hostname), code: $SRC, cgroup: $CG, memory limit: $(memory_limit_gb) GB"
    cgroup_memory
} > $LOG
if [ $(memory_limit_gb) -lt 8 ] && [ -z "${ALLOW_LOW_MEMORY:-}" ]; then
    echo "Only $(memory_limit_gb) GB of memory on $(hostname) (a login node?): mosdepth needs ~4 GB and shapeit ~2 GB," \
        "and are killed below that. Run this on a compute node, or set ALLOW_LOW_MEMORY=1 to run anyway." | tee -a $LOG
    exit 1
fi

cd $W
status=0
PYTHONPATH=$SRC $PY -m pytest -p no:cacheprovider -rA --rootdir=$T ${SELECT:+-k "$SELECT"} $TESTS 2>&1 | tee -a $LOG ||
    status=$?
{
    echo "after the tests:"
    cgroup_memory
} | tee -a $LOG
exit $status
