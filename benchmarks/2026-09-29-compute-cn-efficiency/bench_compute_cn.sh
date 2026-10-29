#!/bin/bash
#SBATCH --job-name=hatchet-cn-bench
#SBATCH --partition=componc_cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=hatchet-cn-bench-%j.log
# bench_compute_cn.sh: runs `hatchet compute-cn` on one BBC/SEG pair with the code before the solver-efficiency changes
# (commit OLD_COMMIT) and with the current working tree, in one job (same node), then prints a comparison of wall/CPU
# time, peak memory, coordinate-descent step times, objectives per ploidy and clone count, and the chosen solution.
#
# Every setting below can be overridden from the environment. Check them first (about 2 minutes, fine on a login node:
# builds both native solvers and validates the exact compute-cn arguments without solving), then submit from this
# directory:
#   PREFLIGHT_ONLY=1 bash bench_compute_cn.sh
#   sbatch bench_compute_cn.sh
# Another dataset, or a quick smoke test:
#   sbatch --export=ALL,BBC=/path/to/bbc/bulk.bbc bench_compute_cn.sh
#   sbatch --export=ALL,SEEDS=4,CLONES=3 bench_compute_cn.sh
# The Pyomo Gurobi backend (the solver in IID_H134708.ini), or both solving modes:
#   sbatch --export=ALL,SOLVER=gurobi bench_compute_cn.sh
#   sbatch --export=ALL,MODES="2 0" bench_compute_cn.sh
# To summarize finished runs again: $PY compare_compute_cn.py $WORK $SOLVER
set -eo pipefail
# lmod's init scripts reference unset variables, so load the module before `set -u`
source /etc/profile.d/modules.sh
module load gurobi
set -u

REPO=/data1/shahs3/users/myersm2/repos/hatchet
HERE=$REPO/benchmarks/2026-09-29-compute-cn-efficiency
PY=/home/myersm2/miniforge3/envs/hatchet-test/bin/python   # hatchet-test env: python, pandas, pyomo, psutil
MEASURE=$REPO/benchmarks/common/measure.py   # wall/CPU time and peak memory of a command
OLD_COMMIT=4baaf666   # solver code before the efficiency changes (same solver as 1d6ab0b)

# ---- settings
# The clustered bins; the SEG file with the same prefix (bulk.seg next to bulk.bbc) must exist too
BBC=${BBC:-/data1/shahs3/users/myersm2/repos/osteo/results/hatchet2/IID_H134708/bbc/bulk.bbc}
# One directory per dataset; each run gets a fresh <WORK>/mode<M>/<version>
WORK=${WORK:-/data1/shahs3/users/myersm2/hatchet-cn-bench/$(basename "$(dirname "$(dirname "$BBC")")")}
SOLVER=${SOLVER:-cpp}             # cpp (the packaged default) or a Pyomo solver such as gurobi
CLONES=${CLONES:-3,5}             # IID_H134708.ini: clones=3,5
EXTRA_ARGS=${EXTRA_ARGS:-"-eD 9 -eT 15"}   # IID_H134708.ini: diploidcmax=9, tetraploidcmax=15
SEEDS=${SEEDS:-100}               # the packaged default is 400; fewer seeds keep the old code's run manageable
JOBS=${JOBS:-${SLURM_CPUS_PER_TASK:-8}}
RANDOMSEED=${RANDOMSEED:-6700}
MODES=${MODES:-2}                 # 2 = coordinate descent only (the default), 0 = CD + exact ILP, 1 = exact ILP only
TIMELIMIT=${TIMELIMIT:-}          # seconds per ILP solve (and per CD step), empty for none; set it when MODES has 0/1
VERSIONS=${VERSIONS:-"old new"}
PREFLIGHT_ONLY=${PREFLIGHT_ONLY:-}

# ---- checks
fail() {
    echo "ERROR: $*" >&2
    exit 1
}
case $BBC in
    *.bbc) INPUT=${BBC%.bbc} ;;
    *) fail "BBC must name a .bbc file: $BBC" ;;
esac
[ -f "$BBC" ] || fail "BBC file not found: $BBC"
[ -f "$INPUT.seg" ] || fail "SEG file not found: $INPUT.seg (compute-cn reads <prefix>.seg and <prefix>.bbc)"
[ -n "${GRB_LICENSE_FILE:-}" ] || fail "GRB_LICENSE_FILE is not set after 'module load gurobi'"
if [ "$SOLVER" == gurobi ]; then
    command -v gurobi_cl > /dev/null || fail "gurobi_cl is not on PATH (the Pyomo gurobi solver needs it)"
fi
$PY -c 'import pandas, pyomo, psutil' || fail "$PY lacks pandas, pyomo or psutil"
mkdir -p "$WORK"
echo "host: $(hostname), BBC: $BBC, WORK: $WORK, SOLVER: $SOLVER, CLONES: $CLONES, SEEDS: $SEEDS, JOBS: $JOBS," \
    "MODES: $MODES, TIMELIMIT: ${TIMELIMIT:-none}, EXTRA_ARGS: $EXTRA_ARGS"

# prepare <version>: sets SRC (the Python sources to put on PYTHONPATH) and SOLVE_BIN (the native solver, cpp only).
# The old code is exported from git into <WORK>/code_old; the new code is the working tree. Native solvers are built
# out of source into <WORK>/build_<version> (incremental after the first build).
prepare() {
    local version=$1 root
    if [ "$version" == old ]; then
        root=$WORK/code_old
        if [ ! -d "$root/src" ]; then
            mkdir -p "$root"
            git -C $REPO archive $OLD_COMMIT CMakeLists.txt FindGUROBI.cmake src | tar -x -C "$root"
        fi
    else
        root=$REPO
    fi
    SRC=$root/src
    SOLVE_BIN=
    if [ "$SOLVER" == cpp ]; then
        local build=$WORK/build_$version
        {
            cmake -S "$root" -B "$build" -DCMAKE_BUILD_TYPE=Release -DGUROBI_HOME="$(dirname "$GUROBI_HOME")" &&
                cmake --build "$build" -j "$JOBS"
        } > "$build.log" 2>&1 || fail "building the $version solver failed, see $build.log"
        SOLVE_BIN=$build/solve
    fi
}

# compute_cn_args <mode> <rundir>: the compute-cn arguments of one run
compute_cn_args() {
    local args=($SOLVE_BIN -i "$INPUT" -x "$2" -n "$CLONES" -p "$SEEDS" -j "$JOBS" -r "$RANDOMSEED" --mode "$1" -v 2)
    [ -n "$TIMELIMIT" ] && args+=(-s "$TIMELIMIT")
    echo "${args[@]} $EXTRA_ARGS"
}

# Runs happen inside the run directory: compute-cn reads a hatchet.ini from its working directory, and the native
# solver appends its per-step times to ./coordinate_descence.tsv. HATCHET_COMPUTE_CN_SOLVER overrides any ini value.
export HATCHET_COMPUTE_CN_SOLVER=$SOLVER
export PYTHONHASHSEED=0

for version in $VERSIONS; do
    prepare $version
    # validate the exact arguments with compute-cn's own parser (files, running directory, solver path, ranges)
    check=$WORK/preflight/$version
    rm -rf "$check"
    mkdir -p "$check"
    (cd "$check" && PYTHONPATH=$SRC $PY -c 'import sys; from hatchet.bin.HATCHet import parsing_arguments
parsing_arguments(sys.argv[1:])' $(compute_cn_args 2 "$check")) > "$check/parse.log" 2>&1 ||
        fail "compute-cn rejected the arguments for $version, see $check/parse.log"
    echo "$version: code $SRC${SOLVE_BIN:+, solver $SOLVE_BIN}; arguments OK"
done
[ -n "$PREFLIGHT_ONLY" ] && exit 0

status=0
for mode in $MODES; do
    for version in $VERSIONS; do
        prepare $version
        run=$WORK/mode$mode/$version
        rm -rf "$run"
        mkdir -p "$run"
        echo "== mode $mode, $version: $PY -m hatchet compute-cn $(compute_cn_args $mode "$run")"
        (cd "$run" && PYTHONPATH=$SRC $PY $MEASURE "$run/time.json" "$run/log.out" "$run/log.err" -- \
            $PY -m hatchet compute-cn $(compute_cn_args $mode "$run")) || {
            status=$?
            echo "mode $mode, $version failed (exit $status); last lines of $run/log.err:"
            tail -20 "$run/log.err"
        }
    done
done

$PY $HERE/compare_compute_cn.py "$WORK" "$SOLVER" | tee "$WORK/summary.md"
exit $status
