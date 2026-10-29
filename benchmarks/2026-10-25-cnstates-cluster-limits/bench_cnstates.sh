#!/bin/bash
#SBATCH --job-name=hatchet-cnstates-bench
#SBATCH --partition=componc_cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=hatchet-cnstates-bench-%A_%a.log
# bench_cnstates.sh: synthetic benchmark of the per-cluster limits on tumor copy-number states (-d/--cnstates with
# --cnstates-exempt or --cluster-cnstates). For each scenario in scenarios.json it simulates inputs with a known truth
# (simulate.py), runs `hatchet compute-cn` once per configuration (free, -d D, -d D with the true or half the hot
# clusters exempt, oracle per-cluster limits, the old pairwise formulation, the pre-change code), and summarizes time,
# memory, model sizes, model selection and accuracy against the truth (summarize_cnstates.py).
#
# One SLURM array task per scenario (index into SCENARIOS). Every setting below can be overridden from the environment.
# Check first (minutes, fine on a login node: simulates, builds both native solvers, validates every argument set with
# compute-cn's parser, no solving), then submit from this directory:
#   PREFLIGHT_ONLY=1 bash bench_cnstates.sh
#   SMOKE=1 sbatch --array=0-5 bench_cnstates.sh       # tiny inputs, 4 seeds: end-to-end check (minutes per task)
#   sbatch --array=0-5 bench_cnstates.sh               # full benchmark
#   sbatch --array=5 --mem=64G bench_cnstates.sh       # 'scale' alone with more memory
# One scenario or configuration, or the Pyomo solver:
#   sbatch --array=0 --export=ALL,CONFIGS="d d_exempt" bench_cnstates.sh
#   sbatch --array=0 --export=ALL,SOLVER=gurobi,MODES=2 bench_cnstates.sh
# To summarize finished runs again: $PY summarize_cnstates.py $WORK > $WORK/summary.md
set -eo pipefail
# lmod's init scripts reference unset variables, so load the module before `set -u`
source /etc/profile.d/modules.sh
module load gurobi
set -u

REPO=/data1/shahs3/users/myersm2/repos/hatchet
HERE=$REPO/benchmarks/2026-10-25-cnstates-cluster-limits
PY=/home/myersm2/miniforge3/envs/hatchet-test/bin/python   # hatchet-test env: python, numpy, pandas, scipy, pyomo, psutil
MEASURE=$REPO/benchmarks/common/measure.py                 # wall/CPU time and peak memory of a command
OLD_COMMIT=9083d8be   # last commit before the per-cluster limits (and the cnstates fixes)

# ---- settings
SMOKE=${SMOKE:-}                   # non-empty: tiny inputs (scenarios.json "smoke") and 4 seeds
SCENARIOS=${SCENARIOS:-"lowhet lowhet_noisy moderate nohot wgd scale"}
WORK=${WORK:-/data1/shahs3/users/myersm2/hatchet-cnstates-bench${SMOKE:+-smoke}}
CONFIGS=${CONFIGS:-}              # empty: each scenario's list in scenarios.json
MODES=${MODES:-}                  # empty: each scenario's modes (2 = coordinate descent, 0 = CD + exact ILP)
TIMELIMIT=${TIMELIMIT:-}          # seconds per ILP solve in modes 0/1; empty: the scenario's timelimit
SEEDS=${SEEDS:-100}               # coordinate-descent seeds (the packaged default is 400)
JOBS=${JOBS:-${SLURM_CPUS_PER_TASK:-8}}
RANDOMSEED=${RANDOMSEED:-6700}
SOLVER=${SOLVER:-cpp}             # cpp (native) or a Pyomo solver such as gurobi
PREFLIGHT_ONLY=${PREFLIGHT_ONLY:-}
[ -n "$SMOKE" ] && SEEDS=4

# ---- checks
fail() {
    echo "ERROR: $*" >&2
    exit 1
}
read -r -a ALL_SCENARIOS <<< "$SCENARIOS"
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    [ "$SLURM_ARRAY_TASK_ID" -lt "${#ALL_SCENARIOS[@]}" ] || fail "array index $SLURM_ARRAY_TASK_ID has no scenario in: $SCENARIOS"
    RUN_SCENARIOS=("${ALL_SCENARIOS[$SLURM_ARRAY_TASK_ID]}")
else
    RUN_SCENARIOS=("${ALL_SCENARIOS[@]}")
fi
[ -n "${GRB_LICENSE_FILE:-}" ] || fail "GRB_LICENSE_FILE is not set after 'module load gurobi'"
$PY -c 'import numpy, pandas, scipy, pyomo, psutil' || fail "$PY lacks numpy, pandas, scipy, pyomo or psutil"
mkdir -p "$WORK"
echo "host: $(hostname), scenarios: ${RUN_SCENARIOS[*]}, WORK: $WORK, SOLVER: $SOLVER, SEEDS: $SEEDS, JOBS: $JOBS," \
    "SMOKE: ${SMOKE:-no}, CONFIGS: ${CONFIGS:-per scenario}, MODES: ${MODES:-per scenario}"

# prepare <version>: sets SRC (the Python sources for PYTHONPATH) and SOLVE_BIN (the native solver; cpp only).
# new = this working tree, old = OLD_COMMIT exported from git, pairwise = the new solver with the old pairwise
# formulation of the state limit (a wrapper adding the hidden --cnstates-formulation option). Builds are shared by the
# array tasks and serialized with a lock.
prepare() {
    local version=$1 root build
    [ "$version" == old ] && build=old || build=new
    if [ "$build" == old ]; then
        root=$WORK/code_old
        (
            flock 9
            if [ ! -d "$root/src" ]; then
                mkdir -p "$root"
                git -C $REPO archive $OLD_COMMIT CMakeLists.txt FindGUROBI.cmake src | tar -x -C "$root"
            fi
        ) 9> "$WORK/.code_old.lock"
    else
        root=$REPO
    fi
    SRC=$root/src
    SOLVE_BIN=
    if [ "$SOLVER" == cpp ]; then
        (
            flock 9
            {
                cmake -S "$root" -B "$WORK/build_$build" -DCMAKE_BUILD_TYPE=Release -DGUROBI_HOME="$(dirname "$GUROBI_HOME")" &&
                    cmake --build "$WORK/build_$build" -j "$JOBS"
            } > "$WORK/build_$build.log" 2>&1
        ) 9> "$WORK/.build_$build.lock" || fail "building the $build solver failed, see $WORK/build_$build.log"
        SOLVE_BIN=$WORK/build_$build/solve
        if [ "$version" == pairwise ]; then
            SOLVE_BIN=$WORK/build_new/solve_pairwise
            printf '#!/bin/bash\nexec %q "$@" --cnstates-formulation pairwise\n' "$WORK/build_new/solve" > "$SOLVE_BIN.tmp.$$"
            chmod +x "$SOLVE_BIN.tmp.$$"
            mv -f "$SOLVE_BIN.tmp.$$" "$SOLVE_BIN"
        fi
    elif [ "$version" == pairwise ]; then
        return 1   # the pairwise formulation exists only in the native solver
    fi
}

# Runs happen inside the run directory: compute-cn reads a hatchet.ini from its working directory, and the native
# solver appends its per-step times and model sizes to ./coordinate_descence.tsv.
export HATCHET_COMPUTE_CN_SOLVER=$SOLVER
export PYTHONHASHSEED=0

status=0
for scenario in "${RUN_SCENARIOS[@]}"; do
    input=$WORK/$scenario/input
    $PY -I $HERE/simulate.py $HERE/scenarios.json "$scenario" "$input" ${SMOKE:+--smoke} || fail "simulating $scenario failed"
    read -r -a common < "$input/common_args.txt"
    { read -r scenario_modes; read -r scenario_timelimit || true; } < "$input/modes.txt"
    modes=${MODES:-$scenario_modes}
    timelimit=${TIMELIMIT:-$scenario_timelimit}

    # runs.tsv: config, version, extra arguments (shell-quoted)
    runs=()
    while IFS=$'\t' read -r config version extra; do
        if [ -n "$CONFIGS" ] && [[ " $CONFIGS " != *" $config "* ]]; then
            continue
        fi
        if [ "$version" == pairwise ] && [ "$SOLVER" != cpp ]; then
            echo "skipping $config: the pairwise formulation exists only in the native solver"
            continue
        fi
        runs+=("$config"$'\t'"$version"$'\t'"$extra")
    done < "$input/runs.tsv"

    for mode in $modes; do
        for run in "${runs[@]}"; do
            IFS=$'\t' read -r config version extra <<< "$run"
            eval "extra_args=($extra)"
            prepare "$version"
            rundir=$WORK/$scenario/mode$mode/$config
            rm -rf "$rundir"
            mkdir -p "$rundir"
            args=($SOLVE_BIN -i "$input/bulk" -x "$rundir" "${common[@]}" -p "$SEEDS" -j "$JOBS" -r "$RANDOMSEED"
                --mode "$mode" -v 2 "${extra_args[@]}")
            if [ "$mode" != 2 ] && [ -n "$timelimit" ]; then
                args+=(-s "$timelimit")
            fi
            printf '%q ' "${args[@]}" > "$rundir/args.txt"
            if [ -n "$PREFLIGHT_ONLY" ]; then
                (cd "$rundir" && PYTHONPATH=$SRC $PY -c 'import sys; from hatchet.bin.HATCHet import parsing_arguments
parsing_arguments(sys.argv[1:])' "${args[@]}") > "$rundir/parse.log" 2>&1 ||
                    fail "compute-cn rejected the arguments of $scenario/mode$mode/$config, see $rundir/parse.log"
                echo "$scenario mode $mode $config ($version): arguments OK"
                continue
            fi
            echo "== $scenario, mode $mode, $config ($version): compute-cn $(cat "$rundir/args.txt")"
            (cd "$rundir" && PYTHONPATH=$SRC $PY $MEASURE "$rundir/time.json" "$rundir/log.out" "$rundir/log.err" -- \
                $PY -m hatchet compute-cn "${args[@]}") || {
                status=$?
                echo "$scenario mode $mode $config failed (exit $status); last lines of $rundir/log.err:"
                tail -20 "$rundir/log.err"
            }
        done
    done
done
if [ -n "$PREFLIGHT_ONLY" ]; then
    echo "pre-flight OK"
    exit 0
fi

# every array task rewrites the combined summary, so the last one to finish covers all scenarios
(
    flock 9
    $PY -I $HERE/summarize_cnstates.py "$WORK" > "$WORK/summary.md.tmp" && mv -f "$WORK/summary.md.tmp" "$WORK/summary.md"
) 9> "$WORK/.summary.lock"
cat "$WORK/summary.md"
exit $status
