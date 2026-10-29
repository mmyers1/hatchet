#!/bin/bash
# compare.sh <old outdir> <new outdir>
# Compares every output file of two runs: VCF records and #CHROM line (the ## header lines hold commands and dates),
# other gzipped files after decompression, everything else byte for byte. Skips logs, indexes and shapeit's
# *_alignments.log (paths, date, random seed of its -check step). Then prints wall/CPU time and peak memory per step.
# Exit status 1 if any file differs or is missing.
source /data1/shahs3/users/myersm2/repos/hatchet/benchmarks/2026-09-24-preprocessing-efficiency/env.sh
A=$1
B=$2
status=0
while IFS= read -r f; do
    rel=${f#$A/}
    case $rel in
        log_* | */log_* | time_*.json | *.tbi | *.csi | old_code/* | */cwd/* | *_alignments.log) continue ;;
    esac
    if [ ! -e "$B/$rel" ]; then
        echo "MISSING    $rel"
        status=1
        continue
    fi
    case $rel in
        *.vcf.gz) ha=$(zcat "$A/$rel" | grep -v '^##' | md5sum) hb=$(zcat "$B/$rel" | grep -v '^##' | md5sum) ;;
        *.gz) ha=$(zcat "$A/$rel" | md5sum) hb=$(zcat "$B/$rel" | md5sum) ;;
        *) ha=$(md5sum < "$A/$rel") hb=$(md5sum < "$B/$rel") ;;
    esac
    if [ "$ha" == "$hb" ]; then echo "IDENTICAL  $rel"; else echo "DIFFERENT  $rel"; status=1; fi
done < <(find "$A" -type f | sort)

$PY - "$A" "$B" <<'EOF'
import glob, json, os, sys
a, b = sys.argv[1:3]
print('\n| step | wall old -> new | CPU old -> new | peak RSS, process tree | peak RSS, python | largest process |')
print('|---|---|---|---|---|---|')
for fa in sorted(glob.glob(os.path.join(a, 'time_*.json'))):
    step = os.path.basename(fa)[5:-5]
    fb = os.path.join(b, os.path.basename(fa))
    if not os.path.exists(fb):
        continue
    x, y = json.load(open(fa)), json.load(open(fb))
    gb = lambda d, k: d.get(k, float('nan')) / 1024
    top = lambda d: max(d.get('peak_single_rss_mb_by_program', {'?': 0}).items(), key=lambda kv: kv[1])
    print(
        f"| {step} | {x['wall_s']:.0f} -> {y['wall_s']:.0f} s | {x['cpu_s']:.0f} -> {y['cpu_s']:.0f} s "
        f"| {gb(x, 'peak_tree_rss_mb'):.2f} -> {gb(y, 'peak_tree_rss_mb'):.2f} GB "
        f"| {gb(x, 'peak_python_rss_mb'):.2f} -> {gb(y, 'peak_python_rss_mb'):.2f} GB "
        f"| {top(x)[0]} {top(x)[1] / 1024:.2f} -> {top(y)[0]} {top(y)[1] / 1024:.2f} GB |"
    )
EOF
exit $status
