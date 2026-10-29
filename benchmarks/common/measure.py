"""Run a command; report wall time, peak summed RSS of its process tree, and max single-process RSS.

usage: measure.py <report.json> <stdout file> <stderr file> -- cmd args...
Peak RSS is sampled every 0.1 s; summed RSS counts pages shared between forked processes once per process.
"""
import json
import resource
import subprocess
import sys
import time

import psutil

report, out, err = sys.argv[1:4]
assert sys.argv[4] == '--'
cmd = sys.argv[5:]

t0 = time.time()
with open(out, 'w') as fo, open(err, 'w') as fe:
    p = subprocess.Popen(cmd, stdout=fo, stderr=fe)
    root = psutil.Process(p.pid)
    peak_tree = 0
    peak_python = 0   # summed over the python processes of the tree (the hatchet code itself)
    peak_by_name = {}   # largest single process, by program name
    while p.poll() is None:
        try:
            procs = [root] + root.children(recursive=True)
            rss = 0
            rss_python = 0
            for q in procs:
                try:
                    r = q.memory_info().rss
                    name = q.name()
                except psutil.Error:
                    continue
                rss += r
                if name.startswith('python'):
                    rss_python += r
                peak_by_name[name] = max(peak_by_name.get(name, 0), r)
            peak_tree = max(peak_tree, rss)
            peak_python = max(peak_python, rss_python)
        except psutil.Error:
            pass
        time.sleep(0.1)
wall = time.time() - t0
ru = resource.getrusage(resource.RUSAGE_CHILDREN)
res = {
    'cmd': cmd,
    'returncode': p.returncode,
    'wall_s': round(wall, 2),
    'cpu_s': round(ru.ru_utime + ru.ru_stime, 2),
    'peak_tree_rss_mb': round(peak_tree / 2**20, 1),
    'max_single_rss_mb': round(ru.ru_maxrss / 1024, 1),
    'peak_python_rss_mb': round(peak_python / 2**20, 1),
    'peak_single_rss_mb_by_program': {
        k: round(v / 2**20, 1) for k, v in sorted(peak_by_name.items(), key=lambda x: -x[1])
    },
}
with open(report, 'w') as f:
    json.dump(res, f, indent=1)
print(json.dumps({k: v for k, v in res.items() if k not in ('cmd', 'peak_single_rss_mb_by_program')}))
sys.exit(p.returncode)
