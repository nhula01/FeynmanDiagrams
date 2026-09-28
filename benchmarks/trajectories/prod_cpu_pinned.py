"""Pinned CPU contribution: one process per L3 block, each writing cpubatch_XXXX.npz (restartable).
usage: prod_cpu_pinned.py config.json outdir seed_offset nbatch B"""
import os, sys, json, time
from multiprocessing import Process
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
cfg = json.load(open(sys.argv[1])); outdir = sys.argv[2]; seed0 = int(sys.argv[3]); nb = int(sys.argv[4]); B = int(sys.argv[5])
os.makedirs(outdir, exist_ok=True)
groups = {}
for cpu in sorted(os.sched_getaffinity(0)):
    key = open(f"/sys/devices/system/cpu/cpu{cpu}/cache/index3/shared_cpu_list").read().strip()
    groups.setdefault(key, []).append(cpu)
picks = [g[0] for g in groups.values()]
print("L3 blocks", len(picks), flush=True)
def worker(k, cpu):
    os.sched_setaffinity(0, {cpu})
    from traj_qjump import run_batch
    for i in range(k, nb, len(picks)):
        final = os.path.join(outdir, f"cpubatch_{i:04d}.npz")
        if os.path.exists(final): continue
        tmp = os.path.join(outdir, f"tmpcpu_{i:04d}.npz"); t0 = time.time()
        run_batch(cfg['params'], cfg['cut'], cfg['T_total'], cfg['dt'], B, seed0 + i, backend='numpy', out=tmp)
        os.replace(tmp, final); print(f"cpubatch {i} done on cpu {cpu} in {time.time()-t0:.0f}s", flush=True)
ps = [Process(target=worker, args=(k, c)) for k, c in enumerate(picks)]
[p.start() for p in ps]; [p.join() for p in ps]; print("ALL DONE", flush=True)

