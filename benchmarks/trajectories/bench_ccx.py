"""Throughput of the quantum-jump solver with one process pinned per L3 block (EPYC 7642: 16 blocks x 3 cores)."""
import os, sys, json, time
from multiprocessing import Process, Queue
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
cfg = json.load(open(sys.argv[1])); T = float(sys.argv[2]); B = int(sys.argv[3])
allowed = sorted(os.sched_getaffinity(0))
groups = {}
for cpu in allowed:
    key = open(f"/sys/devices/system/cpu/cpu{cpu}/cache/index3/shared_cpu_list").read().strip()
    groups.setdefault(key, []).append(cpu)
picks = [g[0] for g in groups.values()]
print("allowed", len(allowed), "L3 blocks in allocation", len(picks), flush=True)
def work(cpu, seed, q):
    os.sched_setaffinity(0, {cpu})
    from traj_qjump import run_batch
    t0 = time.time()
    res = run_batch(cfg['params'], cfg['cut'], T, cfg['dt'], B, seed, backend='numpy', t_burn=5.0)
    q.put((cpu, time.time() - t0, B * (T - 5.0)))
def run(cpus):
    q = Queue(); ps = [Process(target=work, args=(c, 900 + i, q)) for i, c in enumerate(cpus)]
    t0 = time.time(); [p.start() for p in ps]; out = [q.get() for _ in ps]; [p.join() for p in ps]
    wall = time.time() - t0; tt = sum(o[2] for o in out)
    return dict(n=len(cpus), wall=wall, traj_time=tt, rate_total=tt / wall, rate_per_proc=[round(o[2]/o[1], 3) for o in out])
res = {"single": run(picks[:1]), "one_per_L3": run(picks)}
if len(picks) >= 2:
    two = [g[:2] for g in groups.values() if len(g) >= 2]
    res["two_per_L3_first4blocks"] = run([c for g in two[:4] for c in g])
print(json.dumps(res, indent=1), flush=True)
json.dump(res, open("bench_ccx.json", "w"), indent=1)

