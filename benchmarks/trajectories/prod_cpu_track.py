"""Restartable CPU driver: independent batches of quantum-jump trajectories spread over a process pool.
usage: prod_cpu_track.py config.json outdir nproc
Each batch is written atomically to outdir/batch_XXXX.npz; batches already on disk are skipped, so the
job can be cancelled and resubmitted at any time without losing finished work."""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "1"
import sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multiprocessing import Pool
from traj_qjump import run_batch

def work(args):
    cfg, outdir, i = args
    final = os.path.join(outdir, f"batch_{i:04d}.npz"); tmp = os.path.join(outdir, f"tmp_{i:04d}.npz")
    res = run_batch(cfg['params'], cfg['cut'], cfg['T_total'], cfg['dt'], cfg['B'], cfg['seed0'] + i,
                    backend='numpy', out=tmp, t_burn=cfg.get('t_burn', 40.0))
    os.replace(tmp, final)
    return i, res['wall'], res['wall_jumps'], res['D']

if __name__ == "__main__":
    cfg = json.load(open(sys.argv[1])); outdir = sys.argv[2]; nproc = int(sys.argv[3])
    os.makedirs(outdir, exist_ok=True)
    todo = [i for i in range(cfg['nbatch']) if not os.path.exists(os.path.join(outdir, f"batch_{i:04d}.npz"))]
    print(f"{len(todo)} of {cfg['nbatch']} batches to run on {nproc} processes", flush=True)
    t0 = time.time()
    with Pool(nproc) as pool:
        for k, (i, w, wj, D) in enumerate(pool.imap_unordered(work, [(cfg, outdir, i) for i in todo])):
            print(f"batch {i} done ({k+1}/{len(todo)}) wall={time.time()-t0:.0f}s per-batch={w:.0f}s jumps={wj:.0f}s D={D}", flush=True)
    print("ALL DONE", time.time() - t0, flush=True)

