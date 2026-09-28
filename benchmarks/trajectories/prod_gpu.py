"""Restartable GPU production driver for the eight-site counting benchmark.
usage: prod_gpu.py config.json outdir
Each batch is written atomically to outdir/batch_XXXX.npz; batches already on disk are skipped."""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from traj_qjump import run_batch
cfg = json.load(open(sys.argv[1])); outdir = sys.argv[2]; os.makedirs(outdir, exist_ok=True)
p = cfg['params']; t0 = time.time()
for i in range(cfg['nbatch']):
    final = os.path.join(outdir, f"batch_{i:04d}.npz")
    if os.path.exists(final):
        print(f"skip {final}", flush=True); continue
    tmp = os.path.join(outdir, f"tmp_{i:04d}.npz")
    res = run_batch(p, cfg['cut'], cfg['T_total'], cfg['dt'], cfg['B'], cfg['seed0'] + i,
                    backend='cupy', out=tmp, verbose=(i == 0))
    os.replace(tmp, final)
    print(f"batch {i+1}/{cfg['nbatch']} done wall={time.time()-t0:.0f}s per-batch={res['wall']:.0f}s D={res['D']}", flush=True)
print("ALL DONE", time.time() - t0, flush=True)

