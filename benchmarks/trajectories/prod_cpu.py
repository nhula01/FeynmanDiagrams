"""CPU contribution to a GPU production run: writes cpubatch_XXXX.npz (distinct seeds) into the same directory.
usage: prod_cpu.py config.json outdir nproc seed_offset nbatch
Same params/cut/T_total/dt as the GPU config; batch size B from the config's 'B_cpu' (default 8). Restartable (skips existing files)."""
import sys, os, json, time
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from traj_qjump import run_batch
def work(a):
    i, cfg, outdir = a
    final = os.path.join(outdir, f"cpubatch_{i:04d}.npz")
    if os.path.exists(final): return i, 0.0
    tmp = os.path.join(outdir, f"tmpcpu_{i:04d}.npz")
    res = run_batch(cfg['params'], cfg['cut'], cfg['T_total'], cfg['dt'], cfg.get('B_cpu', 8), cfg['seed_cpu'] + i,
                    backend='numpy', out=tmp)
    os.replace(tmp, final); return i, res['wall']
if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    cfg = json.load(open(sys.argv[1])); outdir = sys.argv[2]; nproc = int(sys.argv[3])
    cfg['seed_cpu'] = int(sys.argv[4]); nb = int(sys.argv[5]); os.makedirs(outdir, exist_ok=True); t0 = time.time()
    with Pool(nproc) as pool:
        for k, (i, w) in enumerate(pool.imap_unordered(work, [(i, cfg, outdir) for i in range(nb)])):
            print(f"cpubatch {i} done ({k+1}/{nb}) wall={time.time()-t0:.0f}s per-batch={w:.0f}s", flush=True)
    print("ALL DONE", time.time() - t0, flush=True)

