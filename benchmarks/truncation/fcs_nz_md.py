"""Remaining K=3 maxdeg-convergence points (maxdeg 9 j=2..8, maxdeg 7 all j) in parallel."""
import sys, os, json, time
sys.argv = [sys.argv[0], 'none'] + sys.argv[1:]
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
from multiprocessing import Pool
import fcs_nz as F
if __name__ == '__main__':
    nproc = int(sys.argv[2]); M = 16
    tasks = [("k3", 6, 9, True, 0.1, M, j) for j in range(2, 9)] + [("k3", 6, 7, True, 0.1, M, j) for j in range(9)]
    res = []; t0 = time.time()
    with Pool(nproc) as pool:
        for rr in pool.imap_unordered(F.task, tasks):
            res.append(rr); print(rr[0], f"{rr[2]:.1f}s [{time.time()-t0:.0f}s]", flush=True)
            json.dump([[list(x[0]), x[1], x[2]] for x in res], open('fcs_nz_raw_md.json', 'w'))
    print("done")
