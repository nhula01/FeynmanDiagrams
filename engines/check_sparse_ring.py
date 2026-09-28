"""Validate sparse_ring.SparseRing against multimode.Ring with the pruning of ring_trunc.series_pruned,
and time both.  usage: python3 check_sparse_ring.py [K] [Nmax]"""
import sys, time, numpy as np
import multimode as MM
from sparse_ring import SparseRing, degree

K = int(sys.argv[1]) if len(sys.argv) > 1 else 3
NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 6
P = dict(kappa=1.0, Delta=-1.0, J=0.4, eta=1.0)

def series_pruned(r, O, Nmax, drop=2):
    vals = [O.get(((), ()), 0)]; X = O
    for N in range(1, Nmax+1):
        X = r.V(r.G0(X)); lim = drop*(Nmax-N)
        X = {k: v for k, v in X.items() if degree(k) <= lim and v != 0}
        vals.append(X.get(((), ()), 0))
    return np.array(vals)

R = SparseRing(K=K, **P)
b, bd = R.site_ops(0)
obs = {"n": MM.mul(bd, b), "n2": MM.mul(MM.mul(bd, bd), MM.mul(b, b)), "n3": MM.mul(MM.mul(MM.mul(bd, bd), bd), MM.mul(MM.mul(b, b), b))}
Us = [0.05, 0.15, 0.3]
worst = 0.0
for name, O in obs.items():
    t = time.time(); c, counts = R.coefficients(O, NMAX, return_counts=True); ts = time.time()-t
    Ssp = R.partial_sums(c, Us)
    t = time.time(); Sref = np.array([series_pruned(MM.Ring(K=K, U=U, **P), O, NMAX) for U in Us]).cumsum(axis=1); tr = time.time()-t
    err = np.max(np.abs(Ssp - Sref)/np.maximum(1e-300, np.abs(Sref)))
    worst = max(worst, err)
    print(f"{name:3s} K={K} Nmax={NMAX}: sparse {ts:6.2f} s (all U; {len(c)} orders, labels per layer {counts})   "
          f"multimode {tr:6.2f} s for {len(Us)} U   max rel diff {err:.1e}")
print("worst relative difference:", f"{worst:.1e}", "OK" if worst < 1e-10 else "MISMATCH")
