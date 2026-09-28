"""Timing of the sparse-matrix ring engine (engines/sparse_ring.py) on one core, K = 3: orders 0..2, 0..4, 0..10 of
n, <(a0^dag)^2 a0^2>, <(a0^dag)^3 a0^3> on site 0, valid at every U (the coefficients are independent of U), and the
per-coupling cost of the symbolic recursion (multimode.Ring with degree pruning) for comparison; then the scaling with
K at second and fourth order for g2.  Writes sparse_ring_timing.json.   usage: OMP_NUM_THREADS=1 python3 sparse_ring_timing.py"""
import sys, os, time, json, platform, numpy as np
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
sys.path.insert(0, os.path.join(_ROOT, 'engines'))
import multimode as MM
from sparse_ring import SparseRing, degree
PR = json.load(open(os.path.join(_ROOT, 'notebook', 'data', 'revision_params.json')))["ring3"]
P = dict(kappa=PR["kappa"], Delta=PR["Delta"], J=PR["J"], eta=PR["eta"])

def series_pruned(r, O, Nmax, drop=2):
    vals = [O.get(((), ()), 0)]; X = O
    for N in range(1, Nmax+1):
        X = r.V(r.G0(X)); lim = drop*(Nmax-N)
        X = {k: v for k, v in X.items() if degree(k) <= lim and v != 0}
        vals.append(X.get(((), ()), 0))
    return np.array(vals)

def observables(R):
    b, bd = R.site_ops(0)
    return {"n": MM.mul(bd, b), "n2": MM.mul(MM.mul(bd, bd), MM.mul(b, b)), "n3": MM.mul(MM.mul(MM.mul(bd, bd), bd), MM.mul(MM.mul(b, b), b))}

out = dict(machine=platform.processor() or platform.machine(), python=platform.python_version(), note="one core, wall time; sparse: all couplings at once", K3={}, scaling={})
best = lambda f, n=3: min(f() for _ in range(n))
for Nmax in (2, 4, 10):
    R = SparseRing(K=3, **P); obs = observables(R)
    def run_sparse():
        R._col.clear(); t = time.time(); [R.coefficients(O, Nmax) for O in obs.values()]; return time.time()-t
    t_sparse = best(run_sparse, 3 if Nmax < 10 else 2)
    if Nmax <= 4:
        def run_sym():
            t = time.time(); r = MM.Ring(K=3, U=0.1, **P); [series_pruned(r, O, Nmax) for O in obs.values()]; return time.time()-t
        t_sym = best(run_sym, 2)
    else:
        t_sym = None   # about 54 s per coupling in the stored runs (exact_pool_raw.json)
    counts = {k: R.coefficients(O, Nmax, return_counts=True)[1] for k, O in obs.items()}
    out["K3"][str(Nmax)] = dict(Nmax=Nmax, t_sparse_all_U=t_sparse, t_symbolic_per_U=t_sym, labels_per_layer=counts)
    print(f"K=3 orders 0..{Nmax}: sparse {t_sparse:.2f} s (all U)   symbolic recursion {t_sym if t_sym is None else round(t_sym, 2)} s per U", flush=True)
# effect of the two selection rules (Table S-pruning): n3 on the K=3 ring through order 8
R0 = SparseRing(K=3, **P); b, bd = R0.site_ops(0); O3 = observables(R0)["n3"]
out["pruning"] = {}
for name, drop, mf in (("none", 10**6, False), ("degree_bound", 2, False), ("degree_bound_momentum", 2, True)):
    Rp = SparseRing(K=3, drop=drop, momentum_filter=mf, **P); t = time.time(); c, cnt = Rp.coefficients(O3, 8, return_counts=True)
    out["pruning"][name] = dict(labels_per_layer=cnt, total=int(sum(cnt)), t=time.time()-t, c8=float(c[8].real)); del Rp
    print(f"pruning {name:22s} labels {cnt} total {sum(cnt)} t={out['pruning'][name]['t']:.1f} s c8={c[8].real:.6e}", flush=True)
import gc; gc.collect()
for K in (3, 4, 5, 6, 8, 10, 12):
    R = SparseRing(K=K, **P); b, bd = R.site_ops(0); O2 = MM.mul(MM.mul(bd, bd), MM.mul(b, b)); O1 = MM.mul(bd, b)
    rec = {}
    for Nmax in (2, 4):
        if Nmax == 4 and K > 12: continue
        R._col.clear(); t = time.time(); R.coefficients(O1, Nmax); c, cnt = R.coefficients(O2, Nmax, return_counts=True); rec[str(Nmax)] = dict(t=time.time()-t, labels=cnt)
        print(f"K={K}: g2 through order {Nmax}: {rec[str(Nmax)]['t']:.2f} s, labels {cnt}", flush=True)
    out["scaling"][str(K)] = rec; del R; gc.collect()
Ks = np.array([3, 4, 5, 6, 8, 10, 12], float)
for Nmax in ("2", "4"):
    t = np.array([out["scaling"][str(int(K))][Nmax]["t"] for K in Ks]); sel = Ks >= 6
    out["scaling"][f"exponent_N{Nmax}_K6to12"] = float(np.polyfit(np.log(Ks[sel]), np.log(t[sel]), 1)[0])
    print(f"order {Nmax}: exponent over K=6..12: {out['scaling'][f'exponent_N{Nmax}_K6to12']:.2f}")
json.dump(out, open("sparse_ring_timing.json", "w"), indent=1)
print("wrote sparse_ring_timing.json")
