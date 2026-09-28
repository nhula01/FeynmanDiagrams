"""Combined pool: (i) K=3 uniform ring, momentum-engine series to order NMAX (zero-contraction-free commutator) and exact
steady states at Fock levels 6/7/8 per site; (ii) Table V K=3 exact counting statistics, leading eigenvalue of the tilted
Liouvillian on a circle |chi| = r by power iteration (expm_multiply) started from the chi = 0 steady state, asymmetric
Fock cutoffs [c0, c1, c1] (driven/counted site 0 carries n ~ 1.15, the others much less).  Results stream to JSON.
usage: exact_pool.py NPROC"""
import sys, os, json, time, math, numpy as np
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
import scipy.sparse as sp, scipy.sparse.linalg as spl
from multiprocessing import Pool
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
sys.path.insert(0, os.path.join(_ROOT, 'engines'))
NPROC = int(sys.argv[1]) if len(sys.argv) > 1 else 7
sys.argv = [sys.argv[0], '10', str(NPROC)]
import ring_pool as RPm
import table5_contour as T5
import multimode as MM
from ring_ss_evolve import ring_steady_evolve
P = RPm.P; NMAX = 10
Us13 = [round(float(x), 3) for x in np.linspace(0, 0.3, 13)]
K3CUTS = [(9, 6, 6), (11, 6, 6), (11, 7, 7)]
RADII = {(9, 6, 6): [0.15], (11, 6, 6): [0.15, 0.08], (11, 7, 7): [0.15]}
M = 16

def k3_L(cuts):
    eta = np.zeros(3); eta[T5.P3["eta_site"]] = 1.0
    H, aa = T5.ops(list(cuts), [T5.P3["Delta"]]*3, T5.J3, eta, T5.U3, True)
    L, Jc, D = T5.liouv(H, aa, T5.P3["eta_site"]); return L, Jc, D, aa
def evolve(Lc, D, v, T, steps):
    tri = np.arange(D)*(D+1); hist = []
    for s in range(steps):
        v = spl.expm_multiply(Lc*(T/steps), v); v = v/np.sum(v[tri]); hist.append(np.sum((Lc@v)[tri]))
    return v, hist

def task(a):
    t0 = time.time(); kind = a[0]
    if kind == "series":
        U, out = RPm.task_series(a[1]); return (a, out, time.time()-t0)
    if kind == "ring_exact":
        U, nm = a[1], a[2]; e = ring_steady_evolve(3, U=U, nmax=nm, T=60.0, steps=6, **P)
        return (a, dict(n0=float(e["n0"]), n0n0=float(e["n0n0"]), n0n1=float(e["n0n1"]), n3=float(e["n3"])), time.time()-t0)
    if kind == "k3ss":
        cuts = a[1]; L, Jc, D, aa = k3_L(cuts); v = np.zeros(D*D, complex); v[0] = 1
        v, hist = evolve(L.tocsr(), D, v, 80.0, 8); np.save(f"k3ss_{'_'.join(map(str, cuts))}.npy", v)
        R = v.reshape(D, D, order='F'); nocc = []; tops = []
        for j, x in enumerate(aa):
            nd = (x.conj().T@x).diagonal().real; pop = np.real(np.diag(R))
            nocc.append(float(pop@nd)); tops.append(float(pop[np.isclose(nd, cuts[j]-1)].sum()))
        return (a, dict(n=nocc, top_level_population=tops, theta0=[hist[-1].real, hist[-1].imag], conv=abs(hist[-1]-hist[-2])), time.time()-t0)
    if kind == "k3pt":
        cuts, r, j = a[1], a[2], a[3]; L, Jc, D, aa = k3_L(cuts)
        v0 = np.load(f"k3ss_{'_'.join(map(str, cuts))}.npy"); chi = r*np.exp(2j*np.pi*j/M)
        v, hist = evolve((L + (np.exp(chi)-1)*Jc).tocsr(), D, v0, 40.0, 5)
        return (a, dict(theta=[hist[-1].real, hist[-1].imag], conv=float(abs(hist[-1]-hist[-2])), hist=[[h.real, h.imag] for h in hist]), time.time()-t0)

if __name__ == '__main__':
    t0 = time.time(); res = {}
    def dump():
        json.dump([[list(k), v, dt] for k, (v, dt) in res.items()], open('exact_pool_raw.json', 'w'),
                  default=lambda o: list(o) if isinstance(o, tuple) else str(o))
    stage1 = [("k3ss", c) for c in K3CUTS] + [("series", U) for U in Us13] + [("ring_exact", U, 7) for U in Us13] \
             + [("ring_exact", U, 6) for U in (0.1, 0.15, 0.2, 0.25, 0.3)] + [("ring_exact", 0.3, 8)]
    stage2 = [("k3pt", c, r, j) for c in K3CUTS for r in RADII[c] for j in range(M//2+1)]
    with Pool(NPROC) as pool:
        pending = list(stage1); done_ss = set(); queued2 = False; its = []
        it = pool.imap_unordered(task, stage1)
        results2 = None
        for n in range(len(stage1)):
            a, v, dt = it.next(); res[tuple(a) if not isinstance(a[1], tuple) else (a[0],) + tuple(a[1:])] = (v, dt); dump()
            print(a, f"{dt:.0f}s [{time.time()-t0:.0f}s]", flush=True)
            if a[0] == "k3ss": done_ss.add(a[1])
            if len(done_ss) == len(K3CUTS) and results2 is None:
                results2 = pool.imap_unordered(task, stage2)
        if results2 is None: results2 = pool.imap_unordered(task, stage2)
        for n in range(len(stage2)):
            a, v, dt = results2.next(); res[(a[0],) + tuple(a[1:])] = (v, dt); dump()
            print(a, v["theta"], f"conv {v['conv']:.1e} {dt:.0f}s [{time.time()-t0:.0f}s]", flush=True)
    print("done")
