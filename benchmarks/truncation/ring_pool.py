"""K=3 Kerr ring: g2 and g3 partial sums to order NMAX (momentum engine, degree pruning), exact steady
states (sparse time evolution, cutoffs 7/8/9 levels per site), optimal truncation.  Parallel over U.
usage: ring_pool.py NMAX NPROC"""
import sys, json, os, time, numpy as np
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
from multiprocessing import Pool
sys.path.insert(0, os.path.join(_ROOT, 'engines'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import multimode as MM, tensor_engine as TE
from truncation import truncation_record
from ring_trunc import series_pruned, series_unpruned
from ring_ss_evolve import ring_steady_evolve

# ---- commutator without the zero-contraction terms.  For normal-ordered monomials the uncontracted parts of AB and BA
# are the same operator, so they cancel exactly in [A,B]; MM.comm forms both and subtracts, which leaves roundoff labels
# (|v| ~ 1e-17 relative) of degree d+4 that the degree pruning keeps and that multiply at every later vertex.
from math import comb, factorial
from itertools import product
def mul_nz(A, B):
    out = {}
    for (mA, nA), x in A.items():
        for (mB, nB), y in B.items():
            nAd, mBd = MM.todict(nA), MM.todict(mB)
            modes = [q for q in nAd if q in mBd]
            if not modes: continue
            for ks in product(*[range(min(nAd[q], mBd[q])+1) for q in modes]):
                if not any(ks): continue
                w = 1.0; m = MM.todict(mA); n = MM.todict(nB)
                for q, c in mBd.items(): m[q] = m.get(q, 0) + c
                for q, c in nAd.items(): n[q] = n.get(q, 0) + c
                for q, k in zip(modes, ks):
                    w *= factorial(k)*comb(nAd[q], k)*comb(mBd[q], k); m[q] -= k; n[q] -= k
                kk = MM.key(m, n); out[kk] = out.get(kk, 0) + w*x*y
    return out
def V_nz(r, A):
    return MM.scale(1j, MM.add(mul_nz(r.Hint, A), MM.scale(-1, mul_nz(A, r.Hint))))
def deg(k): return sum(c for _, c in k[0]) + sum(c for _, c in k[1])
def series_nz(r, O, Nmax, counts=None):
    vals = [O.get(((), ()), 0)]; X = O
    if counts is not None: counts.append(len(X))
    for N in range(1, Nmax+1):
        X = V_nz(r, r.G0(X)); lim = 2*(Nmax-N)
        X = {k: v for k, v in X.items() if deg(k) <= lim and v != 0}
        if counts is not None: counts.append(len(X))
        vals.append(X.get(((), ()), 0))
    return np.array(vals)

PR = json.load(open(os.path.join(_ROOT, 'notebook', 'data', 'revision_params.json')))["ring3"]
P = dict(kappa=PR["kappa"], Delta=PR["Delta"], J=PR["J"], eta=PR["eta"])
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10
NPROC = int(sys.argv[2]) if len(sys.argv) > 2 else 8
Us13 = [round(float(x), 3) for x in np.linspace(0, 0.3, 13)]
Us11 = [round(float(u), 3) for u in PR["U_values_fig11"]]

def task_exact(a):
    U, nm = a; t = time.time()
    e = ring_steady_evolve(3, U=U, nmax=nm, T=80.0, steps=8, **P)
    return (U, nm, dict(n0=float(e["n0"]), n0n0=float(e["n0n0"]), n0n1=float(e["n0n1"]), n3=float(e["n3"]), time=time.time()-t))

def task_series(U):
    t = time.time(); r = MM.Ring(K=3, U=U, **P); a0, ad0 = r.site_ops(0)
    On = MM.mul(ad0, a0); O2 = MM.mul(MM.mul(ad0, ad0), MM.mul(a0, a0)); O3 = MM.mul(MM.mul(MM.mul(ad0, ad0), ad0), MM.mul(MM.mul(a0, a0), a0))
    out = {}
    for name, O in (("n", On), ("n2", O2), ("n3", O3)):
        cnt = []; s = series_nz(r, O, NMAX, counts=cnt)
        out[name] = dict(terms_re=[float(np.real(x)) for x in s], terms_im=[float(np.imag(x)) for x in s], labels=cnt)
    out["time"] = time.time()-t
    return (U, out)

def task_cross(_):
    r = MM.Ring(K=3, U=0.15, **P); a0, ad0 = r.site_ops(0)
    O2 = MM.mul(MM.mul(ad0, ad0), MM.mul(a0, a0))
    t = time.time(); sp_ = series_pruned(r, O2, 4); tp = time.time()-t
    t = time.time(); su_ = series_unpruned(r, O2, 4); tu = time.time()-t
    ch = TE.Chain(3, P["kappa"], [P["Delta"]]*3, P["J"], P["eta"], 0.15, ring=True)
    a3, ad3 = ch.site_a(0), ch.site_ad(0)
    Nt = min(NMAX, 8)
    t = time.time(); st_ = ch.series(TE.mul(TE.mul(ad3, ad3), TE.mul(a3, a3)), Nt); tt = time.time()-t
    s6 = series_pruned(r, O2, Nt)
    c_old, c_new = [], []
    t = time.time(); so = series_pruned(r, O2, 6, counts=c_old); to = time.time()-t
    t = time.time(); sn = series_nz(r, O2, 6, counts=c_new); tn = time.time()-t
    return dict(pruned_vs_unpruned_maxdiff_order4=float(np.max(np.abs(sp_-su_))), t_pruned=tp, t_unpruned=tu,
                momentum_vs_tensor_maxdiff=float(np.max(np.abs(s6-st_))), tensor_orders=Nt, t_tensor=tt,
                nz_vs_original_maxdiff_order6=float(np.max(np.abs(so-sn))), labels_original=c_old, labels_nz=c_new, t_original=to, t_nz=tn)

if __name__ == '__main__':
    t0 = time.time()
    ex_tasks = [(U, 8) for U in Us13] + [(U, nm) for U in Us11 for nm in (7, 9)]
    with Pool(NPROC) as pool:
        rs = pool.map_async(task_series, Us13, chunksize=1)
        rc = pool.map_async(task_cross, [0])
        re = pool.map_async(task_exact, ex_tasks, chunksize=1)
        ser = dict(rs.get()); cross = rc.get()[0]; exl = re.get()
    print(f"parallel part done [{time.time()-t0:.0f}s]", flush=True)
    exact = {}
    for U, nm, d in exl: exact.setdefault(f"{U:.3f}", {})[f"nmax{nm}"] = d
    cached = json.load(open(os.path.join(_ROOT, 'notebook', 'data', 'ring_exactU.json')))
    out = dict(parameters=dict(K=3, **P), Nmax=NMAX, note="exact from ring_steady_evolve (T=80, 8 steps); nmaxX = X Fock levels per site; "
               "g2_N = n2_N/n_N^2 and g3_N = n3_N/n_N^3 with partial sums through order N", exact={}, series={}, results={}, results_Nmax6={}, crosscheck=cross)
    for key, d in exact.items():
        e = d["nmax8"]; rec = dict(d); rec["g2"] = e["n0n0"]/e["n0"]**2; rec["g3"] = e["n3"]/e["n0"]**3
        if key in cached:
            rec["cached_ring_exactU"] = dict(n0=cached[key]["n0"], n0n0=cached[key]["n0n0"], g2=cached[key]["n0n0"]/cached[key]["n0"]**2)
        for nm in (7, 9):
            if f"nmax{nm}" in d:
                f = d[f"nmax{nm}"]; rec[f"g2_nmax{nm}"] = f["n0n0"]/f["n0"]**2; rec[f"g3_nmax{nm}"] = f["n3"]/f["n0"]**3
        out["exact"][key] = rec
        print(f"exact U={key}: n0={e['n0']:.6f} g2={rec['g2']:.6f} g3={rec['g3']:.6f} "
              + (f"cached g2={rec['cached_ring_exactU']['g2']:.6f} " if 'cached_ring_exactU' in rec else "")
              + (f"g2(7,9)={rec.get('g2_nmax7', np.nan):.6f},{rec.get('g2_nmax9', np.nan):.6f}" ), flush=True)
    for U in Us13:
        key = f"{U:.3f}"; s = ser[U]; out["series"][key] = s
        S = {k: np.cumsum(np.array(s[k]["terms_re"]) + 1j*np.array(s[k]["terms_im"])).real for k in ("n", "n2", "n3")}
        g2 = S["n2"]/S["n"]**2; g3 = S["n3"]/S["n"]**3; ex = out["exact"][key]
        for tag, M in (("results", NMAX), ("results_Nmax6", 6)):
            out[tag][key] = dict(U=U, n=truncation_record(S["n"][:M+1], exact=ex["nmax8"]["n0"]),
                                 g2=truncation_record(g2[:M+1], exact=ex["g2"]), g3=truncation_record(g3[:M+1], exact=ex["g3"]))
        a = out["results"][key]; b = out["results_Nmax6"][key]
        print(f"U={U}: g2 N*={a['g2']['N_star']} {a['g2']['value']:.6f}+-{a['g2']['delta']:.1e} true {a['g2']['true_error']:.1e} | "
              f"g3 N*={a['g3']['N_star']} {a['g3']['value']:.6f}+-{a['g3']['delta']:.1e} true {a['g3']['true_error']:.1e} | "
              f"(Nmax6: g2 N*={b['g2']['N_star']} true {b['g2']['true_error']:.1e}, g3 N*={b['g3']['N_star']} true {b['g3']['true_error']:.1e}) [{s['time']:.0f}s]", flush=True)
    print("crosscheck", cross, flush=True)
    json.dump(out, open('ring_truncation.json', 'w'), indent=1)
    print(f"done [{time.time()-t0:.0f}s]")
