"""Displaced-frame exact references: K=3 uniform ring moments (13 U, b-levels 4/5/6) and Table V K=3 / K=1 counting
cumulants by contour integration (b-levels [c0,c1,c1]).  usage: exact_disp_pool.py NPROC"""
import sys, os, json, time, math, numpy as np
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
from multiprocessing import Pool
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import exact_displaced as ED
import scipy.sparse.linalg as spl
RP = json.load(open(os.path.join(_ROOT, 'notebook', 'data', 'revision_params.json')))
p3 = RP["table5_fcs"]["K3"]
Us13 = [round(float(x), 3) for x in np.linspace(0, 0.3, 13)]
M = 16
def circle(r): return r*np.exp(2j*np.pi*np.arange(M)/M)

def fcs_model(which, cuts, U):
    if which == "K3":
        eta = np.zeros(3); eta[p3["eta_site"]] = 1.0
        return ED.build(list(cuts), 1.0, p3["Delta"], p3["J"], eta, p3["U"], True, count=p3["eta_site"])
    return ED.build(list(cuts), 1.0, -1.0, 0.0, [1.0], U, False, count=0)

def task(a):
    t0 = time.time()
    if a[0] == "ring":
        _, U, lv = a; m = ED.ring_moments(3, U, lv, method="evolve"); return (a, m, time.time()-t0)
    if a[0] == "fcs":
        _, which, cuts, U, r = a
        L, Jc, D, A, Ad = fcs_model(which, cuts, U)
        R = ED.steady(L, D, "evolve"); v0 = R.reshape(-1, order='F')
        n = [float(np.real(np.sum((Ad[j]@A[j]).multiply(R.T)))) for j in range(len(cuts))]
        z = circle(r); vals = []; conv = []
        for j in range(M//2+1):
            th, cv = ED.theta_lead(L, Jc, D, z[j], v0=v0); vals.append(th); conv.append(float(cv))
        full = vals + [np.conj(vals[M-j]) for j in range(M//2+1, M)]
        co = np.fft.fft(np.array(full))/M; cn = [math.factorial(k)*co[k]/r**k for k in range(6)]
        return (a, dict(c=[[x.real, x.imag] for x in cn], c1=cn[1].real, fano=(cn[2]/cn[1]).real, c3c1=(cn[3]/cn[1]).real,
                        c0_abs=abs(cn[0]), maxconv=max(conv), n=n, D=D), time.time()-t0)

if __name__ == '__main__':
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    tasks = [("fcs", "K3", c, 0.05, r) for c in ((6, 4, 4), (7, 5, 5), (8, 5, 5), (8, 6, 6), (9, 6, 6)) for r in (0.15,)] + [("fcs", "K3", (8, 5, 5), 0.05, 0.08)]
    tasks += [("fcs", "K1", (c,), U, 0.05) for U in (0.05, 0.10) for c in (8, 10, 12, 14)]
    tasks += [("ring", U, lv) for lv in (4, 5, 6) for U in Us13]
    res = []; t0 = time.time()
    with Pool(nproc) as pool:
        for rr in pool.imap_unordered(task, tasks):
            res.append(rr); v = rr[1]
            msg = f"g2 {v['g2']:.10f} g3 {v['g3']:.10f}" if rr[0][0] == "ring" else f"c1 {v['c1']:.8f} F {v['fano']:.8f} c3/c1 {v['c3c1']:.8f} |c0| {v['c0_abs']:.1e} conv {v['maxconv']:.1e} n {np.round(v['n'], 4)}"
            print(rr[0], msg, f"{rr[2]:.0f}s [{time.time()-t0:.0f}s]", flush=True)
            json.dump([[list(x[0]), x[1], x[2]] for x in res], open('exact_disp_raw.json', 'w'))
    print("done")
