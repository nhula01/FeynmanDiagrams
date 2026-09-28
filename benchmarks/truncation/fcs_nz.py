"""Faster tilted-diagram tasks: TiltedChain with (i) the full Rayleigh-Schroedinger recursion and (ii) the commutator with the
Kerr vertex restricted to terms with at least one contraction (the uncontracted terms cancel identically; forming them only
adds roundoff-level blocks).  Tasks: Table V K=3 diagrams on the contour, Gaussian K=1 real-axis Richardson FD, and an RS check
for the eight-site chain.  usage: fcs_nz.py PART NPROC   PART in {smoke, all}"""
import sys, os, json, time, math, numpy as np
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
from multiprocessing import Pool
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
sys.path.insert(0, os.path.join(_ROOT, 'engines'))
import tensor_engine as TE
from fcs_multi import TiltedChain
import table5_contour as T5
from timing_v3 import mul_nz
RP = T5.RP

class TiltedNZ(TiltedChain):
    def Vx(self, A, maxdeg=None):
        return TE.scale(1j, TE.add(mul_nz(self.Hint, A, maxdeg), TE.scale(-1, mul_nz(A, self.Hint, maxdeg))))
    def theta_series(self, Nmax, maxdeg=8, full_rs=True):
        Y = self.Y00(maxdeg); th = [complex(self.theta0)]; psi = [Y]
        for N in range(1, Nmax+1):
            W = self.Vx(psi[N-1], maxdeg=maxdeg)
            c00 = complex(np.asarray(W.get((0, 0), 0))); th.append(c00)
            W = TE.add(W, TE.scale(-c00, Y))
            if full_rs:
                for k in range(1, N): W = TE.add(W, TE.scale(-th[k], psi[N-k]))
            X = TE.prune(self.resolvent(W, maxdeg), maxdeg); psi.append(X)
        return np.array(th)

def k3chain(chi, U=None):
    p = RP["table5_fcs"]["K3"]; eta = np.zeros(3); eta[p["eta_site"]] = 1.0
    ch = TiltedNZ(3, 1.0, [p["Delta"]]*3, p["J"], eta, p["U"] if U is None else U, ring=True); ch.set_tilt(p["eta_site"], chi, 1.0); return ch
def c8chain(chi, U):
    p = RP["chain8_fcs"]; ch = TiltedNZ(p["K"], 1.0, p["detunings"], p["J"], np.array(p["eta"]), U, ring=p["ring"]); ch.set_tilt(p["count_site"], chi, 1.0); return ch

def task(a):
    t = time.time(); kind = a[0]
    if kind == "k3":
        _, N, md, rs, r, M, j = a; th = np.cumsum(k3chain(T5.circle(r, M)[j]).theta_series(N, md, full_rs=rs))
        return (a, [[x.real, x.imag] for x in th], time.time()-t)
    if kind == "c8":
        _, U, N, md, rs, r, M, j = a; th = np.cumsum(c8chain(T5.circle(r, M)[j], U).theta_series(N, md, full_rs=rs))
        return (a, [[x.real, x.imag] for x in th], time.time()-t)
    if kind == "gfd":
        _, U, Nc, h = a
        from gaussian_fcs import gaussian_model
        import scipy.sparse as sp
        H, A, Ad, info = gaussian_model(1.0, -1.0, 1.0, U, Nc=Nc)
        L, Jc, D = T5.liouv(sp.csr_matrix(H), [sp.csr_matrix(A)], 0); L = L.toarray(); Jd = Jc.toarray()
        th = [T5.lead_dense(L + (np.exp(x)-1)*Jd).real for x in np.array([-2, -1, 0, 1, 2])*h]
        return (a, th, time.time()-t)

if __name__ == '__main__':
    part = sys.argv[1]; nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    if part == "smoke":
        from fcs_multi import TiltedChain as TC
        chi = T5.circle(0.1, 16)[3]
        p = RP["table5_fcs"]["K3"]; eta = np.zeros(3); eta[0] = 1.0
        o = TC(3, 1.0, [-1.0]*3, 0.4, eta, 0.05, ring=True); o.set_tilt(0, chi, 1.0)
        t = time.time(); a = o.theta_series(3, 5); t1 = time.time()-t
        t = time.time(); b = k3chain(chi).theta_series(3, 5, full_rs=False); t2 = time.time()-t
        print("nb orig vs nz", np.max(np.abs(a-b)), t1, t2)
        print("rs-nb", np.abs(k3chain(chi).theta_series(4, 6, True) - k3chain(chi).theta_series(4, 6, False)))
        sys.exit()
    M = 16; H5 = list(range(M//2+1)); tasks = []
    tasks += [("k3", 6, 8, True, 0.1, M, j) for j in H5] + [("k3", 6, 8, False, 0.1, M, j) for j in H5]
    tasks += [("k3", N, N+2, rs, 0.1, M, j) for N in (2, 3, 4) for rs in (True, False) for j in H5]
    c8U = [u for u in RP["chain8_fcs"]["U_values"] if u > 0]
    tasks += [("c8", U, 6, 5, rs, 0.1, M, j) for U in c8U for rs in (True, False) for j in H5]
    tasks += [("gfd", U, Nc, h) for U in (0.05, 0.10) for Nc in (30,) for h in (0.05, 0.025, 0.0125, 0.00625)]
    tasks += [("k3", 6, 8, True, 0.2, M, j) for j in H5] + [("k3", 6, 9, True, 0.1, M, j) for j in H5]
    res = []; t0 = time.time()
    with Pool(nproc) as pool:
        for rr in pool.imap_unordered(task, tasks):
            res.append(rr); print(rr[0], f"{rr[2]:.1f}s [{time.time()-t0:.0f}s]", flush=True)
            json.dump([[list(x[0]), x[1], x[2]] for x in res], open('fcs_nz_raw.json', 'w'))
    print("done")
