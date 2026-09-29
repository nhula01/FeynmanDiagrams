"""Counting cumulants with the full Rayleigh-Schroedinger recursion for the tilted eigenvalue.
The notebook's fcs.theta_series / fcs_multi.TiltedChain.theta_series omit the terms -sum_{k=1}^{N-1} theta_k psi_{N-k};
TiltedChainRS restores them.  usage: fcs_rs.py NPROC"""
import sys, json, os, time, numpy as np
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
from multiprocessing import Pool
sys.path.insert(0, os.path.join(_ROOT, 'engines'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tensor_engine as TE
from fcs_multi import TiltedChain
from series import series_record
RP = json.load(open(os.path.join(_ROOT, 'notebook', 'data', 'revision_params.json')))
H = 0.05; CHIS = [-2*H, -H, 0.0, H, 2*H]
NPROC = int(sys.argv[1]) if len(sys.argv) > 1 else 4

class TiltedChainRS(TiltedChain):
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

def cps(table):
    th = np.array(table).real
    c1 = (th[3]-th[1])/(2*H); c2 = (th[3]-2*th[2]+th[1])/H**2; c3 = (th[4]-2*th[3]+2*th[1]-th[0])/(2*H**3)
    return c1, c2/c1, c3/c1

def build(kind, U, chi):
    if kind == "k1":
        ch = TiltedChainRS(1, 1.0, [-1.0], 0.0, [1.0], U, ring=False); ch.set_tilt(0, chi, 1.0); return ch
    if kind == "k3":
        p = RP["table5_fcs"]["K3"]; eta = np.zeros(3); eta[p["eta_site"]] = 1.0
        ch = TiltedChainRS(3, 1.0, [p["Delta"]]*3, p["J"], eta, U, ring=True); ch.set_tilt(p["eta_site"], chi, 1.0); return ch
    p = RP["chain8_fcs"]
    ch = TiltedChainRS(p["K"], 1.0, p["detunings"], p["J"], np.array(p["eta"]), U, ring=p["ring"]); ch.set_tilt(p["count_site"], chi, 1.0); return ch

def task(a):
    kind, U, chi, Nmax, md, rs = a; t = time.time()
    th = build(kind, U, chi).theta_series(Nmax, md, full_rs=rs)
    return (a, [complex(x) for x in th], time.time()-t)

if __name__ == '__main__':
    t0 = time.time()
    tasks = []
    for U in (0.05, 0.10):
        for md in (14, 18):
            tasks += [("k1", U, x, 8, md, True) for x in CHIS] + [("k1", U, x, 8, md, False) for x in CHIS]
    tasks += [("k3", 0.05, x, 6, 8, True) for x in CHIS]
    tasks += [("k3", 0.05, x, N, N+2, True) for N in (2, 3, 4) for x in CHIS]
    ch8U = [u for u in RP["chain8_fcs"]["U_values"] if u > 0]
    tasks += [("chain8", U, x, 6, 5, True) for U in ch8U for x in CHIS]
    tasks += [("chain8", U, x, 6, 6, True) for U in ch8U for x in CHIS]
    tasks += [("k3", 0.05, x, 6, 9, True) for x in CHIS]
    tasks += [("chain8", U, x, 4, 7, True) for U in (0.10, 0.06) for x in CHIS]
    res = {}
    def save():
        json.dump({"|".join(map(str, k)): dict(theta=[[z.real, z.imag] for z in v[0]], time=v[1]) for k, v in res.items()}, open('fcs_rs_raw.json', 'w'))
    with Pool(NPROC) as pool:
        for a, th, dt in pool.imap_unordered(task, tasks):
            res[a] = (th, dt); print(f"{a} {dt:.0f}s [{time.time()-t0:.0f}s]", flush=True); save()
    print("all tasks done", flush=True)
