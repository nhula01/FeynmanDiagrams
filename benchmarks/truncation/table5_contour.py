"""Table V with contour (Cauchy) derivatives of the tilted eigenvalue theta(chi).
c_n = n! (1/M) sum_j theta(r e^{i phi_j}) e^{-i n phi_j} / r^n,  phi_j = 2 pi j/M.  theta(conj chi) = conj theta(chi),
so only j = 0..M/2 are evaluated.  Two radii (0.1, 0.2) are used as an internal check.
Diagrams: full Rayleigh-Schroedinger recursion ('rs'); the notebook recursion ('nb', omits -sum_k theta_k psi_{N-k}) for reference.
usage: table5_contour.py PART NPROC      PART in {k1, k3diag, k3exact, smoke}
"""
import sys, os, json, math, time, numpy as np
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
import scipy.sparse as sp, scipy.sparse.linalg as spl
from multiprocessing import Pool
sys.path.insert(0, os.path.join(_ROOT, 'engines'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
RP = json.load(open(os.path.join(_ROOT, 'notebook', 'data', 'revision_params.json')))
KAP = 1.0; DC = -1.0; ETA = 1.0
P3 = RP["table5_fcs"]["K3"]; U3 = P3["U"]; J3 = P3["J"]

def circle(r, M):
    return r*np.exp(2j*np.pi*np.arange(M)/M)
def half(M): return list(range(M//2+1))
def fill(vals_half, M):
    """vals_half[j] for j=0..M/2 (arrays) -> full list using theta(conj chi)=conj theta(chi)"""
    full = [None]*M
    for j in range(M//2+1): full[j] = np.asarray(vals_half[j])
    for j in range(M//2+1, M): full[j] = np.conj(full[M-j])
    return np.array(full)
def taylor(full, r, nmax=5):
    M = full.shape[0]; co = np.fft.fft(full, axis=0)/M
    return np.array([math.factorial(n)*co[n]/r**n for n in range(nmax+1)])   # (nmax+1, ...)

# ---------------------------------------------------------------- K=1 diagrams (engine.py algebra)
def theta_k1(U, chi, Nmax, maxdeg, full_rs=True):
    from fcs import ladder, Y00
    from engine import add, scale, mul, cav, comm
    z = KAP/2 + 1j*DC; alpha = ETA/z; eps = KAP*(np.exp(chi)-1)
    I1 = np.eye(1); b, bd, one = cav(0, 1, I1), cav(1, 0, I1), cav(0, 0, I1)
    a_ = add(b, scale(alpha, one)); ad_ = add(bd, scale(np.conj(alpha), one))
    HK = scale(U/2, mul(mul(ad_, ad_), mul(a_, a_))); V = lambda X: scale(1j, comm(HK, X))
    Y = Y00(alpha, eps, z, maxdeg); th = [eps*abs(alpha)**2]; psi = [Y]
    for N in range(1, Nmax+1):
        W = V(psi[N-1]); c00 = complex(W.get((0, 0), np.zeros((1, 1)))[0, 0]); th.append(c00)
        W = add(W, scale(-c00, Y))
        if full_rs:
            for k in range(1, N): W = add(W, scale(-th[k], psi[N-k]))
        X = ladder(W, alpha, eps, z, maxdeg)
        psi.append({k: v for k, v in X.items() if sum(k) <= maxdeg and np.max(np.abs(v)) > 1e-16})
    return np.array(th, complex)


def _rs_class():
    import tensor_engine as TE
    from fcs_multi import TiltedChain
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
    return TiltedChainRS

# ---------------------------------------------------------------- exact / Gaussian tilted Liouvillians
def ops(cut, det, J, eta, U, ring):
    cut = list(cut); K = len(cut)
    a1 = [sp.diags(np.sqrt(np.arange(1, n)), 1, shape=(n, n), format='csr') for n in cut]
    I = [sp.identity(n, format='csr') for n in cut]
    def emb(op, j):
        m = sp.identity(1, format='csr')
        for k in range(K): m = sp.kron(m, op if k == j else I[k], format='csr')
        return m
    a = [emb(a1[j], j) for j in range(K)]
    H = 0
    for j in range(K):
        H = H + det[j]*(a[j].T@a[j]) + 0.5*U*(a[j].T@a[j].T@a[j]@a[j]) + 1j*(eta[j]*a[j].T - np.conj(eta[j])*a[j])
    bonds = [(j, j+1) for j in range(K-1)] + ([(K-1, 0)] if (ring and K > 2) else [])
    for j, k in bonds:
        h = a[j].T@a[k]; H = H - J*(h + h.T)
    return sp.csr_matrix(H), a
def liouv(H, a, count, kap=KAP):
    D = H.shape[0]; Id = sp.identity(D, format='csr')
    L = -1j*(sp.kron(Id, H) - sp.kron(H.T, Id))
    for c in a:
        c = np.sqrt(kap)*c; cdc = (c.conj().T@c)
        L = L + sp.kron(c.conj(), c) - 0.5*(sp.kron(Id, cdc) + sp.kron(cdc.T, Id))
    cc = np.sqrt(kap)*a[count]
    return sp.csr_matrix(L), sp.csr_matrix(sp.kron(cc.conj(), cc)), D
def lead_dense(L):
    w = np.linalg.eigvals(L); return w[np.argmax(w.real)]
def theta_evolve(L, Jc, D, chi, T=60.0, steps=6):
    Lc = (L + (np.exp(chi)-1)*Jc).tocsr(); tri = np.arange(D)*(D+1)   # vec(F) diagonal indices
    v = np.zeros(D*D, complex); v[0] = 1.0; hist = []
    for s in range(steps):
        v = spl.expm_multiply(Lc*(T/steps), v); v = v/np.sum(v[tri])
        hist.append(np.sum((Lc@v)[tri]))
    return hist[-1], abs(hist[-1]-hist[-2])

# ---------------------------------------------------------------- tasks
def t_k1_exact(a):
    U, cut, r, M = a; t = time.time()
    H, aa = ops([cut], [DC], 0.0, [ETA], U, False); L, Jc, D = liouv(H, aa, 0)
    L = L.toarray(); Jd = Jc.toarray(); z = circle(r, M)
    vals = fill([lead_dense(L + (np.exp(z[j])-1)*Jd) for j in half(M)], M)
    return ("k1_exact", a, taylor(vals, r).tolist(), time.time()-t)
def t_k1_gauss(a):
    U, Nc, r, M = a; t = time.time()
    from gaussian_fcs import gaussian_model
    H, A, Ad, info = gaussian_model(KAP, DC, ETA, U, Nc=Nc)
    Hs = sp.csr_matrix(H); As = sp.csr_matrix(A); L, Jc, D = liouv(Hs, [As], 0); L = L.toarray(); Jd = Jc.toarray(); z = circle(r, M)
    vals = fill([lead_dense(L + (np.exp(z[j])-1)*Jd) for j in half(M)], M)
    return ("k1_gauss", a, taylor(vals, r).tolist(), time.time()-t)
def t_k1_diag(a):
    U, md, rs, r, M = a; t = time.time(); z = circle(r, M)
    vals = fill([np.cumsum(theta_k1(U, z[j], 8, md, rs)) for j in half(M)], M)     # (M, 9)
    return ("k1_diag", a, taylor(vals, r).tolist(), time.time()-t)
def t_k3_diag_point(a):
    N, md, rs, r, M, j = a; t = time.time()
    TiltedChainRS = _rs_class()
    eta = np.zeros(3); eta[P3["eta_site"]] = 1.0
    ch = TiltedChainRS(3, KAP, [P3["Delta"]]*3, J3, eta, U3, ring=True); ch.set_tilt(P3["eta_site"], circle(r, M)[j], KAP)
    th = np.cumsum(ch.theta_series(N, md, full_rs=rs))
    return ("k3_diag", a, [[x.real, x.imag] for x in th], time.time()-t)
def t_k3_exact_point(a):
    cut, r, M, j = a; t = time.time()
    eta = np.zeros(3); eta[P3["eta_site"]] = 1.0
    H, aa = ops([cut]*3, [P3["Delta"]]*3, J3, eta, U3, True); L, Jc, D = liouv(H, aa, P3["eta_site"])
    th, conv = theta_evolve(L, Jc, D, circle(r, M)[j])
    return ("k3_exact", a, [th.real, th.imag, conv], time.time()-t)

def run(tasks, fn, nproc, outname):
    res = []; t0 = time.time()
    with Pool(nproc) as pool:
        for rr in pool.imap_unordered(fn, tasks):
            res.append(rr); print(rr[0], rr[1], f"{rr[3]:.1f}s [{time.time()-t0:.0f}s]", flush=True)
            json.dump([[x[0], list(x[1]), x[2], x[3]] for x in res], open(outname, 'w'), default=lambda o: [o.real, o.imag] if isinstance(o, complex) else str(o))
    return res

def cplx(o):
    return [o.real, o.imag] if isinstance(o, complex) else o

if __name__ == '__main__':
    part = sys.argv[1]; nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    if part == "smoke":
        z = 0.1*np.exp(0.7j)
        H, aa = ops([20], [DC], 0.0, [ETA], 0.05, False); L, Jc, D = liouv(H, aa, 0)
        te = lead_dense((L + (np.exp(z)-1)*Jc).toarray())
        td = np.cumsum(theta_k1(0.05, z, 8, 14, True)); tn = np.cumsum(theta_k1(0.05, z, 8, 14, False))
        print("K1 exact", te, "diag rs N8", td[-1], "diff", abs(td[-1]-te), "nb diff", abs(tn[-1]-te))
        th, conv = theta_evolve(L, Jc, D, z); print("K1 evolve", th, "diff to dense", abs(th-te), "conv", conv)
        r = t_k1_exact((0.05, 20, 0.1, 16)); print("K1 exact c:", np.round(np.real(r[2][1:4]), 6))
        t = time.time(); r = t_k3_diag_point((2, 4, True, 0.1, 16, 3)); print("K3 diag N2 md4", r[2][-1], f"{time.time()-t:.1f}s")
        t = time.time(); r = t_k3_exact_point((5, 0.1, 16, 3)); print("K3 exact cut5", r[2], f"{time.time()-t:.1f}s")
        sys.exit()
    if part == "k1":
        tasks_e = [(U, cut, r, 32) for U in (0.05, 0.10) for cut in (20, 26, 32) for r in (0.1, 0.2)]
        tasks_g = [(U, Nc, r, 32) for U in (0.05, 0.10) for Nc in (30, 40) for r in (0.1, 0.2)]
        tasks_d = [(U, md, rs, r, 16) for U in (0.05, 0.10) for md in (14, 18) for rs in (True, False) for r in (0.1, 0.2)]
        res = run(tasks_e, t_k1_exact, nproc, 'table5_k1_raw_exact.json') + run(tasks_g, t_k1_gauss, nproc, 'table5_k1_raw_gauss.json') \
            + run(tasks_d, t_k1_diag, nproc, 'table5_k1_raw_diag.json')
        pass
    if part == "k3diag":
        M = 16; tasks = []
        for (N, md, rs, r) in [(6, 8, True, 0.1), (6, 8, True, 0.2), (6, 8, False, 0.1), (2, 4, True, 0.1), (3, 5, True, 0.1), (4, 6, True, 0.1),
                               (2, 4, False, 0.1), (3, 5, False, 0.1), (4, 6, False, 0.1), (6, 9, True, 0.1)]:
            tasks += [(N, md, rs, r, M, j) for j in half(M)]
        res = run(tasks, t_k3_diag_point, nproc, 'table5_k3diag_raw.json')
    if part == "k3exact":
        M = 16; cuts = [int(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else [7, 8, 9, 10]
        tasks = [(cut, r, M, j) for cut in cuts for r in (0.1, 0.2) for j in half(M)]
        res = run(tasks, t_k3_exact_point, nproc, f'table5_k3exact_raw_{"_".join(map(str, cuts))}.json')
    print("done")
