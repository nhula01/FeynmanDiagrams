"""Exact counting statistics for small Kerr arrays from the tilted Liouvillian.

Long-time cumulants c_n: Taylor coefficients of the leading eigenvalue theta(chi) of
L_chi = L_0 + (e^chi - 1) J_0, obtained by a Cauchy contour integral (FFT of theta on a circle
|chi| = r), which is exponentially accurate and avoids finite-difference step errors.
Finite-window cumulants kappa_n(T) of the counts in a window of length T started from the steady
state: Taylor coefficients of log Tr[exp(L_chi T) rho_ss], same contour method.
Usage: exact_fcs.py K1|K3 U cut1,cut2,... [Ts]
"""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, sys, json, math, time

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
        n = a[j].T @ a[j]
        H = H + det[j]*n + 0.5*U*(a[j].T @ a[j].T @ a[j] @ a[j]) + 1j*(eta[j]*a[j].T - np.conj(eta[j])*a[j])
    bonds = [(j, j+1) for j in range(K-1)] + ([(K-1, 0)] if (ring and K > 2) else [])
    for j, k in bonds:
        h = a[j].T @ a[k]; H = H - J*(h + h.T)
    return sp.csr_matrix(H), a

def liouv(H, a, count):
    D = H.shape[0]; Id = sp.identity(D, format='csr')
    L = -1j*(sp.kron(Id, H) - sp.kron(H.T, Id))
    for c in a:
        cdc = (c.conj().T @ c)
        L = L + sp.kron(c.conj(), c) - 0.5*(sp.kron(Id, cdc) + sp.kron(cdc.T, Id))
    Jc = sp.kron(a[count].conj(), a[count])
    return sp.csr_matrix(L), sp.csr_matrix(Jc), D

def steady(L, D):
    # replace one row by trace condition
    tr = sp.csr_matrix(np.eye(D).reshape(1, -1, order='F'))
    A = sp.vstack([L[1:], tr]).tocsc(); b = np.zeros(D*D, complex); b[-1] = 1
    rho = spl.spsolve(A, b)
    return rho

def theta_lead(Lc, dense, sigma=0.0):
    if dense:
        w = np.linalg.eigvals(Lc.toarray()); return w[np.argmax(w.real)]
    w = spl.eigs(Lc.tocsc(), k=1, sigma=sigma, which='LM', return_eigenvectors=False, tol=1e-13)
    return w[0]

def taylor(f, r, N=32, nmax=5):
    z = r*np.exp(2j*np.pi*np.arange(N)/N)
    vals = np.array([f(x) for x in z])
    coef = np.fft.fft(vals)/N
    return np.array([ (coef[n]/r**n) * math.factorial(n) for n in range(nmax+1)])

if __name__ == "__main__":
    which = sys.argv[1]; U = float(sys.argv[2]); cuts = [x if ':' in x else int(x) for x in sys.argv[3].split(',')]
    Ts = ([] if sys.argv[4] == 'none' else [float(x) for x in sys.argv[4].split(',')]) if len(sys.argv) > 4 else [2.5, 5, 7.5, 10, 15, 20, 30, 40, 80]
    radii = [float(x) for x in sys.argv[5].split(',')] if len(sys.argv) > 5 else [0.15, 0.08]
    out = dict(model=which, U=U, results=[])
    for cut in cuts:
        t0 = time.time()
        if which == 'K1':
            cutl = [cut]; det = [-1.0]; J = 0.0; eta = [1.0]; ring = False
        elif which == 'K2':
            # first two sites of the eight-site chain (proxy for the site-1 cutoff), cut given as 'c0:c1'
            cutl = [int(v) for v in cut.split(':')]; det = [-1.0390314686649147, -0.7154882840444468]; J = 0.4; eta = [1.0, 0.0]; ring = False
        else:
            cutl = [cut, cut-3, cut-3]; det = [-1.0]*3; J = 0.4; eta = [1.0, 0, 0]; ring = True
        H, a = ops(cutl, det, J, eta, U, ring)
        L, Jc, D = liouv(H, a, 0)
        dense = D*D <= 400
        rho = steady(L, D); rho = rho/np.trace(rho.reshape(D, D, order='F'))
        R = rho.reshape(D, D, order='F')
        nocc = [float(np.real(np.trace(R @ (x.conj().T @ x).toarray()))) for x in a]
        # top-level population of each site
        tops = []
        for j, x in enumerate(a):
            n = (x.conj().T @ x).diagonal().real
            tops.append(float(np.real(np.diag(R))[np.isclose(n, cutl[j]-1)].sum()))
        th = lambda chi: theta_lead(L + (np.exp(chi)-1)*Jc, dense, sigma=nocc[0]*(np.exp(chi)-1))
        res = dict(cut=cutl, D=D, n=nocc, top=tops)
        for rad in radii:
            cn = taylor(th, rad, N=24 if not dense else 32).real
            res[f'theta_r{rad}'] = dict(c=cn[1:].tolist(), ratios=[cn[1], cn[2]/cn[1], cn[3]/cn[1], cn[4]/cn[1], cn[5]/cn[1]])
        # finite-difference convention of tilted_exact.py (h=0.05)
        h = 0.05; tv = [th(x).real for x in np.array([-2, -1, 0, 1, 2])*h]
        c1 = (tv[3]-tv[1])/(2*h); c2 = (tv[3]-2*tv[2]+tv[1])/h**2; c3 = (tv[4]-2*tv[3]+2*tv[1]-tv[0])/(2*h**3)
        res['fd_h0.05'] = [c1, c2/c1, c3/c1]
        # spectrum of L0 near zero
        if dense:
            w = np.linalg.eigvals(L.toarray())
        else:
            w = spl.eigs(L.tocsc(), k=8, sigma=0.05, which='LM', return_eigenvectors=False, tol=1e-10)
        w = w[np.argsort(-w.real)][:8]; res['L0_eigs'] = [[float(x.real), float(x.imag)] for x in w]
        # finite-window cumulants from steady state
        if not Ts:
            res['finite_window'] = {}; res['wall'] = time.time()-t0; out['results'].append(res); print(json.dumps(res), flush=True); continue
        Tmax = max(Ts); nT = int(round(Tmax/2.5))
        idx = [int(round(T/2.5)) for T in Ts]
        trv = np.eye(D).reshape(-1, order='F')
        def logZ_all(chi):
            Lc = (L + (np.exp(chi)-1)*Jc).tocsc()
            V = spl.expm_multiply(Lc, rho, start=0, stop=Tmax, num=nT+1, endpoint=True)
            return V @ trv
        rad = 0.1; N = 32 if dense else 16
        z = rad*np.exp(2j*np.pi*np.arange(N)/N)
        Zv = np.array([logZ_all(x) for x in z])         # (N, nT+1)
        # continuous logarithm along the closed contour (Z is zero-free inside for these T)
        ph = np.unwrap(np.angle(Zv), axis=0)
        vals = np.log(np.abs(Zv)) + 1j*ph
        res['winding_max'] = float(np.max(np.abs(ph[-1] - ph[0] + np.angle(Zv[0]/Zv[-1]))/(2*np.pi)))   # 0 if Z has no zeros inside
        coef = np.fft.fft(vals, axis=0)/N
        kap = np.array([(coef[n]/rad**n)*math.factorial(n) for n in range(4)]).real   # (4, nT+1)
        fw = {}
        for T, i in zip(Ts, idx):
            k1, k2, k3 = kap[1, i], kap[2, i], kap[3, i]
            fw[str(T)] = [k1/T, k2/k1, k3/k1]
        res['finite_window'] = fw
        res['wall'] = time.time()-t0
        out['results'].append(res)
        print(json.dumps(res), flush=True)
    json.dump(out, open(f"exact_{which}_U{U:.2f}.json", 'w'), indent=1)

