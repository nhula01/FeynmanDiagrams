"""Finite-time diagrammatic expansion (Kerr cavity switched on from vacuum).
Order-N contribution = sum over outgoing labels (p,q) of the order-N amplitude
times the initial-state weight Tr[O_pq rho0] (i.e. the r_qp component)."""
import numpy as np
from scipy.linalg import expm
from engine import *

def kerr_finite_time(kappa, Dc, eta, U, ts, Nmax, O_fn, by_degree=False):
    z = kappa/2 + 1j*Dc; alpha = eta/z; I1 = np.eye(1)
    b = cav(0, 1, I1); bd = cav(1, 0, I1); one = cav(0, 0, I1)
    a = add(b, scale(alpha, one)); ad = add(bd, scale(np.conj(alpha), one))
    HK = scale(U/2, mul(mul(ad, ad), mul(a, a)))
    O = O_fn(a, ad)
    d0 = max(m+n for (m, n) in O)
    Dmax = d0 + 2*Nmax
    labels = [(m, s-m) for s in range(Dmax+1) for m in range(s+1)]
    idx = {l: i for i, l in enumerate(labels)}; L = len(labels)
    A = np.diag([-(z*n + np.conj(z)*m) for (m, n) in labels])
    V = np.zeros((L, L), complex)
    for j, (m, n) in enumerate(labels):
        res = scale(1j, comm(HK, cav(m, n, I1)))
        for k, c in res.items():
            if k in idx:
                V[idx[k], j] += c[0, 0]
    x0 = np.zeros(L, complex)
    for k, c in O.items(): x0[idx[k]] = c[0, 0]
    beta0 = -alpha      # cavity vacuum = coherent state -alpha in the displaced frame
    w = np.array([np.conj(beta0)**m * beta0**n for (m, n) in labels])   # Tr[O_mn rho0]
    B = np.zeros(((Nmax+1)*L, (Nmax+1)*L), complex)
    for k in range(Nmax+1):
        B[k*L:(k+1)*L, k*L:(k+1)*L] = A
        if k < Nmax:
            B[k*L:(k+1)*L, (k+1)*L:(k+2)*L] = V
    from scipy.sparse.linalg import expm_multiply
    import scipy.sparse as sps
    Bs = sps.csr_matrix(B)
    out = np.zeros((len(ts), Nmax+1), complex); parts = {}
    for N in range(Nmax+1):
        v = np.zeros((Nmax+1)*L, complex); v[N*L:(N+1)*L] = x0
        traj = expm_multiply(Bs, v, start=ts[0], stop=ts[-1], num=len(ts), endpoint=True)
        out[:, N] = traj[:, 0:L] @ w
        if by_degree:
            for dgr in range(Dmax+1):
                mask = np.array([(m+n) == dgr for (m, n) in labels])
                parts.setdefault(N, {})[dgr] = traj[:, 0:L] @ (w*mask)
    return (out, parts) if by_degree else out

def kerr_exact_time(kappa, Dc, eta, U, ts, Nc=30):
    A_ = destroy(Nc); Ad = A_.conj().T
    H = Dc*Ad@A_ + 1j*(eta*Ad - np.conj(eta)*A_) + U/2*Ad@Ad@A_@A_
    Lv = liouvillian(H, [np.sqrt(kappa)*A_])
    rho0 = np.zeros((Nc, Nc), complex); rho0[0, 0] = 1
    v0 = rho0.reshape(-1, order='F')
    res = []
    for t in ts:
        r = (expm(Lv*t) @ v0).reshape(Nc, Nc, order='F')
        res.append(np.trace(Ad@A_@r).real)
    return np.array(res)

if __name__ == "__main__":
    ts = np.linspace(0, 6, 13)
    kap, Dc, eta, U = 1.0, -1.0, 1.0, 0.15
    ser = kerr_finite_time(kap, Dc, eta, U, ts, 4, lambda a, ad: mul(ad, a))
    ex = kerr_exact_time(kap, Dc, eta, U, ts)
    ps = np.cumsum(ser.real, axis=1)
    for i, t in enumerate(ts):
        print(round(t, 2), round(ex[i], 5), np.round(ps[i], 5))
