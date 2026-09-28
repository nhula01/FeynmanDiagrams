"""Exact steady states and tilted leading eigenvalues in the displaced frame a_j = alpha_j + b_j, alpha_j the linear (U=0)
displacement.  The displacement is exact; only the b-Fock space is truncated, and the b-occupations are small, so the cutoff
converges much faster than in the laboratory Fock basis (where the error decays only like the amplitude of the top level)."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl

def lin_alpha(K, kappa, det, J, eta, ring):
    Om = np.diag(np.asarray(det, float)*np.ones(K)).astype(complex)
    bonds = ([(k, (k+1) % K) for k in range(K)] if K > 2 else [(0, 1)]) if ring else [(k, k+1) for k in range(K-1)]
    for (k, l) in bonds: Om[k, l] = Om[l, k] = -J
    return np.linalg.solve(kappa/2*np.eye(K) + 1j*Om, np.asarray(eta, complex)*np.ones(K)), bonds

def mf_alpha(K, kappa, det, J, eta, U, ring):
    from scipy.optimize import fsolve
    a0, bonds = lin_alpha(K, kappa, det, J, eta, ring)
    Om = np.diag(np.asarray(det, float)*np.ones(K)).astype(complex)
    for (k, l) in bonds: Om[k, l] = Om[l, k] = -J
    eta = np.asarray(eta, complex)*np.ones(K)
    def f(x):
        a = x[:K] + 1j*x[K:]; r = -(kappa/2*np.eye(K) + 1j*Om)@a - 1j*U*np.abs(a)**2*a + eta; return np.concatenate([r.real, r.imag])
    x = fsolve(f, np.concatenate([a0.real, a0.imag]), xtol=1e-14); return x[:K] + 1j*x[K:]

def build(cuts, kappa, det, J, eta, U, ring, count=0, displace=True):
    K = len(cuts); det = np.asarray(det, float)*np.ones(K); eta = np.asarray(eta, complex)*np.ones(K)
    alpha, bonds = lin_alpha(K, kappa, det, J, eta, ring)
    if displace == "mf": alpha = mf_alpha(K, kappa, det, J, eta, U, ring)
    elif not displace: alpha = 0*alpha
    loc = [sp.diags(np.sqrt(np.arange(1, n)), 1, shape=(n, n), format='csr') for n in cuts]
    Is = [sp.identity(n, format='csr') for n in cuts]
    def emb(o, j):
        m = sp.identity(1, format='csr')
        for k in range(K): m = sp.kron(m, o if k == j else Is[k], format='csr')
        return m
    D = int(np.prod(cuts)); Id = sp.identity(D, format='csr')
    A = [(emb(loc[j], j) + alpha[j]*Id).tocsr() for j in range(K)]; Ad = [x.conj().T.tocsr() for x in A]
    H = sp.csr_matrix((D, D), dtype=complex)
    for k in range(K):
        H = H + det[k]*(Ad[k]@A[k]) + (U/2)*(Ad[k]@Ad[k]@A[k]@A[k]) + 1j*(eta[k]*Ad[k] - np.conj(eta[k])*A[k])
    for (k, l) in bonds:
        H = H - J*(Ad[k]@A[l] + Ad[l]@A[k])
    L = -1j*(sp.kron(Id, H) - sp.kron(H.T, Id))
    for c in A:
        c = np.sqrt(kappa)*c; cdc = (c.conj().T@c).tocsr()
        L = L + sp.kron(c.conj(), c) - 0.5*sp.kron(Id, cdc) - 0.5*sp.kron(cdc.T, Id)
    cc = np.sqrt(kappa)*A[count]
    Jc = sp.kron(cc.conj(), cc)
    return sp.csr_matrix(L), sp.csr_matrix(Jc), D, A, Ad

def steady(L, D, method="auto"):
    if method == "spsolve" or (method == "auto" and D*D <= 60000):
        M = L.tolil(); M[0, :] = np.eye(D).reshape(1, -1, order='F'); b = np.zeros(D*D, complex); b[0] = 1
        v = spl.spsolve(sp.csc_matrix(M), b)
    else:
        v = np.zeros(D*D, complex); v[0] = 1
        for _ in range(8): v = spl.expm_multiply(L*(10.0), v); v = v/np.trace(v.reshape(D, D, order='F'))
    return v.reshape(D, D, order='F')

def ring_moments(K, U, levels, kappa=1.0, Delta=-1.0, J=0.4, eta=1.0, displace=True, method="auto"):
    L, Jc, D, A, Ad = build([levels]*K, kappa, Delta, J, eta, U, True, displace=displace)
    R = steady(L, D, method)
    tr = lambda O: complex(np.sum(O.multiply(R.T)))
    n0 = tr(Ad[0]@A[0]).real; n00 = tr(Ad[0]@Ad[0]@A[0]@A[0]).real; n000 = tr(Ad[0]@Ad[0]@Ad[0]@A[0]@A[0]@A[0]).real
    b0 = A[0] - A[0].diagonal().mean()*0  # placeholder
    return dict(n0=n0, n0n0=n00, n3=n000, g2=n00/n0**2, g3=n000/n0**3, levels=levels, D=D)

def theta_lead(L, Jc, D, chi, v0=None, T=40.0, steps=5, dense=None):
    Lc = (L + (np.exp(chi)-1)*Jc)
    if dense or (dense is None and D*D <= 2500):
        w = np.linalg.eigvals(Lc.toarray()); return w[np.argmax(w.real)], 0.0
    tri = np.arange(D)*(D+1); v = v0.copy() if v0 is not None else np.eye(D).reshape(-1, order='F')/D; hist = []
    Lc = Lc.tocsr()
    for s in range(steps):
        v = spl.expm_multiply(Lc*(T/steps), v); v = v/np.sum(v[tri]); hist.append(np.sum((Lc@v)[tri]))
    return hist[-1], abs(hist[-1]-hist[-2])
