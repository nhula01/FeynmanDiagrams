"""Steady state of the ring by long-time evolution of the Liouvillian (sparse)."""
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply
from ring_exact import ring_H
def ring_steady_evolve(K, kappa, Delta, J, eta, U, nmax, T=60.0, steps=6, ring=True):
    H, A, Ad = ring_H(K, Delta, J, eta, U, nmax, ring=ring)
    D = nmax**K; I = sp.identity(D, format='csr')
    L = -1j*(sp.kron(I, H) - sp.kron(H.T, I))
    for c in A:
        c = (np.sqrt(kappa)*c).tocsr(); cd = c.conj().T.tocsr(); cdc = (cd@c).tocsr()
        L = L + sp.kron(c.conj(), c) - 0.5*sp.kron(I, cdc) - 0.5*sp.kron(cdc.T, I)
    L = sp.csr_matrix(L)
    rho = np.zeros((D, D), complex); rho[0, 0] = 1
    v = rho.reshape(-1, order='F')
    for _ in range(steps):
        v = expm_multiply(L*(T/steps), v)
        v /= np.trace(v.reshape(D, D, order='F')).real
    rho = v.reshape(D, D, order='F')
    n0 = (Ad[0]@A[0]); n00 = (Ad[0]@Ad[0]@A[0]@A[0]); n01 = (Ad[0]@Ad[1]@A[1]@A[0]); n000 = (Ad[0]@Ad[0]@Ad[0]@A[0]@A[0]@A[0])
    tr = lambda O: np.sum(O.multiply(rho.T)) if sp.issparse(O) else np.trace(O@rho)
    return dict(n0=tr(n0).real, n0n0=tr(n00).real, n0n1=tr(n01).real, n3=tr(n000).real,
                a0=np.sum(A[0].multiply(rho.T)), nn=np.sum((Ad[0]@A[1]).multiply(rho.T)))
