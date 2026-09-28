import os
"""K=3 validation of counting statistics: exact tilted Liouvillian vs Gaussian FCS vs diagrams."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, json, time, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from engine import destroy
from hfb_chain import hfb_chain

def ops(K, Nc):
    I = sp.identity(Nc, format='csr'); d = sp.csr_matrix(destroy(Nc))
    A = []
    for k in range(K):
        o = sp.identity(1, format='csr')
        for j in range(K): o = sp.kron(o, d if j == k else I, format='csr')
        A.append(o.tocsr())
    return A, [x.conj().T.tocsr() for x in A]

def theta_exact(K, kappa, det, J, eta, U, chi, site, Nc, ring=True):
    A, Ad = ops(K, Nc); D = Nc**K
    H = sp.csr_matrix((D, D), dtype=complex)
    for k in range(K):
        H = H + det[k]*(Ad[k]@A[k]) + (U/2)*(Ad[k]@Ad[k]@A[k]@A[k]) + 1j*(eta[k]*Ad[k] - np.conj(eta[k])*A[k])
    bonds = ([(k, (k+1) % K) for k in range(K)] if K > 2 else [(0, 1)]) if ring else [(k, k+1) for k in range(K-1)]
    for (k, l) in bonds: H = H - J*(Ad[k]@A[l] + Ad[l]@A[k])
    Iq = sp.identity(D, format='csr')
    L = -1j*(sp.kron(Iq, H) - sp.kron(H.T, Iq))
    for k in range(K):
        c = (np.sqrt(kappa)*A[k]).tocsr(); cd = c.conj().T.tocsr(); cdc = (cd@c).tocsr()
        tilt = np.exp(chi) if k == site else 1.0
        L = L + tilt*sp.kron(c.conj(), c) - 0.5*sp.kron(Iq, cdc) - 0.5*sp.kron(cdc.T, Iq)
    w = spl.eigs(sp.csr_matrix(L), k=1, which='LR', return_eigenvectors=False, maxiter=20000, tol=1e-10)
    return float(np.real(w[0]))

def cumulants(f, h=0.05):
    th = np.array([f(x) for x in [-2*h, -h, 0.0, h, 2*h]])
    c1 = (th[3]-th[1])/(2*h); c2 = (th[3]-2*th[2]+th[1])/h**2
    c3 = (th[4]-2*th[3]+2*th[1]-th[0])/(2*h**3)
    return c1, c2/c1, c3/c1
