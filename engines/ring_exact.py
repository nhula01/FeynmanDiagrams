import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
def destroy_s(n): return sp.diags(np.sqrt(np.arange(1, n)), 1, format='csr')
def ring_H(K, Delta, J, eta, U, nmax, ring=True):
    I = sp.identity(nmax, format='csr'); d = destroy_s(nmax)
    def op(o, site):
        out = sp.identity(1, format='csr')
        for k in range(K):
            out = sp.kron(out, o if k == site else I, format='csr')
        return out
    A = [op(d, k) for k in range(K)]
    Ad = [x.conj().T.tocsr() for x in A]
    D = nmax**K
    H = sp.csr_matrix((D, D), dtype=complex)
    bonds = ([(k, (k+1) % K) for k in range(K)] if K > 2 else [(0, 1)]) if ring else [(k, k+1) for k in range(K-1)]
    Delta = np.asarray(Delta, float)*np.ones(K)
    for k in range(K):
        H = H + Delta[k]*(Ad[k]@A[k]) + (U/2)*(Ad[k]@Ad[k]@A[k]@A[k]) + 1j*(eta*Ad[k] - np.conj(eta)*A[k])
    for (k, l) in bonds:
        H = H - J*(Ad[k]@A[l] + Ad[l]@A[k])
    return H, A, Ad
def ring_exact(K, kappa, Delta, J, eta, U, nmax, want=("n0",)):
    H, A, Ad = ring_H(K, Delta, J, eta, U, nmax)
    D = nmax**K; I = sp.identity(D, format='csr')
    L = -1j*(sp.kron(I, H) - sp.kron(H.T, I))
    for c in A:
        c = np.sqrt(kappa)*c; cd = c.conj().T.tocsr(); cdc = (cd@c).tocsr()
        L = L + sp.kron(c.conj(), c) - 0.5*sp.kron(I, cdc) - 0.5*sp.kron(cdc.T, I)
    L = L.tolil(); L[0, :] = np.eye(D).reshape(-1, order='F')
    rhs = np.zeros(D*D, complex); rhs[0] = 1
    v = spl.spsolve(sp.csc_matrix(L), rhs)
    rho = v.reshape(D, D, order='F')
    out = {}
    if "n0" in want: out["n0"] = np.trace((Ad[0]@A[0]).toarray()@rho).real
    if "a0" in want: out["a0"] = np.trace(A[0].toarray()@rho)
    if "nn" in want and K > 1:
        out["nn"] = np.trace((Ad[0]@A[1]).toarray()@rho)
    if "g2" in want:
        out["n0n0"] = np.trace((Ad[0]@Ad[0]@A[0]@A[0]).toarray()@rho).real
        if K > 1:
            out["n0n1"] = np.trace((Ad[0]@Ad[1]@A[1]@A[0]).toarray()@rho).real
            out["n1"] = np.trace((Ad[1]@A[1]).toarray()@rho).real
    return out
