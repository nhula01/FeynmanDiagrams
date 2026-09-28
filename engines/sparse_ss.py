import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
def steady_sparse(H, cops):
    D = H.shape[0]; I = sp.identity(D, format='csr')
    H = sp.csr_matrix(H)
    L = -1j*(sp.kron(I, H) - sp.kron(H.T, I))
    for c in cops:
        c = sp.csr_matrix(c); cd = c.conj().T; cdc = cd @ c
        L = L + sp.kron(c.conj(), c) - 0.5*sp.kron(I, cdc) - 0.5*sp.kron(cdc.T, I)
    L = sp.lil_matrix(L)
    L[0, :] = np.eye(D).reshape(-1, order='F')
    rhs = np.zeros(D*D, complex); rhs[0] = 1
    return spl.spsolve(sp.csc_matrix(L), rhs).reshape(D, D, order='F')
