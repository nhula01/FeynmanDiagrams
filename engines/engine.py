"""
Diagrammatic steady-state perturbation engine for Lindblad master equations.

Heisenberg-side operators on (cavity mode) x (optional d-level emitter) are
stored as dicts  {(m,n): d x d complex matrix}, meaning  sum (b^dag)^m b^n (x) A_mn.
The order-N steady-state correction is
    <O>^(N) = Tr[ ((L_int^dag G)^N O)_{00} rho_a ]
with G = (-L_0^dag)^+ the pseudo-inverse that kills the steady-state direction.
This is exactly the sum of all N-vertex diagrams of the Feynman rules.
"""
import numpy as np
from math import comb, factorial

# ---------------------------------------------------------------- algebra ----
def add(*ops):
    out = {}
    for op in ops:
        for k, v in op.items():
            out[k] = out.get(k, 0) + v
    return out

def scale(c, op):
    return {k: c * v for k, v in op.items()}

def mul(A, B):
    """normal-ordered product of two operators"""
    out = {}
    for (m, n), X in A.items():
        for (p, q), Y in B.items():
            XY = X @ Y
            for k in range(min(n, p) + 1):
                w = factorial(k) * comb(n, k) * comb(p, k)
                key = (m + p - k, n + q - k)
                out[key] = out.get(key, 0) + w * XY
    return out

def comm(A, B):
    return add(mul(A, B), scale(-1, mul(B, A)))

def prune(op, tol=1e-300, maxdeg=None):
    return {k: v for k, v in op.items()
            if np.max(np.abs(v)) > tol and (maxdeg is None or sum(k) <= maxdeg)}

class Model:
    """
    d      : emitter dimension (1 = cavity only, 2 = two-level emitter)
    z      : complex cavity rate  kappa/2 + i Delta_c   (L0^dag O_mn = -(z n + z* m) O_mn)
    La_dag : function X -> L_a^dag X  (emitter Heisenberg generator in L0), d x d
    rho_a  : emitter steady state of L0
    Hint   : interaction Hamiltonian as operator dict
    """
    def __init__(self, d, z, La_dag, rho_a, Hint, vertex_drop):
        self.d, self.z, self.rho_a, self.Hint = d, z, rho_a, Hint
        self.vertex_drop = vertex_drop
        I = np.eye(d)
        # superoperator matrix of La_dag on column-stacked vec
        M = np.zeros((d * d, d * d), complex)
        for j in range(d * d):
            E = np.zeros(d * d, complex); E[j] = 1
            M[:, j] = La_dag(E.reshape(d, d, order='F')).reshape(-1, order='F')
        self.Ma = M
        P = np.outer(I.reshape(-1, order='F'), rho_a.T.reshape(-1, order='F'))  # X -> Tr[X rho] I
        self.P = P
        self.inv00 = -np.linalg.inv(M - P) @ (np.eye(d * d) - P)
        self._cache = {}

    def G(self, op):
        out = {}
        for (m, n), Y in op.items():
            if (m, n) == (0, 0):
                K = self.inv00
            else:
                c = self.z * n + np.conj(self.z) * m
                key = (m, n)
                if key not in self._cache:
                    self._cache[key] = np.linalg.inv(c * np.eye(self.d ** 2) - self.Ma)
                K = self._cache[key]
            X = (K @ Y.reshape(-1, order='F')).reshape(self.d, self.d, order='F')
            out[(m, n)] = X
        return out

    def V(self, op):
        return scale(1j, comm(self.Hint, op))

    def series(self, O, Nmax):
        """return list of <O>^(N), N = 0..Nmax"""
        vals = []
        X = O
        vals.append(np.trace(X.get((0, 0), np.zeros((self.d, self.d))) @ self.rho_a))
        for N in range(1, Nmax + 1):
            X = self.V(self.G(X))
            X = prune(X, maxdeg=self.vertex_drop * (Nmax - N))
            vals.append(np.trace(X.get((0, 0), np.zeros((self.d, self.d))) @ self.rho_a))
        return np.array(vals)

def cav(m, n, A):
    return {(m, n): np.array(A, dtype=complex)}

# ---------------------------------------------------------------- exact ------
def destroy(N):
    return np.diag(np.sqrt(np.arange(1, N)), 1).astype(complex)

def liouvillian(H, cops):
    D = H.shape[0]; I = np.eye(D)
    L = -1j * (np.kron(I, H) - np.kron(H.T, I))
    for c in cops:
        cd = c.conj().T; cdc = cd @ c
        L += np.kron(c.conj(), c) - 0.5 * np.kron(I, cdc) - 0.5 * np.kron(cdc.T, I)
    return L

def steady(H, cops):
    import scipy.sparse as sp, scipy.sparse.linalg as spl
    D = H.shape[0]
    L = liouvillian(H, cops)
    tr = np.eye(D).reshape(-1, order='F')
    A = L.copy(); A[0, :] = tr
    rhs = np.zeros(D * D, complex); rhs[0] = 1
    if D * D <= 4000:
        v = np.linalg.solve(A, rhs)
    else:
        v = spl.spsolve(sp.csc_matrix(A), rhs)
    return v.reshape(D, D, order='F')
