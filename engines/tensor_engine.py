"""Dense-tensor diagram engine for a general Gaussian L0 (normal drift) with on-site Kerr.

Modes: normal modes beta_q of the drift Gamma = kappa/2 + i Omega, Omega real symmetric
(on-site detunings + hopping).  Gamma is normal, so the eigenvectors are orthonormal and
contractions conserve the mode label.  Operators are dicts {(p,q): T} with T of shape
(K,)*p + (K,)*q, meaning sum T[i..,j..] beta^dag_{i1}..beta^dag_{ip} beta_{j1}..beta_{jq}.
"""
import numpy as np
from math import comb, factorial
from itertools import permutations

def sym(T, p, q):
    """symmetrize a (p,q) tensor over its creation and annihilation index groups"""
    if T.ndim == 0: return T
    out = np.zeros_like(T); cnt = 0
    for pa in permutations(range(p)):
        for pb in permutations(range(q)):
            axes = list(pa) + [p + b for b in pb]
            out += np.transpose(T, axes); cnt += 1
    return out/cnt

def mul(A, B, maxdeg=None):
    """normal-ordered product of two tensor operators (sectors above maxdeg skipped)"""
    out = {}
    for (p1, q1), T1 in A.items():
        for (p2, q2), T2 in B.items():
            for k in range(min(q1, p2)+1):
                if maxdeg is not None and (p1+p2-k)+(q1+q2-k) > maxdeg: continue
                w = factorial(k)*comb(q1, k)*comb(p2, k)
                # contract last k annihilation indices of T1 with first k creation indices of T2
                if k == 0:
                    R = np.multiply.outer(T1, T2) if (T1.ndim or T2.ndim) else T1*T2
                else:
                    R = np.tensordot(T1, T2, axes=(list(range(p1+q1-k, p1+q1)), list(range(k))))
                # R indices: [T1 cre (p1)] [T1 ann (q1-k)] [T2 cre (p2-k)] [T2 ann (q2)]
                # reorder to [cre: p1 + (p2-k)] [ann: (q1-k) + q2]
                n1c, n1a, n2c, n2a = p1, q1-k, p2-k, q2
                if R.ndim:
                    order = list(range(n1c)) + list(range(n1c+n1a, n1c+n1a+n2c)) + list(range(n1c, n1c+n1a)) + list(range(n1c+n1a+n2c, R.ndim))
                    R = np.transpose(R, order)
                key = (p1+p2-k, q1+q2-k)
                out[key] = out.get(key, 0) + w*R
    return {k: sym(np.asarray(v), *k) for k, v in out.items()}

def add(*ops):
    out = {}
    for o in ops:
        for k, v in o.items(): out[k] = out.get(k, 0) + v
    return out
def scale(c, A): return {k: c*v for k, v in A.items()}
def comm(A, B, maxdeg=None): return add(mul(A, B, maxdeg), scale(-1, mul(B, A, maxdeg)))
def prune(A, maxdeg, tol=1e-300):
    return {k: v for k, v in A.items() if sum(k) <= maxdeg and np.max(np.abs(v)) > tol}

class Chain:
    def __init__(self, K, kappa, detunings, J, eta, U, ring=True):
        self.K, self.U = K, U
        Om = np.diag(np.asarray(detunings, float)).astype(complex)
        for j in range(K-1): Om[j, j+1] = Om[j+1, j] = -J
        if ring and K > 2: Om[0, K-1] = Om[K-1, 0] = -J
        Gam = kappa/2*np.eye(K) + 1j*Om
        eps, Vm = np.linalg.eigh(Om)              # Om real symmetric -> orthonormal modes
        self.V = Vm.astype(complex); self.zeta = kappa/2 + 1j*eps
        eta = np.asarray(eta, complex)*np.ones(K)
        self.alpha_site = np.linalg.solve(Gam, eta)
        self.alpha = self.V.conj().T @ self.alpha_site           # mode-basis displacement
        self.Hint = self._kerr()
    # site operator a_j = alpha_j + sum_q V[j,q] beta_q
    def site_a(self, j):
        return {(0, 0): np.array(self.alpha_site[j]), (0, 1): self.V[j, :].copy()}
    def site_ad(self, j):
        return {(0, 0): np.array(np.conj(self.alpha_site[j])), (1, 0): self.V[j, :].conj().copy()}
    def _kerr(self):
        K, U = self.K, self.U
        H = {}
        for j in range(K):
            ad, a = self.site_ad(j), self.site_a(j)
            term = mul(mul(ad, ad), mul(a, a))
            H = add(H, term)
        return {k: (U/2)*v for k, v in H.items()}
    def Lam(self, p, q):
        """propagator denominators for a (p,q) tensor: sum_j zeta_j + sum_i zeta_i^*"""
        z = self.zeta; zc = np.conj(z)
        shape = [1]*(p+q); D = np.zeros((self.K,)*(p+q), complex)
        for i in range(p):
            s = [1]*(p+q); s[i] = self.K; D = D + zc.reshape(s)
        for j in range(q):
            s = [1]*(p+q); s[p+j] = self.K; D = D + z.reshape(s)
        return D
    def G0(self, A):
        out = {}
        for (p, q), T in A.items():
            if (p, q) == (0, 0): continue
            out[(p, q)] = T/self.Lam(p, q)
        return out
    def Vx(self, A, maxdeg=None): return scale(1j, comm(self.Hint, A, maxdeg))
    def series(self, O, Nmax, drop=2):
        vals = [complex(O.get((0, 0), 0))]
        X = O
        for N in range(1, Nmax+1):
            X = prune(self.Vx(self.G0(X), maxdeg=drop*(Nmax-N)), drop*(Nmax-N))
            vals.append(complex(X.get((0, 0), 0)))
        return np.array(vals)
