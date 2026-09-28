"""Diagrammatic expansion for a driven-dissipative Bose-Hubbard ring.

L0 : K coupled linear cavities, uniform drive, uniform loss  (Gaussian)
Lint: on-site Kerr term, expanded in the displaced momentum modes.

Operators are sparse normal-ordered monomials in the momentum modes:
   key = (m, n) with m, n dicts {mode: power}, meaning prod (b_q^dag)^m_q prod b_q^n_q
"""
import numpy as np
from math import comb, factorial
from itertools import product

# ------------------------------------------------------------------ algebra --
def key(m, n):
    return (tuple(sorted((q, c) for q, c in m.items() if c)),
            tuple(sorted((q, c) for q, c in n.items() if c)))
def todict(t):
    return {q: c for q, c in t}
def add(A, B):
    out = dict(A)
    for k, v in B.items():
        out[k] = out.get(k, 0) + v
        if out[k] == 0: del out[k]
    return out
def scale(c, A):
    return {k: c*v for k, v in A.items()} if c != 0 else {}
def mul(A, B):
    """normal-ordered product: contractions between n of A and m of B, mode by mode"""
    out = {}
    for (mA, nA), x in A.items():
        for (mB, nB), y in B.items():
            nAd, mBd = todict(nA), todict(mB)
            modes = [q for q in nAd if q in mBd]
            ranges = [range(min(nAd[q], mBd[q])+1) for q in modes]
            for ks in product(*ranges) if modes else [()]:
                w = 1.0
                m = todict(mA); n = todict(nB)
                for q, c in mBd.items(): m[q] = m.get(q, 0) + c
                for q, c in nAd.items(): n[q] = n.get(q, 0) + c
                for q, k in zip(modes, ks):
                    w *= factorial(k)*comb(nAd[q], k)*comb(mBd[q], k)
                    m[q] -= k; n[q] -= k
                kk = key(m, n)
                out[kk] = out.get(kk, 0) + w*x*y
    return {k: v for k, v in out.items() if v != 0}
def comm(A, B):
    return add(mul(A, B), scale(-1, mul(B, A)))

# ------------------------------------------------------------------- model ---
class Ring:
    """H = sum_j [Delta n_j + (U/2) a_j^dag^2 a_j^2] - J sum_j (a_j^dag a_{j+1} + h.c.)
           + i sum_j (eta a_j^dag - h.c.),   kappa D[a_j]"""
    def __init__(self, K, kappa, Delta, J, eta, U, eps=None):
        self.K, self.U = K, U
        qs = 2*np.pi*np.arange(K)/K
        self.eps = (Delta - 2*J*np.cos(qs)) if eps is None else np.asarray(eps, float)  # mode detunings
        self.z = kappa/2 + 1j*self.eps                 # complex rates, per mode
        self.alpha = eta/(kappa/2 + 1j*self.eps[0])   # uniform displacement (q=0)
        self.Hint = self._kerr()

    def _kerr(self):
        """(U/2) sum_j (a_j^dag)^2 a_j^2 in displaced momentum modes."""
        K, U, al = self.K, self.U, self.alpha
        alc = np.conj(al); H = {}
        def cnt(ms):
            d = {}
            for q in ms: d[q] = d.get(q, 0) + 1
            return d
        def put(ms, ns, c):
            k = key(cnt(ms), cnt(ns)); H[k] = H.get(k, 0) + c
        put([], [], U/2 * K*abs(al)**4)                                   # constant
        put([0], [], U/2 * 2*alc*al**2*np.sqrt(K))                        # b0^dag
        put([], [0], U/2 * 2*alc**2*al*np.sqrt(K))                        # b0
        for q in range(K):
            put([q], [q], U/2 * 4*abs(al)**2)                             # b^dag b
            put([q, (-q) % K], [], U/2 * al**2)                           # b^dag b^dag
            put([], [q, (-q) % K], U/2 * alc**2)                          # b b
        for q1 in range(K):
            for q2 in range(K):
                q3 = (q1+q2) % K
                put([q1, q2], [q3], U/2 * 2*al/np.sqrt(K))
                put([q3], [q1, q2], U/2 * 2*alc/np.sqrt(K))
        for q1 in range(K):
            for q2 in range(K):
                for q3 in range(K):
                    q4 = (q1+q2-q3) % K
                    put([q1, q2], [q3, q4], U/2 / K)
        return {k: v for k, v in H.items() if v != 0}

    def Lam(self, k):
        m, n = todict(k[0]), todict(k[1])
        return -(sum(self.z[q]*c for q, c in n.items())
                 + sum(np.conj(self.z[q])*c for q, c in m.items()))

    def G0(self, A):
        out = {}
        for k, v in A.items():
            lam = self.Lam(k)
            if abs(lam) > 1e-12:
                out[k] = -v/lam          # (0 - L0)^-1 ; kills the (0,0) label
        return out

    def V(self, A):
        return scale(1j, comm(self.Hint, A))

    def series(self, O, Nmax):
        vals = [O.get(((), ()), 0)]
        X = O
        for _ in range(Nmax):
            X = self.V(self.G0(X))
            vals.append(X.get(((), ()), 0))
        return np.array(vals)

    # ---- observables on site j, in displaced momentum modes -----------------
    def site_ops(self, j=0):
        K, al = self.K, self.alpha
        ph = np.exp(2j*np.pi*np.arange(K)*j/K)/np.sqrt(K)
        a = {((), ()): al}
        for q in range(K):
            a[key({}, {q: 1})] = ph[q]
        ad = {((), ()): np.conj(al)}
        for q in range(K):
            ad[key({q: 1}, {})] = np.conj(ph[q])
        return a, ad
