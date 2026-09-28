"""Counting statistics from the diagrammatic rules (single driven Kerr mode).

Tilted generator, displaced frame:
    L_chi^dag X = L_0^dag X + eps [ b^dag X b + a* X b + a b^dag X + |a|^2 X ],  eps = kappa(e^chi - 1)
which only raises the label, so
  * the spectrum is Lambda^chi_mn = -(z n + z* m) + eps|alpha|^2,
  * theta_0(chi) = eps |alpha|^2 (Poissonian),
  * the left eigenoperator Y_00 is the ladder resummation of the raising part,
  * (theta_0 - Lambda^chi_mn)^{-1} = (z n + z* m)^{-1}: the untilted propagator.
Because the ladder only raises, the Y_00 component of any operator equals its O_00
component, so the interacting corrections are

    theta^(N) = [ (L_int^dag G_chi)^N Y_00 ]_00 ,   G_chi = ladder-resummed resolvent.
"""
import numpy as np
from engine import add, scale, mul, cav, comm

def ladder(X, alpha, eps, z, maxdeg, tol=1e-14):
    """apply (theta_0 - L_chi^dag)^{-1} = sum_j [(diag)^{-1} P]^j (diag)^{-1}"""
    I1 = np.eye(1)
    def diag_inv(A):
        out = {}
        for (m, n), v in A.items():
            d = z*n + np.conj(z)*m
            if abs(d) > 1e-14: out[(m, n)] = v/d
        return out
    def P(A):                      # raising part of the tilt
        out = {}
        for (m, n), v in A.items():
            for key, w in [((m+1, n+1), 1.0), ((m, n+1), np.conj(alpha)), ((m+1, n), alpha)]:
                if sum(key) <= maxdeg:
                    out[key] = out.get(key, 0) + eps*w*v
        return out
    out = diag_inv(X); term = out
    for _ in range(200):
        term = diag_inv(P(term))
        term = {k: v for k, v in term.items() if np.max(np.abs(v)) > tol}
        if not term: break
        out = add(out, term)
    return out

def Y00(alpha, eps, z, maxdeg, tol=1e-14):
    """left eigenoperator of the tilted Gaussian generator for the leading eigenvalue"""
    I1 = np.eye(1)
    Y = {(0, 0): np.array([[1.0+0j]])}
    term = Y
    for _ in range(200):
        new = {}
        for (m, n), v in term.items():
            for key, w in [((m+1, n+1), 1.0), ((m, n+1), np.conj(alpha)), ((m+1, n), alpha)]:
                if sum(key) <= maxdeg:
                    d = z*key[1] + np.conj(z)*key[0]
                    new[key] = new.get(key, 0) + eps*w*v/d
        new = {k: v for k, v in new.items() if np.max(np.abs(v)) > tol}
        if not new: break
        Y = add(Y, new); term = new
    return Y

def theta_series(kappa, Dc, eta, U, chi, Nmax, maxdeg=14):
    z = kappa/2 + 1j*Dc; alpha = eta/z; eps = kappa*(np.exp(chi)-1)
    I1 = np.eye(1)
    b, bd, one = cav(0, 1, I1), cav(1, 0, I1), cav(0, 0, I1)
    a_ = add(b, scale(alpha, one)); ad_ = add(bd, scale(np.conj(alpha), one))
    HK = scale(U/2, mul(mul(ad_, ad_), mul(a_, a_)))
    V = lambda X: scale(1j, comm(HK, X))
    Y = Y00(alpha, eps, z, maxdeg)
    th = [eps*abs(alpha)**2]                       # theta_0
    psi = [Y]                                      # psi_0 = Y_00
    for N in range(1, Nmax + 1):
        X = V(psi[N-1])
        c00 = complex(X.get((0, 0), np.zeros((1, 1)))[0, 0])
        th.append(c00)
        # Rayleigh-Schroedinger: remove the Y_00 direction and the renormalization terms theta_k psi_{N-k}
        X = add(X, scale(-c00, Y))
        for k in range(1, N):
            X = add(X, scale(-th[k], psi[N-k]))
        X = ladder(X, alpha, eps, z, maxdeg)
        X = {k: v for k, v in X.items() if sum(k) <= maxdeg and np.max(np.abs(v)) > 1e-16}
        psi.append(X)
    return np.array(th)
