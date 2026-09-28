"""Two-time correlators and spectra from the diagrammatic rules.

  G(s) = int_0^infty dtau e^{-s tau} <A(tau) B(0)>_ss
       = Tr[ ((s - L^dag)^{-1} A) B rho_ss ]

Both factors are diagram sums:
  * (s - L^dag)^{-1} A  is the Laplace-space Heisenberg expansion of the
    observable A, i.e. the finite-time rules with propagator 1/(s + z n + z* m)
    and no projection on the (0,0) label;
  * Tr[O_pq B rho_ss] is a steady-state moment, i.e. the usual diagram sum with
    external operator O_pq B.
The product of the two series, truncated at total order N, is the correlator.
"""
import numpy as np
from engine import add, scale, mul, cav, Model

def resolvent_series(mod, A, Nmax, s, project=False):
    """[ (s-L0^dag - Lint^dag)^{-1} A ]  as a list of order-N label dictionaries."""
    def G(X):
        out = {}
        for (m, n), v in X.items():
            if project and (m, n) == (0, 0):
                continue                                   # traceless complement
            lam = s + mod.z*n + np.conj(mod.z)*m          # (s - Lambda_mn)
            K = np.linalg.inv(lam*np.eye(mod.d**2) - mod.Ma) if mod.d > 1 else None
            if mod.d == 1:
                out[(m, n)] = v/lam
            else:
                out[(m, n)] = (K @ v.reshape(-1, order='F')).reshape(mod.d, mod.d, order='F')
        return out
    terms = []
    X = G(A)
    terms.append(X)
    for _ in range(Nmax):
        X = G(mod.V(X))
        X = {k: v for k, v in X.items() if np.max(np.abs(v)) > 0}
        terms.append(X)
    return terms

def correlator(mod, A, B, Nmax, s, Bmoment, project=False):
    """G(s) to total order Nmax.  Bmoment(p, q) returns the series (length Nmax+1)
    of the steady-state moment Tr[O_pq B rho_ss]."""
    terms = resolvent_series(mod, A, Nmax, s, project=project)
    tot = np.zeros(Nmax+1, complex)
    for NA, X in enumerate(terms):
        if NA > Nmax: break
        for (p, q), v in X.items():
            mom = Bmoment(p, q, Nmax-NA)
            c = v[0, 0] if mod.d == 1 else None
            for NB, mv in enumerate(mom):
                tot[NA+NB] += c*mv
    return tot
