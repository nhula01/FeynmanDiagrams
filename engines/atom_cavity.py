"""Atom (two-level) vs classical dipole (harmonic) coupled to a coherently driven cavity.
Diagrammatic series in g with L0 = displaced cavity (x) undriven emitter; the emitter algebra
is either Pauli (d=2) or a truncated boson (classical dipole, d=dc)."""
import numpy as np
from engine import *

def emitter_ops(kind, dc=12):
    if kind == "atom":
        s = np.array([[0, 0], [1, 0]], complex)      # basis (e, g): sigma^- = |g><e|
        return s, 2
    c = destroy(dc)
    return c, dc

def series_cavity(kind, kap, gam, Dc, Da, eta, g, Nmax, observables, dc=12):
    lo, d = emitter_ops(kind, dc)
    lod = lo.conj().T; Id = np.eye(d)
    z = kap/2 + 1j*Dc; alpha = eta/z
    Ha = Da*lod@lo
    La = lambda X: 1j*(Ha@X - X@Ha) + gam*(lod@X@lo - 0.5*(lod@lo@X + X@lod@lo))
    rho_a = np.zeros((d, d), complex)
    rho_a[(1 if kind == "atom" else 0), (1 if kind == "atom" else 0)] = 1   # ground state / vacuum
    b = cav(0, 1, Id); bd = cav(1, 0, Id); one = cav(0, 0, Id)
    a = add(b, scale(alpha, one)); ad = add(bd, scale(np.conj(alpha), one))
    S = cav(0, 0, lo); Sd = cav(0, 0, lod)
    Hint = scale(g, add(mul(ad, S), mul(a, Sd)))
    mod = Model(d, z, La, rho_a, Hint, vertex_drop=2)
    out = {}
    for name, fn in observables.items():
        out[name] = mod.series(fn(a, ad), Nmax)
    return out, alpha

def exact_cavity(kind, kap, gam, Dc, Da, eta, g, Nc=10, dc=10):
    lo, d = emitter_ops(kind, dc)
    A = np.kron(destroy(Nc), np.eye(d)); Ad = A.conj().T
    Lo = np.kron(np.eye(Nc), lo); Lod = Lo.conj().T
    H = Dc*Ad@A + Da*Lod@Lo + g*(Ad@Lo + A@Lod) + 1j*(eta*Ad - np.conj(eta)*A)
    from sparse_ss import steady_sparse
    rho = steady_sparse(H, [np.sqrt(kap)*A, np.sqrt(gam)*Lo])
    n = np.trace(Ad@A@rho).real; n2 = np.trace(Ad@Ad@A@A@rho).real
    return n, n2, np.trace(A@rho)

if __name__ == "__main__":
    kap, gam, g, eta = 1.0, 1.0, 0.35, 0.02
    obs = {"n": lambda a, ad: mul(ad, a), "n2": lambda a, ad: mul(mul(ad, ad), mul(a, a))}
    for D in [0.0, 0.5, 1.0]:
        sa, al = series_cavity("atom", kap, gam, D, D, eta, g, 10, obs)
        sc, _ = series_cavity("cl", kap, gam, D, D, eta, g, 10, obs)
        ea = exact_cavity("atom", kap, gam, D, D, eta, g)
        ec = exact_cavity("cl", kap, gam, D, D, eta, g)
        pa_n = np.cumsum(sa["n"].real); pa_n2 = np.cumsum(sa["n2"].real)
        pc_n = np.cumsum(sc["n"].real); pc_n2 = np.cumsum(sc["n2"].real)
        print("D", D, "exact atom g2", ea[1]/ea[0]**2, " classical g2", ec[1]/ec[0]**2)
        print("   series atom g2 by order", np.round((pa_n2/pa_n**2)[[2,4,6,8,10]], 5))
        print("   series cl g2 by order", np.round((pc_n2/pc_n**2)[[2,4,6,8,10]], 5))
        print("   n atom exact", ea[0], "series", np.round(pa_n[[2,4,8,10]],7), " cl exact", ec[0], np.round(pc_n[[2,4,8,10]],7))
        print("   coeffs n atom-cl", np.round((sa["n"]-sc["n"]).real[:7], 9))
