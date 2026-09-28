import math
"""Counting statistics of an interacting multimode array (tensor engine + counting field).

Tilted generator in the displaced frame, counting the photons leaving site c:
    L_chi^dag X = L_0^dag X + eps [ a_c^dag X a_c ]_displaced ,   eps = kappa_c (e^chi - 1)
In the normal-mode (b) basis, a_c = alpha_c + sum_q V[c,q] beta_q, so the tilt piece is
    eps [ (B_c^dag) X (B_c) + alpha_c^* X B_c + alpha_c B_c^dag X + |alpha_c|^2 X ],
with B_c = sum_q V[c,q] beta_q.  Every term raises the total degree or leaves it alone,
so the generator is upper triangular:  theta_0 = eps |alpha_c|^2,
the propagator is the untilted one, and the left eigenoperator Y_00 is the ladder sum.
"""
import numpy as np
import tensor_engine as TE

class TiltedChain(TE.Chain):
    def set_tilt(self, site, chi, kappa):
        self.eps = kappa*(np.exp(chi)-1.0)
        self.c = site
        self.theta0 = self.eps*abs(self.alpha_site[site])**2
        self.Bc = {(0, 1): self.V[site, :].copy()}                 # B_c (annihilation part)
        self.Bcd = {(1, 0): self.V[site, :].conj().copy()}         # B_c^dag
        self.alc = self.alpha_site[site]
    def raise_op(self, X, maxdeg):
        """the raising part P of the tilt"""
        out = {}
        for (p, q), T in X.items():
            for A, B, w in [(self.Bcd, self.Bc, 1.0), (None, self.Bc, np.conj(self.alc)), (self.Bcd, None, self.alc)]:
                Y = {(p, q): T}
                if A is not None: Y = TE.mul(A, Y, maxdeg)
                if B is not None: Y = TE.mul(Y, B, maxdeg)
                for k, v in Y.items():
                    if sum(k) <= maxdeg: out[k] = out.get(k, 0) + self.eps*w*v
        return out
    def diag_inv(self, X):
        out = {}
        for (p, q), T in X.items():
            if (p, q) == (0, 0): continue
            out[(p, q)] = T/self.Lam(p, q)
        return out
    def resolvent(self, X, maxdeg, tol=1e-13, nmax=60):
        out = self.diag_inv(X); term = out
        for _ in range(nmax):
            term = self.diag_inv(self.raise_op(term, maxdeg))
            term = {k: v for k, v in term.items() if np.max(np.abs(v)) > tol}
            if not term: break
            out = TE.add(out, term)
        return out
    def Y00(self, maxdeg, tol=1e-13, nmax=60):
        Y = {(0, 0): np.array(1.0+0j)}; term = Y
        for _ in range(nmax):
            new = self.diag_inv(self.raise_op(term, maxdeg))
            new = {k: v for k, v in new.items() if np.max(np.abs(v)) > tol}
            if not new: break
            Y = TE.add(Y, new); term = new
        return Y
    def theta_series(self, Nmax, maxdeg=8):
        Y = self.Y00(maxdeg)
        th = [complex(self.theta0)]
        psi = [Y]
        for N in range(1, Nmax + 1):
            X = self.Vx(psi[N-1], maxdeg=maxdeg)
            c00 = complex(np.asarray(X.get((0, 0), 0)))
            th.append(c00)
            # Rayleigh-Schroedinger: remove the Y_00 direction and the renormalization terms theta_k psi_{N-k}
            X = TE.add(X, TE.scale(-c00, Y))
            for k in range(1, N):
                X = TE.add(X, TE.scale(-th[k], psi[N-k]))
            X = TE.prune(self.resolvent(X, maxdeg), maxdeg)
            psi.append(X)
        return np.array(th)

def cumulants(build, chis, Nmax, maxdeg=8, h=0.06):
    """returns (c1, Fano, c3/c1) from theta(chi) by central differences"""
    th = []
    for x in [-2*h, -h, 0.0, h, 2*h]:
        ch = build(x); th.append(np.cumsum(ch.theta_series(Nmax, maxdeg)).real[Nmax])
    th = np.array(th)
    c1 = (th[3]-th[1])/(2*h); c2 = (th[3]-2*th[2]+th[1])/h**2
    c3 = (th[4]-2*th[3]+2*th[1]-th[0])/(2*h**3)
    return c1, c2/c1, c3/c1


def contour_cumulants(theta_of_chi, r=0.08, M=32):
    """(c1, c2/c1, c3/c1) from an analytic theta(chi) by a discrete Cauchy integral on |chi| = r.
    theta_of_chi must accept complex chi. Keep r below the distance to the nearest eigenvalue crossing."""
    phi = 2*np.pi*np.arange(M)/M
    th = np.array([theta_of_chi(r*np.exp(1j*p)) for p in phi])
    c = [np.real(np.mean(th*np.exp(-1j*n*phi)))*math.factorial(n)/r**n for n in (1, 2, 3)]
    return c[0], c[1]/c[0], c[2]/c[0]


def cumulants_contour(build, Nmax, maxdeg=8, r=0.08, M=32):
    """contour version of cumulants(); build(chi) must accept complex chi"""
    return contour_cumulants(lambda x: np.sum(build(x).theta_series(Nmax, maxdeg)), r, M)
