"""Counting statistics of a Gaussian (HFB) model of a Kerr array, by tilted moment equations.

The tilted generator is quadratic, so the unnormalized state stays Gaussian and the
normalized tilted moments obey a closed (Riccati) system.  With xi = delta a - m the
fluctuation about the tilted mean, mu = alpha_HFB + m the total mean, and
    n_ij = <xi_i^dag xi_j>,   q_ij = <xi_i xi_j>,   eps = kappa_c (e^chi - 1),
    dm   = -Gamma m - i B m*            + eps [ mu_c* q_.c + mu_c n_c. ]
    dn   = -Gamma* n - n Gamma + i B* q - i q* B + eps [ q*_.c q_.c^T + n_.c n_c. ]
    dq   = -Gamma q - q Gamma^T - i B n - i (n^T + 1) B + eps [ n_c. (x) q_.c + q_.c (x) n_c. ]
    theta(chi) = eps ( |mu_c|^2 + n_cc ).
Gamma = kappa/2 + i(Omega + diag(A)) and B are the HFB coefficients.
"""
import numpy as np
from scipy.integrate import solve_ivp
from hfb_lattice import hfb_lattice

def gaussian_theta(Om, kappa, eta, U, chi, site=0, T=800.0):
    K = Om.shape[0]
    al, N0, M0, res = hfb_lattice(Om, kappa, eta, U)
    A = 2*U*(np.abs(al)**2 + np.real(np.diag(N0)))
    B = U*(al**2 + np.diag(M0))
    Gam = kappa/2*np.eye(K) + 1j*(Om + np.diag(A))
    eps = kappa*(np.exp(chi) - 1.0)
    Bm = np.diag(B)
    def pack(m, n, q):
        return np.concatenate([m.real, m.imag, n.real.ravel(), n.imag.ravel(), q.real.ravel(), q.imag.ravel()])
    def unpack(y):
        m = y[:K] + 1j*y[K:2*K]
        n = (y[2*K:2*K+K*K] + 1j*y[2*K+K*K:2*K+2*K*K]).reshape(K, K)
        q = (y[2*K+2*K*K:2*K+3*K*K] + 1j*y[2*K+3*K*K:]).reshape(K, K)
        return m, n, q
    c = site
    def rhs(t, y):
        m, n, q = unpack(y)
        mu = al + m
        dm = -Gam@m - 1j*B*np.conj(m) + eps*(np.conj(mu[c])*q[:, c] + mu[c]*n[c, :])
        dn = (-np.conj(Gam)@n - n@Gam + 1j*np.conj(Bm)@q - 1j*np.conj(q)@Bm
              + eps*(np.outer(np.conj(q[:, c]), q[:, c]) + np.outer(n[:, c], n[c, :])))
        dq = (-Gam@q - q@Gam.T - 1j*Bm@n - 1j*(n.T + np.eye(K))@Bm
              + eps*(np.outer(n[c, :], q[:, c]) + np.outer(q[:, c], n[c, :])))
        return pack(dm, dn, dq)
    y0 = pack(np.zeros(K, complex), N0.copy(), M0.copy())
    sol = solve_ivp(rhs, (0, T), y0, method='LSODA', rtol=1e-11, atol=1e-13)
    m, n, q = unpack(sol.y[:, -1])
    mu = al + m
    resid = np.max(np.abs(rhs(0, sol.y[:, -1])))
    return float(np.real(eps*(abs(mu[c])**2 + n[c, c]))), resid

def cumulants(theta_fn, h=0.05):
    th = np.array([theta_fn(x) for x in [-2*h, -h, 0.0, h, 2*h]])
    c1 = (th[3]-th[1])/(2*h); c2 = (th[3]-2*th[2]+th[1])/h**2
    c3 = (th[4]-2*th[3]+2*th[1]-th[0])/(2*h**3)
    return c1, c2/c1, c3/c1
