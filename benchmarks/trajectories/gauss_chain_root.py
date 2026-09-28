"""Gaussian (HFB) counting cumulants of the eight-site chain with exact chi-derivatives.

theta_G(chi) from the tilted Riccati system of gaussian_fcs_riccati.py at real chi.  The Riccati equations
contain complex conjugates of the tilted moments, so theta_G is not evaluated off the real axis (a complex-chi
contour would not be an analytic continuation).  Instead the five-point central differences are taken at
h = 0.05, 0.025, 0.0125 and Richardson-extrapolated in h^2 twice, which removes the O(h^2) and O(h^4) errors.
The h = 0.05 values (the paper's convention) are stored alongside.  The stationary point is found by Powell's hybrid root finder instead of time integration (identical fixed point, much faster).
usage: gauss_chain_fd.py U1,U2,... -> gauss_chain8.json
"""
import sys, os, json, math, time, numpy as np
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
sys.path.insert(0, os.path.join(_ROOT, 'engines'))
from scipy.integrate import solve_ivp
from scipy.optimize import root
from hfb_lattice import hfb_lattice
RP = json.load(open(os.path.join(_ROOT, 'notebook', 'data', 'revision_params.json')))['chain8_fcs']
K, det, J, eta, c = RP['K'], RP['detunings'], RP['J'], np.array(RP['eta'], complex), RP['count_site']
Om = np.diag(np.asarray(det, float)).astype(complex)
for j in range(K-1): Om[j, j+1] = Om[j+1, j] = -J

def theta_G(U, chi, kappa=1.0, T=800.0, _cache={}):
    if U not in _cache: _cache[U] = hfb_lattice(Om, kappa, eta, U)
    al, N0, M0, res = _cache[U]
    A = 2*U*(np.abs(al)**2 + np.real(np.diag(N0))); B = U*(al**2 + np.diag(M0))
    Gam = kappa/2*np.eye(K) + 1j*(Om + np.diag(A)); eps = kappa*(np.exp(chi) - 1.0); Bm = np.diag(B)
    pack = lambda m, n, q: np.concatenate([m.real, m.imag, n.real.ravel(), n.imag.ravel(), q.real.ravel(), q.imag.ravel()])
    def unpack(y):
        m = y[:K] + 1j*y[K:2*K]
        n = (y[2*K:2*K+K*K] + 1j*y[2*K+K*K:2*K+2*K*K]).reshape(K, K)
        q = (y[2*K+2*K*K:2*K+3*K*K] + 1j*y[2*K+3*K*K:]).reshape(K, K)
        return m, n, q
    def rhs(t, y):
        m, n, q = unpack(y); mu = al + m
        dm = -Gam@m - 1j*B*np.conj(m) + eps*(np.conj(mu[c])*q[:, c] + mu[c]*n[c, :])
        dn = (-np.conj(Gam)@n - n@Gam + 1j*np.conj(Bm)@q - 1j*np.conj(q)@Bm
              + eps*(np.outer(np.conj(q[:, c]), q[:, c]) + np.outer(n[:, c], n[c, :])))
        dq = (-Gam@q - q@Gam.T - 1j*Bm@n - 1j*(n.T + np.eye(K))@Bm
              + eps*(np.outer(n[c, :], q[:, c]) + np.outer(q[:, c], n[c, :])))
        return pack(dm, dn, dq)
    # stationary tilted moments: root of the Riccati right-hand side, continued from the untilted HFB state
    sol = root(lambda y: rhs(0, y), pack(np.zeros(K, complex), N0.copy(), M0.copy()), method='hybr', tol=1e-15)
    m, n, q = unpack(sol.x); mu = al + m
    return complex(eps*(abs(mu[c])**2 + n[c, c])), float(np.max(np.abs(rhs(0, sol.x)))), al

if __name__ == "__main__":
    Us = [float(x) for x in sys.argv[1].split(',')]
    out = dict(description=__doc__, params=RP, results={})
    def fd(th0, thf, h):
        tm2, tm1, tp1, tp2 = thf
        c1 = (tp1-tm1)/(2*h); c2 = (tp1-2*th0+tm1)/h**2; c3 = (tp2-2*tp1+2*tm1-tm2)/(2*h**3)
        return np.array([c1, c2, c3])
    for U in Us:
        t0 = time.time(); rec = dict(U=U); resid = 0.0
        th0, rs, al = theta_G(U, 0.0); th0 = th0.real
        hs = [0.05, 0.025, 0.0125]; C = []
        for h in hs:
            vals = []
            for x in np.array([-2, -1, 1, 2])*h:
                t, rs, _ = theta_G(U, x); vals.append(t.real); resid = max(resid, rs)
            C.append(fd(th0, vals, h))
        C = np.array(C)
        R1 = [(4*C[1]-C[0])/3, (4*C[2]-C[1])/3]; R2 = (16*R1[1]-R1[0])/15
        rat = lambda c: [float(c[0]), float(c[1]/c[0]), float(c[2]/c[0])]
        rec['fd_h0.05'] = rat(C[0]); rec['fd_h0.025'] = rat(C[1]); rec['fd_h0.0125'] = rat(C[2])
        rec['richardson'] = rat(R2); rec['richardson_change_last'] = [abs(a-b) for a, b in zip(rat(R2), rat(R1[1]))]
        rec['max_residual'] = resid; rec['alpha_hfb_site0_sq'] = float(abs(al[c])**2); rec['time'] = time.time()-t0
        out['results'][f'{U:.2f}'] = rec
        print(U, 'richardson', np.round(rec['richardson'], 6), 'h.05', np.round(rec['fd_h0.05'], 6), f"resid {resid:.1e} {rec['time']:.0f}s", flush=True)
        json.dump(out, open('gauss_chain8_root.json', 'w'), indent=1)

