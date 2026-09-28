"""Second-order cumulant (Gaussian) closure for the driven Kerr cavity and the uniform Kerr ring."""
import numpy as np
from scipy.optimize import fsolve

def kerr_cumulant_ss(kappa, Delta, eta, U, J=0.0, x0=None):
    """Uniform solution: single mode with effective detuning Delta-2J (ring, q=0 mean field).
    Variables: alpha (complex), n_c (real), m_c (complex)."""
    z = kappa/2 + 1j*(Delta - 2*J)
    def rhs(v):
        al = v[0]+1j*v[1]; nc = v[2]; mc = v[3]+1j*v[4]
        alc = np.conj(al)
        A3 = alc*al**2 + alc*mc + 2*al*nc
        A4 = alc*al**3 + 3*abs(al)**2*mc + 3*al**2*nc + 3*nc*mc
        n = nc + abs(al)**2; m = mc + al**2
        dal = -z*al + eta - 1j*U*A3
        dn = -kappa*n + eta*alc + np.conj(eta)*al
        dm = -(kappa + 2j*(Delta-2*J))*m + 2*eta*al - 1j*U*(2*A4 + m)
        dnc = dn - 2*np.real(alc*dal)
        dmc = dm - 2*al*dal
        return [dal.real, dal.imag, dnc.real, dmc.real, dmc.imag]
    if x0 is None:
        a0 = eta/z; x0 = [a0.real, a0.imag, 0.0, 0.0, 0.0]
    sol, info, ier, msg = fsolve(rhs, x0, full_output=True, xtol=1e-12)
    al = sol[0]+1j*sol[1]; nc = sol[2]; mc = sol[3]+1j*sol[4]
    n = nc + abs(al)**2
    n2 = abs(al)**4 + 4*abs(al)**2*nc + np.conj(al)**2*mc + al**2*np.conj(mc) + 2*nc**2 + abs(mc)**2
    return dict(alpha=al, n=n, g2=(n2/n**2).real, ok=(ier == 1), x=sol)

def kerr_meanfield(kappa, Delta, eta, U, J=0.0, eps0=None):
    e0 = (Delta-2*J) if eps0 is None else eps0
    z = kappa/2 + 1j*e0
    r = np.roots([U**2, 2*e0*U, abs(z)**2, -abs(eta)**2]) if U > 0 else [abs(eta/z)**2]
    r = np.array(r); return float(np.min(r[np.abs(np.imag(r)) < 1e-9].real))

if __name__ == "__main__":
    import json
    d = json.load(open("numbers.json"))
    for U in [0.05, 0.1, 0.15, 0.2, 0.25, 0.3]:
        c = kerr_cumulant_ss(1.0, -1.0, 1.0, U)
        print(U, "cumulant n", round(c['n'],4), "g2", round(c['g2'],4), c['ok'], " meanfield", round(kerr_meanfield(1.0,-1.0,1.0,U),4))

def ring_cumulant_ss(K, kappa, Delta, J, eta, U, T=400.0, eps=None):
    """Uniform second-order cumulant (Gaussian / HFB) closure for the Kerr ring,
    in momentum space: alpha (q=0 mean), n_q = <b_q^dag b_q>, m_q = <b_q b_{-q}>."""
    from scipy.integrate import solve_ivp
    qs = 2*np.pi*np.arange(K)/K; eps = (Delta - 2*J*np.cos(qs)) if eps is None else np.asarray(eps, float); z0 = kappa/2 + 1j*eps[0]
    def unpack(y):
        al = y[0]+1j*y[1]; n = y[2:2+K]; m = y[2+K:2+2*K] + 1j*y[2+2*K:2+3*K]
        return al, n, m
    def rhs(t, y):
        al, n, m = unpack(y)
        nl = n.mean(); ml = m.mean()
        A = 2*U*(abs(al)**2 + nl); B = U*(al**2 + ml)
        dal = -z0*al + eta - 1j*U*(np.conj(al)*al**2 + np.conj(al)*ml + 2*al*nl)
        dn = -kappa*n + 1j*(np.conj(B)*m - B*np.conj(m))
        dm = -(kappa + 2j*eps)*m - 2j*A*m - 1j*B*(2*n + 1)
        return np.concatenate([[dal.real, dal.imag], dn.real, dm.real, dm.imag])
    y0 = np.zeros(2+3*K); a0 = eta/z0; y0[0], y0[1] = a0.real, a0.imag
    sol = solve_ivp(rhs, (0, T), y0, method='LSODA', rtol=1e-10, atol=1e-12)
    al, n, m = unpack(sol.y[:, -1])
    nl = n.mean(); ml = m.mean()
    ntot = nl + abs(al)**2
    n2 = abs(al)**4 + 4*abs(al)**2*nl + np.conj(al)**2*ml + al**2*np.conj(ml) + 2*nl**2 + abs(ml)**2
    n01 = None
    return dict(alpha=al, n=ntot.real, g2=(n2/ntot**2).real, n_q=n, m_q=m,
                converged=np.max(np.abs(rhs(0, sol.y[:, -1]))) < 1e-7)
