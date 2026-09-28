"""Second-order cumulant (HFB) closure for a disordered Kerr chain, site basis."""
import numpy as np
from scipy.integrate import solve_ivp

def hfb_chain(K, kappa, detunings, J, eta, U, ring=True, T=400.0):
    Om = np.diag(np.asarray(detunings, float)).astype(complex)
    for j in range(K-1): Om[j, j+1] = Om[j+1, j] = -J
    if ring and K > 2: Om[0, K-1] = Om[K-1, 0] = -J
    Gam = kappa/2*np.eye(K) + 1j*Om
    eta = np.asarray(eta, complex)*np.ones(K)
    def unpack(y):
        al = y[:K] + 1j*y[K:2*K]
        N = (y[2*K:2*K+K*K] + 1j*y[2*K+K*K:2*K+2*K*K]).reshape(K, K)
        M = (y[2*K+2*K*K:2*K+3*K*K] + 1j*y[2*K+3*K*K:]).reshape(K, K)
        return al, N, M
    def rhs(t, y):
        al, N, M = unpack(y)
        nl = np.real(np.diag(N)); ml = np.diag(M)
        A = 2*U*(np.abs(al)**2 + nl); B = U*(al**2 + ml)
        Ge = Gam + 1j*np.diag(A); Bm = np.diag(B)
        dal = -Gam@al + eta - 1j*U*(np.conj(al)*al**2 + np.conj(al)*ml + 2*al*nl)
        dN = -np.conj(Ge)@N - N@Ge.T + 1j*np.conj(Bm)@M - 1j*np.conj(M)@Bm
        dM = -Ge@M - M@Ge.T - 1j*Bm@N - 1j*(N.T + np.eye(K))@Bm
        return np.concatenate([dal.real, dal.imag, dN.real.ravel(), dN.imag.ravel(), dM.real.ravel(), dM.imag.ravel()])
    a0 = np.linalg.solve(Gam, eta)
    y0 = np.concatenate([a0.real, a0.imag, np.zeros(4*K*K)])
    sol = solve_ivp(rhs, (0, T), y0, method='LSODA', rtol=1e-10, atol=1e-12)
    al, N, M = unpack(sol.y[:, -1])
    out = {}
    for j in range(K):
        A2 = abs(al[j])**2; n = N[j, j].real; m = M[j, j]
        ntot = A2 + n
        n2 = A2**2 + 4*A2*n + np.conj(al[j])**2*m + al[j]**2*np.conj(m) + 2*n**2 + abs(m)**2
        n3 = (abs(al[j])**6 + 9*A2**2*n + 3*A2*(np.conj(al[j])**2*m + al[j]**2*np.conj(m)) + 18*A2*n**2 + 9*A2*abs(m)**2
              + 9*n*(np.conj(al[j])**2*m + al[j]**2*np.conj(m)) + 6*n**3 + 9*n*abs(m)**2)
        out[j] = dict(n=ntot.real, g2=(n2/ntot**2).real, g3=(n3/ntot**3).real)
    return out, np.max(np.abs(rhs(0, sol.y[:, -1])))
