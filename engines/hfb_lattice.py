"""Second-order cumulant (HFB) closure for a Kerr lattice with complex hoppings."""
import numpy as np
from scipy.integrate import solve_ivp

def hfb_lattice(Om, kappa, eta, U, T=600.0):
    K = Om.shape[0]
    Gam = kappa/2*np.eye(K) + 1j*Om
    eta = np.asarray(eta, complex)*np.ones(K)
    def unpack(y):
        al = y[:K] + 1j*y[K:2*K]
        N = (y[2*K:2*K+K*K] + 1j*y[2*K+K*K:2*K+2*K*K]).reshape(K, K)   # N_ij = <da_i^dag da_j>
        M = (y[2*K+2*K*K:2*K+3*K*K] + 1j*y[2*K+3*K*K:]).reshape(K, K)  # M_ij = <da_i da_j>
        return al, N, M
    def rhs(t, y):
        al, N, M = unpack(y)
        nl = np.real(np.diag(N)); ml = np.diag(M)
        A = 2*U*(np.abs(al)**2 + nl)          # site-diagonal
        B = U*(al**2 + ml)
        Ge = Gam + 1j*np.diag(A); Bm = np.diag(B)
        dal = -Gam@al + eta - 1j*U*(np.conj(al)*al**2 + np.conj(al)*ml + 2*al*nl)
        dN = -np.conj(Ge)@N - N@Ge.T + 1j*np.conj(Bm)@M - 1j*np.conj(M)@Bm
        dM = -Ge@M - M@Ge.T - 1j*Bm@N - 1j*(N.T + np.eye(K))@Bm
        return np.concatenate([dal.real, dal.imag, dN.real.ravel(), dN.imag.ravel(),
                               dM.real.ravel(), dM.imag.ravel()])
    a0 = np.linalg.solve(Gam, eta)
    y0 = np.concatenate([a0.real, a0.imag, np.zeros(4*K*K)])
    sol = solve_ivp(rhs, (0, T), y0, method='LSODA', rtol=1e-10, atol=1e-12)
    al, N, M = unpack(sol.y[:, -1])
    return al, N, M, np.max(np.abs(rhs(0, sol.y[:, -1])))

def g2_cross_gaussian(al, N, M, i, j):
    """<a_i^dag a_j^dag a_j a_i>/(n_i n_j) by Wick, for a Gaussian state with mean alpha."""
    ni = abs(al[i])**2 + N[i, i].real
    nj = abs(al[j])**2 + N[j, j].real
    aic, ajc = np.conj(al[i]), np.conj(al[j])
    # <a_i^dag a_j^dag a_j a_i> with a = alpha + da, Wick on the fluctuations
    t  = abs(al[i])**2*abs(al[j])**2
    t += abs(al[j])**2*N[i, i] + abs(al[i])**2*N[j, j]
    t += (aic*al[j]*N[j, i] + np.conj(aic*al[j]*N[j, i]))          # <da_j^dag da_i> cross terms
    t += aic*ajc*M[j, i] + np.conj(aic*ajc*M[j, i])
    t += N[i, i]*N[j, j] + abs(N[i, j])**2 + abs(M[i, j])**2
    return (t/(ni*nj)).real, ni, nj
