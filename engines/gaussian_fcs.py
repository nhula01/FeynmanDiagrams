"""Proper Gaussian (HFB) counting statistics for the driven Kerr cavity.

Step 1: solve the untilted HFB equations for (alpha, n, m).
Step 2: freeze them into a quadratic fluctuation Hamiltonian
        H_G = Delta a^dag a + i(eta a^dag - h.c.) + A da^dag da + (B/2)(da^dag^2 + h.c.),
        A = 2U(|alpha|^2 + n),  B = U(alpha^2 + m),  da = a - alpha,
        which is quadratic + linear in a, so the model is Gaussian.
Step 3: compute the FCS of that Gaussian model exactly (tilted leading eigenvalue).
This is the "strong scalable competitor": it is squeezed, hence NOT Poissonian.
"""
import numpy as np
from engine import destroy
from cumulant import kerr_cumulant_ss

def gaussian_model(kappa, Dc, eta, U, Nc=30):
    c = kerr_cumulant_ss(kappa, Dc, eta, U)
    al = c['alpha']; n = c['x'][2]; m = c['x'][3] + 1j*c['x'][4]
    A = 2*U*(abs(al)**2 + n); B = U*(al**2 + m)
    a = destroy(Nc); ad = a.conj().T; I = np.eye(Nc)
    da = a - al*I; dad = da.conj().T
    # effective drive that makes the Gaussian model's mean exactly the HFB mean:
    # stationarity of <a> requires eta_eff = (kappa/2 + i Delta) alpha
    eta_eff = (kappa/2 + 1j*Dc)*al
    H = Dc*ad@a + 1j*(eta_eff*ad - np.conj(eta_eff)*a) + A*(dad@da) + 0.5*(B*(dad@dad) + np.conj(B)*(da@da))
    return H, a, ad, dict(alpha=al, n=n, m=m, A=A, B=B)

def theta_lead(H, a, kappa, chi):
    ad = a.conj().T; D = a.shape[0]; I = np.eye(D)
    L = np.zeros((D*D, D*D), complex); E = np.zeros((D, D), complex)
    for i in range(D):
        for j in range(D):
            E[:] = 0; E[i, j] = 1
            col = -1j*(H@E - E@H) + kappa*(np.exp(chi)*a@E@ad - 0.5*(ad@a@E + E@ad@a))
            L[:, i + D*j] = col.reshape(-1, order='F')
    w = np.linalg.eigvals(L); return w[np.argmax(w.real)].real

def cumulant_ratios(theta_of_chi, h=0.05):
    xs = np.array([-2, -1, 0, 1, 2])*h
    th = np.array([theta_of_chi(x) for x in xs])
    c1 = (th[3]-th[1])/(2*h); c2 = (th[3]-2*th[2]+th[1])/h**2
    c3 = (th[4]-2*th[3]+2*th[1]-th[0])/(2*h**3)
    return c1, c2/c1, c3/c1
