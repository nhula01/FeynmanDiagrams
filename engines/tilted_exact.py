"""Exact leading eigenvalue of the counting-field-tilted Liouvillian (correct sign)."""
import numpy as np
def theta_lead(H, jumps, kappa_list, chi_list, sparse=False, complex_out=False):
    """L_chi(rho) = -i[H,rho] + sum_k ( e^{chi_k} c_k rho c_k^dag - 1/2{c_k^dag c_k, rho} ),
    with c_k = sqrt(kappa_k) * jumps[k]."""
    D = H.shape[0]
    cs = [np.sqrt(k)*j for k, j in zip(kappa_list, jumps)]
    L = np.zeros((D*D, D*D), complex); E = np.zeros((D, D), complex)
    for i in range(D):
        for j in range(D):
            E[:] = 0; E[i, j] = 1
            col = -1j*(H@E - E@H)
            for c, chi in zip(cs, chi_list):
                cd = c.conj().T
                col = col + np.exp(chi)*(c@E@cd) - 0.5*(cd@c@E + E@cd@c)
            L[:, i + D*j] = col.reshape(-1, order='F')
    w = np.linalg.eigvals(L)
    lead = w[np.argmax(w.real)]
    return lead if complex_out else lead.real
def cumulants(theta_of_chi, h=0.05):
    xs = np.array([-2, -1, 0, 1, 2])*h
    th = np.array([theta_of_chi(x) for x in xs])
    c1 = (th[3]-th[1])/(2*h); c2 = (th[3]-2*th[2]+th[1])/h**2
    c3 = (th[4]-2*th[3]+2*th[1]-th[0])/(2*h**3)
    return c1, c2/c1, c3/c1
