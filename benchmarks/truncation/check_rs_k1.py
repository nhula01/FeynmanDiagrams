
import sys, os, numpy as np
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
sys.path.insert(0, os.path.join(_ROOT, 'engines'))
from fcs import ladder, Y00
from engine import add, scale, mul, cav, comm
from tilted_exact import theta_lead
from engine import destroy

def theta_series_rs(kappa, Dc, eta, U, chi, Nmax, maxdeg=14, full_rs=True):
    """Rayleigh-Schrodinger recursion for the leading eigenvalue of the tilted generator.
    full_rs=True includes the terms -sum_{k=1}^{N-1} theta_k psi_{N-k}; False reproduces fcs.theta_series."""
    z = kappa/2 + 1j*Dc; alpha = eta/z; eps = kappa*(np.exp(chi)-1)
    I1 = np.eye(1)
    b, bd, one = cav(0, 1, I1), cav(1, 0, I1), cav(0, 0, I1)
    a_ = add(b, scale(alpha, one)); ad_ = add(bd, scale(np.conj(alpha), one))
    HK = scale(U/2, mul(mul(ad_, ad_), mul(a_, a_)))
    V = lambda X: scale(1j, comm(HK, X))
    Y = Y00(alpha, eps, z, maxdeg)
    th = [eps*abs(alpha)**2]; psi = [Y]
    for N in range(1, Nmax+1):
        W = V(psi[N-1])
        c00 = complex(W.get((0, 0), np.zeros((1, 1)))[0, 0]); th.append(c00)
        W = add(W, scale(-c00, Y))
        if full_rs:
            for k in range(1, N):
                W = add(W, scale(-th[k], psi[N-k]))
        X = ladder(W, alpha, eps, z, maxdeg)
        X = {k: v for k, v in X.items() if sum(k) <= maxdeg and np.max(np.abs(v)) > 1e-16}
        psi.append(X)
    return np.array(th)

kappa, Dc, eta = 1.0, -1.0, 1.0
A = destroy(30); Ad = A.conj().T
for U in (0.05, 0.1):
    Hx = Dc*Ad@A + 1j*(eta*Ad - np.conj(eta)*A) + U/2*Ad@Ad@A@A
    for chi in (0.05, 0.1):
        te = theta_lead(Hx, [A], [kappa], [chi])
        old = np.cumsum(theta_series_rs(kappa, Dc, eta, U, chi, 8, full_rs=False)).real
        new = np.cumsum(theta_series_rs(kappa, Dc, eta, U, chi, 8, full_rs=True)).real
        print(f"U={U} chi={chi}: exact {te:.10f}")
        print("   old (paper) diffs by order:", np.array2string(old-te, formatter={'float_kind':lambda x: f'{x:+.1e}'}))
        print("   full RS     diffs by order:", np.array2string(new-te, formatter={'float_kind':lambda x: f'{x:+.1e}'}))
def cums(f, h=0.05):
    th = np.array([f(x) for x in [-2*h,-h,0,h,2*h]])
    c1=(th[3]-th[1])/(2*h); c2=(th[3]-2*th[2]+th[1])/h**2; c3=(th[4]-2*th[3]+2*th[1]-th[0])/(2*h**3)
    return c1, c2/c1, c3/c1
for U in (0.05, 0.1):
    Hx = Dc*Ad@A + 1j*(eta*Ad - np.conj(eta)*A) + U/2*Ad@Ad@A@A
    ce = cums(lambda x: theta_lead(Hx,[A],[kappa],[x]))
    cn = cums(lambda x: np.cumsum(theta_series_rs(kappa,Dc,eta,U,x,8)).real[-1])
    co = cums(lambda x: np.cumsum(theta_series_rs(kappa,Dc,eta,U,x,8,full_rs=False)).real[-1])
    print(f"U={U}: exact F {ce[1]:.5f} c3/c1 {ce[2]:.5f} | full RS N8 F {cn[1]:.5f} c3/c1 {cn[2]:.5f} | paper N8 F {co[1]:.5f} c3/c1 {co[2]:.5f}")
