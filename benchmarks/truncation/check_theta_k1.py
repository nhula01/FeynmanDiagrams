
import sys, os, numpy as np
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
sys.path.insert(0, os.path.join(_ROOT, 'engines'))
from fcs import theta_series
from tilted_exact import theta_lead
from engine import destroy, steady
kappa, Dc, eta, U = 1.0, -1.0, 1.0, 0.05
A = destroy(30); Ad = A.conj().T
Hx = Dc*Ad@A + 1j*(eta*Ad - np.conj(eta)*A) + U/2*Ad@Ad@A@A
nss = np.trace(Ad@A@steady(Hx,[A])).real
print("kappa n_ss =", nss)
for chi in [0.005, 0.01, 0.02, 0.05, 0.1, 0.2]:
    te = theta_lead(Hx, [A], [kappa], [chi])
    td = np.cumsum(theta_series(kappa, Dc, eta, U, chi, 8, maxdeg=14)).real
    print(f"chi={chi}: exact {te:.10f}  diag N8 {td[-1]:.10f}  diff {td[-1]-te:+.2e}   diag N7 diff {td[-2]-te:+.2e}")
def cums(f, h):
    th = np.array([f(x) for x in [-2*h,-h,0,h,2*h]])
    c1=(th[3]-th[1])/(2*h); c2=(th[3]-2*th[2]+th[1])/h**2; c3=(th[4]-2*th[3]+2*th[1]-th[0])/(2*h**3)
    return c1, c2/c1, c3/c1
for h in [0.05, 0.02, 0.01]:
    ce = cums(lambda x: theta_lead(Hx,[A],[kappa],[x]), h)
    cd = cums(lambda x: np.cumsum(theta_series(kappa,Dc,eta,U,x,8,maxdeg=14)).real[-1], h)
    print(f"h={h}: exact c1 {ce[0]:.6f} F {ce[1]:.6f} c3/c1 {ce[2]:.6f} | diag c1 {cd[0]:.6f} F {cd[1]:.6f} c3/c1 {cd[2]:.6f} | dF {cd[1]-ce[1]:+.1e} dc3 {cd[2]-ce[2]:+.1e}")
