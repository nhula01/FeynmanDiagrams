"""Rayleigh-Schroedinger counting cumulants (as exact_rs.py) with the linear solves replaced by relaxation:
for traceless r, x(t) from dx/dt = L_0 x + r, x(0)=0, converges to x = -L_0^{-1} r at rate |Re lambda_1|;
the steady state comes from d rho/dt = L_0 rho.  Fixed-step RK4 on the sparse Liouvillian; no factorization,
so memory stays at a few vectors.  usage: exact_rs_ode.py K3 U cut [Tmax dt]"""
import sys, json, math, time, numpy as np
sys.path.insert(0, '..'); sys.path.insert(0, '.')
from exact_fcs import ops, liouv

def relax(L, x0, r, T, dt):
    f = (lambda x: L @ x + r) if r is not None else (lambda x: L @ x)
    x = x0.copy()
    for _ in range(int(round(T/dt))):
        k1 = f(x); k2 = f(x + 0.5*dt*k1); k3 = f(x + 0.5*dt*k2); k4 = f(x + dt*k3)
        x = x + dt/6*(k1 + 2*k2 + 2*k3 + k4)
    return x

if __name__ == "__main__":
    which, U, cut = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
    Tm = float(sys.argv[4]) if len(sys.argv) > 4 else 80.0; dt = float(sys.argv[5]) if len(sys.argv) > 5 else 0.02
    t0 = time.time()
    cutl, det, J, eta, ring = ([cut], [-1.0], 0.0, [1.0], False) if which == 'K1' else ([cut, cut-3, cut-3], [-1.0]*3, 0.4, [1.0, 0, 0], True)
    H, a = ops(cutl, det, J, eta, U, ring); L, Jc, D = liouv(H, a, 0); L = L.tocsr(); Jc = Jc.tocsr()
    tr = np.eye(D).reshape(-1, order='F')
    rho = np.zeros(D*D, complex); rho[0] = 1.0
    rho = relax(L, rho, None, Tm, dt); rho = rho/(tr @ rho)
    res0 = float(np.linalg.norm(L @ rho))
    rhos = [rho]; th = [0.0]; resid = [res0]
    for k in range(1, 4):
        s = sum((Jc @ rhos[k-j])/math.factorial(j) for j in range(1, k+1))
        thk = tr @ s; th.append(thk)
        r = s - sum(th[j]*rhos[k-j] for j in range(1, k+1))
        x = relax(L, np.zeros_like(rho), r, Tm, dt)     # L x = -r
        x = x - (tr @ x)*rho                              # remove any trace component
        resid.append(float(np.linalg.norm(L @ x + r)/np.linalg.norm(r))); rhos.append(x)
    c = [math.factorial(n)*th[n].real for n in range(1, 4)]
    R = rho.reshape(D, D, order='F')
    n = [float(np.real(np.trace(R @ (y.conj().T @ y).toarray()))) for y in a]
    rec = dict(model=which, U=U, cut=cutl, D=D, c=c, ratios=[c[0], c[1]/c[0], c[2]/c[0]], n=n, residuals=resid, Tmax=Tm, dt=dt, wall=time.time()-t0, method='rayleigh-schroedinger (relaxation)')
    print(json.dumps(rec), flush=True)

