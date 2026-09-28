"""Exact long-time counting cumulants by Rayleigh-Schroedinger expansion of the leading eigenvalue of
L_chi = L_0 + (e^chi - 1) J (one sparse LU of the bordered L_0, no step error, no eigensolver).
theta_k = sum_{j=1..k} Tr[L_j rho_{k-j}],  L_j = J/j!,  rho_k from L_0 rho_k = -(sum_j L_j rho_{k-j} - sum_j theta_j rho_{k-j}), Tr rho_k = 0.
c_n = n! theta_n.  usage: exact_rs.py K1|K3 U cut [cut...]   (K3 cut c -> [c, c-3, c-3])"""
import sys, json, math, time, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
sys.path.insert(0, '..'); sys.path.insert(0, '.')
from exact_fcs import ops, liouv

def rs_cumulants(L, Jc, D, nmax=5):
    tr = np.eye(D).reshape(-1, order='F')
    A = sp.vstack([sp.csr_matrix(tr.reshape(1, -1)), L[1:]]).tocsc()
    lu = spl.splu(A)
    def solve(b):                    # L0 x = b on the traceless subspace, Tr x = 0
        bb = b.copy(); bb[0] = 0.0; return lu.solve(bb)
    b = np.zeros(D*D, complex); b[0] = 1.0; rho0 = lu.solve(b)
    rhos = [rho0]; th = [0.0]
    for k in range(1, nmax+1):
        s = sum((Jc @ rhos[k-j])/math.factorial(j) for j in range(1, k+1))
        thk = tr @ s; th.append(thk)
        rhs = s - sum(th[j]*rhos[k-j] for j in range(1, k+1))
        rhos.append(solve(-rhs))
    c = [math.factorial(n)*th[n].real for n in range(1, nmax+1)]
    n0 = None
    return c, rho0

if __name__ == "__main__":
    which, U = sys.argv[1], float(sys.argv[2]); out = []
    for cut in [int(x) for x in sys.argv[3:]]:
        t0 = time.time()
        if which == 'K1': cutl, det, J, eta, ring = [cut], [-1.0], 0.0, [1.0], False
        else: cutl, det, J, eta, ring = [cut, cut-3, cut-3], [-1.0]*3, 0.4, [1.0, 0, 0], True
        H, a = ops(cutl, det, J, eta, U, ring); L, Jc, D = liouv(H, a, 0)
        c, rho = rs_cumulants(L, Jc, D)
        R = rho.reshape(D, D, order='F'); R = R/np.trace(R)
        n = [float(np.real(np.trace(R @ (x.conj().T @ x).toarray()))) for x in a]
        rec = dict(model=which, U=U, cut=cutl, D=D, c=c, ratios=[c[0], c[1]/c[0], c[2]/c[0]], n=n, wall=time.time()-t0)
        out.append(rec); print(json.dumps(rec), flush=True)
    json.dump(out, open(f'exact_rs_{which}_U{U:.2f}.json', 'w'), indent=1)

