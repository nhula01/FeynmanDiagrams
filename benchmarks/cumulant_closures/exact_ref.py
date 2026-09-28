"""Driver: exact reference (with nmax convergence), cumulant closures N=2,3,4, diagrams N=2,4,
validation checks, timings.  Writes cumulant_closure_results.json."""
import sys, os, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(HERE, '..', '..', 'engines')
sys.path.insert(0, HERE); sys.path.insert(0, CODE)
import scipy.sparse as sp, scipy.sparse.linalg as spl
from cumulant_closure import ring_closure, linear_coherent_alpha, unit
import json
from cumulant import ring_cumulant_ss
from hfb_chain import hfb_chain
from ring_exact import ring_H
import multimode as MM

P = json.load(open(os.path.join(HERE, '..', '..', 'notebook', 'data', 'revision_params.json')))['ring3']
K, kappa, Delta, J, eta = P['K'], P['kappa'], P['Delta'], P['J'], P['eta']
PP = dict(kappa=kappa, Delta=Delta, J=J, eta=eta)
U_fig11 = P['U_values_fig11']
U_grid = [round(0.025*i, 3) for i in range(13)]          # 0 .. 0.30, matches ring_exactU.json keys
U_check = [0.0, 0.02]
U_all = sorted(set(U_grid) | set(U_fig11) | set(U_check))
out = dict(params=P, U_grid=U_grid, U_fig11=U_fig11)
log = lambda *a: (print(*a), sys.stdout.flush())

# ------------------------------------------------------------- exact reference
def build_L(U, nmax):
    H, A, Ad = ring_H(K, Delta, J, eta, U, nmax)
    D = nmax**K; I = sp.identity(D, format='csr')
    L = -1j*(sp.kron(I, H) - sp.kron(H.T, I))
    for cop in A:
        cop = (np.sqrt(kappa)*cop).tocsr(); cd = cop.conj().T.tocsr(); cdc = (cd@cop).tocsr()
        L = L + sp.kron(cop.conj(), cop) - 0.5*sp.kron(I, cdc) - 0.5*sp.kron(cdc.T, I)
    return sp.csr_matrix(L), A, Ad, D

def coherent_vec(alpha, nmax):
    n = np.arange(nmax); from math import factorial
    v = np.exp(-abs(alpha)**2/2)*alpha**n/np.sqrt([factorial(int(k)) for k in n])
    return v/np.linalg.norm(v)

def exact_ss(U, nmax, dt=10.0, Tmax=400.0, tol=1e-12):
    """Steady state by Krylov time evolution of the vectorised Liouvillian (sparse LU fills in to dense
    for this problem and is impractical beyond nmax=5).  Start from the product coherent state of the linear
    problem; evolve in steps of dt until the residual max|L rho| / max|rho| < tol."""
    from scipy.sparse.linalg import expm_multiply
    t0 = time.time()
    L, A, Ad, D = build_L(U, nmax)
    al = linear_coherent_alpha(K, kappa, Delta, J, eta)
    psi = np.array([1.0+0j])
    for j in range(K): psi = np.kron(psi, coherent_vec(al[j], nmax))
    rho = np.outer(psi, psi.conj()); v = rho.reshape(-1, order='F')
    T = 0.0; res = np.inf; hist = []
    while T < Tmax:
        v = expm_multiply(L*dt, v); T += dt
        v /= np.trace(v.reshape(D, D, order='F')).real
        res = float(np.max(np.abs(L@v))/np.max(np.abs(v))); hist.append((T, res))
        if res < tol: break
    rho = v.reshape(D, D, order='F')
    tr = lambda O: complex(np.sum(O.multiply(rho.T)))
    n1 = Ad[0]@A[0]; n2 = Ad[0]@Ad[0]@A[0]@A[0]; n3 = Ad[0]@Ad[0]@Ad[0]@A[0]@A[0]@A[0]
    pop_top = float(np.real(sum(rho[i, i] for i in range(D) if max(np.unravel_index(i, (nmax,)*K)) == nmax-1)))
    return dict(n=tr(n1).real, n2=tr(n2).real, n3=tr(n3).real, t=time.time()-t0, nmax=nmax, T=T, residual=res,
                herm=float(np.max(np.abs(rho-rho.conj().T))), pop_top_level=pop_top, hist=hist)

NMAXES = (8, 9, 10)
part, nparts = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (0, 1)
U_all = U_all[part::nparts]
exact = {}
old = json.load(open(os.path.join(HERE, '..', '..', 'notebook', 'data', 'ring_exactU.json')))
for U in U_all:
    runs = {}
    for nmax in NMAXES:
        runs[nmax] = exact_ss(U, nmax)
        log(f"exact U={U} nmax={nmax}: n={runs[nmax]['n']:.10f} g2={runs[nmax]['n2']/runs[nmax]['n']**2:.10f} g3={runs[nmax]['n3']/runs[nmax]['n']**3:.10f} ({runs[nmax]['t']:.1f}s, T={runs[nmax]['T']}, res={runs[nmax]['residual']:.1e}, top-level pop {runs[nmax]['pop_top_level']:.1e})")
    r = runs[NMAXES[-1]]
    g3s = [runs[m]['n3']/runs[m]['n']**3 for m in NMAXES]
    ex = dict(n=r['n'], n2=r['n2'], n3=r['n3'], g2=r['n2']/r['n']**2, g3=r['n3']/r['n']**3,
              nmax=NMAXES[-1], t_solve=r['t'], T=r['T'], residual=r['residual'], pop_top_level=r['pop_top_level'], herm=r['herm'], g3_change_last=abs(g3s[-1]-g3s[-2]), g3_change_prev=abs(g3s[-2]-g3s[-3]), g3_by_nmax={str(m): runs[m]['n3']/runs[m]['n']**3 for m in NMAXES}, n_by_nmax={str(m): runs[m]['n'] for m in NMAXES}, t_by_nmax={str(m): runs[m]['t'] for m in NMAXES})
    k = f"{U:.3f}"
    if k in old:
        ex['n_vs_ring_exactU'] = abs(r['n'] - old[k]['n0']); ex['n2_vs_ring_exactU'] = abs(r['n2'] - old[k]['n0n0'])
    exact[k] = ex
    json.dump(exact, open(os.path.join(HERE, f'exact_ref_part{part}of{nparts}.partial.json'), 'w'), indent=1, default=float)
json.dump(exact, open(os.path.join(HERE, f'exact_ref_part{part}of{nparts}.json'), 'w'), indent=1, default=float)
log(f'wrote exact_ref_part{part}of{nparts}.json')
parts = [os.path.join(HERE, f'exact_ref_part{i}of{nparts}.json') for i in range(nparts)]
if all(os.path.exists(p) for p in parts):
    time.sleep(3); merged = {}
    for p in parts: merged.update(json.load(open(p)))
    merged = dict(sorted(merged.items()))
    tmp = os.path.join(HERE, 'exact_ref.json.tmp'); json.dump(merged, open(tmp, 'w'), indent=1, default=float); os.replace(tmp, os.path.join(HERE, 'exact_ref.json'))
    log('wrote exact_ref.json')

