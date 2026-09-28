"""Exact K=3 ring steady state for one U by Krylov time evolution of the vectorised Liouvillian,
with a Fock-cutoff convergence sequence.  Usage: python exact_ref2.py U nmax_start nmax_cap [tol]
The cutoff is raised from nmax_start until |g3(nmax)-g3(nmax-1)| < tol (default 1e-6) or nmax_cap is reached;
each run is warm-started from the steady state at the previous cutoff.
Writes exact_U{U:.3f}.json.  (Same method as exact_ref.py, which is kept unchanged; this file
exists so that each U can run as its own single-core process after the node reboot.)"""
import sys, os, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(HERE, '..', '..', 'engines')
sys.path.insert(0, HERE); sys.path.insert(0, CODE)
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply
from cumulant_closure import linear_coherent_alpha
from ring_exact import ring_H
from math import factorial

P = json.load(open(os.path.join(HERE, '..', '..', 'notebook', 'data', 'revision_params.json')))['ring3']
K, kappa, Delta, J, eta = P['K'], P['kappa'], P['Delta'], P['J'], P['eta']
log = lambda *a: (print(*a), sys.stdout.flush())

def build_L(U, nmax):
    H, A, Ad = ring_H(K, Delta, J, eta, U, nmax)
    D = nmax**K; I = sp.identity(D, format='csr')
    L = -1j*(sp.kron(I, H) - sp.kron(H.T, I))
    for cop in A:
        cop = (np.sqrt(kappa)*cop).tocsr(); cd = cop.conj().T.tocsr(); cdc = (cd@cop).tocsr()
        L = L + sp.kron(cop.conj(), cop) - 0.5*sp.kron(I, cdc) - 0.5*sp.kron(cdc.T, I)
    return sp.csr_matrix(L), A, Ad, D

def coherent_vec(alpha, nmax):
    n = np.arange(nmax)
    v = np.exp(-abs(alpha)**2/2)*alpha**n/np.sqrt([factorial(int(k)) for k in n])
    return v/np.linalg.norm(v)

def exact_ss(U, nmax, rho0=None, dt=10.0, Tmax=400.0, tol=1e-12):
    t0 = time.time()
    L, A, Ad, D = build_L(U, nmax)
    t_build = time.time() - t0
    if rho0 is None:
        al = linear_coherent_alpha(K, kappa, Delta, J, eta)
        psi = np.array([1.0+0j])
        for j in range(K): psi = np.kron(psi, coherent_vec(al[j], nmax))
        rho_init = np.outer(psi, psi.conj())
    else:   # embed the lower-cutoff steady state (tensor index order of ring_H: site 0 slowest)
        m = rho0.shape[0]; mo = round(m**(1/K))
        R = np.zeros((nmax,)*(2*K), complex)
        R[(slice(0, mo),)*(2*K)] = rho0.reshape((mo,)*(2*K))
        rho_init = R.reshape(D, D)
    v = rho_init.reshape(-1, order='F')
    T = 0.0; res = np.inf
    while T < Tmax:
        v = expm_multiply(L*dt, v); T += dt
        v /= np.trace(v.reshape(D, D, order='F')).real
        res = float(np.max(np.abs(L@v))/np.max(np.abs(v)))
        if res < tol: break
    rho = v.reshape(D, D, order='F')
    tr = lambda O: complex(np.sum(O.multiply(rho.T)))
    n1 = Ad[0]@A[0]; n2 = Ad[0]@Ad[0]@A[0]@A[0]; n3 = Ad[0]@Ad[0]@Ad[0]@A[0]@A[0]@A[0]
    n1b = Ad[1]@A[1]; n3b = Ad[2]@Ad[2]@Ad[2]@A[2]@A[2]@A[2]
    idx = np.array([max(np.unravel_index(i, (nmax,)*K)) for i in range(D)])
    pop_top = float(np.real(np.sum(np.diag(rho)[idx == nmax-1])))
    out = dict(n=tr(n1).real, n2=tr(n2).real, n3=tr(n3).real, t=time.time()-t0, t_build=t_build, nmax=nmax, T=T, residual=res,
               herm=float(np.max(np.abs(rho-rho.conj().T))), pop_top_level=pop_top, warm_start=rho0 is not None,
               site_asym=max(abs(tr(n1b).real-tr(n1).real), abs(tr(n3b).real-tr(n3).real)))
    return out, rho

U = float(sys.argv[1]); n_start, n_cap = int(sys.argv[2]), int(sys.argv[3])
TOL = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-6
k = f"{U:.3f}"
old = json.load(open(os.path.join(HERE, '..', '..', 'notebook', 'data', 'ring_exactU.json')))
runs = {}; rho = None; nmax = n_start
while True:
    r, rho = exact_ss(U, nmax, rho0=rho)
    runs[nmax] = r
    ms = sorted(runs); g3s = [runs[m]['n3']/runs[m]['n']**3 for m in ms]; g2s = [runs[m]['n2']/runs[m]['n']**2 for m in ms]
    log(f"exact U={U} nmax={nmax}: n={r['n']:.10f} g2={g2s[-1]:.10f} g3={g3s[-1]:.10f} ({r['t']:.1f}s, build {r['t_build']:.1f}s, T={r['T']}, "
        f"res={r['residual']:.1e}, top-level pop {r['pop_top_level']:.1e}, site asym {r['site_asym']:.1e})")
    ex = dict(n=r['n'], n2=r['n2'], n3=r['n3'], g2=g2s[-1], g3=g3s[-1], nmax=nmax, t_solve=r['t'], T=r['T'], residual=r['residual'],
              pop_top_level=r['pop_top_level'], herm=r['herm'], site_asym=r['site_asym'],
              g3_by_nmax={str(m): g3s[i] for i, m in enumerate(ms)}, g2_by_nmax={str(m): g2s[i] for i, m in enumerate(ms)},
              n_by_nmax={str(m): runs[m]['n'] for m in ms}, t_by_nmax={str(m): runs[m]['t'] for m in ms},
              T_by_nmax={str(m): runs[m]['T'] for m in ms}, pop_top_by_nmax={str(m): runs[m]['pop_top_level'] for m in ms})
    if len(ms) >= 2:
        ex['g3_change_last'] = abs(g3s[-1]-g3s[-2]); ex['g2_change_last'] = abs(g2s[-1]-g2s[-2]); ex['n_change_last'] = abs(runs[ms[-1]]['n']-runs[ms[-2]]['n'])
    if len(ms) >= 3:
        # geometric (Aitken) estimate of the remaining cutoff error
        for q, gs in (('g3', g3s), ('g2', g2s)):
            d1, d2 = gs[-2]-gs[-3], gs[-1]-gs[-2]
            rr = d2/d1 if d1 != 0 else 0.0
            ex[q+'_ratio'] = rr
            ex[q+'_extrap'] = gs[-1] + (d2*rr/(1-rr) if abs(rr) < 0.9 else np.nan)
            ex[q+'_cutoff_err_est'] = abs(d2*rr/(1-rr)) if abs(rr) < 0.9 else abs(d2)
    if k in old:
        ex['n_vs_ring_exactU'] = abs(r['n'] - old[k]['n0']); ex['n2_vs_ring_exactU'] = abs(r['n2'] - (old[k]['n0n0'] - old[k]['n0']))
    ex['converged_in_cutoff'] = bool(len(ms) >= 2 and ex['g3_change_last'] < TOL)
    json.dump({k: ex}, open(os.path.join(HERE, f'exact_U{k}.json'), 'w'), indent=1, default=float)
    if ex['converged_in_cutoff'] or nmax >= n_cap: break
    nmax += 1
log(f'wrote exact_U{k}.json (final nmax {nmax}, converged={ex["converged_in_cutoff"]})')
