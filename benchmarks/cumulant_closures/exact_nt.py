"""Exact K=3 ring steady state in the total-photon-number truncated Fock space {n : sum_j n_j <= Nt}.
Every configuration with a single site holding up to Nt photons is kept, which is what the on-site g3 needs
(the per-site cube n_j < nmax of exact_ref.py wastes most of its D^2 on multiply occupied configurations
with negligible weight).  H_trunc = P H P exactly (a_j maps the space into itself; a_j^dag is projected).
Steady state by Krylov time evolution (expm_multiply) of the vectorised Liouvillian, warm-started from the
previous cutoff; the cutoff is raised until |g3(Nt)-g3(Nt-1)| < tol or Nt_cap.
Usage: python exact_nt.py U Nt_start Nt_cap [tol]   ->  exactNt_U{U:.3f}.json"""
import sys, os, json, time, itertools
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = _ROOT
P = json.load(open(os.path.join(ROOT, 'notebook', 'data', 'revision_params.json')))['ring3']
K, kappa, Delta, J, eta = P['K'], P['kappa'], P['Delta'], P['J'], P['eta']
log = lambda *a: (print(*a), sys.stdout.flush())

def basis(Nt):
    st = [s for s in itertools.product(range(Nt+1), repeat=K) if sum(s) <= Nt]
    return st, {s: i for i, s in enumerate(st)}

def ops(Nt):
    st, ix = basis(Nt); D = len(st); A = []
    for j in range(K):
        r, c, v = [], [], []
        for i, s in enumerate(st):
            if s[j] > 0:
                t = list(s); t[j] -= 1; r.append(ix[tuple(t)]); c.append(i); v.append(np.sqrt(s[j]))
        A.append(sp.csr_matrix((np.array(v, complex), (r, c)), shape=(D, D)))
    return st, ix, A

def build(U, Nt):
    st, ix, A = ops(Nt); D = len(st); Ad = [a.conj().T.tocsr() for a in A]
    H = sp.csr_matrix((D, D), dtype=complex)
    for j in range(K):
        H = H + Delta*(Ad[j]@A[j]) + U/2*(Ad[j]@Ad[j]@A[j]@A[j]) + 1j*(eta*Ad[j] - np.conj(eta)*A[j])
    for j in range(K):
        l = (j+1) % K
        H = H - J*(Ad[j]@A[l] + Ad[l]@A[j])
    I = sp.identity(D, format='csr', dtype=complex)
    L = -1j*(sp.kron(I, H) - sp.kron(H.T, I))
    for a in A:
        c = np.sqrt(kappa)*a; cdc = (c.conj().T@c).tocsr()
        L = L + sp.kron(c.conj(), c) - 0.5*sp.kron(I, cdc) - 0.5*sp.kron(cdc.T, I)
    return sp.csr_matrix(L), st, ix, A, Ad, D

def linear_alpha():
    Om = Delta*np.eye(K, dtype=complex)
    for j in range(K): Om[j, (j+1) % K] = Om[(j+1) % K, j] = -J
    return np.linalg.solve(kappa/2*np.eye(K) + 1j*Om, eta*np.ones(K, complex))

def exact_ss(U, Nt, prev=None, dt=10.0, Tmax=400.0, tol=1e-12):
    t0 = time.time()
    L, st, ix, A, Ad, D = build(U, Nt); t_build = time.time()-t0
    if prev is None:
        al = linear_alpha(); from math import factorial
        psi = np.array([np.prod([np.exp(-abs(al[j])**2/2)*al[j]**s[j]/np.sqrt(factorial(s[j])) for j in range(K)]) for s in st])
        psi /= np.linalg.norm(psi); rho = np.outer(psi, psi.conj())
    else:
        rho_p, st_p = prev; m = np.array([ix[s] for s in st_p])
        rho = np.zeros((D, D), complex); rho[np.ix_(m, m)] = rho_p
    v = rho.reshape(-1, order='F'); T = 0.0; res = np.inf
    while T < Tmax:
        v = expm_multiply(L*dt, v); T += dt
        v /= np.trace(v.reshape(D, D, order='F')).real
        res = float(np.max(np.abs(L@v))/np.max(np.abs(v)))
        if res < tol: break
    rho = v.reshape(D, D, order='F')
    tr = lambda O: complex(np.sum(O.multiply(rho.T)))
    n1 = Ad[0]@A[0]; n2 = Ad[0]@Ad[0]@A[0]@A[0]; n3 = Ad[0]@Ad[0]@Ad[0]@A[0]@A[0]@A[0]
    n1b = Ad[1]@A[1]; n3b = Ad[2]@Ad[2]@Ad[2]@A[2]@A[2]@A[2]
    tot = np.array([sum(s) for s in st]); pd = np.real(np.diag(rho))
    return dict(n=tr(n1).real, n2=tr(n2).real, n3=tr(n3).real, t=time.time()-t0, t_build=t_build, Nt=Nt, D=D, T=T, residual=res,
                herm=float(np.max(np.abs(rho-rho.conj().T))), pop_top_shell=float(pd[tot == Nt].sum()),
                site_asym=max(abs(tr(n1b).real-tr(n1).real), abs(tr(n3b).real-tr(n3).real)), warm=prev is not None), (rho, st)

if __name__ == "__main__":
    U = float(sys.argv[1]); n0, ncap = int(sys.argv[2]), int(sys.argv[3]); TOL = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-6
    k = f"{U:.3f}"; runs = {}; prev = None; Nt = n0
    while True:
        r, prev = exact_ss(U, Nt, prev); runs[Nt] = r
        ms = sorted(runs); g3s = [runs[m]['n3']/runs[m]['n']**3 for m in ms]; g2s = [runs[m]['n2']/runs[m]['n']**2 for m in ms]
        log(f"exactNt U={U} Nt={Nt} D={r['D']}: n={r['n']:.10f} g2={g2s[-1]:.10f} g3={g3s[-1]:.10f} ({r['t']:.1f}s, build {r['t_build']:.1f}s, "
            f"T={r['T']}, res={r['residual']:.1e}, top-shell pop {r['pop_top_shell']:.1e}, site asym {r['site_asym']:.1e})")
        ex = dict(n=r['n'], n2=r['n2'], n3=r['n3'], g2=g2s[-1], g3=g3s[-1], truncation='total photon number sum_j n_j <= Nt', Nt=Nt, D=r['D'],
                  t_solve=r['t'], T=r['T'], residual=r['residual'], pop_top_shell=r['pop_top_shell'], herm=r['herm'], site_asym=r['site_asym'],
                  g3_by_Nt={str(m): g3s[i] for i, m in enumerate(ms)}, g2_by_Nt={str(m): g2s[i] for i, m in enumerate(ms)},
                  n_by_Nt={str(m): runs[m]['n'] for m in ms}, t_by_Nt={str(m): runs[m]['t'] for m in ms},
                  pop_top_by_Nt={str(m): runs[m]['pop_top_shell'] for m in ms})
        if len(ms) >= 2:
            ex['g3_change_last'] = abs(g3s[-1]-g3s[-2]); ex['g2_change_last'] = abs(g2s[-1]-g2s[-2]); ex['n_change_last'] = abs(runs[ms[-1]]['n']-runs[ms[-2]]['n'])
        if len(ms) >= 3:
            for q, gs in (('g3', g3s), ('g2', g2s)):
                d1, d2 = gs[-2]-gs[-3], gs[-1]-gs[-2]; rr = d2/d1 if d1 != 0 else 0.0
                ex[q+'_ratio'] = rr
                ex[q+'_cutoff_err_est'] = abs(d2*rr/(1-rr)) if abs(rr) < 0.9 else abs(d2)
        ex['converged_in_cutoff'] = bool(len(ms) >= 2 and ex['g3_change_last'] < TOL)
        json.dump({k: ex}, open(os.path.join(HERE, f'exactNt_U{k}.json'), 'w'), indent=1, default=float)
        if ex['converged_in_cutoff'] or Nt >= ncap: break
        Nt += 1
    log(f"wrote exactNt_U{k}.json (final Nt {Nt}, converged={ex['converged_in_cutoff']})")
