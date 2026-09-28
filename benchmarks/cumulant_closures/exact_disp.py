"""Exact K=3 ring steady state in a displaced, total-photon-number-truncated Fock basis.

a_j = alpha_j + b_j with alpha the mean-field (classical) steady state, and the fluctuation Fock space truncated
at sum_j n_j(b) <= Nt.  The displacement is a unitary change of frame, so the Liouvillian is unchanged; only the
truncation is different, and because the coherent part of the state no longer has to be represented the
truncation error falls much faster with Nt.  With P the projector on the truncated space, every normal-ordered
product satisfies P a^dag^k a^l P = (P a^dag P)^k (P a P)^l (a = alpha + b lowers or keeps the fluctuation number),
so H, the jump terms and all observables are the exact projections.

Steady state: fixed-step RK4 on d rho/dt = -i(Heff rho - rho Heff^dag) + sum_j c_j rho c_j^dag with dense rho and
sparse operators (matrix-free).  The RK4 map is a polynomial p(hL) with p(0)=1, so its fixed point is exactly the
null vector of L whatever the step; h is set from a power-iteration estimate of the spectral radius.  Iterate until
max|L rho| / max|rho| < tol.  Warm start from the previous cutoff.
Usage: python exact_disp.py U Nt_start Nt_cap [tol_g3] [time_budget_s]  ->  exactD_U{U:.3f}.json"""
import sys, os, json, time, itertools
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
import numpy as np
import scipy.sparse as sp
from scipy.optimize import fsolve
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = _ROOT
P = json.load(open(os.path.join(ROOT, 'notebook', 'data', 'revision_params.json')))['ring3']
K, kappa, Delta, J, eta = P['K'], P['kappa'], P['Delta'], P['J'], P['eta']
log = lambda *a: (print(*a), sys.stdout.flush())

def mean_field(U):
    def f(y):
        a = y[:K] + 1j*y[K:]
        d = -(kappa/2 + 1j*Delta)*a - 1j*U*np.abs(a)**2*a + 1j*J*(np.roll(a, 1) + np.roll(a, -1)) + eta
        return np.concatenate([d.real, d.imag])
    Om = Delta*np.eye(K, dtype=complex)
    for j in range(K): Om[j, (j+1) % K] = Om[(j+1) % K, j] = -J
    a0 = np.linalg.solve(kappa/2*np.eye(K) + 1j*Om, eta*np.ones(K, complex))
    y = fsolve(f, np.concatenate([a0.real, a0.imag]), xtol=1e-14)
    assert np.max(np.abs(f(y))) < 1e-10
    return y[:K] + 1j*y[K:]

def basis(Nt):
    st = [s for s in itertools.product(range(Nt+1), repeat=K) if sum(s) <= Nt]
    return st, {s: i for i, s in enumerate(st)}

def operators(U, Nt, alpha):
    st, ix = basis(Nt); D = len(st); B = []
    for j in range(K):
        r, c, v = [], [], []
        for i, s in enumerate(st):
            if s[j] > 0:
                t = list(s); t[j] -= 1; r.append(ix[tuple(t)]); c.append(i); v.append(np.sqrt(s[j]))
        B.append(sp.csr_matrix((np.array(v, complex), (r, c)), shape=(D, D)))
    I = sp.identity(D, dtype=complex, format='csr')
    A = [(alpha[j]*I + B[j]).tocsr() for j in range(K)]
    Ad = [a.conj().T.tocsr() for a in A]
    H = sp.csr_matrix((D, D), dtype=complex)
    for j in range(K):
        H = H + Delta*(Ad[j]@A[j]) + U/2*(Ad[j]@Ad[j]@A[j]@A[j]) + 1j*(eta*Ad[j] - np.conj(eta)*A[j])
    for j in range(K):
        l = (j+1) % K; H = H - J*(Ad[j]@A[l] + Ad[l]@A[j])
    Heff = (H - 0.5j*kappa*sum(Ad[j]@A[j] for j in range(K))).tocsr()
    C = [(np.sqrt(kappa)*a).tocsr() for a in A]
    return st, ix, D, A, Ad, Heff, C

def make_L(Heff, C):
    def L(r):                       # valid for Hermitian r (all RK4 stages are Hermitian)
        X = -1j*(Heff @ r)
        out = X + X.conj().T
        for c in C:
            Y = c @ r
            out += c @ Y.conj().T   # c r c^dag = c (c r)^dag for Hermitian r
        return out
    return L

def spectral_radius(L, D, it=40, seed=1):
    rng = np.random.default_rng(seed); r = rng.standard_normal((D, D)) + 1j*rng.standard_normal((D, D)); r = r + r.conj().T
    lam = 0.0
    for _ in range(it):
        w = L(r); nw = np.linalg.norm(w); lam = max(lam, nw/np.linalg.norm(r)); r = w/nw
    return lam

def steady(U, Nt, alpha, prev=None, tol=1e-11, maxT=2000.0):
    t0 = time.time()
    st, ix, D, A, Ad, Heff, C = operators(U, Nt, alpha); L = make_L(Heff, C)
    lam = spectral_radius(L, D); h = 2.2/(1.15*lam)
    if prev is None:
        r = np.zeros((D, D), complex); r[0, 0] = 1.0
    else:
        rp, stp = prev; m = np.array([ix[s] for s in stp]); r = np.zeros((D, D), complex); r[np.ix_(m, m)] = rp
    T = 0.0; nstep = 0; res = np.inf; res_hist = []
    while T < maxT:
        for _ in range(50):
            k1 = L(r); k2 = L(r + 0.5*h*k1); k3 = L(r + 0.5*h*k2); k4 = L(r + h*k3)
            r = r + (h/6)*(k1 + 2*k2 + 2*k3 + k4); nstep += 1
        T += 50*h
        r = 0.5*(r + r.conj().T); r /= np.trace(r).real
        res = float(np.max(np.abs(L(r)))/np.max(np.abs(r))); res_hist.append(res)
        if not np.isfinite(res) or res > 1e3:
            raise RuntimeError(f"RK4 unstable at h={h}")
        if res < tol: break
    tr = lambda O: complex(np.sum(O.multiply(r.T)))
    obs = {}
    for j in range(K):
        n1 = Ad[j]@A[j]; n2 = Ad[j]@n1@A[j]; n3 = Ad[j]@n2@A[j]
        obs[j] = (tr(n1).real, tr(n2).real, tr(n3).real)
    tot = np.array([sum(s) for s in st]); pd = np.real(np.diag(r))
    n, n2, n3 = obs[0]
    out = dict(n=n, n2=n2, n3=n3, g2=n2/n**2, g3=n3/n**3, Nt=Nt, D=D, t=time.time()-t0, T=T, steps=nstep, h=h, spectral_radius=lam,
               residual=res, pop_top_shell=float(pd[tot == Nt].sum()), mean_b_number=float(np.sum(pd*tot)),
               site_asym=max(abs(obs[j][q]-obs[0][q]) for j in range(K) for q in range(3)), warm=prev is not None,
               a0=complex(tr(A[0])))
    return out, (r, st)

if __name__ == "__main__":
    U = float(sys.argv[1]); n0, ncap = int(sys.argv[2]), int(sys.argv[3]); TOL = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-8
    TBUDGET = float(sys.argv[5]) if len(sys.argv) > 5 else 1e9   # stop raising Nt once one solve exceeds this (s)
    k = f"{U:.3f}"; alpha = mean_field(U); runs = {}; prev = None; Nt = n0
    while True:
        r, prev = steady(U, Nt, alpha, prev); runs[Nt] = r
        ms = sorted(runs); g3s = [runs[m]['g3'] for m in ms]; g2s = [runs[m]['g2'] for m in ms]
        log(f"exactD U={U} Nt={Nt} D={r['D']}: n={r['n']:.12f} g2={r['g2']:.12f} g3={r['g3']:.12f} ({r['t']:.1f}s, T={r['T']:.1f}, "
            f"res={r['residual']:.1e}, top-shell {r['pop_top_shell']:.1e}, <N_b>={r['mean_b_number']:.3e}, asym {r['site_asym']:.1e})")
        ex = dict(n=r['n'], n2=r['n2'], n3=r['n3'], g2=r['g2'], g3=r['g3'], truncation='displaced frame a=alpha_mf+b, sum_j n_j(b) <= Nt',
                  alpha_mf=[complex(x) for x in alpha], a0=r['a0'], Nt=Nt, D=r['D'], t_solve=r['t'], T=r['T'], steps=r['steps'], residual=r['residual'],
                  pop_top_shell=r['pop_top_shell'], mean_b_number=r['mean_b_number'], site_asym=r['site_asym'],
                  g3_by_Nt={str(m): g3s[i] for i, m in enumerate(ms)}, g2_by_Nt={str(m): g2s[i] for i, m in enumerate(ms)},
                  n_by_Nt={str(m): runs[m]['n'] for m in ms}, t_by_Nt={str(m): runs[m]['t'] for m in ms},
                  pop_top_by_Nt={str(m): runs[m]['pop_top_shell'] for m in ms})
        if len(ms) >= 2:
            ex['g3_change_last'] = abs(g3s[-1]-g3s[-2]); ex['g2_change_last'] = abs(g2s[-1]-g2s[-2]); ex['n_change_last'] = abs(runs[ms[-1]]['n']-runs[ms[-2]]['n'])
        if len(ms) >= 3:
            for q, gs in (('g3', g3s), ('g2', g2s)):
                d1, d2 = gs[-2]-gs[-3], gs[-1]-gs[-2]; rr = d2/d1 if d1 != 0 else 0.0
                ex[q+'_ratio'] = rr; ex[q+'_cutoff_err_est'] = abs(d2*rr/(1-rr)) if abs(rr) < 0.9 else abs(d2)
        ex['converged_in_cutoff'] = bool(len(ms) >= 2 and ex['g3_change_last'] < TOL)
        json.dump({k: ex}, open(os.path.join(HERE, f'exactD_U{k}.json'), 'w'), indent=1, default=lambda o: [o.real, o.imag] if isinstance(o, complex) else float(o))
        ex['stopped_by'] = 'converged' if ex['converged_in_cutoff'] else ('cap' if Nt >= ncap else ('time' if r['t'] > TBUDGET else None))
        json.dump({k: ex}, open(os.path.join(HERE, f'exactD_U{k}.json'), 'w'), indent=1, default=lambda o: [o.real, o.imag] if isinstance(o, complex) else float(o))
        if ex['stopped_by']: break
        Nt += 1
    log(f"wrote exactD_U{k}.json (final Nt {Nt}, converged={ex['converged_in_cutoff']})")
