"""Quantum-jump (Monte Carlo wave-function) solver for a driven-dissipative Kerr array with
photon counting on one site.  Runs on numpy (CPU) or cupy (GPU); many trajectories are the
columns of one dense state matrix.

Model (kappa_j = 1):
    H = sum_j [Delta_j n_j + (U/2) n_j(n_j-1) + i(eta_j a_j^dag - eta_j^* a_j)] - J sum_<jk> (a_j^dag a_k + h.c.)
    loss D[a_j] on every site.
Unravelling: waiting-time (norm-decay) method.  The no-jump state is propagated with a fixed-step
classical RK4 for  d psi/dt = -i H_eff psi,  H_eff = H - (i/2) sum_j n_j.  When the squared norm of
a column falls below its uniform random threshold r, the crossing time inside the step is located by
log-linear interpolation of the norm (two secant refinements), the pre-step state is re-propagated
to the crossing time, the jump site is drawn with probability ~ <n_j>, a_j is applied, a new r is
drawn and the step is completed.  Several jumps per step are handled by iterating.  Jump times on
the counting site are stored; all-site jump counts, time-averaged occupations and the population of
the top Fock level of every site are accumulated as truncation diagnostics.
Per-site Fock cutoffs may differ (list `cut`, site j keeps levels 0..cut[j]-1).
"""
import numpy as np, scipy.sparse as sp, time, sys, json, os


def build(cut, det, J, eta, U, ring, kappa=1.0):
    cut = np.asarray(cut, int); K = len(cut); D = int(np.prod(cut))
    strides = np.ones(K, int)
    for j in range(K-2, -1, -1): strides[j] = strides[j+1]*cut[j+1]
    idx = np.arange(D)
    occ = np.stack([(idx // strides[j]) % cut[j] for j in range(K)])          # (K, D)
    a = []
    for j in range(K):
        m = occ[j] > 0
        rows = idx[m] - strides[j]; cols = idx[m]; vals = np.sqrt(occ[j][m]).astype(complex)
        a.append(sp.csr_matrix((vals, (rows, cols)), shape=(D, D)))
    diag = np.zeros(D, complex)
    for j in range(K):
        n = occ[j].astype(float)
        diag += det[j]*n + 0.5*U*n*(n-1) - 0.5j*kappa*n
    H = sp.diags(diag).tocsr()
    for j in range(K):
        if eta[j] != 0:
            H = H + 1j*(eta[j]*a[j].conj().T - np.conj(eta[j])*a[j])
    bonds = [(j, j+1) for j in range(K-1)] + ([(K-1, 0)] if (ring and K > 2) else [])
    for (j, k) in bonds:
        hop = a[j].conj().T @ a[k]
        H = H - J*(hop + hop.conj().T)
    H = H.tocsr(); H.sort_indices()
    top = np.stack([(occ[j] == cut[j]-1).astype(float) for j in range(K)])
    return dict(K=K, D=D, cut=cut, occ=occ.astype(float), top=top, a=a, Heff=H, kappa=kappa)


class Engine:
    def __init__(self, ops, backend='numpy'):
        self.backend = backend
        Hi = (-1j*ops['Heff']).tocsr(); Hi.sort_indices()          # -i H_eff, so that dpsi/dt = Hi psi
        if backend == 'cupy':
            import cupy as cp, cupyx.scipy.sparse as csp
            self.xp = cp
            self.Hi = csp.csr_matrix(Hi); self.a = [csp.csr_matrix(A) for A in ops['a']]
            self.occ = cp.asarray(ops['occ']); self.top = cp.asarray(ops['top'])
            self._axpy = cp.ElementwiseKernel('complex128 x, complex128 k, float64 c', 'complex128 y', 'y = x + c*k', 'axpy_c')
            self._fin = cp.ElementwiseKernel('complex128 x, complex128 k1, complex128 k2, complex128 k3, complex128 k4, float64 c',
                                             'complex128 y', 'y = x + c*(k1 + 2.0*k2 + 2.0*k3 + k4)', 'rk4_final')
            self._nrm = cp.ReductionKernel('complex128 x', 'float64 y', 'x.real()*x.real() + x.imag()*x.imag()', 'a + b', 'y = a', '0', 'norm2_cols')
        else:
            self.xp = np
            self.Hi = Hi; self.a = ops['a']; self.occ = ops['occ']; self.top = ops['top']
            self._axpy = lambda x, k, c: x + c*k
            self._fin = lambda x, k1, k2, k3, k4, c: x + c*(k1 + 2*k2 + 2*k3 + k4)
            self._nrm = None
        self.K = ops['K']; self.D = ops['D']; self.kappa = ops['kappa']

    def rk4(self, psi, dt):
        """classical RK4 step of dpsi/dt = Hi psi; dt scalar or (1, nb) array"""
        H = self.Hi; xp = self.xp
        if not np.isscalar(dt): dt = xp.asarray(dt, dtype=float)
        k1 = H @ psi
        k2 = H @ self._axpy(psi, k1, 0.5*dt)
        k3 = H @ self._axpy(psi, k2, 0.5*dt)
        k4 = H @ self._axpy(psi, k3, dt)
        return self._fin(psi, k1, k2, k3, k4, dt/6.0)

    def norms(self, psi):
        if self._nrm is not None: return self._nrm(psi, axis=0)
        return np.einsum('ij,ij->j', psi.conj(), psi).real

    def advance(self, psi, t0, tau, r, rng, rec, cols, nrefine=2, maxiter=12):
        """batched propagation of columns `psi` (D, nb) over [t0, t0+tau] with jumps.
        tau: scalar or (nb,) array; r: (nb,) thresholds (norm^2 of psi >= r).  Jumps are appended to
        rec as (time, site, column-id) using the ids in `cols`.  Returns psi, r at the end."""
        xp = self.xp
        nb = psi.shape[1]
        rem = xp.full(nb, float(tau)) if np.isscalar(tau) else xp.asarray(tau)
        tcur = xp.full(nb, float(t0))
        out = xp.empty_like(psi); rout = xp.empty(nb); active = xp.arange(nb)
        for it in range(maxiter):
            N0 = self.norms(psi)
            psi1 = self.rk4(psi, rem[None, :])
            N1 = self.norms(psi1)
            ok = N1 >= r
            if it == maxiter-1: ok = xp.ones_like(ok)   # safety: accept remaining columns
            iok = xp.nonzero(ok)[0]
            out[:, active[iok]] = psi1[:, iok]; rout[active[iok]] = r[iok]
            ibad = xp.nonzero(~ok)[0]
            if ibad.size == 0: break
            psi, N0, N1, rem, tcur, r, active = psi[:, ibad], N0[ibad], N1[ibad], rem[ibad], tcur[ibad], r[ibad], active[ibad]
            # crossing time by log-linear interpolation + secant refinements on the bracket
            lo = xp.zeros_like(rem); hi = rem.copy(); Nlo = N0.copy(); Nhi = N1.copy()
            delta = rem*xp.log(N0/r)/xp.log(N0/N1)
            for _ in range(nrefine):
                ps = self.rk4(psi, delta[None, :]); Ns = self.norms(ps)
                up = Ns > r
                lo = xp.where(up, delta, lo); Nlo = xp.where(up, Ns, Nlo)
                hi = xp.where(up, hi, delta); Nhi = xp.where(up, Nhi, Ns)
                delta = lo + (hi-lo)*xp.log(Nlo/r)/xp.log(Nlo/Nhi)
                delta = xp.clip(delta, lo, hi)
            ps = self.rk4(psi, delta[None, :])
            # jump site
            P = self.occ @ (xp.abs(ps)**2)                       # (K, nb_bad), unnormalized
            cdf = xp.cumsum(P, axis=0)/P.sum(0)[None, :]
            u = xp.asarray(rng.uniform(size=ps.shape[1]))
            site = (cdf < u[None, :]).sum(0)
            site = xp.minimum(site, self.K-1)
            new = xp.empty_like(ps)
            for j in range(self.K):
                cj = xp.nonzero(site == j)[0]
                if cj.size: new[:, cj] = self.a[j] @ ps[:, cj]
            new = new/xp.sqrt(self.norms(new))[None, :]
            tj = tcur + delta
            rec.append((self._np(tj), self._np(site), np.asarray(cols)[self._np(active)]))
            r = xp.asarray(rng.uniform(size=ps.shape[1]))
            tcur = tj; rem = rem - delta; psi = new
        return out, rout

    def _np(self, x):
        return x.get() if self.backend == 'cupy' else np.asarray(x)


def run_batch(params, cut, T_total, dt, B, seed, backend='numpy', t_burn=40.0, sample_every=10,
              out=None, verbose=False, start='vacuum'):
    rng = np.random.default_rng(seed)
    ops = build(cut, params['det'], params['J'], params['eta'], params['U'], params['ring'])
    eng = Engine(ops, backend); xp = eng.xp
    D = ops['D']; K = ops['K']
    Psi = xp.zeros((D, B), complex); Psi[0, :] = 1.0
    r = xp.asarray(rng.uniform(size=B))
    nsteps = int(round(T_total/dt))
    rec = []; cols = np.arange(B)
    occ_sum = np.zeros(K); top_sum = np.zeros(K); nsamp = 0
    occ0_traj = np.zeros(B); nsamp_t = 0          # per-trajectory time average of <n_c>, for the paired c1 check
    t = 0.0; t0 = time.time(); tj = 0.0
    for s in range(nsteps):
        Psi_prev = Psi
        Psi = eng.rk4(Psi, dt)
        N = eng.norms(Psi)
        bad = xp.nonzero(N < r)[0]
        if bad.size:
            tb = time.time()
            Psi[:, bad], r[bad] = eng.advance(Psi_prev[:, bad], t, dt, r[bad], rng, rec, cols[eng._np(bad)])
            tj += time.time()-tb
        t += dt
        if t > t_burn and (s % sample_every == 0):
            N = eng.norms(Psi)
            P = xp.abs(Psi)**2 / N[None, :]
            occ_sum += eng._np((eng.occ @ P).sum(1)); top_sum += eng._np((eng.top @ P).sum(1)); nsamp += B
            occ0_traj += eng._np(eng.occ[params['count_site']] @ P); nsamp_t += 1
        if verbose and s % max(1, nsteps//10) == 0:
            print(f"  step {s}/{nsteps} t={t:.1f} wall={time.time()-t0:.0f}s jumps-wall={tj:.0f}s", flush=True)
    wall = time.time()-t0
    times = np.concatenate([x[0] for x in rec]) if rec else np.zeros(0)
    sites = np.concatenate([x[1] for x in rec]).astype(int) if rec else np.zeros(0, int)
    colid = np.concatenate([x[2] for x in rec]).astype(int) if rec else np.zeros(0, int)
    counts_all = np.bincount(sites, minlength=K); counts_all_after = np.bincount(sites[times > t_burn], minlength=K)
    m = sites == params['count_site']
    order = np.lexsort((times[m], colid[m]))
    jt = times[m][order]; cid = colid[m][order]
    offsets = np.concatenate([[0], np.cumsum(np.bincount(cid, minlength=B))])
    res = dict(jump_times=jt, offsets=offsets, T_total=T_total, t_burn=t_burn, dt=dt, B=B, seed=seed,
               cut=np.array(cut), occ_mean=occ_sum/max(nsamp, 1), top_mean=top_sum/max(nsamp, 1),
               occ0_traj=occ0_traj/max(nsamp_t, 1),
               counts_all=counts_all, counts_all_after=counts_all_after, wall=wall, wall_jumps=tj, D=D,
               backend=backend)
    if out: np.savez(out, **res)
    return res


def _worker(args):
    params, cut, T_total, dt, B, seed, out = args
    return run_batch(params, cut, T_total, dt, B, seed, out=out)


if __name__ == "__main__":
    # usage: traj_qjump.py config.json outdir nproc [backend]
    cfg = json.load(open(sys.argv[1])); outdir = sys.argv[2]; nproc = int(sys.argv[3])
    backend = sys.argv[4] if len(sys.argv) > 4 else 'numpy'
    os.makedirs(outdir, exist_ok=True)
    params = cfg['params']; cut = cfg['cut']; T_total = cfg['T_total']; dt = cfg['dt']; B = cfg['B']
    nbatch = cfg['nbatch']; seed0 = cfg.get('seed0', 1000)
    t0 = time.time()
    if backend == 'cupy':
        for i in range(nbatch):
            res = run_batch(params, cut, T_total, dt, B, seed0+i, backend='cupy',
                            out=os.path.join(outdir, f"batch_{i:04d}.npz"), verbose=(i == 0))
            print(f"batch done {i+1}/{nbatch} wall={time.time()-t0:.0f}s per-batch={res['wall']:.0f}s (jumps {res['wall_jumps']:.0f}s) D={res['D']}", flush=True)
    else:
        tasks = [(params, cut, T_total, dt, B, seed0 + i, os.path.join(outdir, f"batch_{i:04d}.npz")) for i in range(nbatch)]
        from multiprocessing import Pool
        with Pool(nproc) as pool:
            for k, res in enumerate(pool.imap_unordered(_worker, tasks)):
                print(f"batch done {k+1}/{nbatch} wall={time.time()-t0:.0f}s per-batch={res['wall']:.0f}s (jumps {res['wall_jumps']:.0f}s) D={res['D']}", flush=True)
    print("ALL DONE", time.time()-t0)


