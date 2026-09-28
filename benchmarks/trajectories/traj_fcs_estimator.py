"""Long-time counting cumulants from quantum-jump records, corrected for the finite counting window.

For every trajectory the counting-site jump times after the burn-in define the counting process N(t).
Counts in overlapping windows, n_T(s) = N(s+T) - N(s) with window starts s on a grid of spacing ds,
give pooled estimates of the window cumulants kappa_m(T) (m = 1, 2, 3) of the stationary process.
For a Markovian process kappa_m(T) = c_m T + b_m + O(e^{-g T}), where g is the Liouvillian gap, so the
plain rate kappa_m(T)/T carries a bias b_m/T.  The corrected estimator is the two-window difference
    c_m[T] = [kappa_m(2T) - kappa_m(T)] / T,
which removes b_m exactly and leaves only O(e^{-gT}) terms.  A least-squares slope of kappa_m(T)
over a range of T is reported as a cross-check.  Standard errors come from a nonparametric bootstrap
over trajectories (trajectories are independent; windows inside a trajectory are not).
usage: traj_fcs_estimator.py rundir [rundir2 ...] --out result.json [--T0 10] [--ds 0.5]
"""
import numpy as np, glob, sys, json, argparse, os

def load(dirs):
    trajs = []; meta = dict(dirs=dirs, nbatch=0, ntraj=0, occ=[], top=[], wall=0.0, files=[])
    occ_traj = []; batch_pairs = []
    for d in dirs:
        for f in sorted(glob.glob(os.path.join(d, "batch_*.npz")) + glob.glob(os.path.join(d, "cpubatch_*.npz"))):   # GPU and CPU batches pooled
            z = np.load(f)
            jt, off, tb, T = z['jump_times'], z['offsets'], float(z['t_burn']), float(z['T_total'])
            nb_count = 0
            for b in range(len(off)-1):
                t = jt[off[b]:off[b+1]]; trajs.append((t[t > tb] - tb, T - tb)); nb_count += int((t > tb).sum())
                occ_traj.append(float(z['occ0_traj'][b]) if 'occ0_traj' in z.files else None)
            batch_pairs.append((nb_count/((len(off)-1)*(T - tb)), float(z['occ_mean'][0]), len(off)-1))
            meta['files'].append(os.path.basename(f)); meta['nbatch'] += 1; meta['occ'].append(z['occ_mean']*int(z['B'])); meta['top'].append(z['top_mean']*int(z['B']))
            meta['wall'] += float(z['wall']); meta['dt'] = float(z['dt']); meta['cut'] = z['cut'].tolist()
            meta['D'] = int(z['D']); meta['t_burn'] = tb
    meta['ntraj'] = len(trajs)
    meta['occ'] = (np.sum(meta['occ'], 0)/len(trajs)).tolist(); meta['top'] = (np.sum(meta['top'], 0)/len(trajs)).tolist()
    meta['total_time'] = float(sum(Te for _, Te in trajs))
    meta['occ_traj'] = occ_traj; meta['batch_pairs'] = batch_pairs
    return trajs, meta

def c1_check(trajs, meta, nboot=1000, seed=3):
    """c1 from jump counts versus kappa <n_c> from the same trajectories (kappa = 1).
    Paired per trajectory when the batches carry occ0_traj; otherwise paired per batch.  The martingale
    N_T - int I dt has variance E[N_T] = c1 T, which gives the expected SE of the difference."""
    n = np.array([len(t) for t, _ in trajs], float); Te = np.array([T for _, T in trajs])
    c1 = n.sum()/Te.sum(); out = dict(c1_counts=float(c1), se_martingale=float(np.sqrt(c1/Te.sum())))
    bp = np.array(meta['batch_pairs'])
    out['occupation_batchmean'] = float(np.average(bp[:, 1], weights=bp[:, 2]))
    if all(o is not None for o in meta['occ_traj']) and len(trajs):
        o = np.array(meta['occ_traj']); rng = np.random.default_rng(seed); diffs = []
        for _ in range(nboot):
            i = rng.integers(0, len(n), len(n)); diffs.append(n[i].sum()/Te[i].sum() - (o[i]*Te[i]).sum()/Te[i].sum())
        occ = (o*Te).sum()/Te.sum()
        out.update(pairing='per trajectory', occupation=float(occ), diff=float(c1-occ), se_diff_bootstrap=float(np.std(diffs, ddof=1)))
    else:
        d = bp[:, 0] - bp[:, 1]
        out.update(pairing='per batch', occupation=out['occupation_batchmean'], diff=float(c1 - out['occupation_batchmean']),
                   se_diff_batches=float(d.std(ddof=1)/np.sqrt(len(d))) if len(d) >= 3 else None, nbatch=len(d))
    out['z_martingale'] = out['diff']/out['se_martingale']
    return out

def suff_stats(trajs, Ts, ds, overlap=True):
    """per-trajectory (nw, S1, S2, S3) for each window T -> array (ntraj, nT, 4)"""
    S = np.zeros((len(trajs), len(Ts), 4))
    for i, (t, Te) in enumerate(trajs):
        ng = int(np.floor(Te/ds))
        Nc = np.concatenate([[0], np.cumsum(np.bincount(np.minimum((t/ds).astype(int), ng), minlength=ng+1)[:ng])]).astype(float)
        for k, T in enumerate(Ts):
            m = int(round(T/ds))
            if overlap:
                n = Nc[m:] - Nc[:-m]
            else:
                g = Nc[::m]; n = g[1:] - g[:-1]
            S[i, k] = [len(n), n.sum(), (n*n).sum(), (n**3).sum()]
    return S

def kappas(Ssum):
    """Ssum (nT, 4) -> kappa1..3 (nT, 3)"""
    N = Ssum[:, 0]; m1, m2, m3 = Ssum[:, 1]/N, Ssum[:, 2]/N, Ssum[:, 3]/N
    return np.stack([m1, m2 - m1**2, m3 - 3*m1*m2 + 2*m1**3], 1)

def estimators(K, Ts, T0, fitT):
    Ts = list(Ts); out = {}
    for k, T in enumerate(Ts):
        c = K[k]/T; out[f"plain T={T:g}"] = np.array([c[0], c[1]/c[0], c[2]/c[0]])
    for T in Ts:
        if 2*T in Ts:
            c = (K[Ts.index(2*T)] - K[Ts.index(T)])/T; out[f"diff T={T:g}"] = np.array([c[0], c[1]/c[0], c[2]/c[0]])
    sel = [Ts.index(T) for T in fitT]; x = np.array(fitT)
    A = np.stack([x, np.ones_like(x)], 1); coef = np.linalg.lstsq(A, K[sel], rcond=None)[0][0]
    out[f"fit T={fitT[0]:g}-{fitT[-1]:g}"] = np.array([coef[0], coef[1]/coef[0], coef[2]/coef[0]])
    out["primary"] = out[f"diff T={T0:g}"]
    return out

def analyze(trajs, Ts, T0, fitT, ds, nboot=1000, seed=1, overlap=True):
    S = suff_stats(trajs, Ts, ds, overlap)
    est = estimators(kappas(S.sum(0)), Ts, T0, fitT)
    rng = np.random.default_rng(seed); n = len(trajs); boots = {k: [] for k in est}
    for _ in range(nboot):
        e = estimators(kappas(S[rng.integers(0, n, n)].sum(0)), Ts, T0, fitT)
        for k in e: boots[k].append(e[k])
    return {k: dict(est=est[k].tolist(), se=np.std(boots[k], 0, ddof=1).tolist()) for k in est}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument('dirs', nargs='+'); ap.add_argument('--out', required=True)
    ap.add_argument('--T0', type=float, default=10.0); ap.add_argument('--ds', type=float, default=0.5)
    ap.add_argument('--Ts', default="2.5,5,7.5,10,15,20,30,40,80"); ap.add_argument('--fit', default="5,7.5,10,15,20")
    ap.add_argument('--nboot', type=int, default=1000); ap.add_argument('--nooverlap', action='store_true')
    a = ap.parse_args()
    Ts = [float(x) for x in a.Ts.split(',')]; fitT = [float(x) for x in a.fit.split(',')]
    trajs, meta = load(a.dirs)
    res = analyze(trajs, Ts, a.T0, fitT, a.ds, a.nboot, overlap=not a.nooverlap)
    res['c1_check'] = c1_check(trajs, meta)
    meta.pop('occ_traj')
    meta.update(T0=a.T0, ds=a.ds, overlap=not a.nooverlap, nboot=a.nboot)
    json.dump(dict(meta=meta, results=res), open(a.out, 'w'), indent=1)
    print('c1_check', json.dumps(res['c1_check']))
    p = res['primary']; print(a.out, meta['ntraj'], f"{meta['total_time']:.3g}", [round(x, 4) for x in p['est']], [round(x, 4) for x in p['se']])

