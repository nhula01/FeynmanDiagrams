"""Validation of the trajectory counting estimator against exact tilted-Liouvillian cumulants.
For each validation run: the primary estimate (overlapping-window two-window difference, T0 = 10) is compared
with the exact long-time cumulants at the same Fock cutoff; the plain window rates kappa_n(T)/T are compared
with the exact finite-window values; and the bootstrap standard error is checked against the scatter of
independent sub-ensembles (12 groups of trajectories).  usage: traj_validate.py  -> traj_validation.json"""
import json, os, numpy as np, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traj_fcs_estimator as E

def exact_from_log(path, cut):
    for l in open(path):
        d = json.loads(l)
        if d['cut'] == cut:
            th = d.get('theta_r0.08') or d.get('theta_r0.15')
            if th is None:      # exact_rs.py output (Rayleigh-Schroedinger, long-time values only)
                return dict(longtime=d['ratios'][:3], finite_window={}, n=d['n'], top=None, cut=cut, method='rayleigh-schroedinger')
            return dict(longtime=th['ratios'][:3], finite_window=d['finite_window'], n=d['n'], top=d['top'], cut=cut, method='contour')
    return None

Ts = [2.5, 5, 7.5, 10, 15, 20, 30, 40]; fitT = [5, 7.5, 10, 15, 20]
cases = [('K1_U0.05_dt0.05', 'runs/cfg_K1_U0.05_dt0.05', 'exact/K1_U0.05.log', [14]),
         ('K1_U0.05_dt0.1', 'runs/cfg_K1_U0.05_dt0.1', 'exact/K1_U0.05.log', [14]),
         ('K1_U0.10_dt0.05', 'runs/cfg_K1_U0.10_dt0.05', 'exact/K1_U0.10.log', [14]),
         ('K1_U0.10_dt0.1', 'runs/cfg_K1_U0.10_dt0.1', 'exact/K1_U0.10.log', [14]),
         ('K3ring_U0.05', 'runs/cpu_K3_U0.05', 'exact/rsode_K3.log', [10, 7, 7])]
out = dict(description=__doc__, table5_published=dict(K1_U0_05=[0.8600, 1.0747, 1.2523], K1_U0_10=[0.9469, 1.2710, 2.4310], K3_U0_05=[1.1457, 1.1587, 1.6116]),
           note_table5="Published Table V values use five-point differences with h=0.05; the exact references here use a Cauchy contour (K=1) or the Rayleigh-Schroedinger expansion (K=3), both free of step error, at the same Fock cutoff as the trajectories.",
           K3_cutoff_note="K=3 exact at the trajectory cutoff (10,7,7): Fano 1.163562, c3/c1 1.627539; converged value (truncation track, displaced frame): 1.165260, 1.639251. The lab-basis cutoff (10,7,7) therefore lowers Fano by 0.0017 and c3/c1 by 0.012; the trajectories are compared with the same-cutoff value.",
           results={})
rng = np.random.default_rng(7)
for name, d, logf, cut in cases:
    if not os.path.exists(d) or not os.path.exists(logf): continue
    ex = exact_from_log(logf, cut)
    if ex is None: continue
    trajs, meta = E.load([d])
    res = E.analyze(trajs, Ts, 10.0, fitT, 0.5, nboot=1000)
    p = res['primary']; z = [(a-b)/s for a, b, s in zip(p['est'], ex['longtime'], p['se'])]
    fw = {}
    for T in (5, 10, 20):
        key = next((k for k in ex['finite_window'] if float(k) == T), None)
        if key is None: continue
        e_T = ex['finite_window'][key]; t = res[f'plain T={T:g}']
        fw[f'T={T}'] = dict(traj=t['est'], se=t['se'], exact=e_T, z=[(a-b)/s for a, b, s in zip(t['est'], e_T, t['se'])])
    # sub-ensemble check of the bootstrap error
    idx = rng.permutation(len(trajs)); groups = np.array_split(idx, 12); ests = []
    S = E.suff_stats(trajs, Ts, 0.5)
    for g in groups:
        ests.append(E.estimators(E.kappas(S[g].sum(0)), Ts, 10.0, fitT)['primary'])
    ests = np.array(ests); se_groups = (ests.std(0, ddof=1)/np.sqrt(len(groups))).tolist()
    out['results'][name] = dict(ntraj=meta['ntraj'], total_time=meta['total_time'], cut=meta['cut'], dt=meta['dt'],
        exact_longtime=ex['longtime'], exact_method=ex['method'], estimate=p['est'], se_bootstrap=p['se'], se_subensembles=se_groups, z=z,
        finite_window=fw, c1_check=E.c1_check(trajs, meta), occupations_traj=meta['occ'], occupations_exact=ex['n'], cross_checks={k: v for k, v in res.items() if not k.startswith('plain')})
    print(name, [round(x, 4) for x in p['est']], [round(x, 4) for x in p['se']], 'exact', [round(x, 4) for x in ex['longtime']], 'z', [round(x, 2) for x in z], 'se_grp', [round(x, 4) for x in se_groups], flush=True)
json.dump(out, open('traj_validation.json', 'w'), indent=1)

