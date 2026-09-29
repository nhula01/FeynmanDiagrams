"""Assemble the eight-site counting benchmark: trajectory estimates, diagrams (contour derivatives),
Gaussian FCS (Richardson differences), mean field; write chain8_contour.json, traj_fcs_results.json
and fig_fcs.pdf.  usage: traj_fcs_final.py  (run in benchmarks/trajectories)"""
import json, os, numpy as np
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
from scipy.integrate import solve_ivp
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

RP = json.load(open(os.path.join(_ROOT, 'notebook', 'data', 'revision_params.json')))['chain8_fcs']
J = lambda f: json.load(open(f)) if os.path.exists(f) else None
import glob
def merge(pattern):
    fs = sorted(glob.glob(pattern))
    if not fs: return None
    out = json.load(open(fs[0])); out['results'] = {}
    for f in fs: out['results'].update(json.load(open(f))['results'])
    return out
d5, d6 = merge('diag/diagRS_chain8_maxdeg5_U*.json'), merge('diag/diagRS_chain8_maxdeg6_U*.json')   # full RS recursion
d5nb = J('diag/diag_chain8_maxdeg5.json')          # legacy no-renormalization notebook snapshot; kept for historical reference
gs = J('diag/gauss_chain8_root.json')

def mean_field(U):
    K = RP['K']; Om = np.diag(RP['detunings']).astype(complex)
    for j in range(K-1): Om[j, j+1] = Om[j+1, j] = -RP['J']
    G = 0.5*np.eye(K) + 1j*Om; eta = np.array(RP['eta'], complex)
    f = lambda t, y: (lambda a: np.concatenate([(-G@a + eta - 1j*U*abs(a)**2*a).real, (-G@a + eta - 1j*U*abs(a)**2*a).imag]))(y[:K]+1j*y[K:])
    a0 = np.linalg.solve(G, eta); s = solve_ivp(f, (0, 400), np.concatenate([a0.real, a0.imag]), rtol=1e-12, atol=1e-14)
    a = s.y[:K, -1] + 1j*s.y[K:, -1]; return float(abs(a[0])**2)

chain = dict(description="Eight-site chain (revision_params.json chain8_fcs), counting on site 0. Diagrams: partial sums of "
             "theta_N(chi) through order N with the full Rayleigh-Schroedinger recursion (TiltedNZ, benchmarks/truncation/fcs_nz.py; "
             "'diagrams_maxdeg5_notebook_recursion' = a legacy checked-in no-renormalization snapshot; current fcs_multi.TiltedChain includes -sum_k theta_k psi_(N-k)), chi-derivatives from a Cauchy contour |chi|=0.1 (16 points; "
             "|chi|=0.05 check). Gaussian: stationary point of the tilted Riccati system (root finder), five-point differences at h=0.05,0.025,0.0125 Richardson-"
             "extrapolated. Mean field: c1=kappa|alpha_0|^2, Fano=c3/c1=1. 'fd_h0.05' = the paper's finite-difference convention. "
             "Entries are [c1/kappa, c2/c1, c3/c1]; the diagrams are reported at sixth order, the highest computed.",
             results={})
for key in (d5 or {}).get('results', {}):
    U = float(key); r5 = d5['results'][key]; rec = dict(U=U)
    rec['diagrams_maxdeg5'] = {'orders': r5['contour_r0.1'], 'orders_r0.05': r5['contour_r0.05'], 'orders_fd_h0.05': r5['fd_h0.05']}
    if d5nb and key in d5nb['results']:
        rec['diagrams_maxdeg5_notebook_recursion'] = {'orders': d5nb['results'][key]['contour_r0.1'], 'orders_fd_h0.05': d5nb['results'][key]['fd_h0.05']}
    if d6 and key in d6['results']:
        r6 = d6['results'][key]; rec['diagrams_maxdeg6'] = {'orders': r6['contour_r0.1'], 'orders_r0.05': r6['contour_r0.05'], 'orders_fd_h0.05': r6['fd_h0.05']}
    if gs and key in gs['results']:
        g = gs['results'][key]; rec['gaussian'] = {'value': g['richardson'], 'fd_h0.05': g['fd_h0.05'], 'richardson_change_last': g['richardson_change_last']}
    mf = mean_field(U); rec['mean_field'] = [mf, 1.0, 1.0]
    chain['results'][key] = rec
json.dump(chain, open('chain8_contour.json', 'w'), indent=1)

# ---- trajectory estimates
MIN_TIME = 1.0e6          # trajectory time below which a point is reported as provisional and left out of the figure
runs = [(0.04, 'est/est_gpu_U0.04_c10.json'), (0.08, 'est/est_gpu_U0.08_c12.json'), (0.10, 'est/est_gpu_U0.10_c14.json')]
traj = dict(description="Quantum-jump (waiting-time MCWF) counting cumulants of site 0, eight-site chain. Estimator: overlapping "
            "windows (start spacing 0.5/kappa), two-window difference c_n=[kappa_n(2T0)-kappa_n(T0)]/T0 with T0=10/kappa; "
            "bootstrap over trajectories (1000 resamples). Values [c1/kappa, c2/c1, c3/c1].", results={})
for U, f in runs + [(0.10, 'est/est_gpu_U0.10_cut8.json')]:
    e = J(f)
    if e is None: continue
    m = e['meta']; tag = f"{U:.2f}" + ('_cut8' if 'cut8' in f else '')
    traj['results'][tag] = dict(U=U, estimate=e['results']['primary']['est'], se=e['results']['primary']['se'],
        cross_checks={k: v for k, v in e['results'].items() if (k.startswith('diff') or k.startswith('fit'))},
        plain_windows={k: v for k, v in e['results'].items() if k.startswith('plain')},
        ntraj=m['ntraj'], total_time=m['total_time'], cut=m['cut'], D=m['D'], dt=m['dt'], t_burn=m['t_burn'],
        status='final' if m['total_time'] >= MIN_TIME else f'provisional (total time {m["total_time"]:.2g} < {MIN_TIME:.0e}; not plotted)',
        nbatch=m['nbatch'], c1_check=e['results'].get('c1_check'), occupations=m['occ'], top_level_population=m['top'], wall_s=m['wall'])
for tag, v in traj['results'].items():
    key = f"{v['U']:.2f}"
    if key not in chain['results']: continue
    ch = chain['results'][key]; est, se = np.array(v['estimate']), np.array(v['se'])
    refs = {f'diagrams_N{N}': ch['diagrams_maxdeg5']['orders'][N] for N in range(1, 7)}
    refs['diagrams_N3_paper_fd_h0.05'] = ch['diagrams_maxdeg5']['orders_fd_h0.05'][3]
    if 'diagrams_maxdeg5_notebook_recursion' in ch:
        refs['diagrams_N3_notebook_recursion_as_published'] = ch['diagrams_maxdeg5_notebook_recursion']['orders_fd_h0.05'][3]
    if 'diagrams_maxdeg6' in ch:
        refs.update({f'diagrams_maxdeg6_N{N}': ch['diagrams_maxdeg6']['orders'][N] for N in range(3, len(ch['diagrams_maxdeg6']['orders']))})
    refs['gaussian'] = ch['gaussian']['value']; refs['mean_field'] = ch['mean_field']
    v['comparison'] = dict(note="z = (trajectory - reference)/SE_trajectory for [c1, c2/c1, c3/c1]",
        references={k: dict(value=r, z=((est - np.array(r))/se).round(2).tolist()) for k, r in refs.items()})
json.dump(traj, open('traj_fcs_results.json', 'w'), indent=1)

# ---- figure
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9, "legend.frameon": False,
                     "xtick.direction": "in", "ytick.direction": "in"})
COL = ["#d95f02", "#1b9e77", "#7570b3", "#e7298a", "#66a61e", "#a6761d"]
keys = sorted(chain['results'], key=float); Us = np.array([float(k) for k in keys])
N6 = np.array([chain['results'][k]['diagrams_maxdeg5']['orders'][-1] for k in keys])      # diagrams through sixth order
Gs = np.array([chain['results'][k]['gaussian']['value'] for k in keys])
fig, ax = plt.subplots(figsize=(3.375, 2.95))
ax.axhline(1.0, color='0.5', ls=(0, (4, 2)), lw=1.0, label='mean field')
for i, (q, lab) in enumerate([(1, r'$c_2/c_1$'), (2, r'$c_3/c_1$')]):
    ax.plot(Us, Gs[:, q], color=COL[i], lw=0.9, ls=(0, (1, 1.2)))
    ax.plot(Us, N6[:, q], 'o-', color=COL[i], ms=3, lw=1.1, label=lab + ' diagrams')
first = True
for k, v in sorted(traj['results'].items()):
    if k.endswith('cut8') or v['status'] != 'final': continue
    for i, q in enumerate([1, 2]):
        ax.errorbar(v['U'], v['estimate'][q], yerr=v['se'][q], fmt='D', ms=3.2, mfc='white', mec='k', ecolor='k', elinewidth=0.8, capsize=1.5,
                    label='trajectories' if first else None, zorder=5); first = False
ax.plot([], [], color='0.3', lw=0.9, ls=(0, (1, 1.2)), label='Gaussian')
ax.set_xlabel(r'$U/\kappa$'); ax.set_ylabel('counting cumulant ratio'); ax.set_xlim(-0.004, 0.104)
ymax = max(2.0, max([v['estimate'][2] + v['se'][2] for k, v in traj['results'].items() if not k.endswith('cut8') and v['status'] == 'final'] + [0]) + 0.08)
ymax = max(ymax, float(np.max(N6[:, 2])) + 0.05)
ax.set_ylim(0.95, ymax)
ax.legend(fontsize=6.5, loc='lower left', bbox_to_anchor=(0.0, 1.01), ncol=3, handlelength=1.6, columnspacing=0.9, borderaxespad=0.0)
ins = ax.inset_axes([0.17, 0.52, 0.38, 0.32])
for i, q in enumerate([1, 2]):
    ins.plot(Us, N6[:, q] - Gs[:, q], 'o-', color=COL[i], ms=2, lw=0.9)
    for k, v in traj['results'].items():
        if k.endswith('cut8') or v['status'] != 'final': continue
        g = Gs[list(Us).index(v['U'])] if v['U'] in list(Us) else None
        if g is not None:
            ins.errorbar(v['U'] + (0.0015 if q == 2 else -0.0015), v['estimate'][q] - g[q], yerr=v['se'][q], fmt='D', ms=2, mfc='white', mec=COL[i], ecolor=COL[i], elinewidth=0.6, capsize=1)
ins.axhline(0, color='0.6', lw=0.5); ins.tick_params(labelsize=6); ins.set_title('minus Gaussian', fontsize=6.5, pad=2); ins.set_xticks([0, 0.05, 0.1])
fig.tight_layout(); fig.savefig('fig_fcs.pdf'); fig.savefig('fig_fcs_preview.png', dpi=200)
print(json.dumps({k: [round(x, 4) for x in v['estimate']] + [round(x, 4) for x in v['se']] for k, v in traj['results'].items()}))

