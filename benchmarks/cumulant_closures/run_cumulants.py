"""Driver: exact reference (with nmax convergence), cumulant closures N=2,3,4, diagrams N=2,4,
validation checks, timings.  Writes cumulant_closure_results.json."""
import sys, os, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(HERE, '..', '..', 'engines')
sys.path.insert(0, HERE); sys.path.insert(0, CODE)
import scipy.sparse as sp, scipy.sparse.linalg as spl
from cumulant_closure import ring_closure, linear_coherent_alpha, unit
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

# ------------------------------------------------------------- closures
def run_closure(N, U, x0=None, T=300.0):
    cc = ring_closure(K, N, U=U, **PP)
    if x0 is None:
        x0 = cc.coherent_state(linear_coherent_alpha(K, **PP))
    res = cc.steady_state(x0, T=T)
    rec = dict(N=N, U=U, Nv=cc.Nv, Nh=cc.Nh, Nterms=cc.Nterms, t_generate=cc.t_generate, t_closure=cc.t_closure,
               t_ivp=res['t_ivp'], t_polish=res['t_polish'], t_solve=res['t_ivp']+res['t_polish'],
               converged=res['converged'], diverged=res['diverged'], residual=res.get('residual', None),
               ivp_residual=res['ivp_residual'], t_end=res['t_end'], nfev=res['ivp_nfev'], hermiticity=res.get('hermiticity'))
    if res['x'] is not None:
        st = cc.site_stats(res['x'], 0)
        rec.update(n=st['n'].real, g2=st['g2'], g3=st['g3'], imag_parts=max(st['imag_n'], st['imag_n2'], st['imag_n3']))
        # translation symmetry check (not assumed): compare site 1 and 2 with site 0
        st1 = cc.site_stats(res['x'], 1); st2 = cc.site_stats(res['x'], 2)
        rec['translation_asymmetry'] = max(abs(st1['g3']-st['g3']), abs(st2['g3']-st['g3']), abs(st1['n'].real-st['n'].real))
        # linear stability of the fixed point
        Jm = cc.jacobian(res['x']); ev = np.linalg.eigvals(Jm)
        rec['max_re_eig'] = float(np.max(ev.real)); rec['slowest_rate'] = float(-np.max(ev.real[ev.real < 0])) if np.any(ev.real < 0) else None
    return cc, res, rec

closure = {2: {}, 3: {}, 4: {}}
for N in (2, 3, 4):
    prev = None
    for U in U_all:
        cc, res, rec = run_closure(N, U)
        rec['init'] = 'coherent'
        if not rec['converged'] and prev is not None:
            log(f"  N={N} U={U}: coherent-state start did not converge (diverged={rec['diverged']}); retrying by continuation from U={prev[0]}")
            cc, res, rec2 = run_closure(N, U, x0=prev[1], T=600.0)
            rec2['init'] = f'continuation from U={prev[0]}'
            rec2['coherent_start'] = {k: rec[k] for k in ('converged', 'diverged', 'residual', 'ivp_residual', 't_end')}
            rec = rec2
        if rec['converged']: prev = (U, res['x'])
        closure[N][f"{U:.3f}"] = rec
        log(f"N={N} U={U}: conv={rec['converged']} res={rec.get('residual')} n={rec.get('n')} g2={rec.get('g2')} g3={rec.get('g3')} "
            f"t={rec['t_solve']:.2f}s maxReEig={rec.get('max_re_eig')} init={rec['init']}")
json.dump({str(N): closure[N] for N in closure}, open(os.path.join(HERE, 'closure_partial.json'), 'w'), indent=1, default=float)

# ------------------------------------------------------------- diagrams (multimode Ring), orders 2 and 4
diag = {}
for U in U_all:
    t0 = time.time()
    r = MM.Ring(K=K, U=U, **PP); b0, bd0 = r.site_ops(0)
    n_ser = np.cumsum(r.series(MM.mul(bd0, b0), 4)).real
    n2_ser = np.cumsum(r.series(MM.mul(MM.mul(bd0, bd0), MM.mul(b0, b0)), 4)).real
    t_after_g2 = time.time() - t0
    t1 = time.time()
    op3 = MM.mul(MM.mul(MM.mul(bd0, bd0), bd0), MM.mul(MM.mul(b0, b0), b0))
    n3_ser = np.cumsum(r.series(op3, 4)).real
    t_g3 = time.time() - t1
    rec = {}
    for N in (2, 4):
        g2 = n2_ser[N]/n_ser[N]**2; g3 = n3_ser[N]/n_ser[N]**3
        rec[str(N)] = dict(n=float(n_ser[N]), g2=float(g2), g3=float(g3))
    rec['series_n'] = [float(v) for v in n_ser]; rec['series_n2'] = [float(v) for v in n2_ser]; rec['series_n3'] = [float(v) for v in n3_ser]
    rec['t_all_orders_to_4_n_g2'] = t_after_g2; rec['t_all_orders_to_4_g3'] = t_g3; rec['t_total'] = time.time()-t0
    # time for order 2 only
    t2 = time.time(); r2 = MM.Ring(K=K, U=U, **PP); b0, bd0 = r2.site_ops(0)
    r2.series(MM.mul(bd0, b0), 2); r2.series(MM.mul(MM.mul(bd0, bd0), MM.mul(b0, b0)), 2); r2.series(op3, 2)
    rec['t_order2_only'] = time.time()-t2
    diag[f"{U:.3f}"] = rec
    log(f"diagrams U={U}: N2 g3={rec['2']['g3']:.6f} N4 g3={rec['4']['g3']:.6f}  t={rec['t_total']:.1f}s (order2 only {rec['t_order2_only']:.1f}s)")
out['diagrams'] = diag

# ------------------------------------------------------------- exact reference (from exact_ref.py, possibly still running)
fx = os.path.join(HERE, 'exact_ref.json')
while not os.path.exists(fx):
    log('waiting for exact_ref.json ...'); time.sleep(60)
time.sleep(5); exact = json.load(open(fx)); out['exact'] = exact
for N in (2, 3, 4):
    for k, rec in closure[N].items():
        ex = exact[k]
        if rec.get('g3') is not None:
            rec['err_n'] = abs(rec['n']/ex['n']-1); rec['err_g2'] = abs(rec['g2']/ex['g2']-1); rec['err_g3'] = abs(rec['g3']/ex['g3']-1)
        log(f"N={N} U={k}: errs n={rec.get('err_n')} g2={rec.get('err_g2')} g3={rec.get('err_g3')}")
for k, rec in diag.items():
    ex = exact[k]
    for N in ('2', '4'):
        rec[N]['err_n'] = abs(rec[N]['n']/ex['n']-1); rec[N]['err_g2'] = abs(rec[N]['g2']/ex['g2']-1); rec[N]['err_g3'] = abs(rec[N]['g3']/ex['g3']-1)
    log(f"diagrams U={k}: errs g3 N2={rec['2']['err_g3']:.2e} N4={rec['4']['err_g3']:.2e}; g2 errs N2={rec['2']['err_g2']:.2e} N4={rec['4']['err_g2']:.2e}")
out['closure'] = {str(N): closure[N] for N in closure}

# ------------------------------------------------------------- validation
val = {}
# (i) second order vs cumulant.py (momentum-space HFB) and hfb_chain.py (site basis, gives g3 by Wick)
v2 = {}
for U in U_all:
    rc = ring_cumulant_ss(K, kappa, Delta, J, eta, U)
    hc, resid = hfb_chain(K, kappa, [Delta]*K, J, eta, U, ring=True)
    mine = closure[2][f"{U:.3f}"]
    v2[f"{U:.3f}"] = dict(cumulant_py=dict(n=rc['n'], g2=rc['g2'], converged=bool(rc['converged'])),
                          hfb_chain=dict(n=hc[0]['n'], g2=hc[0]['g2'], g3=hc[0]['g3'], residual=float(resid)),
                          diff_n_cumulant_py=abs(mine['n']-rc['n']), diff_g2_cumulant_py=abs(mine['g2']-rc['g2']),
                          diff_n_hfb_chain=abs(mine['n']-hc[0]['n']), diff_g2_hfb_chain=abs(mine['g2']-hc[0]['g2']),
                          diff_g3_hfb_chain=abs(mine['g3']-hc[0]['g3']))
    log(f"validate N=2 U={U}: |dn|={v2[f'{U:.3f}']['diff_n_cumulant_py']:.2e} |dg2|={v2[f'{U:.3f}']['diff_g2_cumulant_py']:.2e} |dg3 vs hfb_chain|={v2[f'{U:.3f}']['diff_g3_hfb_chain']:.2e}")
val['second_order_vs_existing'] = v2
# (ii) U=0: every order must reproduce the coherent state exactly (g2=g3=1, n=|alpha|^2)
al = linear_coherent_alpha(K, **PP)
val['U0_coherent'] = {str(N): dict(n_minus_alpha2=abs(closure[N]['0.000']['n']-abs(al[0])**2), g2_minus_1=abs(closure[N]['0.000']['g2']-1),
                                   g3_minus_1=abs(closure[N]['0.000']['g3']-1), exact_g2_minus_1=abs(exact['0.000']['g2']-1),
                                   exact_g3_minus_1=abs(exact['0.000']['g3']-1)) for N in (2, 3, 4)}
# (iii) small U: errors must decrease with order
val['small_U_0.02'] = {str(N): dict(err_n=closure[N]['0.020']['err_n'], err_g2=closure[N]['0.020']['err_g2'], err_g3=closure[N]['0.020']['err_g3']) for N in (2, 3, 4)}
out['validation'] = val

# ------------------------------------------------------------- summary table at fig11 points
tab = []
for U in U_fig11:
    k = f"{U:.3f}"; row = dict(U=U, exact_g2=exact[k]['g2'], exact_g3=exact[k]['g3'])
    for N in (2, 3, 4):
        r = closure[N][k]
        row[f'closure{N}_g3'] = r.get('g3'); row[f'closure{N}_err_g2'] = r.get('err_g2'); row[f'closure{N}_err_g3'] = r.get('err_g3'); row[f'closure{N}_converged'] = r['converged']; row[f'closure{N}_t'] = r['t_solve']
    for N in (2, 4):
        row[f'diag{N}_g3'] = diag[k][str(N)]['g3']; row[f'diag{N}_err_g2'] = diag[k][str(N)]['err_g2']; row[f'diag{N}_err_g3'] = diag[k][str(N)]['err_g3']
    row['diag_t_total_to_order4'] = diag[k]['t_total']; row['diag_t_order2'] = diag[k]['t_order2_only']
    tab.append(row)
out['table_fig11'] = tab
out['variables_per_order'] = {str(N): closure[N]['0.100']['Nv'] for N in (2, 3, 4)}
json.dump(out, open(os.path.join(HERE, 'cumulant_closure_results.json'), 'w'), indent=1, default=float)
log("wrote cumulant_closure_results.json")
