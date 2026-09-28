"""Assemble cumulant_closure_results.json for the K=3 ring (revision track "Cumulants").
Inputs in this directory: closure_partial.json (closure steady states, orders 2-4), closure_timing_1core.json,
exact_U*.json (+ exact_ref_part*of2.partial.json for U=0, 0.02, 0.025), diag_U*.json."""
import sys, os, json, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', '..', 'engines'))
from cumulant import ring_cumulant_ss
from hfb_chain import hfb_chain
from cumulant_closure import linear_coherent_alpha
J_ = lambda f: json.load(open(os.path.join(HERE, f)))
P = json.load(open(os.path.join(HERE, '..', '..', 'notebook', 'data', 'revision_params.json')))['ring3']
K, kappa, Delta, Jh, eta = P['K'], P['kappa'], P['Delta'], P['J'], P['eta']
U11 = P['U_values_fig11']
closure = J_('closure_partial.json')
timing = J_('closure_timing_1core.json') if os.path.exists(os.path.join(HERE, 'closure_timing_1core.json')) else {}
exact, percube = {}, {}
for f in ('exact_ref_part0of2.partial.json', 'exact_ref_part1of2.partial.json'):
    for k, v in J_(f).items():
        v['source'] = f + ' (per-site cutoff nmax=8,9,10, run before the node reboot)'; percube[k] = v
for f in sorted(glob.glob(os.path.join(HERE, 'exactD_U*.json'))):
    for k, v in json.load(open(f)).items():
        v['source'] = os.path.basename(f) + ' (displaced frame a=alpha_mf+b, fluctuation number cutoff Nt, matrix-free RK4)'; exact[k] = v
for k, v in percube.items():          # fall back to the per-site-cutoff runs only where no Nt run exists
    if k not in exact: exact[k] = v
diag = {}
for f in sorted(glob.glob(os.path.join(HERE, 'diag_U*.json'))):
    diag.update(json.load(open(f)))
import re
diag_t0 = {}
for line in open(os.path.join(HERE, 'run_cumulants.log'), errors='ignore'):
    m = re.match(r"diagrams U=([0-9.]+): .* t=([0-9.]+)s \(order2 only ([0-9.]+)s\)", line)
    if m: diag_t0[f"{float(m.group(1)):.3f}"] = dict(t4=float(m.group(2)), t2=float(m.group(3)))
keys = sorted(set(closure['2']) & set(exact))
rel = lambda a, b: abs(a/b - 1)
def cut_err(ex, q):
    """conservative cutoff uncertainty: the larger of the last two successive changes in the cutoff sequence"""
    seq = ex.get(q+'_by_Nt', ex.get(q+'_by_nmax'))
    if not seq or len(seq) < 2: return None
    v = [seq[k] for k in sorted(seq, key=int)]
    d = [abs(v[i+1]-v[i]) for i in range(len(v)-1)]
    return max(d[-2:])
out = dict(params=P, description="K=3 uniform Kerr ring, on-site n, g2(0), g3(0) of site 0. Cumulant closures of order 2, 3, 4 "
           "(all joint cumulants above the order set to zero) from cumulant_closure.py; diagrams through order 2 and 4 from "
           "multimode.Ring; exact steady state in a displaced frame a=alpha_mf+b with the fluctuation Fock space truncated at total number Nt, converged in Nt (per-site-cutoff runs at U<=0.025 kept as a cross-check).",
           variables_per_order={N: closure[N][keys[0]]['Nv'] for N in ('2', '3', '4')},
           closure_polynomial_terms_at_U0p1={N: closure[N]['0.100']['Nterms'] for N in ('2', '3', '4')})
res = {}
for k in keys:
    ex = exact[k]; row = dict(U=float(k), exact=dict(n=ex['n'], g2=ex['g2'], g3=ex['g3'], cutoff=ex.get('Nt', ex.get('nmax')), truncation=ex.get('truncation', 'per-site n_j < nmax'), D=ex.get('D'),
               g3_change_last=ex.get('g3_change_last'), g2_change_last=ex.get('g2_change_last'),
               g3_cutoff_err_est=cut_err(ex, 'g3'), g2_cutoff_err_est=cut_err(ex, 'g2'),
               g3_by_cutoff=ex.get('g3_by_Nt', ex.get('g3_by_nmax')), t_final_cutoff=ex['t_solve'], t_all_cutoffs=float(sum(ex.get('t_by_Nt', ex.get('t_by_nmax', {})).values())),
               converged_in_cutoff=ex.get('converged_in_cutoff'), source=ex['source']))
    for N in ('2', '3', '4'):
        c = closure[N][k]
        d = dict(n=c['n'], g2=c['g2'], g3=c['g3'], converged=c['converged'], residual=c['residual'], max_re_eig=c.get('max_re_eig'),
                 translation_asymmetry=c.get('translation_asymmetry'), imag_parts=c.get('imag_parts'), init=c.get('init'),
                 err_n=rel(c['n'], ex['n']), err_g2=rel(c['g2'], ex['g2']), err_g3=rel(c['g3'], ex['g3']),
                 t_solve_original_run=c['t_solve'], t_build_original_run=c['t_generate'] + c['t_closure'])
        if N in timing and k in timing[N]:
            d.update(t_build_1core=timing[N][k]['t_build'], t_solve_1core=timing[N][k]['t_solve'], t_total_1core=timing[N][k]['t_total'])
        row[f'closure{N}'] = d
    if k in diag:
        dg = diag[k]
        for N in ('2', '4'):
            row[f'diagrams{N}'] = dict(n=dg[N]['n'], g2=dg[N]['g2'], g3=dg[N]['g3'], err_n=rel(dg[N]['n'], ex['n']),
                                        t_same_run_as_closures=diag_t0.get(k, {}).get('t4' if N == '4' else 't2'),
                                        err_g2=rel(dg[N]['g2'], ex['g2']), err_g3=rel(dg[N]['g3'], ex['g3']),
                                        t_total_1core=dg['t_order4_g2_g3'] if N == '4' else dg['t_order2_g2_g3'],
                                        t_g2_only_1core=dg['t_order4_g2'] if N == '4' else dg['t_order2_g2'])
        row['diagram_series'] = dict(n=dg['series_n'], n2=dg['series_n2'], n3=dg['series_n3'])
    res[k] = row
out['results'] = res
# validation
al = linear_coherent_alpha(K, kappa, Delta, Jh, eta)
v2 = {}
for k in keys:
    U = float(k); rc = ring_cumulant_ss(K, kappa, Delta, Jh, eta, U); hc, r_ = hfb_chain(K, kappa, [Delta]*K, Jh, eta, U, ring=True)
    c = closure['2'][k]
    v2[k] = dict(diff_n_vs_cumulant_py=abs(c['n']-float(rc['n'])), diff_g2_vs_cumulant_py=abs(c['g2']-float(rc['g2'])),
                 diff_n_vs_hfb_chain=abs(c['n']-float(hc[0]['n'])), diff_g2_vs_hfb_chain=abs(c['g2']-float(hc[0]['g2'])),
                 g3_hfb_chain=float(hc[0]['g3']), g3_closure2=c['g3'], diff_g3_vs_hfb_chain=abs(c['g3']-float(hc[0]['g3'])))
xc = {k: dict(g3_Nt=exact[k]['g3'], g3_percube_nmax10=percube[k]['g3'], diff_g3=abs(exact[k]['g3']-percube[k]['g3']),
              g2_Nt=exact[k]['g2'], g2_percube_nmax10=percube[k]['g2'], diff_g2=abs(exact[k]['g2']-percube[k]['g2']),
              n_diff=abs(exact[k]['n']-percube[k]['n']), percube_g3_change_last=percube[k]['g3_change_last'])
      for k in percube if 'Nt' in exact[k]}
out['validation'] = dict(
    exact_truncation_crosscheck=xc,
    second_order_vs_existing=v2,
    second_order_max_abs_diff=dict(n=max(v['diff_n_vs_cumulant_py'] for v in v2.values()), g2=max(v['diff_g2_vs_cumulant_py'] for v in v2.values())),
    hfb_chain_g3_note="hfb_chain.py evaluates the Gaussian <a^dag^3 a^3> with 18 n|m|^2 in the cubic fluctuation term; Wick's theorem gives "
                      "<b^dag^3 b^3> = 6n^3 + 9n|m|^2 (checked on a displaced squeezed thermal state in QuTiP). The g3 difference in "
                      "second_order_vs_existing is entirely this term; n and g2 agree to machine precision.",
    U0_coherent={N: dict(n_minus_alpha2=abs(closure[N]['0.000']['n']-abs(al[0])**2), g2_minus_1=abs(closure[N]['0.000']['g2']-1),
                         g3_minus_1=abs(closure[N]['0.000']['g3']-1)) for N in ('2', '3', '4')},
    small_U_0_02={N: dict(err_n=res['0.020'][f'closure{N}']['err_n'], err_g2=res['0.020'][f'closure{N}']['err_g2'],
                          err_g3=res['0.020'][f'closure{N}']['err_g3']) for N in ('2', '3', '4')},
    small_U_note="the exact reference at U=0.02 has |g3(10)-g3(9)| = %.1e, so errors below ~1e-6 are not resolved" %
                 exact['0.020'].get('g3_change_last', float('nan')))
out['exact_per_site_cutoff_runs'] = percube
# Fig. 11 table + verdict
tab = []
for U in U11:
    k = f"{U:.3f}"; r = res[k]; row = dict(U=U, exact_g2=r['exact']['g2'], exact_g3=r['exact']['g3'], exact_cutoff=r['exact']['cutoff'],
                                          exact_g3_cutoff_err=r['exact']['g3_cutoff_err_est'], exact_t=r['exact']['t_final_cutoff'])
    for m in ('closure2', 'closure3', 'closure4', 'diagrams2', 'diagrams4'):
        if m in r:
            t0 = (r[m]['t_build_original_run'] + r[m]['t_solve_original_run']) if m.startswith('closure') else r[m].get('t_same_run_as_closures')
            row[m] = dict(g2=r[m]['g2'], g3=r[m]['g3'], err_g2=r[m]['err_g2'], err_g3=r[m]['err_g3'], t=t0, t_rerun_contended=r[m].get('t_total_1core'))
    best = {q: min(('closure2', 'closure3', 'closure4'), key=lambda m: r[m][f'err_{q}']) for q in ('g2', 'g3')}
    row['best_closure'] = best
    if 'diagrams4' in r:
        row['ratio_best_closure_err_over_diag4_err'] = {q: r[best[q]][f'err_{q}']/r['diagrams4'][f'err_{q}'] for q in ('g2', 'g3')}
        row['ratio_closure4_err_over_diag4_err'] = {q: r['closure4'][f'err_{q}']/r['diagrams4'][f'err_{q}'] for q in ('g2', 'g3')}
    tab.append(row)
out['table_fig11'] = tab
out['timing_note'] = ("'t' in table_fig11 and t_*_original_run / t_same_run_as_closures: wall-clock per evaluation (build + steady state for the closures; "
    "n, g2 and g3 series through the stated order for the diagrams) measured in one process, one after the other, in the run of 2026-09-24 15:33-16:13 "
    "(run_cumulants.py). t_*_1core: single-thread re-measurement on 2026-09-24 17:30-18:40 while the node was oversubscribed (load 80-150 on 96 threads); "
    "those numbers scatter by up to a factor 10 between identical calls and are kept only as a record.")
json.dump(out, open(os.path.join(HERE, 'cumulant_closure_results.json'), 'w'), indent=1, default=float)
for row in tab:
    print(f"U={row['U']:.2f} ex g3={row['exact_g3']:.6f}(cut {row['exact_cutoff']}, ±{row['exact_g3_cutoff_err'] or 0:.1e}) | " +
          " ".join(f"{m}: g3={row[m]['g3']:.5f} e3={row[m]['err_g3']:.2e} e2={row[m]['err_g2']:.2e}" for m in ('closure2','closure3','closure4','diagrams2','diagrams4') if m in row)
          + f" | best {row['best_closure']} ratio {row.get('ratio_best_closure_err_over_diag4_err')}")
print('validation max diffs', out['validation']['second_order_max_abs_diff'], out['validation']['small_U_0_02'])

# optimally truncated diagrams from the truncation track, re-scored against the converged exact reference
ebf = os.path.join(HERE, '..', 'truncation', 'fig_g3_errorbars.json')
if os.path.exists(ebf):
    E = json.load(open(ebf)); ot = {}
    for i, u in enumerate(E['U']):
        k = f"{u:.3f}"
        if k not in res: continue
        ex = res[k]['exact']; rec = dict(U=u)
        for q in ('g2', 'g3'):
            v, d = E[q][i], E[q+'_delta'][i]; err = abs(v/ex[q]-1)
            best = min(res[k][m]['err_'+q] for m in ('closure2', 'closure3', 'closure4'))
            rec[q] = dict(value=v, delta=d, N_star=E[q+'_Nstar'][i], their_exact=E[q+'_exact'][i], err_vs_converged_exact=err,
                          true_abs_error_over_delta=(abs(v-ex[q])/d if d > 0 else None), best_closure_err=best,
                          best_closure_err_over_optimal_diagram_err=(best/err if err > 0 else None),
                          their_exact_rel_diff_from_converged=abs(E[q+'_exact'][i]/ex[q]-1))
        ot[k] = rec
    out['optimal_truncation_diagrams'] = dict(source='../truncation/fig_g3_errorbars.json', results=ot)
    for k, rec in ot.items():
        print('opt', k, 'g3 err %.2e delta_rel %.2e true/delta %.1f bestclosure/opt %.2f' % (rec['g3']['err_vs_converged_exact'], rec['g3']['delta']/res[k]['exact']['g3'],
              rec['g3']['true_abs_error_over_delta'] or 0, rec['g3']['best_closure_err_over_optimal_diagram_err'] or 0),
              '| g2 err %.2e bestclosure/opt %.2f' % (rec['g2']['err_vs_converged_exact'], rec['g2']['best_closure_err_over_optimal_diagram_err'] or 0))
json.dump(out, open(os.path.join(HERE, 'cumulant_closure_results.json'), 'w'), indent=1, default=float)
