"""Single-core re-timing of the cumulant closures (orders 2,3,4) at the Fig. 11 couplings.
The steady states themselves are in closure_partial.json; this pass only measures wall-clock
per evaluation with OMP/BLAS threads = 1 and checks that the same fixed point is reached."""
import sys, os, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', '..', 'engines'))
from cumulant_closure import ring_closure, linear_coherent_alpha
P = json.load(open(os.path.join(HERE, '..', '..', 'notebook', 'data', 'revision_params.json')))['ring3']
K = P['K']; PP = dict(kappa=P['kappa'], Delta=P['Delta'], J=P['J'], eta=P['eta'])
part = json.load(open(os.path.join(HERE, 'closure_partial.json')))
out = {}
for N in (2, 3, 4):
    out[str(N)] = {}
    for U in P['U_values_fig11']:
        reps = []
        for rep in range(3):
            t0 = time.time()
            cc = ring_closure(K, N, U=U, **PP)
            t_build = time.time() - t0
            res = cc.steady_state(cc.coherent_state(linear_coherent_alpha(K, **PP)), T=300.0)
            t_tot = time.time() - t0
            st = cc.site_stats(res['x'], 0)
            reps.append(dict(t_build=t_build, t_solve=res['t_ivp']+res['t_polish'], t_total=t_tot, g3=st['g3'], g2=st['g2']))
        k = f"{U:.3f}"; ref = part[str(N)][k]
        med = lambda q: float(np.median([r[q] for r in reps]))
        out[str(N)][k] = dict(t_build=med('t_build'), t_solve=med('t_solve'), t_total=med('t_total'), Nv=cc.Nv,
                              dg3_vs_partial=abs(reps[0]['g3']-ref['g3']), dg2_vs_partial=abs(reps[0]['g2']-ref['g2']))
        print(N, U, out[str(N)][k], flush=True)
json.dump(out, open(os.path.join(HERE, 'closure_timing_1core.json'), 'w'), indent=1)
