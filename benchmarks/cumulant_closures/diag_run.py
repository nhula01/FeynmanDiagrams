"""Diagrams (multimode.Ring) through order 4 for the K=3 ring: n, g2, g3 of site 0 with single-core timings.
Usage: python diag_run.py U1,U2,...   Writes diag_U{U:.3f}.json per U."""
import sys, os, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(HERE, '..', '..', 'engines')
sys.path.insert(0, CODE)
import multimode as MM
P = json.load(open(os.path.join(HERE, '..', '..', 'notebook', 'data', 'revision_params.json')))['ring3']
K, kappa, Delta, J, eta = P['K'], P['kappa'], P['Delta'], P['J'], P['eta']
PP = dict(kappa=kappa, Delta=Delta, J=J, eta=eta)
for U in [float(x) for x in sys.argv[1].split(',')]:
    rec = {}
    # order-2 only (separate object so no cache is shared)
    t2 = time.time(); r2 = MM.Ring(K=K, U=U, **PP); b0, bd0 = r2.site_ops(0)
    op1 = MM.mul(bd0, b0); op2 = MM.mul(MM.mul(bd0, bd0), MM.mul(b0, b0)); op3 = MM.mul(MM.mul(MM.mul(bd0, bd0), bd0), MM.mul(MM.mul(b0, b0), b0))
    s1 = r2.series(op1, 2); s2 = r2.series(op2, 2); t_g2_o2 = time.time()-t2; s3 = r2.series(op3, 2)
    rec['t_order2_g2'] = t_g2_o2; rec['t_order2_g2_g3'] = time.time()-t2
    # through order 4
    t0 = time.time(); r = MM.Ring(K=K, U=U, **PP); b0, bd0 = r.site_ops(0)
    n_ser = np.cumsum(r.series(op1, 4)).real; n2_ser = np.cumsum(r.series(op2, 4)).real
    rec['t_order4_g2'] = time.time()-t0
    n3_ser = np.cumsum(r.series(op3, 4)).real
    rec['t_order4_g2_g3'] = time.time()-t0
    for N in (2, 4):
        rec[str(N)] = dict(n=float(n_ser[N]), g2=float(n2_ser[N]/n_ser[N]**2), g3=float(n3_ser[N]/n_ser[N]**3))
    rec['series_n'] = [float(v) for v in n_ser]; rec['series_n2'] = [float(v) for v in n2_ser]; rec['series_n3'] = [float(v) for v in n3_ser]
    rec['order2_check'] = dict(n=float(np.sum(s1).real), n2=float(np.sum(s2).real), n3=float(np.sum(s3).real))
    json.dump({f"{U:.3f}": rec}, open(os.path.join(HERE, f'diag_U{U:.3f}.json'), 'w'), indent=1)
    print(f"diagrams U={U}: N2 g3={rec['2']['g3']:.6f} N4 g3={rec['4']['g3']:.6f} t4={rec['t_order4_g2_g3']:.1f}s t2={rec['t_order2_g2_g3']:.1f}s", flush=True)
