"""Single-core runtime measurements for fig_ring(b).  usage: timing.py tensor2|tensor4|momentum|exact"""
import sys, json, os, time, numpy as np
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
sys.path.insert(0, os.path.join(_ROOT, 'engines'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
P = dict(kappa=1.0, Delta=-1.0, J=0.4, eta=1.0); U = 0.15
what = sys.argv[1]
env = {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "SLURM_CPUS_PER_TASK")}
out = dict(what=what, env=env, U=U, **P, times={}, load_start=list(os.getloadavg()), load_end=None)
def save():
    out['load_end'] = list(os.getloadavg())
    json.dump(out, open(f'timing_{what}.json', 'w'), indent=1)

if what.startswith("tensor"):
    import tensor_engine as TE
    Nord = int(what[-1])
    Ks = [4, 6, 8, 12, 16, 24, 32] if Nord == 2 else [4, 6, 8, 10, 12, 14]
    for K in Ks:
        t0 = time.perf_counter()
        ch = TE.Chain(K, P["kappa"], [P["Delta"]]*K, P["J"], P["eta"], U, ring=True)
        a0, ad0 = ch.site_a(0), ch.site_ad(0)
        t1 = time.perf_counter()
        n = np.cumsum(ch.series(TE.mul(ad0, a0), Nord)).real
        t2 = time.perf_counter()
        n2 = np.cumsum(ch.series(TE.mul(TE.mul(ad0, ad0), TE.mul(a0, a0)), Nord)).real
        t3 = time.perf_counter()
        out["times"][str(K)] = dict(setup=t1-t0, n=t2-t1, n2=t3-t2, g2_total=t3-t0, g2=float(n2[-1]/n[-1]**2), n_value=float(n[-1]))
        print(f"tensor N<={Nord} K={K}: setup {t1-t0:.3f}s  n {t2-t1:.3f}s  n2 {t3-t2:.3f}s  total {t3-t0:.3f}s  g2={n2[-1]/n[-1]**2:.6f}", flush=True)
        save()

elif what == "momentum":
    import multimode as MM
    from ring_trunc import series_pruned  # only the function; module-level code guarded below
    for K in [3, 4, 5, 6, 7, 8]:
        rec = {}
        r = MM.Ring(K=K, U=U, **P); b0, bd0 = r.site_ops(0)
        On = MM.mul(bd0, b0); Og = MM.mul(MM.mul(bd0, bd0), MM.mul(b0, b0))
        if K <= 6:
            t = time.perf_counter(); r.series(On, 2); rec["n_unpruned_N2"] = time.perf_counter()-t
        if K <= 5:
            t = time.perf_counter(); r.series(Og, 2); rec["g2_unpruned_N2"] = time.perf_counter()-t
        t = time.perf_counter(); series_pruned(r, On, 2); rec["n_pruned_N2"] = time.perf_counter()-t
        t = time.perf_counter(); series_pruned(r, Og, 2); series_pruned(r, On, 2); rec["g2_pruned_N2"] = time.perf_counter()-t
        if K <= 6:
            t = time.perf_counter(); series_pruned(r, Og, 4); series_pruned(r, On, 4); rec["g2_pruned_N4"] = time.perf_counter()-t
        out["times"][str(K)] = rec
        print(f"momentum K={K}: {json.dumps({k: round(v, 3) for k, v in rec.items()})}", flush=True)
        save()

elif what == "exact":
    from ring_ss_evolve import ring_steady_evolve
    from ring_exact import ring_exact
    K = int(sys.argv[2]); nmax = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    t = time.perf_counter(); e = ring_steady_evolve(K, U=U, nmax=nmax, **P); te = time.perf_counter()-t
    out["times"][str(K)] = dict(evolve=te, nmax=nmax, g2=float(e["n0n0"]/e["n0"]**2), n0=float(e["n0"]))
    print(f"exact evolve K={K} nmax={nmax}: {te:.2f}s g2={e['n0n0']/e['n0']**2:.6f}", flush=True)
    save()
    if K <= 3:
        t = time.perf_counter(); e2 = ring_exact(K, P["kappa"], P["Delta"], P["J"], P["eta"], U, nmax, want=("n0", "g2")); ts = time.perf_counter()-t
        out["times"][str(K)]["spsolve"] = ts; out["times"][str(K)]["g2_spsolve"] = float(e2["n0n0"]/e2["n0"]**2)
        print(f"exact spsolve K={K} nmax={nmax}: {ts:.2f}s g2={e2['n0n0']/e2['n0']**2:.6f}", flush=True)
        save()
