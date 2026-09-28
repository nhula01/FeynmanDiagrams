"""Single-core runtimes for fig_ring(b).  usage: timing_v2.py tensor2|tensor4|momentum|exact K nmax
Writes timing2_<tag>.json with the load average before/after every point."""
import sys, json, os, time, numpy as np
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
sys.path.insert(0, os.path.join(_ROOT, 'engines'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
P = dict(kappa=1.0, Delta=-1.0, J=0.4, eta=1.0); U = 0.15
what = sys.argv[1]
tag = what if what != "exact" else f"exact_K{sys.argv[2]}_nmax{sys.argv[3]}"
env = {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "SLURM_CPUS_PER_TASK", "SLURM_JOB_ID")}
aff = sorted(os.sched_getaffinity(0))
out = dict(what=what, env=env, cpu_affinity=aff, U=U, **P, times={}, load_start=list(os.getloadavg()))
def save():
    out['load_end'] = list(os.getloadavg()); json.dump(out, open(f'timing2_{tag}.json', 'w'), indent=1)
def rec(key, d):
    d["load1"] = os.getloadavg()[0]; out["times"][key] = d; save()

if what.startswith("tensor"):
    import tensor_engine as TE
    Nord = int(what[-1])
    Ks = [4, 6, 8, 12, 16, 24, 32] if Nord == 2 else [4, 6, 8, 10, 12, 14]
    for K in Ks:
        best = None
        for rep in range(3 if K <= 16 else 1):
            t0 = time.perf_counter()
            ch = TE.Chain(K, P["kappa"], [P["Delta"]]*K, P["J"], P["eta"], U, ring=True)
            a0, ad0 = ch.site_a(0), ch.site_ad(0)
            t1 = time.perf_counter()
            n = np.cumsum(ch.series(TE.mul(ad0, a0), Nord)).real
            t2 = time.perf_counter()
            n2 = np.cumsum(ch.series(TE.mul(TE.mul(ad0, ad0), TE.mul(a0, a0)), Nord)).real
            t3 = time.perf_counter()
            d = dict(setup=t1-t0, n=t2-t1, n2=t3-t2, g2_total=t3-t0, g2=float(n2[-1]/n[-1]**2), n_value=float(n[-1]))
            if best is None or d["g2_total"] < best["g2_total"]: best = d
        best["repeats"] = 3 if K <= 16 else 1
        rec(str(K), best)
        print(f"tensor N<={Nord} K={K}: setup {best['setup']:.3f}s n {best['n']:.3f}s n2 {best['n2']:.3f}s total {best['g2_total']:.3f}s g2={best['g2']:.6f}", flush=True)

elif what == "momentum":
    import multimode as MM
    from ring_trunc import series_pruned
    for K in [3, 4, 5, 6, 7, 8]:
        d = {}
        r = MM.Ring(K=K, U=U, **P); b0, bd0 = r.site_ops(0)
        On = MM.mul(bd0, b0); Og = MM.mul(MM.mul(bd0, bd0), MM.mul(b0, b0))
        if K <= 6:   # the notebook curve: <n> series to N<=2, unpruned (ring_timing.json "diagrams")
            t = time.perf_counter(); r2 = MM.Ring(K=K, U=U, **P); c0, cd0 = r2.site_ops(0); r2.series(MM.mul(cd0, c0), 2); d["n_unpruned_N2_incl_setup"] = time.perf_counter()-t
        t = time.perf_counter(); series_pruned(r, On, 2); d["n_pruned_N2"] = time.perf_counter()-t
        t = time.perf_counter(); series_pruned(r, Og, 2); series_pruned(r, On, 2); d["g2_pruned_N2"] = time.perf_counter()-t
        if K <= 6:
            t = time.perf_counter(); series_pruned(r, Og, 4); series_pruned(r, On, 4); d["g2_pruned_N4"] = time.perf_counter()-t
        rec(str(K), d)
        print(f"momentum K={K}: {json.dumps({k: round(v, 3) for k, v in d.items()})}", flush=True)

elif what == "exact":
    from ring_ss_evolve import ring_steady_evolve
    from ring_exact import ring_exact
    K = int(sys.argv[2]); nmax = int(sys.argv[3])
    t = time.perf_counter(); e = ring_steady_evolve(K, U=U, nmax=nmax, **P); te = time.perf_counter()-t
    d = dict(evolve=te, nmax_levels=nmax, g2=float(e["n0n0"]/e["n0"]**2), n0=float(e["n0"]))
    rec(str(K), d); print(f"exact evolve K={K} nmax={nmax}: {te:.2f}s g2={d['g2']:.6f}", flush=True)
    if K <= 3 and nmax <= 6:
        t = time.perf_counter(); e2 = ring_exact(K, P["kappa"], P["Delta"], P["J"], P["eta"], U, nmax, want=("n0", "g2")); ts = time.perf_counter()-t
        d["spsolve"] = ts; d["g2_spsolve"] = float(e2["n0n0"]/e2["n0"]**2); rec(str(K), d)
        print(f"exact spsolve K={K} nmax={nmax}: {ts:.2f}s", flush=True)
