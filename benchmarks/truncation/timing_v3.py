"""Second single-core timing pass: wall and CPU time, tensor engine as in the paper ('te') and with the commutator
restricted to terms with at least one contraction ('te_nz'; the uncontracted terms cancel identically), momentum engine
original and zero-contraction-free.  usage: timing_v3.py tensor|momentum"""
import sys, json, os, time, numpy as np
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
from math import comb, factorial
sys.path.insert(0, os.path.join(_ROOT, 'engines'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tensor_engine as TE
P = dict(kappa=1.0, Delta=-1.0, J=0.4, eta=1.0); U = 0.15
what = sys.argv[1]
out = dict(what=what, env={k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "SLURM_CPUS_PER_TASK", "SLURM_JOB_ID")},
           cpu_affinity=sorted(os.sched_getaffinity(0)), times={}, load_start=list(os.getloadavg()))
def save():
    out['load_end'] = list(os.getloadavg()); json.dump(out, open(f'timing3_{what}.json', 'w'), indent=1)

def mul_nz(A, B, maxdeg=None):
    o = {}
    for (p1, q1), T1 in A.items():
        for (p2, q2), T2 in B.items():
            for k in range(1, min(q1, p2)+1):
                if maxdeg is not None and (p1+p2-k)+(q1+q2-k) > maxdeg: continue
                w = factorial(k)*comb(q1, k)*comb(p2, k)
                R = np.tensordot(T1, T2, axes=(list(range(p1+q1-k, p1+q1)), list(range(k))))
                n1c, n1a, n2c = p1, q1-k, p2-k
                if R.ndim:
                    order = list(range(n1c)) + list(range(n1c+n1a, n1c+n1a+n2c)) + list(range(n1c, n1c+n1a)) + list(range(n1c+n1a+n2c, R.ndim))
                    R = np.transpose(R, order)
                key = (p1+p2-k, q1+q2-k); o[key] = o.get(key, 0) + w*R
    return {k: TE.sym(np.asarray(v), *k) for k, v in o.items()}
class ChainNZ(TE.Chain):
    def Vx(self, A, maxdeg=None):
        return TE.scale(1j, TE.add(mul_nz(self.Hint, A, maxdeg), TE.scale(-1, mul_nz(A, self.Hint, maxdeg))))

def timed(f):
    w0, c0 = time.perf_counter(), time.process_time(); r = f(); return r, time.perf_counter()-w0, time.process_time()-c0

if what == "tensor":
    for Nord, Ks in ((2, [4, 6, 8, 12, 16, 24, 32]), (4, [4, 6, 8, 10, 12, 14])):
        for K in Ks:
            for name, cls in (("te", TE.Chain), ("te_nz", ChainNZ)):
                best = None
                for rep in range(3 if K <= 12 else 1):
                    def run():
                        ch = cls(K, P["kappa"], [P["Delta"]]*K, P["J"], P["eta"], U, ring=True); a0, ad0 = ch.site_a(0), ch.site_ad(0)
                        n = np.cumsum(ch.series(TE.mul(ad0, a0), Nord)).real
                        n2 = np.cumsum(ch.series(TE.mul(TE.mul(ad0, ad0), TE.mul(a0, a0)), Nord)).real
                        return float(n2[-1]/n[-1]**2)
                    g2, tw, tc = timed(run)
                    if best is None or tw < best["wall"]: best = dict(wall=tw, cpu=tc, g2=g2, load1=os.getloadavg()[0])
                out["times"][f"{name}_N{Nord}_K{K}"] = best; save()
                print(f"{name} N<={Nord} K={K}: wall {best['wall']:.3f}s cpu {best['cpu']:.3f}s g2={best['g2']:.8f} load {best['load1']:.0f}", flush=True)
elif what == "momentum":
    import multimode as MM
    sys.argv = ['x']; import ring_pool as RPm
    for K in [3, 4, 5, 6, 7, 8]:
        r = MM.Ring(K=K, U=U, **P); b0, bd0 = r.site_ops(0)
        On = MM.mul(bd0, b0); Og = MM.mul(MM.mul(bd0, bd0), MM.mul(b0, b0))
        d = {}
        if K <= 6:
            _, d["n_unpruned_N2_wall"], d["n_unpruned_N2_cpu"] = timed(lambda: MM.Ring(K=K, U=U, **P).series(On, 2))
        _, d["g2_pruned_N2_wall"], d["g2_pruned_N2_cpu"] = timed(lambda: (RPm.series_pruned(r, Og, 2), RPm.series_pruned(r, On, 2)))
        _, d["g2_nz_N2_wall"], d["g2_nz_N2_cpu"] = timed(lambda: (RPm.series_nz(r, Og, 2), RPm.series_nz(r, On, 2)))
        if K <= 6:
            _, d["g2_nz_N4_wall"], d["g2_nz_N4_cpu"] = timed(lambda: (RPm.series_nz(r, Og, 4), RPm.series_nz(r, On, 4)))
        d["load1"] = os.getloadavg()[0]; out["times"][str(K)] = d; save()
        print(f"momentum K={K}: {json.dumps({k: round(v, 3) for k, v in d.items()})}", flush=True)
