"""Part (C): labels visited by the engines versus K, and the analytic bound.
Degree bound: a Kerr vertex (quartic in the displaced modes) raises the degree of a label by at most 2 and lowers it by at most 2.
After s of N vertices a label can still reach the (0,0) label only if its degree is <= 2(N-s); it cannot exceed d0+2s,
d0 = degree of the observable.  Hence D_s = min(d0+2s, 2(N-s)) (s>=1), D_0 = d0, and the number of labels at step s is at most
the number of normal-ordered monomials of degree <= D_s in 2K operators, C(2K+D_s, D_s) ~ (2K)^D_s / D_s!.
On the translation-invariant ring the vertices conserve the lattice momentum of a label, so only the zero-momentum
sector (about 1/K of all labels) can reach (0,0).
usage: label_counts.py momentum|tensor"""
import sys, json, os, time, numpy as np
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
from math import comb
sys.path.insert(0, os.path.join(_ROOT, 'engines'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
P = dict(kappa=1.0, Delta=-1.0, J=0.4, eta=1.0); U = 0.15
what = sys.argv[1]
BUDGET = float(sys.argv[2]) if len(sys.argv) > 2 else 900.0

def Dsteps(d0, N):
    return [d0] + [min(d0+2*s, 2*(N-s)) for s in range(1, N+1)]

def bound_counts(K, D):
    """number of monomials in b_q^dag (momentum +q), b_q (momentum -q), q=0..K-1, of degree <= D: total and zero-momentum"""
    # dp[d][P] over 2K variables
    dp = np.zeros((D+1, K), dtype=object); dp[0, 0] = 1
    for var in range(2*K):
        mom = (var % K) if var < K else (-(var - K)) % K
        new = dp.copy()
        # unbounded multiplicity: iterate degree ascending
        for d in range(1, D+1):
            for Pm in range(K):
                new[d, Pm] += new[d-1, (Pm - mom) % K]
        dp = new
    tot = int(sum(dp[d, p] for d in range(D+1) for p in range(K))); zero = int(sum(dp[d, 0] for d in range(D+1)))
    assert tot == comb(2*K+D, D)
    return tot, zero

def netmom(k, K):
    return (sum(q*c for q, c in k[0]) - sum(q*c for q, c in k[1])) % K

out = dict(what=what, U=U, **P, rule="D_s = min(d0+2s, 2(N-s)); bound = C(2K+D_s, D_s)", runs={})
if what == "momentum":
    import multimode as MM
    for obs, d0 in (("n", 2), ("g2num", 4)):
        for N in (2, 3, 4):
            last = 0.0
            for K in [3, 4, 5, 6, 7, 8, 10, 12, 14, 16]:
                if last > BUDGET/20 and K > 8: break
                if last > BUDGET/4: break
                r = MM.Ring(K=K, U=U, **P); a0, ad0 = r.site_ops(0)
                O = MM.mul(ad0, a0) if obs == "n" else MM.mul(MM.mul(ad0, ad0), MM.mul(a0, a0))
                X = O; tot = [len(X)]; zero = [sum(1 for k in X if netmom(k, K) == 0)]; dmax = [max(sum(c for _, c in k[0])+sum(c for _, c in k[1]) for k in X)]
                t = time.time()
                for s in range(1, N+1):
                    X = r.V(r.G0(X)); lim = 2*(N-s)
                    X = {k: v for k, v in X.items() if sum(c for _, c in k[0])+sum(c for _, c in k[1]) <= lim and v != 0}
                    tot.append(len(X)); zero.append(sum(1 for k in X if netmom(k, K) == 0))
                    dmax.append(max([sum(c for _, c in k[0])+sum(c for _, c in k[1]) for k in X] or [0]))
                last = time.time()-t
                Ds = Dsteps(d0, N); b = [bound_counts(K, D) for D in Ds]
                out["runs"][f"{obs}_N{N}_K{K}"] = dict(obs=obs, d0=d0, N=N, K=K, labels=tot, labels_zero_momentum=zero, max_degree=dmax,
                                                       D_bound=Ds, bound_total=[x[0] for x in b], bound_zero_momentum=[x[1] for x in b], time=last)
                print(f"momentum {obs} N={N} K={K}: labels {tot} P=0 {zero} deg {dmax} | bound {[x[0] for x in b]} P=0 bound {[x[1] for x in b]} [{last:.1f}s]", flush=True)
                json.dump(out, open('label_counts_momentum.json', 'w'), indent=1)
else:
    import tensor_engine as TE
    for obs, d0 in (("n", 2), ("g2num", 4)):
        for N, Ks in ((2, [3, 4, 6, 8, 12, 16, 24]), (3, [3, 4, 6, 8, 10, 12]), (4, [3, 4, 6, 8, 10, 12, 14])):
            for K in Ks:
                ch = TE.Chain(K, P["kappa"], [P["Delta"]]*K, P["J"], P["eta"], U, ring=True)
                a0, ad0 = ch.site_a(0), ch.site_ad(0)
                O = TE.mul(ad0, a0) if obs == "n" else TE.mul(TE.mul(ad0, ad0), TE.mul(a0, a0))
                X = O; t = time.time()
                def stats(X):
                    ent = sum(int(np.asarray(T).size) for T in X.values())
                    nz = sum(int(np.count_nonzero(np.abs(np.asarray(T)) > 1e-14*max(1e-300, np.max(np.abs(np.asarray(T)))))) for T in X.values())
                    dist = sum(comb(K+p-1, p)*comb(K+q-1, q) for (p, q) in X)
                    return ent, nz, dist, sorted(X.keys())
                st = [stats(X)]
                for s in range(1, N+1):
                    X = TE.prune(ch.Vx(ch.G0(X), maxdeg=2*(N-s)), 2*(N-s)); st.append(stats(X))
                dt = time.time()-t
                out["runs"][f"{obs}_N{N}_K{K}"] = dict(obs=obs, d0=d0, N=N, K=K, dense_entries=[x[0] for x in st], nonzero_entries=[x[1] for x in st],
                                                       distinct_labels_in_blocks=[x[2] for x in st], blocks=[[list(b) for b in x[3]] for x in st],
                                                       D_bound=Dsteps(d0, N), time=dt)
                print(f"tensor {obs} N={N} K={K}: dense {[x[0] for x in st]} distinct {[x[2] for x in st]} [{dt:.2f}s]", flush=True)
                json.dump(out, open('label_counts_tensor.json', 'w'), indent=1)
