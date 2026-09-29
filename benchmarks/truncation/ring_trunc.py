"""K=3 Kerr ring: g2 and g3 partial sums to order 6 (momentum engine with degree pruning),
tensor-engine cross-check, exact reference (time evolution, cutoff 8), error of every partial sum,
and label-count instrumentation versus K for part (C)."""
import sys, json, os, time, numpy as np
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
sys.path.insert(0, os.path.join(_ROOT, 'engines'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import multimode as MM, tensor_engine as TE
from series import series_record
from ring_ss_evolve import ring_steady_evolve
PR = json.load(open(os.path.join(_ROOT, 'notebook', 'data', 'revision_params.json')))["ring3"]
P = dict(kappa=PR["kappa"], Delta=PR["Delta"], J=PR["J"], eta=PR["eta"])

def deg(k):
    return sum(c for _, c in k[0]) + sum(c for _, c in k[1])

def series_pruned(r, O, Nmax, drop=2, counts=None):
    """same as Ring.series but discarding, after vertex N, every label whose degree exceeds
    drop*(Nmax-N): a Kerr vertex lowers the degree by at most two, so those labels cannot
    return to the steady state.  counts (optional list) receives the number of labels kept."""
    vals = [O.get(((), ()), 0)]
    X = O
    if counts is not None: counts.append(len(X))
    for N in range(1, Nmax+1):
        X = r.V(r.G0(X))
        lim = drop*(Nmax-N)
        X = {k: v for k, v in X.items() if deg(k) <= lim and v != 0}
        if counts is not None: counts.append(len(X))
        vals.append(X.get(((), ()), 0))
    return np.array(vals)

def series_unpruned(r, O, Nmax, counts=None):
    vals = [O.get(((), ()), 0)]; X = O
    if counts is not None: counts.append(len(X))
    for N in range(1, Nmax+1):
        X = r.V(r.G0(X)); X = {k: v for k, v in X.items() if v != 0}
        if counts is not None: counts.append(len(X))
        vals.append(X.get(((), ()), 0))
    return np.array(vals)

if __name__ == '__main__':
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    out = dict(parameters=dict(K=3, **P), Nmax=NMAX, exact={}, results={}, fig_ring_a={}, label_counts={})
    # ---------------------------------------------------------------- exact reference (g3 needs n3)
    Us_fig11 = PR["U_values_fig11"]
    Us_fig10 = [round(x, 3) for x in np.linspace(0, 0.3, 13)]
    cached = json.load(open(os.path.join(_ROOT, 'notebook', 'data', 'ring_exactU.json')))
    for U in sorted(set(Us_fig11)):
        t = time.time()
        e8 = ring_steady_evolve(3, U=U, nmax=8, T=80.0, steps=8, **P)
        e7 = ring_steady_evolve(3, U=U, nmax=7, T=80.0, steps=8, **P)
        g2_8 = e8["n0n0"]/e8["n0"]**2; g3_8 = e8["n3"]/e8["n0"]**3
        g2_7 = e7["n0n0"]/e7["n0"]**2; g3_7 = e7["n3"]/e7["n0"]**3
        out["exact"][f"{U:.3f}"] = dict(n0=e8["n0"], n0n0=e8["n0n0"], n3=e8["n3"], g2=g2_8, g3=g3_8,
                                        cutoff_change_g2=abs(g2_8-g2_7), cutoff_change_g3=abs(g3_8-g3_7),
                                        cached_n0n0=cached[f"{U:.3f}"]["n0n0"], cached_diff_n0n0=abs(cached[f"{U:.3f}"]["n0n0"]-e8["n0n0"]))
        print(f"exact U={U}: g2={g2_8:.6f} (d7-8 {abs(g2_8-g2_7):.1e}) g3={g3_8:.6f} (d7-8 {abs(g3_8-g3_7):.1e}) "
              f"cached n0n0 diff {abs(cached[f'{U:.3f}']['n0n0']-e8['n0n0']):.1e}  [{time.time()-t:.0f}s]", flush=True)

    # ---------------------------------------------------------------- diagrams K=3 to order NMAX
    def ring_series(U, Nmax):
        r = MM.Ring(K=3, U=U, **P); a0, ad0 = r.site_ops(0)
        n = series_pruned(r, MM.mul(ad0, a0), Nmax)
        n2 = series_pruned(r, MM.mul(MM.mul(ad0, ad0), MM.mul(a0, a0)), Nmax)
        n3 = series_pruned(r, MM.mul(MM.mul(MM.mul(ad0, ad0), ad0), MM.mul(MM.mul(a0, a0), a0)), Nmax)
        return np.cumsum(n).real, np.cumsum(n2).real, np.cumsum(n3).real

    # cross-checks at U=0.15: pruned vs unpruned momentum engine (order 4) and vs tensor engine (order NMAX)
    r = MM.Ring(K=3, U=0.15, **P); a0, ad0 = r.site_ops(0)
    O2 = MM.mul(MM.mul(ad0, ad0), MM.mul(a0, a0))
    t = time.time(); sp_ = series_pruned(r, O2, 4); tp = time.time()-t
    t = time.time(); su_ = series_unpruned(r, O2, 4); tu = time.time()-t
    ch = TE.Chain(3, P["kappa"], [P["Delta"]]*3, P["J"], P["eta"], 0.15, ring=True)
    a3, ad3 = ch.site_a(0), ch.site_ad(0)
    t = time.time(); st_ = ch.series(TE.mul(TE.mul(ad3, ad3), TE.mul(a3, a3)), NMAX); tt = time.time()-t
    s6 = series_pruned(r, O2, NMAX)
    out["crosscheck"] = dict(pruned_vs_unpruned_maxdiff=float(np.max(np.abs(sp_-su_))), t_pruned=tp, t_unpruned=tu,
                             momentum_vs_tensor_maxdiff=float(np.max(np.abs(s6-st_))), t_tensor=tt, orders=NMAX)
    print("crosscheck:", out["crosscheck"], flush=True)

    for U in Us_fig11:
        t = time.time()
        n, n2, n3 = ring_series(U, NMAX)
        g2 = n2/n**2; g3 = n3/n**3
        ex = out["exact"][f"{U:.3f}"]
        out["results"][f"{U:.3f}"] = dict(U=U, n=series_record(n, exact=ex["n0"]),
                                          g2=series_record(g2, exact=ex["g2"]),
                                          g3=series_record(g3, exact=ex["g3"]))
        rg2 = out["results"][f"{U:.3f}"]["g2"]; rg3 = out["results"][f"{U:.3f}"]["g3"]
        print(f"U={U}: g2 N={rg2['N']} {rg2['value']:.5f} (err {rg2['error']:.1e}, exact {ex['g2']:.5f}) | "
              f"g3 N={rg3['N']} {rg3['value']:.5f} (err {rg3['error']:.1e}, exact {ex['g3']:.5f}) [{time.time()-t:.0f}s]", flush=True)

    # fig_ring panel (a): g2 on the 13-point U grid, orders to NMAX, exact from cache (n0, n0n0)
    for U in Us_fig10:
        n, n2, n3 = ring_series(U, NMAX)
        g2 = n2/n**2
        exg2 = cached[f"{U:.3f}"]["n0n0"]/cached[f"{U:.3f}"]["n0"]**2
        out["fig_ring_a"][f"{U:.3f}"] = dict(U=U, g2=series_record(g2, exact=exg2), n=series_record(n, exact=cached[f"{U:.3f}"]["n0"]))
        rg = out["fig_ring_a"][f"{U:.3f}"]["g2"]
        print(f"fig10a U={U}: g2 N={rg['N']} {rg['value']:.5f} err {rg['error']:.1e}", flush=True)
    json.dump(out, open('ring_truncation.json', 'w'), indent=1)

    # ---------------------------------------------------------------- label counts versus K (part C)
    def dense_count(ch, O, Nmax, drop=2):
        """tensor engine: number of tensor entries formed per step (dense labels incl. permutation copies)
        and number of distinct labels (multisets) they represent"""
        from math import comb
        X = O; entries = [sum(np.asarray(T).size for T in X.values())]
        distinct = [sum(comb(ch.K+p-1, p)*comb(ch.K+q-1, q) for (p, q) in X)]
        for N in range(1, Nmax+1):
            X = TE.prune(ch.Vx(ch.G0(X), maxdeg=drop*(Nmax-N)), drop*(Nmax-N))
            entries.append(sum(np.asarray(T).size for T in X.values()))
            distinct.append(sum(comb(ch.K+p-1, p)*comb(ch.K+q-1, q) for (p, q) in X))
        return entries, distinct
    LC = {}
    for Nord, Ks in [(2, [3, 4, 5, 6, 7, 8, 10, 12]), (4, [3, 4, 5, 6, 7, 8])]:
        for K in Ks:
            rK = MM.Ring(K=K, U=0.15, **P); aK, adK = rK.site_ops(0)
            OK = MM.mul(MM.mul(adK, adK), MM.mul(aK, aK))     # on-site g2 observable
            cp, cu = [], []
            t = time.time(); series_pruned(rK, OK, Nord, counts=cp); tp = time.time()-t
            if (Nord == 2 and K <= 6) or (Nord == 4 and K <= 4):
                t = time.time(); series_unpruned(rK, OK, Nord, counts=cu); tu = time.time()-t
            else:
                tu = None
            chK = TE.Chain(K, P["kappa"], [P["Delta"]]*K, P["J"], P["eta"], 0.15, ring=True)
            aT, adT = chK.site_a(0), chK.site_ad(0)
            ent, dis = dense_count(chK, TE.mul(TE.mul(adT, adT), TE.mul(aT, aT)), Nord)
            LC[f"N{Nord}_K{K}"] = dict(N=Nord, K=K, momentum_pruned_labels=cp, momentum_unpruned_labels=cu,
                                       tensor_entries=ent, tensor_distinct_labels=dis, t_pruned=tp, t_unpruned=tu)
            print(f"labels N={Nord} K={K}: pruned {cp} (max {max(cp)}), unpruned {cu if cu else '-'}, tensor entries {ent}", flush=True)
    out["label_counts"] = LC
    # fitted exponents of the maximal label count versus K
    def fit(Ks, ys):
        Ks, ys = np.array(Ks, float), np.array(ys, float); m = ys > 0
        return float(np.polyfit(np.log(Ks[m]), np.log(ys[m]), 1)[0])
    fits = {}
    for Nord in (2, 4):
        keys = [k for k in LC if LC[k]["N"] == Nord]
        Ks = [LC[k]["K"] for k in keys]
        fits[f"N{Nord}"] = dict(K=Ks,
            pruned_max_labels=fit(Ks, [max(LC[k]["momentum_pruned_labels"]) for k in keys]),
            pruned_labels_after_first_vertex=fit(Ks, [LC[k]["momentum_pruned_labels"][1] for k in keys]),
            tensor_entries_max=fit(Ks, [max(LC[k]["tensor_entries"]) for k in keys]),
            tensor_distinct_max=fit(Ks, [max(LC[k]["tensor_distinct_labels"]) for k in keys]),
            pruned_total_labels=fit(Ks, [sum(LC[k]["momentum_pruned_labels"]) for k in keys]),
            pruned_time=fit(Ks, [LC[k]["t_pruned"] for k in keys]))
        uk = [k for k in keys if LC[k]["momentum_unpruned_labels"]]
        if len(uk) >= 3:
            fits[f"N{Nord}"]["unpruned_max_labels"] = fit([LC[k]["K"] for k in uk], [max(LC[k]["momentum_unpruned_labels"]) for k in uk])
            fits[f"N{Nord}"]["unpruned_time"] = fit([LC[k]["K"] for k in uk], [LC[k]["t_unpruned"] for k in uk])
    out["label_count_fits"] = fits
    print(json.dumps(fits, indent=1))
    json.dump(out, open('ring_truncation.json', 'w'), indent=1)
