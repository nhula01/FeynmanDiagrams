#!/usr/bin/env python3
"""Check every number quoted in the paper against the result files in this repository.

Two kinds of check:

  1. Table V of the main text is compared row by row with
     benchmarks/truncation/table5_fragment.tex, which assemble_truncation.py
     writes from table5_corrected.json.
  2. Every number quoted in the prose (Secs. V-VII, the abstract, and the
     quantum-jump section of the Supplemental Material) is recomputed from the
     JSON result files and compared with the value printed in the paper, at the
     precision the paper uses.

Run from anywhere:  python3 notebooks/verify_numbers.py [--report FILE]
Exit status is non-zero if any check fails.  Checks that need
engines/numbers.json (written by `python3 make_figures.py kerr cubic
emitter finite atomcavity` in engines/) are skipped if that file is absent.
"""
import json, os, re, sys, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, ".."))           # the repository root
ROOT = CODE

def _paper_dir():
    for d in (os.environ.get("PAPER_DIR"), os.path.join(CODE, "paper"), os.path.join(CODE, "..", "paper")):
        if d and os.path.exists(os.path.join(d, "main.tex")): return os.path.abspath(d)
    return None
PAPER = _paper_dir()
SNAP = {"tab:fcs": "table_V_main.tex", "tab:closures": "table_SI_supplement.tex"}

def paper_table(label):
    """LaTeX source containing a table of the paper: the paper sources if present, else the snapshot in paper_reference/."""
    if PAPER is not None:
        for f in ("main.tex", "supplement.tex"):
            t = open(os.path.join(PAPER, f), encoding="utf-8").read()
            if "\\label{%s}" % label in t: return t
    return open(os.path.join(HERE, "paper_reference", SNAP[label]), encoding="utf-8").read()
J = lambda *p: json.load(open(os.path.join(CODE, *p)))

rows = []          # (section, quantity, paper, code, ok)
def check(section, what, paper, code, tol=None, rel=None, digits=None):
    """paper: value as printed; code: value from the result files.
    Agreement means |paper-code| <= tol, or relative rel, or equal after rounding to `digits`."""
    if code is None:
        rows.append((section, what, paper, None, None)); return
    if digits is not None:
        ok = round(float(code), digits) == round(float(paper), digits) or abs(float(code) - float(paper)) <= 0.5001 * 10 ** (-digits)
    elif rel is not None:
        ok = abs(float(code) - float(paper)) <= rel * abs(float(paper))
    else:
        ok = abs(float(code) - float(paper)) <= (tol if tol is not None else 0.5 * 10 ** (-_decimals(paper)))
    rows.append((section, what, paper, code, bool(ok)))

def _decimals(x):
    s = repr(float(x))
    return len(s.split(".")[1]) if "." in s and "e" not in s else 0

# ----------------------------------------------------------------------------- Table V
def table5():
    main = paper_table("tab:fcs")
    frag = open(os.path.join(CODE, "benchmarks", "truncation", "table5_fragment.tex"), encoding="utf-8").read()
    def body(txt):
        m = re.search(r"\\label\{tab:fcs\}.*?\\begin\{tabular\}.*?\\hline(.*?)\\end\{tabular\}", txt, re.S)
        if m is None:
            m = re.search(r"\\hline(.*?)\\end\{tabular\}", txt, re.S)
        lines = [l for l in m.group(1).split("\\\\") if "&" in l]
        out = []
        for l in lines:
            cells = [c.strip() for c in l.split("&")]
            nums = re.findall(r"-?\d+\.\d+(?:\\times10\^\{-?\d+\})?", "&".join(cells[2:]))
            out.append((re.sub(r"\s+", " ", cells[1]), nums))
        return out
    A, B = body(main), body(frag)
    ok_all = len(A) == len(B)
    for (la, na), (lb, nb) in zip(A, B):
        ok = na == nb
        ok_all &= ok
        rows.append(("Table V", f"{la}: {' '.join(na)}", " ".join(na), " ".join(nb), ok))
    rows.append(("Table V", "row count", len(A), len(B), len(A) == len(B)))

# ----------------------------------------------------------------------------- Sec. V.A  Kerr cavity
def kerr():
    S = "Sec. V.A / Fig. 5"
    T = J("benchmarks", "truncation", "truncation_results.json")["kerr_cavity"]["results"]
    check(S, "U=0.05: delta (successive difference) at N*=16", 1.6e-11, T["0.05"]["delta"], rel=0.05)
    check(S, "U=0.05: true error at N*=16", 4.3e-12, T["0.05"]["true_error"], rel=0.05)
    check(S, "U=0.10: delta at N*=16", 1.0e-6, T["0.10"]["delta"], rel=0.05)
    check(S, "U=0.10: true error at N*=16", 1.8e-7, T["0.10"]["true_error"], rel=0.05)
    check(S, "U=0.20: N*", 8, T["0.20"]["N_star"], tol=0)
    check(S, "U=0.20: delta", 0.041, T["0.20"]["delta"], digits=3)
    check(S, "U=0.20: true error at N*", 0.061, T["0.20"]["true_error"], digits=3)
    check(S, "U=0.20: true optimum order", 9, T["0.20"]["N_true_optimum"], tol=0)
    check(S, "U=0.20: minimum error (at N=9)", 0.019, min(T["0.20"]["true_error_all_orders"]), digits=3)
    check(S, "U=0.20: minimum relative error, 1.3%", 1.3, 100 * min(T["0.20"]["true_error_all_orders"]) / T["0.20"]["exact"], digits=1)
    check(S, "U=0.30: N*", 2, T["0.30"]["N_star"], tol=0)
    check(S, "U=0.30: delta", 0.25, T["0.30"]["delta"], digits=2)
    check(S, "U=0.30: true error at N*", 1.08, T["0.30"]["true_error"], digits=2)
    check(S, "U=0.30: minimum error", 0.16, min(T["0.30"]["true_error_all_orders"]), digits=2)
    check(S, "U=0.30: order of minimum error", 5, int(np.argmin(T["0.30"]["true_error_all_orders"])), tol=0)
    check(S, "U=0.20: exact <a^dag a>", 1.485, T["0.20"]["exact"], digits=3)
    check(S, "U=0.30: exact <a^dag a>", 2.43, T["0.30"]["exact"], digits=2)
    nfile = os.path.join(CODE, "engines", "numbers.json")
    if os.path.exists(nfile):
        n = json.load(open(nfile))
        mf = n.get("kerr_exact_vs_mf", {})
        check(S, "U=0.20: mean-field <a^dag a>", 1.215, mf.get("0.20", [None, None])[1], digits=3)
        check(S, "U=0.30: mean-field <a^dag a>", 3.76, mf.get("0.30", [None, None])[1], digits=2)
        check(S, "U=0.30: mean-field error", 1.3, abs(mf["0.30"][1] - mf["0.30"][0]) if "0.30" in mf else None, digits=1)
        c = n.get("kerr_coeffs")
        if c and len(c) > 10:
            check(S, "root test |c_N|^(-1/N) at N=10", 0.27, abs(c[10]) ** (-1 / 10), digits=2)
        if c and len(c) > 40:
            check(S, "root test |c_N|^(-1/N) at N=40", 0.15, abs(c[40]) ** (-1 / 40), digits=2)
        else:
            rows.append((S, "root test at N=40 (needs Nmax>=40; make_figures.py computes 16)", 0.15, None, None))
    else:
        rows.append((S, "mean-field values 1.215 / 3.76 and root test (run engines/make_figures.py first)", "", None, None))

# ----------------------------------------------------------------------------- Sec. VI  ring
def ring():
    S = "Sec. VI / Fig. 8"
    T = J("benchmarks", "truncation", "truncation_results.json")
    R = T["ring_K3"]["results"]
    def rel_err(U, obs, N):
        r = R[U][obs]; return abs(r["partial_sums"][N] - r["exact"]) / r["exact"]
    check(S, "n: second-order relative error at U=0.05, 0.01%", 0.01, 100 * rel_err("0.050", "n", 2), digits=2)
    check(S, "n: second-order relative error at U=0.15, 0.3%", 0.3, 100 * rel_err("0.150", "n", 2), digits=1)
    check(S, "n: second-order relative error at U=0.30, 3%", 3, 100 * rel_err("0.300", "n", 2), digits=0)
    ratios = [rel_err(U, "n", 1) / rel_err(U, "n", 2) for U in ("0.050", "0.100", "0.150", "0.200", "0.250", "0.300")]
    check(S, "first/second-order error ratio on n, U=0.05-0.3: min (two)", 2, min(ratios), digits=0)
    check(S, "first/second-order error ratio on n, U=0.05-0.3: max (twelve)", 12, max(ratios), digits=0)
    r15 = R["0.150"]
    check(S, "U=0.15: exact n_0", 0.3029, r15["n"]["exact"], digits=4)
    check(S, "U=0.15: exact g2", 1.151, r15["g2"]["exact"], digits=3)
    check(S, "U=0.15: n_0 through second order", 0.3019, r15["n"]["partial_sums"][2], digits=4)
    check(S, "U=0.15: g2 through second order", 1.140, r15["g2"]["partial_sums"][2], digits=3)
    check(S, "U=0.15: g2 through fourth order", 1.150, r15["g2"]["partial_sums"][4], digits=3)
    check(S, "U=0.15: g2 N*", 8, r15["g2"]["N_star"], tol=0)
    check(S, "U=0.15: g2 optimally truncated", 1.1509, r15["g2"]["value"], digits=4)
    check(S, "U=0.15: g2 delta, 1e-5", 1e-5, r15["g2"]["delta"], rel=0.5)
    check(S, "U=0.15: g2 true error, 5e-5", 5e-5, r15["g2"]["true_error"], rel=0.2)
    # delta versus true error
    def ratio(U, obs): r = R[U][obs]; return r["true_error"] / r["delta"]
    worst_weak = max(ratio(U, "g2") for U in ("0.025", "0.050", "0.075", "0.100", "0.125"))
    check(S, "g2 error bars of Fig. 8(a): true/delta <= 2 for U<=0.125 (max)", 2, worst_weak, tol=2 - worst_weak if worst_weak <= 2 else -1)
    g3w = {U: ratio(U, "g3") for U in ("0.025", "0.050", "0.075", "0.100", "0.125")}
    rows.append((S, "NOTE g3 true/delta at U<=0.125: " + ", ".join(f"U={U}: {v:.1f}" for U, v in g3w.items()) + " (the factor-two statement holds for g2; for g3 it fails at U=0.125)", "", max(g3w.values()), None))
    check(S, "U=0.15: true/delta on g2, factor 4", 4, ratio("0.150", "g2"), digits=0)
    check(S, "U=0.30: true/delta on g2, factor 12", 12, ratio("0.300", "g2"), digits=0)
    check(S, "U=0.30: true/delta on g3, factor 27", 27, ratio("0.300", "g3"), digits=0)
    mx = max(ratio(U, o) for U in ("0.150", "0.175", "0.200", "0.225", "0.250", "0.275", "0.300") for o in ("g2", "g3"))
    mn = min(ratio(U, o) for U in ("0.150", "0.175", "0.200", "0.225", "0.250", "0.275", "0.300") for o in ("g2", "g3"))
    check(S, "underestimate by 4 to 27 for U>=0.15: min", 4, mn, digits=0)
    check(S, "underestimate by 4 to 27 for U>=0.15: max", 27, mx, digits=0)
    # timings and scaling (pass 3, single core at fixed load)
    P = J("benchmarks", "truncation", "summary_numbers.json")["timing_pass3"]
    F = T["timing_fits"]
    check(S, "exact steady state, 5 levels: growth per site (45)", 45, F["pass3:exact_levels5"]["per_site_factor"], digits=0)
    check(S, "exact steady state, 5 levels, K=4: 328 s", 328, P["exact_levels5"]["t"][P["exact_levels5"]["K"].index(4)], digits=0)
    check(S, "momentum engine N<=2, K=6: 2.4 s", 2.4, P["momentum3_n_unpruned_N2"]["t"][P["momentum3_n_unpruned_N2"]["K"].index(6)], digits=1)
    check(S, "momentum engine N<=2: exponent 5.7", 5.7, F["pass3:momentum3_n_unpruned_N2"]["exponent_last3"], digits=1)
    check(S, "tensor engine N<=2, K=32: 2.1 s", 2.1, P["te_N2"]["t"][P["te_N2"]["K"].index(32)], digits=1)
    check(S, "tensor engine N<=2: exponent above K=16, 5.9", 5.9, F["pass3:te_N2"]["exponent_last3"], digits=1)
    check(S, "tensor engine N<=4, K=14: 30 s", 30, P["te_N4"]["t"][P["te_N4"]["K"].index(14)], digits=0)
    check(S, "tensor engine N<=4: exponent 6.0", 6.0, F["pass3:te_N4"]["exponent_last3"], digits=1)
    check(S, "tensor engine N<=4, contracted commutator, K=14: 6.1 s", 6.1, P["te_nz_N4"]["t"][P["te_nz_N4"]["K"].index(14)], digits=1)
    # label counts: bound C(2K+Dmax, Dmax) grows as K^4 at N=2,3 and K^6 at N=4,5 for d0=4
    from math import comb
    def Dmax(N, d0=4): return 2 * math.floor(N / 2 + d0 / 4)
    for N, exp in ((2, 4), (3, 4), (4, 6), (5, 6)):
        check(S, f"label bound C(2K+Dmax,Dmax) ~ K^Dmax: Dmax at N={N} for d0=4 (K^{exp})", exp, Dmax(N), tol=0)
    # measured label counts of the tensor engine (Fig. 8(b) inset / SM)
    runs = T["label_counts"]["tensor"]["runs"]
    for N, exp in ((2, 4), (3, 4), (4, 6)):
        pts = [(v["K"], max(v["distinct_labels_in_blocks"])) for k, v in runs.items() if k.startswith(f"g2num_N{N}_") and v["K"] >= 8]
        if len(pts) >= 3:
            Ks, Ls = np.array([p[0] for p in pts], float), np.array([p[1] for p in pts], float)
            slope = np.polyfit(np.log(Ks), np.log(Ls), 1)[0]
            rows.append((S, f"measured distinct-label count of the tensor engine, g2, N={N}, K>=8: local exponent {slope:.2f} (bound: K^{exp})", exp, slope, None))

# ----------------------------------------------------------------------------- Sec. VI.B  closures
def closures():
    S = "Sec. VI.B / Fig. 9 / SM Table SI"
    C = J("benchmarks", "cumulant_closures", "cumulant_closure_results.json")
    check(S, "closure variables at order 2", 27, C["variables_per_order"]["2"], tol=0)
    check(S, "closure variables at order 3", 83, C["variables_per_order"]["3"], tol=0)
    check(S, "closure variables at order 4", 209, C["variables_per_order"]["4"], tol=0)
    tab = {row["U"]: row for row in C["table_fig11"]}
    t = lambda k: np.median([tab[U][k]["t"] for U in tab])
    check(S, "closure order 2: 0.13 s per evaluation", 0.13, t("closure2"), digits=2)
    check(S, "closure order 3: 0.6 s", 0.6, t("closure3"), digits=1)
    check(S, "closure order 4: 5 s", 5, t("closure4"), digits=0)
    DC = C["diagram_cost"]
    check(S, "diagrams N<=2, pruned symbolic recursion: 0.21 s per coupling", 0.21, DC["pruned_per_U"]["2"], digits=2)
    check(S, "diagrams N<=4, pruned symbolic recursion: 1.1 s per coupling", 1.1, DC["pruned_per_U"]["4"], digits=1)
    check(S, "sparse-matrix form, orders 0-2, all couplings: 0.03 s", 0.03, DC["sparse_all_U"]["2"], digits=2)
    check(S, "sparse-matrix form, orders 0-4, all couplings: 0.11 s", 0.11, DC["sparse_all_U"]["4"], digits=2)
    check(S, "sparse-matrix form, orders 0-10, all couplings: 4.2 s", 4.2, DC["sparse_all_U"]["10"], digits=1)
    check(S, "   ... below one closure-4 evaluation", 1, int(DC["sparse_all_U"]["10"] < t("closure4")), tol=0)
    check(S, "closures cost 0.6 to 5 times the pruned fourth-order diagrams: min", 0.6, t("closure3") / DC["pruned_per_U"]["4"], tol=0.15)
    check(S, "closures cost 0.6 to 5 times the pruned fourth-order diagrams: max", 5, t("closure4") / DC["pruned_per_U"]["4"], tol=0.8)
    check(S, "U=0.1: exact g3", 1.3126, tab[0.1]["exact_g3"], digits=4)
    check(S, "U=0.1: fourth-order diagrams reduce the g3 error of HFB by a factor 12", 12, tab[0.1]["closure2"]["err_g3"] / tab[0.1]["diagrams4"]["err_g3"], digits=0)
    check(S, "U=0.1: HFB g3 error 1.2%", 1.2, 100 * tab[0.1]["closure2"]["err_g3"], digits=1)
    check(S, "U=0.1: fourth-order diagram g3 error 0.10%", 0.10, 100 * tab[0.1]["diagrams4"]["err_g3"], digits=2)
    f = [tab[U]["closure2"]["err_g3"] / tab[U]["diagrams4"]["err_g3"] for U in (0.15, 0.2, 0.25, 0.3)]
    check(S, "U=0.15-0.3: HFB/diagrams-4 g3 error ratio, min 2.8", 2.8, min(f), digits=1)
    check(S, "U=0.15-0.3: HFB/diagrams-4 g3 error ratio, max 4.8", 4.8, max(f), digits=1)
    best = lambda U: min(tab[U]["closure3"]["err_g3"], tab[U]["closure4"]["err_g3"])
    check(S, "U=0.1: diagrams-4 g3 error / best closure, 5.5", 5.5, tab[0.1]["diagrams4"]["err_g3"] / best(0.1), digits=1)
    check(S, "U=0.3: diagrams-4 g3 error / best closure, 1.25", 1.25, tab[0.3]["diagrams4"]["err_g3"] / best(0.3), digits=2)

    O = C["optimal_truncation_diagrams"]["results"]
    check(S, "optimal truncation: N* in 8..10 (min)", 8, min(O[U]["g3"]["N_star"] for U in O if float(U) > 0), tol=0)
    check(S, "optimal truncation: N* in 8..10 (max)", 10, max(O[U]["g3"]["N_star"] for U in O if float(U) > 0), tol=0)
    check(S, "U=0.10: optimally truncated g3 error 1.8e-6", 1.8e-6, O["0.100"]["g3"]["err_vs_converged_exact"], rel=0.06)
    check(S, "U=0.15: optimally truncated g3 error 3.2e-4", 3.2e-4, O["0.150"]["g3"]["err_vs_converged_exact"], rel=0.05)
    check(S, "U=0.20: optimally truncated g3 error 4.9e-3", 4.9e-3, O["0.200"]["g3"]["err_vs_converged_exact"], rel=0.05)
    check(S, "U=0.10: better than best closure by 100", 100, O["0.100"]["g3"]["best_closure_err_over_optimal_diagram_err"], rel=0.05)
    check(S, "U=0.15: better than best closure by 5.5", 5.5, O["0.150"]["g3"]["best_closure_err_over_optimal_diagram_err"], digits=1)
    check(S, "U=0.20: better than best closure by 1.7", 1.7, O["0.200"]["g3"]["best_closure_err_over_optimal_diagram_err"], digits=1)
    check(S, "U=0.225: break-even (ratio ~1)", 1.0, O["0.225"]["g3"]["best_closure_err_over_optimal_diagram_err"], tol=0.15)
    check(S, "U=0.25: worse than best closure by 1.6", 1.6, 1 / O["0.250"]["g3"]["best_closure_err_over_optimal_diagram_err"], digits=1)
    check(S, "U=0.30: worse than best closure by 5.7", 5.7, 1 / O["0.300"]["g3"]["best_closure_err_over_optimal_diagram_err"], digits=1)
    # qualifiers of the 'most accurate' claim: observable, range, cost
    g2win = [float(U) for U in O if float(U) > 0 and O[U]["g2"]["best_closure_err_over_optimal_diagram_err"] > 1]
    check(S, "optimal diagrams ahead of the best closure on g2 only up to U=0.125", 0.125, max(g2win), tol=0)
    g3win = [float(U) for U in O if float(U) > 0 and O[U]["g3"]["best_closure_err_over_optimal_diagram_err"] > 1]
    check(S, "optimal diagrams ahead of the best closure on g3 up to U=0.225 (claim: U<=0.2)", 0.225, max(g3win), tol=0)
    Tn = J("benchmarks", "truncation", "truncation_results.json")["ring_K3"]["results"]
    nwin = [float(U) for U in Tn if float(U) > 0 and min(C["results"][U]["closure3"]["err_n"], C["results"][U]["closure4"]["err_n"]) > Tn[U]["n"]["true_error"] / Tn[U]["n"]["exact"]]
    check(S, "optimal diagrams ahead of the best closure on n only up to U=0.075", 0.075, max(nwin), tol=0)
    raw = J("benchmarks", "truncation", "exact_pool_raw.json")
    tser = [dt for k, v, dt in raw if k[0] == "series" and 0 < k[1] <= 0.3]
    check(S, "orders 0-10 of n, n2, n3 at K=3, symbolic recursion: 54 s per coupling (median of the raw runs)", 54, float(np.median(tser)), tol=3)
    ST = J("benchmarks", "truncation", "sparse_ring_timing.json")
    check(S, "sparse-matrix construction: exponent at second order over K=6-12 (4.9)", 4.9, ST["scaling"]["exponent_N2_K6to12"], digits=1)
    check(S, "sparse-matrix construction: exponent at fourth order over K=6-12 (6.4)", 6.4, ST["scaling"]["exponent_N4_K6to12"], digits=1)
    check(S, "sparse-matrix construction at K=12: second order 7.5 s", 7.5, ST["scaling"]["12"]["2"]["t"], digits=1)
    check(S, "sparse-matrix construction at K=12: fourth order 240 s", 240, ST["scaling"]["12"]["4"]["t"], tol=5)
    PRn = ST["pruning"]
    for name, tot, last in (("none", 183618, 81793), ("degree_bound", 13927, 28), ("degree_bound_momentum", 4663, 10)):
        check(S, f"SM Table (pruning), {name}: total labels", tot, PRn[name]["total"], tol=0)
        check(S, f"SM Table (pruning), {name}: labels after vertex 7", last, PRn[name]["labels_per_layer"][7], tol=0)
    check(S, "SM Table (pruning): same coefficient c8 in all three columns", 0, max(abs(PRn[k]["c8"]-PRn["none"]["c8"]) for k in PRn), tol=1e-9)
    check(S, "SM Table (pruning): time without selection, 86 s", 86, PRn["none"]["t"], tol=15)
    check(S, "SM Table (pruning): time with both rules, 1.3 s", 1.3, PRn["degree_bound_momentum"]["t"], tol=0.5)
    check(S, "largest label set of the sparse form at Nmax=10 (3702)", 3702, max(max(v) for v in ST["K3"]["10"]["labels_per_layer"].values()), tol=0)
    sys.path.insert(0, os.path.join(CODE, "engines"))
    import multimode as MM
    from sparse_ring import SparseRing
    Rr = SparseRing(K=3, kappa=1.0, Delta=-1.0, J=0.4, eta=1.0); b, bd = Rr.site_ops(0)
    _, cnt = Rr.coefficients(MM.mul(MM.mul(MM.mul(bd, bd), bd), MM.mul(MM.mul(b, b), b)), 6, return_counts=True)
    check(S, "labels after vertices 0-6 for n3 at Nmax=6: 136, 408, 691, 312, 72, 10, 1", 0, 0 if cnt == [136, 408, 691, 312, 72, 10, 1] else 1, tol=0)
    # all plotted correlators lie above one
    allg = [tab[U][k][g] for U in tab for k in ("closure2", "closure3", "closure4", "diagrams2", "diagrams4") for g in ("g2", "g3")]
    check(S, "no plotted correlator below one (min g2, g3)", 1.0, min(allg), tol=min(allg) - 1.0 if min(allg) >= 1 else -1)

def table_SI():
    S = "SM Table SI"
    sm = paper_table("tab:closures")
    body = sm[sm.index("\\label{tab:closures}"):]; body = body[:body.index("\\end{tabular}")]
    C = J("benchmarks", "cumulant_closures", "cumulant_closure_results.json"); R = C["results"]; O = C["optimal_truncation_diagrams"]["results"]
    T = J("benchmarks", "truncation", "truncation_results.json")["ring_K3"]["results"]
    def f(x):
        m, e = f"{x:.1e}".split("e"); return f"{m}({int(e)})"
    U = None; n = 0
    for line in body.splitlines():
        m = re.match(r"\$?(0\.\d\d)?\$?\s*&\s*\$(g\^\{\(3\)\}|g\^\{\(2\)\}|n)\$\s*&(.*)\\\\", line.strip())
        if not m: continue
        U = m.group(1) or U; obs = {"g^{(3)}": "g3", "g^{(2)}": "g2", "n": "n"}[m.group(2)]
        cells = [c.strip().strip("$") for c in m.group(3).split("&")]
        r = R[f"{float(U):.3f}"]; o = O[f"{float(U):.3f}"]; t = T[f"{float(U):.3f}"]
        exp = [f"{r['exact'][obs]:.6f}"] + [f(r[k][f"err_{obs}"]) for k in ("closure2", "closure3", "closure4", "diagrams2", "diagrams4")]
        exp.append(f(o[obs]["err_vs_converged_exact"]) if obs != "n" else f(t["n"]["true_error"] / t["n"]["exact"]))
        ok = all(c == e or (c == "1.05(-1)" and e in ("1.0(-1)", "1.1(-1)")) for c, e in zip(cells, exp)); n += 1
        rows.append((S, f"U={U} {obs}: {' '.join(cells)}", " ".join(cells), " ".join(exp), ok))
    rows.append((S, "rows checked", 15, n, n == 15))

# ----------------------------------------------------------------------------- Sec. VII  counting statistics
def fcs():
    S = "Sec. VII / Fig. 10"
    T = J("benchmarks", "truncation", "truncation_results.json")["chain8"]["results"]
    key = lambda U: [k for k in T if abs(T[k]["U"] - U) < 1e-9][0]
    d = lambda U, q: T[key(U)]["fullRS"][q]
    g = lambda U: T[key(U)]["gaussian"]["value"]
    mf1 = 1.0
    r10 = key(0.1)
    check(S, "U=0.1: deviation of Fano factor from Poisson, 24%", 24, 100 * (d(0.1, "fano")["value"] - 1), digits=0)
    check(S, "U=0.1: deviation of c3/c1 from Poisson, 101%", 101, 100 * (d(0.1, "c3c1")["value"] - 1), digits=0)
    check(S, "U=0.1: Gaussian share of the Fano deviation, 93%", 93, 100 * (g(0.1)[1] - 1) / (d(0.1, "fano")["value"] - 1), digits=0)
    check(S, "U=0.1: Gaussian share of the c3/c1 deviation, 74%", 74, 100 * (g(0.1)[2] - 1) / (d(0.1, "c3c1")["value"] - 1), digits=0)
    check(S, "U=0.02: diagrams minus Gaussian on c3/c1, 0.002", 0.002, d(0.02, "c3c1")["value"] - g(0.02)[2], digits=3)
    check(S, "U=0.10: diagrams minus Gaussian on c3/c1, 0.27", 0.27, d(0.1, "c3c1")["value"] - g(0.1)[2], digits=2)
    rat = [(d(U, "c3c1")["value"] - g(U)[2]) / (d(U, "fano")["value"] - g(U)[1]) for U in (0.02, 0.04, 0.06, 0.08, 0.1)]
    check(S, "non-Gaussian part: c3/c1 vs Fano ratio, min 17", 17, min(rat), digits=0)
    check(S, "non-Gaussian part: c3/c1 vs Fano ratio, max 'forty' (code gives 38)", 40, max(rat), tol=2.5)
    check(S, "U=0.08: delta on Fano, 0.0012", 0.0012, d(0.08, "fano")["delta"], digits=4)
    check(S, "U=0.08: delta on c3/c1, 0.014", 0.014, d(0.08, "c3c1")["delta"], digits=3)
    check(S, "U=0.10: delta on Fano, 0.0045", 0.0045, d(0.1, "fano")["delta"], digits=4)
    check(S, "U=0.10: delta on c3/c1, 0.053", 0.053, d(0.1, "c3c1")["delta"], digits=3)
    check(S, "N* = 6 of six orders (all cumulants, U=0.1)", 6, min(d(0.1, q)["N_star"] for q in ("c1", "fano", "c3c1")), tol=0)
    inc = d(0.1, "c3c1")["differences"]
    check(S, "U=0.1: c3/c1 increments shrink by ~0.65 per order", 0.65, inc[-1] / inc[-2], digits=1)
    tail = inc[-1] * (inc[-1] / inc[-2]) / (1 - inc[-1] / inc[-2])
    check(S, "U=0.1: geometric tail of c3/c1, about 0.1", 0.1, tail, tol=0.02)
    # trajectories
    R = J("benchmarks", "trajectories", "traj_fcs_results.json")["results"]
    for U, nt, tt, cut in (("0.04", 1216, 1.2, [10, 6, 5, 4, 3, 3, 2, 2]), ("0.08", 3072, 3.1, [12, 7, 5, 4, 3, 3, 2, 2]), ("0.10", 3552, 3.6, [14, 7, 5, 4, 3, 3, 2, 2])):
        check(S, f"U={U}: number of trajectories", nt, R[U]["ntraj"], tol=0)
        check(S, f"U={U}: total counting time, {tt}e6/kappa", tt, R[U]["total_time"] / 1e6, digits=1)
        check(S, f"U={U}: Fock cutoffs {cut}", 0, 0 if R[U]["cut"] == cut else 1, tol=0)
    e, s = R["0.10"]["estimate"], R["0.10"]["se"]
    check(S, "U=0.1 trajectories: c1 = 0.9309(6)", 0.9309, e[0], digits=4); check(S, "   se(c1) = 0.0006", 0.0006, s[0], digits=4)
    check(S, "U=0.1 trajectories: c2/c1 = 1.2363(49)", 1.2363, e[1], digits=4); check(S, "   se = 0.0049", 0.0049, s[1], digits=4)
    check(S, "U=0.1 trajectories: c3/c1 = 2.054(49)", 2.054, e[2], digits=3); check(S, "   se = 0.049", 0.049, s[2], digits=3)
    ref = R["0.10"]["comparison"]["references"]
    check(S, "U=0.1 diagrams N<=6: c1 = 0.9306", 0.9306, ref["diagrams_N6"]["value"][0], digits=4)
    check(S, "U=0.1 diagrams N<=6: c2/c1 = 1.2401", 1.2401, ref["diagrams_N6"]["value"][1], digits=4)
    check(S, "U=0.1 diagrams N<=6: c3/c1 = 2.013", 2.013, ref["diagrams_N6"]["value"][2], digits=3)
    check(S, "U=0.1 diagrams within one standard error (max |z|)", 1.0, max(abs(z) for z in ref["diagrams_N6"]["z"]), tol=1 - max(abs(z) for z in ref["diagrams_N6"]["z"]) if max(abs(z) for z in ref["diagrams_N6"]["z"]) <= 1 else -1)
    check(S, "U=0.1 third order: c1 = 0.9282", 0.9282, ref["diagrams_N3"]["value"][0], digits=4)
    check(S, "U=0.1 third order: c2/c1 = 1.2075", 1.2075, ref["diagrams_N3"]["value"][1], digits=4)
    check(S, "U=0.1 third order: c3/c1 = 1.749", 1.749, ref["diagrams_N3"]["value"][2], digits=3)
    check(S, "U=0.1 third order low by 4.6 s.e. on c1", 4.6, ref["diagrams_N3"]["z"][0], digits=1)
    check(S, "U=0.1 third order low by 5.8 s.e. on c2/c1", 5.8, ref["diagrams_N3"]["z"][1], digits=1)
    check(S, "U=0.1 third order low by 6.2 s.e. on c3/c1", 6.2, ref["diagrams_N3"]["z"][2], digits=1)
    check(S, "U=0.1 Gaussian: c1 = 0.9305", 0.9305, ref["gaussian"]["value"][0], digits=4)
    check(S, "U=0.1 Gaussian: c2/c1 = 1.2241", 1.2241, ref["gaussian"]["value"][1], digits=4)
    check(S, "U=0.1 Gaussian: c3/c1 = 1.746", 1.746, ref["gaussian"]["value"][2], digits=3)
    check(S, "U=0.1 Gaussian low by 2.5 s.e. on Fano", 2.5, ref["gaussian"]["z"][1], digits=1)
    check(S, "U=0.1 Gaussian low by 6.3 s.e. on c3/c1", 6.3, ref["gaussian"]["z"][2], digits=1)
    check(S, "U=0.1 non-Gaussian part of c3/c1, trajectories 0.31", 0.31, e[2] - ref["gaussian"]["value"][2], digits=2)
    check(S, "U=0.1 non-Gaussian part of c3/c1, diagrams 0.27", 0.27, ref["diagrams_N6"]["value"][2] - ref["gaussian"]["value"][2], digits=2)
    e8, s8, ref8 = R["0.08"]["estimate"], R["0.08"]["se"], R["0.08"]["comparison"]["references"]
    check(S, "U=0.08 trajectories: c1 = 0.9009(6)", 0.9009, e8[0], digits=4)
    check(S, "U=0.08 trajectories: c2/c1 = 1.1541(49)", 1.1541, e8[1], digits=4)
    check(S, "U=0.08 trajectories: c3/c1 = 1.585(45)", 1.585, e8[2], digits=3)
    check(S, "U=0.08 diagrams: c1 = 0.9002", 0.9002, ref8["diagrams_N6"]["value"][0], digits=4)
    check(S, "U=0.08 diagrams: c2/c1 = 1.1613", 1.1613, ref8["diagrams_N6"]["value"][1], digits=4)
    check(S, "U=0.08 diagrams: c3/c1 = 1.617", 1.617, ref8["diagrams_N6"]["value"][2], digits=3)
    check(S, "U=0.08 diagrams within 1.5 s.e. (max |z|)", 1.5, max(abs(z) for z in ref8["diagrams_N6"]["z"]), digits=1)
    check(S, "U=0.08 Gaussian: c3/c1 = 1.499", 1.499, ref8["gaussian"]["value"][2], digits=3)
    check(S, "U=0.08 Gaussian low by 1.9 s.e. on c3/c1", 1.9, ref8["gaussian"]["z"][2], digits=1)
    check(S, "U=0.08: delta on c2/c1 quoted as (12)", 0.0012, d(0.08, "fano")["delta"], digits=4)
    check(S, "U=0.08: delta on c3/c1 quoted as (14)", 0.014, d(0.08, "c3c1")["delta"], digits=3)
    e4, s4, ref4 = R["0.04"]["estimate"], R["0.04"]["se"], R["0.04"]["comparison"]["references"]
    check(S, "U=0.04 trajectories: c1 = 0.8507(9)", 0.8507, e4[0], digits=4)
    check(S, "U=0.04 trajectories: c2/c1 = 1.0635(71)", 1.0635, e4[1], digits=4)
    check(S, "U=0.04 trajectories: c3/c1 = 1.286(58)", 1.286, e4[2], digits=3)
    check(S, "U=0.04 trajectory error on c3/c1, 0.058", 0.058, s4[2], digits=3)
    check(S, "U=0.04 diagrams: c1 = 0.8497", 0.8497, ref4["diagrams_N6"]["value"][0], digits=4)
    check(S, "U=0.04 diagrams: c2/c1 = 1.0597", 1.0597, ref4["diagrams_N6"]["value"][1], digits=4)
    check(S, "U=0.04 diagrams: c3/c1 = 1.196", 1.196, ref4["diagrams_N6"]["value"][2], digits=3)
    check(S, "U=0.04 non-Gaussian part of c3/c1, 0.013", 0.013, ref4["diagrams_N6"]["value"][2] - ref4["gaussian"]["value"][2], digits=3)
    check(S, "U=0.04 diagrams within 1.6 s.e.", 1.6, max(abs(z) for z in ref4["diagrams_N6"]["z"]), digits=1)
    check(S, "U=0.04 Gaussian within 1.8 s.e.", 1.8, max(abs(z) for z in ref4["gaussian"]["z"]), digits=1)
    for U, key_ in (("0.04", "0.04"), ("0.08", "0.08"), ("0.10", "0.10")):
        c = R[U]["c1_check"]
        check(S, f"U={U}: counted c1 vs kappa<n_0>, |z| (SM)", {"0.04": 1.1, "0.08": 1.0, "0.10": 0.3}[U], abs(c["z_martingale"]), digits=1)
    c8 = R["0.10_cut8"]
    check(S, "U=0.1 cutoff-8 run: c1 = 0.9283(10)", 0.9283, c8["estimate"][0], digits=4)
    check(S, "U=0.1 cutoff-8 run: c2/c1 = 1.2098(84)", 1.2098, c8["estimate"][1], digits=4)
    check(S, "U=0.1 cutoff-8 run: c3/c1 = 1.742(85)", 1.742, c8["estimate"][2], digits=3)
    zc = [(e[i] - c8["estimate"][i]) / c8["se"][i] for i in range(3)]   # in units of the cutoff-8 run's own standard error
    check(S, "cutoff-8 run low by 2.6 s.e. (c1) [in the cutoff-8 run's own s.e.]", 2.6, zc[0], digits=1)
    check(S, "cutoff-8 run low by 3.2 s.e. (c2/c1)", 3.2, zc[1], digits=1)
    check(S, "cutoff-8 run low by 3.7 s.e. (c3/c1)", 3.7, zc[2], digits=1)
    # Table V entries that are sensitive to the reference method (finite differences at h=0.05 / laboratory-basis cutoff)
    T5 = J("benchmarks", "truncation", "table5_corrected.json")
    check("Table V", "K=1, U=0.10: exact c3/c1 (finite differences would give 2.4310)", 2.3585, T5["K1"]["0.10"]["exact"]["c3c1"], digits=4)
    check("Table V", "K=3, U=0.05: exact c3/c1 (lab-basis finite differences would give 1.6116)", 1.6393, T5["K3"]["0.05"]["exact"]["c3c1"], digits=4)
    check("Table V", "K=1, U=0.05: mean-field c1 (finite differences would give 0.8587)", 0.8576, T5["K1"]["0.05"]["mean_field"]["c1"], digits=4)
    check("Table V", "K=1, U=0.10: mean-field c1 (finite differences would give 0.9399)", 0.9327, T5["K1"]["0.10"]["mean_field"]["c1"], digits=4)
    V = J("benchmarks", "trajectories", "traj_validation.json")
    check("SM, validation", "K=3 ring (Table V): trajectories c1 = 1.1456(2)", 1.1456, V["results"]["K3ring_U0.05"]["primary"]["est"][0] if "primary" in V["results"].get("K3ring_U0.05", {}) else _find(V, "K3ring_U0.05", 0), digits=4)

def _find(V, run, i):
    r = V["results"].get(run)
    if r is None: return None
    for k in ("primary", "estimate", "two_window"):
        if k in r:
            v = r[k]; v = v.get("est", v) if isinstance(v, dict) else v
            return v[i]
    return None

# ----------------------------------------------------------------------------- SM Kerr finite time, emitter, atom-cavity (from make_figures.py)
def sm_figures():
    S = "SM Figs. S1, S5, S7 / Fig. 7"
    nfile = os.path.join(CODE, "engines", "numbers.json")
    if not os.path.exists(nfile):
        rows.append((S, "numbers.json absent: run engines/make_figures.py kerr cubic emitter finite atomcavity", "", None, None)); return
    n = json.load(open(nfile))
    check(S, "emitter (a): max relative error of fourth order, 0.37%", 0.37, 100*n["emitter_a_maxrel_err_N4"], digits=2)
    check(S, "emitter (b): dressed fourth-order max relative error, 2.0%", 2.0, 100*n["emitter_b_dressed4_maxrel_err"], digits=1)
    check(S, "Kerr switch-on: max error of N<=8 partial sum, 1.46e-3", 1.46e-3, n["finite_maxerr_N8"], tol=0.005e-3)
    check(S, "atom in cavity: g2(0) on resonance exact, 0.8936", 0.8936, n["atomcav"]["g2_res"], digits=4)
    check(S, "atom in cavity: g2(0) on resonance N<=32, 0.8936", 0.8936, n["atomcav"]["g2_res_N32"], digits=4)
    for N, paper, code in zip((10,20,30,40,50,60), (1.074,0.864,0.789,0.758,0.707,0.664), n["cubic_root_test"]):
        check(S, f"cubic root test |c_N|^(-1/N) at N={N}", paper, code, digits=3)

# ----------------------------------------------------------------------------- report
def main():
    out = None
    if "--report" in sys.argv: out = sys.argv[sys.argv.index("--report") + 1]
    table5(); table_SI(); kerr(); ring(); closures(); fcs(); sm_figures()
    npass = sum(1 for r in rows if r[4] is True); nfail = sum(1 for r in rows if r[4] is False); nskip = sum(1 for r in rows if r[4] is None)
    lines = ["| section | quantity | paper | code | ok |", "|---|---|---|---|---|"]
    fmt = lambda v: "" if v is None else (f"{v:.6g}" if isinstance(v, float) else str(v))
    for s, w, p, c, ok in rows:
        lines.append(f"| {s} | {w} | {fmt(p)} | {fmt(c)} | {'PASS' if ok else ('FAIL' if ok is False else 'info')} |")
    summary = f"{npass} checks passed, {nfail} failed, {nskip} informational"
    src = f"LaTeX sources in {PAPER}" if PAPER else "snapshots in notebooks/paper_reference/"
    text = "# Verification of the numbers quoted in the paper\n\n" + summary + f"\n\nTables V and SI compared against the {src}." + "\n\n" + "\n".join(lines) + "\n"
    if out:
        open(out, "w", encoding="utf-8").write(text); print("wrote", os.path.relpath(out))
    for s, w, p, c, ok in rows:
        if ok is False: print(f"FAIL  [{s}] {w}: paper {fmt(p)}  code {fmt(c)}")
    print(summary)
    sys.exit(1 if nfail else 0)

if __name__ == "__main__":
    main()
