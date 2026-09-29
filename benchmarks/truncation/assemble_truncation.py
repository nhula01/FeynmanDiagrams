"""Assemble the series benchmarks into truncation_results.json, timing_results.json, table5_corrected.json,
table5_fragment.tex; redraw fig_kerr_convergence.pdf and fig_ring.pdf.

Every series is reported at a stated order (the highest order computed unless a lower one is named) with, where an
exact reference exists, the error of every partial sum.  No order is selected from the series itself."""
import sys, os, json, glob, math, numpy as np
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(_ROOT, 'engines'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from series import series_record
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 8, "axes.labelsize": 8, "legend.fontsize": 6.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.6, "lines.linewidth": 1.0, "xtick.direction": "in",
    "ytick.direction": "in", "xtick.top": True, "ytick.right": True, "legend.frameon": False})
COLW = 3.375
ORDER_COLORS = ["#d95f02", "#1b9e77", "#7570b3", "#e7298a", "#66a61e", "#a6761d"]
J = lambda f: json.load(open(f))
ex = lambda f: os.path.exists(f)
RP = J(os.path.join(_ROOT, 'notebook', 'data', 'revision_params.json'))
res = dict(convention="each series: partial sums S_0..S_Nmax, the value S_N at the stated order N (the highest computed unless named), "
                      "and, where an exact reference exists, the error |S_n - exact| of every partial sum")
# =========================================================== Kerr cavity
kerr = J('kerr_truncation.json'); res["kerr_cavity"] = kerr
# =========================================================== ring K=3
raw = J('exact_pool_raw.json')
series = {round(k[1], 3): v for k, v, dt in raw if k[0] == "series"}
rex = {}
for k, v, dt in raw:
    if k[0] == "ring_exact": rex.setdefault(round(k[1], 3), {})[k[2]] = dict(v, time=dt)
ring = dict(parameters=dict(K=3, **{q: RP["ring3"][q] for q in ("kappa", "Delta", "J", "eta")}),
            note="momentum engine with degree pruning and the zero-contraction-free commutator, orders 0..10; g2_N = n2_N/n_N^2, "
                 "g3_N = n3_N/n_N^3 from partial sums through order N; exact = sparse time evolution to T=60 at 7 Fock levels per site "
                 "(checks at 6 and 8 levels)", exact={}, results={}, results_Nmax6={}, labels={})
cached = J(os.path.join(_ROOT, 'notebook', 'data', 'ring_exactU.json'))
DISP = {}
for f in ('exact_disp_raw.json', 'exact_disp_raw2.json'):
    if not ex(f): continue
    for k, v, dt in J(f):
        if k[0] != "ring": continue
        d = k[3] if len(k) > 3 else "lin"; DISP.setdefault(round(k[1], 3), {})[(k[2], d)] = v
ring["note"] = ("momentum engine with degree pruning and the zero-contraction-free commutator, orders 0..10; g2_N = n2_N/n_N^2, g3_N = n3_N/n_N^3 "
                "from partial sums through order N; exact reference = displaced-frame steady state (a_j = alpha_j + b_j, b-Fock space truncated), "
                "highest b-cutoff available, uncertainty = change from the next lower b-cutoff; laboratory-basis values (7 Fock levels) kept for comparison")
for U in sorted(series):
    key = f"{U:.3f}"; e = rex[U][7]
    lab = dict(n0=e["n0"], g2=e["n0n0"]/e["n0"]**2, g3=e["n3"]/e["n0"]**3, levels=7)
    for nm in (6, 8):
        if nm in rex[U]:
            f = rex[U][nm]; lab[f"g2_levels{nm}"] = f["n0n0"]/f["n0"]**2; lab[f"g3_levels{nm}"] = f["n3"]/f["n0"]**3
    dd = DISP.get(U, {}); rec = dict(lab_basis=lab, displaced={f"{lv}_{d}": dict(g2=v["g2"], g3=v["g3"], n0=v["n0"]) for (lv, d), v in sorted(dd.items())})
    order = [(8, "mf"), (8, "lin"), (7, "mf"), (7, "lin"), (6, "mf"), (6, "lin"), (5, "lin"), (4, "lin")]
    avail = [o for o in order if o in dd]
    if avail:
        best = avail[0]; lower = [o for o in avail if o[0] < best[0]]
        rec.update(n0=dd[best]["n0"], g2=dd[best]["g2"], g3=dd[best]["g3"], source=f"displaced {best[1]}, {best[0]} b-levels")
        if lower:
            lo = lower[0]; rec["g2_cutoff_uncertainty"] = abs(dd[best]["g2"]-dd[lo]["g2"]); rec["g3_cutoff_uncertainty"] = abs(dd[best]["g3"]-dd[lo]["g3"])
            rec["cutoff_compared"] = f"{lo[0]}_{lo[1]}"
    else:
        rec.update(n0=lab["n0"], g2=lab["g2"], g3=lab["g3"], source="laboratory basis, 7 levels")
    if key in cached: rec["cached_ring_exactU"] = dict(n0=cached[key]["n0"], g2=cached[key]["n0n0"]/cached[key]["n0"]**2)
    g2e, g3e = rec["g2"], rec["g3"]
    ring["exact"][key] = rec
    s = series[U]; S = {q: np.cumsum(np.array(s[q]["terms_re"]) + 1j*np.array(s[q]["terms_im"])).real for q in ("n", "n2", "n3")}
    g2 = S["n2"]/S["n"]**2; g3 = S["n3"]/S["n"]**3
    for tag, Mx in (("results", 10), ("results_Nmax6", 6)):
        ring[tag][key] = dict(U=U, n=series_record(S["n"][:Mx+1], exact=rec["n0"]), g2=series_record(g2[:Mx+1], exact=g2e),
                              g3=series_record(g3[:Mx+1], exact=g3e))
    ring["labels"][key] = {q: s[q]["labels"] for q in ("n", "n2", "n3")}
res["ring_K3"] = ring
# =========================================================== label counts (part C)
LM = J('label_counts_momentum.json'); LT = J('label_counts_tensor.json')
res["label_counts"] = dict(momentum=LM, tensor=LT, analytic="D_max = max_s min(d0+2s, 2(N-s)) = 2*floor(N/2 + d0/4); labels <= C(2K+D_max, D_max) ~ (2K)^D_max/D_max!; "
                           "zero-momentum labels on the ring ~ K^(D_max-1)")
# =========================================================== Table V (contour derivatives)
def taylor_half(vals_half, r, M, nmax=5):
    full = [None]*M
    for j in range(M//2+1): full[j] = np.asarray(vals_half[j])
    for j in range(M//2+1, M): full[j] = np.conj(full[M-j])
    co = np.fft.fft(np.array(full), axis=0)/M
    return np.array([math.factorial(n)*co[n]/r**n for n in range(nmax+1)])
def C(x):
    x = np.array(x, float); return x[..., 0] + 1j*x[..., 1]
t5 = dict(method="cumulants c_n = n! [chi^n] theta(chi) by the discrete Cauchy integral on |chi| = r (M points, trapezoidal rule = FFT), "
                 "theta(conj chi) = conj theta(chi); radius checked for independence", K1={}, K3={})
nzr = (J('fcs_nz_raw.json') if ex('fcs_nz_raw.json') else []) + (J('fcs_nz_raw_md.json') if ex('fcs_nz_raw_md.json') else [])
groups = {}; c8g = {}; gfd = {}
for a, th, dt in nzr:
    if a[0] == "k3":
        _, N, md, rs, r, M, j = a; groups.setdefault((N, md, rs, r, M), {})[j] = C(th)
    elif a[0] == "c8":
        _, U, N, md, rs, r, M, j = a; c8g.setdefault((round(U, 3), N, md, rs, r, M), {})[j] = C(th)
    elif a[0] == "gfd":
        _, U, Nc, h = a; gfd[(round(U, 3), Nc, h)] = np.array(th)
k1e = [x for x in J('table5_k1_raw_exact.json')]; k1g = J('table5_k1_raw_gauss.json'); k1d = J('table5_k1_raw_diag.json')
k1r = J('k1_radius.json') if ex('k1_radius.json') else {}
def ratios(cn): return dict(c1=float(cn[1].real), fano=float((cn[2]/cn[1]).real), c3c1=float((cn[3]/cn[1]).real))
for U in (0.05, 0.10):
    d = {}
    ee = {tuple(a[1:]): ratios(C(cc)) for kind, a, cc, dt in k1e if a[0] == U}
    d["exact_all"] = {f"cut{a[0]}_r{a[1]}_M{a[2]}": v for a, v in ee.items()}
    rr = {k: v for k, v in k1r.items() if k.startswith(f"U{U:.2f}_r")}
    d["exact_radius_scan_cut26"] = {k: dict(c1=v["c1"], fano=v["fano"], c3c1=v["c3c1"], c5c1=v["c5c1"], c0_abs=v["c0_abs"]) for k, v in rr.items()}
    if f"U{U:.2f}_r0.05_M32" in k1r:
        v = k1r[f"U{U:.2f}_r0.05_M32"]; d["exact"] = dict(c1=v["c1"], fano=v["fano"], c3c1=v["c3c1"], source="cut 26, r=0.05, M=32")
    else:
        d["exact"] = dict(ee[(26, 0.1, 32)], source="cut 26, r=0.1, M=32")
    gg = {tuple(a[1:]): ratios(C(cc)) for kind, a, cc, dt in k1g if a[0] == U}
    d["gaussian_contour"] = dict(gg[(40, 0.1, 32)], source="gaussian_fcs.gaussian_model (frozen HFB alpha, n, m; linear tilted Liouvillian), cut 40, r=0.1, M=32")
    d["gaussian_all"] = {str(k): v for k, v in gg.items()}
    hs = sorted([h for (u, nc, h) in gfd if u == round(U, 3)], reverse=True)
    if hs:
        fd = []
        for h in hs:
            th = gfd[(round(U, 3), 30, h)]
            c1 = (th[3]-th[1])/(2*h); c2 = (th[3]-2*th[2]+th[1])/h**2; c3 = (th[4]-2*th[3]+2*th[1]-th[0])/(2*h**3); fd.append([c1, c2, c3])
        fd = np.array(fd)                                   # rows h, h/2, h/4, ...
        R1 = (4*fd[1:] - fd[:-1])/3; R2 = (16*R1[1:] - R1[:-1])/15       # remove h^2 then h^4
        best = R2[-1]; d["gaussian_richardson"] = dict(c1=float(best[0]), fano=float(best[1]/best[0]), c3c1=float(best[2]/best[0]), h=hs,
            fd_raw=[dict(h=h, c1=float(x[0]), fano=float(x[1]/x[0]), c3c1=float(x[2]/x[0])) for h, x in zip(hs, fd)],
            change_last=[float(abs(R2[-1][i]/R2[-1][0] - R2[-2][i]/R2[-2][0])) if len(R2) > 1 else None for i in range(3)], cut=30)
    d["gaussian"] = d.get("gaussian_richardson", d["gaussian_contour"])
    # mean field
    Dc, kap, eta = -1.0, 1.0, 1.0
    rts = np.roots([U**2, 2*Dc*U, Dc**2 + kap**2/4, -abs(eta)**2]); nmf = float(np.min(rts[np.abs(rts.imag) < 1e-9].real))
    d["mean_field"] = dict(c1=kap*nmf, fano=1.0, c3c1=1.0)
    for kind, a, cc, dt in k1d:
        if a[0] != U: continue
        _, md, rs, r, M = a; cn = C(cc)                   # (6, 9): cumulant index, order
        seq = dict(c1=cn[1].real, fano=(cn[2]/cn[1]).real, c3c1=(cn[3]/cn[1]).real)
        d[f"diagrams_maxdeg{md}_{'fullRS' if rs else 'notebookRS'}_r{r}"] = {q: series_record(seq[q], exact=d["exact"][q]) for q in seq}
    t5["K1"][f"{U:.2f}"] = d
# K=3 diagrams
k3pts = {}
for k, v, dt in raw:
    if k[0] == "k3pt": k3pts.setdefault((tuple(k[1]), k[2]), {})[k[3]] = v["theta"][0] + 1j*v["theta"][1]
k3ss = {tuple(k[1]): v for k, v, dt in raw if k[0] == "k3ss"}
d3 = dict(parameters=RP["table5_fcs"]["K3"], exact_all={}, steady_state_checks={str(k): v for k, v in k3ss.items()})
for (cuts, r), pts in sorted(k3pts.items()):
    if len(pts) == 9:
        cn = taylor_half([pts[j] for j in range(9)], r, 16); d3["exact_all"][f"cut{'-'.join(map(str, cuts))}_r{r}"] = dict(ratios(cn), c0_abs=float(abs(cn[0])))
d3["exact_lab_basis"] = dict(d3["exact_all"]); d3["exact_all"] = {}
if ex('exact_disp_raw.json'):
    for k, v, dt in J('exact_disp_raw.json'):
        if k[0] == "fcs" and k[1] == "K3": d3["exact_all"][f"displaced_cut{'-'.join(map(str, k[2]))}_r{k[4]}"] = dict(c1=v["c1"], fano=v["fano"], c3c1=v["c3c1"], c0_abs=v["c0_abs"])
        if k[0] == "fcs" and k[1] == "K1": t5["K1"][f"{k[3]:.2f}"].setdefault("exact_displaced_frame", {})[f"cut{k[2][0]}"] = dict(c1=v["c1"], fano=v["fano"], c3c1=v["c3c1"])
b = d3["exact_all"].get("displaced_cut9-6-6_r0.15"); b2 = d3["exact_all"].get("displaced_cut8-6-6_r0.15")
d3["exact"] = dict(b, source="displaced frame, b-cutoffs [9,6,6], r=0.15, M=16",
                   cutoff_uncertainty={q: abs(b[q]-b2[q]) for q in ("c1", "fano", "c3c1")}) if b and b2 else None
# K=3 mean field: U=0 displacement (paper) and self-consistent
from scipy.optimize import fsolve
p3 = RP["table5_fcs"]["K3"]; U3 = p3["U"]; eta3 = np.zeros(3); eta3[p3["eta_site"]] = 1.0
Om = np.diag([p3["Delta"]]*3).astype(complex)
for (i, k) in ((0, 1), (1, 2), (2, 0)): Om[i, k] = Om[k, i] = -p3["J"]
al0 = np.linalg.solve(0.5*np.eye(3) + 1j*Om, eta3)
def f(x):
    a = x[:3] + 1j*x[3:]; rhs = -(0.5*np.eye(3) + 1j*Om)@a - 1j*U3*np.abs(a)**2*a + eta3; return np.concatenate([rhs.real, rhs.imag])
xs = fsolve(f, np.concatenate([al0.real, al0.imag]), xtol=1e-14); asc = xs[:3] + 1j*xs[3:]
d3["mean_field_U0_displacement"] = dict(c1=float(abs(al0[p3["eta_site"]])**2), fano=1.0, c3c1=1.0)
d3["mean_field_selfconsistent"] = dict(c1=float(abs(asc[p3["eta_site"]])**2), fano=1.0, c3c1=1.0)
for (N, md, rs, r, M), pts in sorted(groups.items()):
    if len(pts) < M//2+1: continue
    cn = taylor_half([pts[j] for j in range(M//2+1)], r, M)      # (6, N+1)
    seq = dict(c1=cn[1].real, fano=(cn[2]/cn[1]).real, c3c1=(cn[3]/cn[1]).real)
    exd = d3["exact"] or {}
    d3[f"diagrams_N{N}_maxdeg{md}_{'fullRS' if rs else 'notebookRS'}_r{r}"] = {q: series_record(seq[q], exact=exd.get(q)) for q in seq}
t5["K3"]["0.05"] = d3
res["table5"] = t5
json.dump(t5, open('table5_corrected.json', 'w'), indent=1)
# LaTeX fragment
def row(lbl, meth, v):
    fmt = lambda q: f"${v[q]:.4f}$"
    return f"{lbl:<24}& {meth:<34}& {fmt('c1')} & {fmt('fano')} & {fmt('c3c1')} \\\\"
lines = [r"\begin{tabular}{llccc}", r" & method & $c_1/\kappa$ & $c_2/c_1$ & $c_3/c_1$ \\ \hline"]
for U in ("0.05", "0.10"):
    d = t5["K1"][U]; lab = f"$K=1$, $U={float(U):.2f}\\kappa$"
    dg = d["diagrams_maxdeg14_fullRS_r0.1"]; val = {q: dg[q]["value"] for q in dg}; Nd = max(dg[q]["N"] for q in dg)
    lines += [row(lab, "exact", d["exact"]), row("", "mean field", d["mean_field"]), row("", "Gaussian FCS", d["gaussian"]),
              row("", f"diagrams $N\\leq{Nd}$", val) + r"[2pt]"]
if d3.get("exact"):
    lab = r"$K=3$, $U=0.05\kappa$"; lines += [row(lab, "exact", d3["exact"]), row("", "mean field ($U=0$ displacement)", d3["mean_field_U0_displacement"])]
    for N in (2, 3, 4):
        k = f"diagrams_N{N}_maxdeg{N+2}_fullRS_r0.1"
        if k in d3: lines.append(row("", f"diagrams $N\\leq{N}$", {q: d3[k][q]["partial_sums"][-1] for q in ("c1", "fano", "c3c1")}))
    k = "diagrams_N6_maxdeg8_fullRS_r0.1"
    if k in d3:
        lines.append(row("", f"diagrams $N\\leq{d3[k]['c1']['N']}$", {q: d3[k][q]["value"] for q in ("c1", "fano", "c3c1")}))
lines.append(r"\end{tabular}")
open('table5_fragment.tex', 'w').write("\n".join(lines) + "\n")
# =========================================================== eight-site chain: from the trajectories track
TR = os.path.join(_ROOT, 'benchmarks', 'trajectories')
fcs = dict(description="Disordered eight-site Kerr chain (revision_params.json chain8_fcs), counting on site 0. "
                       "chain8_contour.json stores the full Rayleigh-Schroedinger diagram series under diagrams_maxdeg5 and the "
                       "recursion without the renormalization terms under diagrams_maxdeg5_notebook_recursion; both are assembled here "
                       "order by order (orders 0..6, reported at N = 6). The fcs_nz raw outputs in this directory independently reproduce both series.",
           source="benchmarks/trajectories/chain8_contour.json", U=[], results={})
fcs_source = os.path.join(TR, 'chain8_contour.json')
if ex(fcs_source):
    cc8 = J(fcs_source)
    for key, rec in sorted(cc8["results"].items()):
        U = rec["U"]; fcs["U"].append(U)
        full_arr = np.array(rec["diagrams_maxdeg5"]["orders"], float)
        notebook_rec = rec.get("diagrams_maxdeg5_notebook_recursion", rec["diagrams_maxdeg5"])
        notebook_arr = np.array(notebook_rec["orders"], float)
        out = dict(U=U, gaussian=rec.get("gaussian"), mean_field=rec.get("mean_field"),
                   fullRS={q: series_record(full_arr[:, i]) for i, q in enumerate(("c1", "fano", "c3c1"))},
                   notebookRS={q: series_record(notebook_arr[:, i]) for i, q in enumerate(("c1", "fano", "c3c1"))})
        grp = [(k, v) for k, v in c8g.items() if k[0] == round(U, 3)]
        for (u, N, md, rs, r, M), pts in grp:
            if len(pts) < M//2+1: continue
            cn = taylor_half([pts[j] for j in range(M//2+1)], r, M)
            seq = dict(c1=cn[1].real, fano=(cn[2]/cn[1]).real, c3c1=(cn[3]/cn[1]).real)
            here = {q: series_record(seq[q]) for q in seq}
            tag = "fullRS" if rs else "notebookRS"
            diff = float(max(np.max(np.abs(np.array(here[q]["partial_sums"]) - np.array(out[tag][q]["partial_sums"]))) for q in ("c1", "fano", "c3c1")))
            if diff > 1e-10: raise RuntimeError(f"chain8 {tag} cross-check failed at U={U}: maxdiff={diff}")
        fcs["results"][key] = out
res["chain8"] = fcs
# =========================================================== timing
def fitexp(Ks, ts):
    p = np.polyfit(np.log(np.array(Ks, float)), np.log(np.array(ts, float)), 1); return float(p[0])
def curves(T2, T3):
    cv = {}
    for Nord in (2, 4):
        t = T2.get(f"tensor{Nord}", {}).get("times", {})
        if t:
            Ks = sorted(int(k) for k in t); cv[f"tensor_N{Nord}_pass1script"] = dict(K=Ks, t=[t[str(k)]["g2_total"] for k in Ks], load1=[t[str(k)]["load1"] for k in Ks])
    m = T2.get("momentum", {}).get("times", {})
    for q in ("n_unpruned_N2_incl_setup", "g2_pruned_N2", "g2_pruned_N4"):
        Ks = sorted(int(k) for k in m if q in m[k])
        if Ks: cv[f"momentum_{q}"] = dict(K=Ks, t=[m[str(k)][q] for k in Ks], load1=[m[str(k)]["load1"] for k in Ks])
    for nm in (4, 5):
        Ks = sorted(int(k.split("_")[1][1:]) for k in T2 if k.startswith("exact_") and k.endswith(f"nmax{nm}"))
        if Ks: cv[f"exact_levels{nm}"] = dict(K=Ks, t=[T2[f"exact_K{k}_nmax{nm}"]["times"][str(k)]["evolve"] for k in Ks],
                                             load1=[T2[f"exact_K{k}_nmax{nm}"]["times"][str(k)]["load1"] for k in Ks],
                                             g2=[T2[f"exact_K{k}_nmax{nm}"]["times"][str(k)]["g2"] for k in Ks])
    tt = T3.get("tensor", {}).get("times", {})
    for name in ("te", "te_nz"):
        for Nord in (2, 4):
            ks = sorted(int(k.split("_K")[1]) for k in tt if k.startswith(f"{name}_N{Nord}_"))
            if ks: cv[f"{name}_N{Nord}"] = dict(K=ks, t=[tt[f"{name}_N{Nord}_K{k}"]["wall"] for k in ks], cpu=[tt[f"{name}_N{Nord}_K{k}"]["cpu"] for k in ks],
                                                load1=[tt[f"{name}_N{Nord}_K{k}"]["load1"] for k in ks], g2=[tt[f"{name}_N{Nord}_K{k}"]["g2"] for k in ks])
    mm = T3.get("momentum", {}).get("times", {})
    for q in ("n_unpruned_N2", "g2_pruned_N2", "g2_nz_N2", "g2_nz_N4"):
        ks = sorted(int(k) for k in mm if f"{q}_wall" in mm[k])
        if ks: cv[f"momentum3_{q}"] = dict(K=ks, t=[mm[str(k)][f"{q}_wall"] for k in ks], load1=[mm[str(k)]["load1"] for k in ks])
    return cv
def load(d, pat): return {os.path.basename(f)[8:-5]: J(f) for f in glob.glob(os.path.join(d, pat))}
timing = dict(note="single core (SLURM --cpus-per-task=1 = one physical core, BLAS threads 1) on a shared 48-core node; wall time equals CPU time "
                   "(no time sharing), the spread between passes is cache/memory-bandwidth contention; load averages recorded per point",
              parameters=dict(kappa=1.0, Delta=-1.0, J=0.4, eta=1.0, U=0.15, observable="on-site g2 at site 0"),
              passes=dict(pass1_2=curves(load('.', 'timing2_*.json'), load('.', 'timing3_*.json')),
                          pass3=curves(load('pass3', 'timing2_*.json'), load('pass3', 'timing3_*.json'))), fits={})
for ps, cvs in timing["passes"].items():
    for name, cv in cvs.items():
        if len(cv["K"]) < 3: continue
        Ks, y = cv["K"], cv["t"]
        if name.startswith("exact"):
            b = np.polyfit(Ks, np.log(y), 1)[0]; timing["fits"][f"{ps}:{name}"] = dict(per_site_factor=float(np.exp(b)), K_range=[Ks[0], Ks[-1]]); continue
        rec = dict(exponent_all=fitexp(Ks, y), K_range=[Ks[0], Ks[-1]])
        if len(Ks) >= 4: rec["exponent_last3"] = fitexp(Ks[-3:], y[-3:])
        timing["fits"][f"{ps}:{name}"] = rec
timing["paper_ring_timing_json"] = J(os.path.join(_ROOT, 'notebook', 'data', 'ring_timing.json'))
json.dump(timing, open('timing_results.json', 'w'), indent=1)
res["timing_fits"] = timing["fits"]
json.dump(res, open('truncation_results.json', 'w'), indent=1, default=float)
# =========================================================== fig_kerr_convergence.pdf
from engine import destroy, steady
kp = kerr["parameters"]; kappa, Dc, eta = kp["kappa"], kp["Delta"], kp["eta"]
cn = np.array(kerr["coefficients"]); Nmax = kerr["Nmax"]
Us = np.linspace(0, 0.3, 61); Nc = 32; A = destroy(Nc); Ad = A.conj().T
exact = np.array([np.trace(Ad@A@steady(Dc*Ad@A + 1j*(eta*Ad - np.conj(eta)*A) + U/2*Ad@Ad@A@A, [np.sqrt(kappa)*A])).real for U in Us])
mf = []
for U in Us:
    r_ = np.roots([U**2, 2*Dc*U, Dc**2 + kappa**2/4, -abs(eta)**2]) if U > 0 else np.array([kp["alpha2"]])
    mf.append(np.min(r_[np.abs(np.imag(r_)) < 1e-9].real))
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(COLW, 4.3))
ax1.plot(Us, exact, "k-", lw=1.6, label="exact")
ax1.plot(Us, mf, color="gray", ls=(0, (4, 2)), lw=1.0, label="mean field ($L=0$)")
pw = Us[:, None]**np.arange(Nmax+1)[None, :]
for col, N in zip(ORDER_COLORS, [1, 2, 4, 8]):
    ax1.plot(Us, (pw[:, :N+1]*cn[None, :N+1]).sum(1), color=col, lw=0.9, label=f"$N\\leq{N}$")
ax1.set_ylim(0.6, 3.2); ax1.set_xlim(0, Us[-1]); ax1.set_xlabel(r"$U/\kappa$"); ax1.set_ylabel(r"$\langle a^\dagger a\rangle$")
lg = ax1.legend(ncol=1, loc="upper left", bbox_to_anchor=(0.06, 1.0)); [l.set_linewidth(1.4) for l in lg.get_lines()]
ax1.text(0.97, 0.05, "(a)", transform=ax1.transAxes, ha="right", va="bottom")
for col, (key, rec) in zip(["#1b9e77", "#7570b3", "#d95f02", "#e7298a"], sorted(kerr["results"].items())):
    err = np.array(rec["error_all_orders"]); U = rec["U"]
    ax2.semilogy(np.arange(Nmax+1), err, "o-", color=col, ms=2.6, lw=0.8, label=f"$U/\\kappa={U}$")
ax2.set_xlabel(r"order $N$"); ax2.set_ylabel(r"$|\langle a^\dagger a\rangle_N-\langle a^\dagger a\rangle|$")
ax2.set_xlim(-0.3, Nmax+0.3); ax2.set_ylim(1e-12, 1e3); ax2.set_xticks(range(0, Nmax+1, 2)); ax2.legend(loc="lower left", ncol=2)
ax2.text(0.45, 0.95, "(b)", transform=ax2.transAxes, ha="left", va="top")
fig.tight_layout(pad=0.3, h_pad=0.6); fig.savefig("fig_kerr_convergence.pdf")
# =========================================================== fig_ring.pdf
Us13 = np.array(sorted(float(k) for k in ring["results"]))
fa = [ring["results"][f"{U:.3f}"]["g2"] for U in Us13]
PS = np.array([r["partial_sums"] for r in fa]); exG = np.array([r["exact"] for r in fa])
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(COLW, 4.5))
ax1.plot(Us13, exG, "k-", lw=1.6, label="exact")
for col, N, ls in zip(ORDER_COLORS, [1, 2, 4, 6, 10], [(0, (4, 2)), ":", "-.", (0, (6, 1.5, 1, 1.5)), "-"]):
    ax1.plot(Us13, PS[:, N], color=col, lw=1.0, ls=ls, label=f"$N\\leq{N}$")
ax1.set_xlabel(r"$U/\kappa$"); ax1.set_ylabel(r"on-site $g^{(2)}(0)$"); ax1.set_xlim(0, Us13[-1])
lg = ax1.legend(loc="upper left", fontsize=6.2); [l.set_linewidth(1.4) for l in lg.get_lines()]
ax1.text(0.97, 0.06, "(a) $K=3$", transform=ax1.transAxes, ha="right", va="bottom")
PS_ = timing["passes"]["pass3"] if "te_N4" in timing["passes"]["pass3"] else timing["passes"]["pass1_2"]
PX_ = timing["passes"]["pass1_2"]
def curve(ax, cv, mk, col, lab, ls="-"):
    K, t = cv["K"], cv["t"]; e = fitexp(K[-3:], t[-3:]) if len(K) >= 4 else fitexp(K, t)
    ax.loglog(K, t, mk + ls, color=col, ms=3.3, lw=1.0, mfc="white" if mk in "o" else col, label=lab + f" ($K^{{{e:.1f}}}$)")
ex5 = PS_.get("exact_levels5") or PX_.get("exact_levels5")
if ex5:
    ke, te = ex5["K"], ex5["t"]; Kref = np.arange(2, 8)
    ax2.loglog(Kref, te[-1]*(25.0**(Kref-ke[-1])), color="#d95f02", ls=":", lw=0.9, label=r"$\propto 25^{K}$")
    ax2.loglog(ke, te, "s", color="#d95f02", ms=3.5, label="exact steady state (5 levels/site)")
mo = PS_.get("momentum3_n_unpruned_N2") or PX_.get("momentum_n_unpruned_N2_incl_setup")
if mo: curve(ax2, mo, "o", "#1b9e77", "momentum engine, $N\\leq2$")
if "te_N2" in PS_: curve(ax2, PS_["te_N2"], "^", "#7570b3", "tensor engine, $N\\leq2$")
if "te_N4" in PS_: curve(ax2, PS_["te_N4"], "v", "#e7298a", "tensor engine, $N\\leq4$")
if "te_nz_N4" in PS_: curve(ax2, PS_["te_nz_N4"], "D", "#66a61e", "tensor, $N\\leq4$, contracted terms only", ls="--")
ax2.set_xlabel(r"ring size $K$"); ax2.set_ylabel("runtime (s)"); ax2.set_xlim(1.8, 36); ax2.set_ylim(1e-3, 1e6)
ax2.set_xticks([2, 3, 4, 6, 8, 12, 16, 24, 32]); ax2.set_xticklabels([str(k) for k in [2, 3, 4, 6, 8, 12, 16, 24, 32]])
ax2.tick_params(axis="x", which="minor", bottom=False, top=False)
from matplotlib.ticker import NullFormatter, NullLocator
ax2.xaxis.set_minor_formatter(NullFormatter()); ax2.xaxis.set_minor_locator(NullLocator())
lg = ax2.legend(loc="upper left", fontsize=5.6); [l.set_linewidth(1.4) for l in lg.get_lines()]
ax2.text(0.97, 0.06, "(b)", transform=ax2.transAxes, ha="right", va="bottom")
fig.tight_layout(pad=0.3, h_pad=0.5); fig.savefig("fig_ring.pdf")
print(json.dumps(timing["fits"], indent=1)); print("assembled")
