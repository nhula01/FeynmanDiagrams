"""Assemble truncation_results.json / timing_results.json and redraw fig_kerr_convergence.pdf, fig_ring.pdf.
Runs on the host in revision/truncation/ after the SLURM jobs have finished."""
import sys, os, json, glob, numpy as np
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(_ROOT, 'engines'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 8,
    "axes.labelsize": 8, "legend.fontsize": 6.5, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "axes.linewidth": 0.6, "lines.linewidth": 1.0,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True, "legend.frameon": False,
})
COLW = 3.375
ORDER_COLORS = ["#d95f02", "#1b9e77", "#7570b3", "#e7298a", "#66a61e", "#a6761d"]
J = lambda f: json.load(open(f))
kerr = J('kerr_truncation.json'); ring = J('ring_truncation.json')
k1 = J('fcs_k1_truncation.json'); k3 = J('fcs_k3_truncation_maxdeg8.json')
k3b = J('fcs_k3_truncation_maxdeg10.json') if os.path.exists('fcs_k3_truncation_maxdeg10.json') else None
ch8 = J('fcs_chain8_maxdeg5_N6.json')
ch8b = J('fcs_chain8_maxdeg6_N4.json') if os.path.exists('fcs_chain8_maxdeg6_N4.json') else None
ch8c = J('fcs_chain8_maxdeg7_N4.json') if os.path.exists('fcs_chain8_maxdeg7_N4.json') else None

# ------------------------------------------------------------------ truncation_results.json
res = dict(rule="N* = argmin_{1<=N<=Nmax} |S_N - S_{N-1}| (the highest computed order if the differences are still decreasing); "
                "value S_{N*}; error estimate delta = |S_{N*} - S_{N*-1}|",
           kerr_cavity=kerr, ring_K3=dict(parameters=ring["parameters"], exact=ring["exact"], results=ring["results"],
                                          fig_ring_a=ring["fig_ring_a"], crosscheck=ring["crosscheck"]),
           table5_K1=k1, table5_K3=dict(maxdeg8=k3, maxdeg10=k3b), chain8=dict(maxdeg5=ch8, maxdeg6_check=ch8b, maxdeg7_check=ch8c),
           label_counts=dict(counts=ring["label_counts"], fitted_exponents=ring["label_count_fits"]))
json.dump(res, open('truncation_results.json', 'w'), indent=1)

# ------------------------------------------------------------------ fig_g3_errorbars.json (K=3 ring, Fig. 11)
g3 = dict(description="K=3 uniform Kerr ring (Delta=-1, J=0.4, eta=1, kappa=1): optimally truncated diagrammatic g2 and g3 with "
                      "successive-difference error bars, exact reference from time evolution at Fock cutoff 8. "
                      "partial_sums[N] is the value through order N (N=0..6).",
          instruction="Plot value with yerr=delta at each U; annotate N_star. The exact g3 column replaces any earlier reference.",
          U=[], g2=[], g2_delta=[], g2_Nstar=[], g2_exact=[], g3=[], g3_delta=[], g3_Nstar=[], g3_exact=[],
          g2_partial_sums=[], g3_partial_sums=[], g2_true_error=[], g3_true_error=[], n_partial_sums=[])
for key in sorted(ring["results"]):
    r = ring["results"][key]; g3["U"].append(r["U"])
    for q in ("g2", "g3"):
        g3[q].append(r[q]["value"]); g3[q+"_delta"].append(r[q]["delta"]); g3[q+"_Nstar"].append(r[q]["N_star"])
        g3[q+"_exact"].append(r[q]["exact"]); g3[q+"_partial_sums"].append(r[q]["partial_sums"]); g3[q+"_true_error"].append(r[q]["true_error"])
    g3["n_partial_sums"].append(r["n"]["partial_sums"])
json.dump(g3, open('fig_g3_errorbars.json', 'w'), indent=1)

# ------------------------------------------------------------------ fig_fcs_errorbars.json (eight-site chain, Fig. 12)
fcs = dict(description="Disordered eight-site Kerr chain (revision_params.json 'chain8_fcs'): counting cumulants of site 0 from "
                       "the tilted diagrams, orders 0..6 at maxdeg 5, optimally truncated with successive-difference error bars.",
           instruction="Plot fano and c3c1 versus U with yerr = *_delta; the mean-field line is 1. Order actually used is *_Nstar. "
                       "The partial sums (orders 0..6) are given so any fixed order can be drawn as well; maxdeg sensitivity is in maxdeg_check.",
           parameters=ch8["parameters"], U=[], c1=[], c1_delta=[], fano=[], fano_delta=[], fano_Nstar=[], c3c1=[], c3c1_delta=[], c3c1_Nstar=[],
           fano_partial_sums=[], c3c1_partial_sums=[], c1_partial_sums=[], maxdeg_check={})
for key in sorted(ch8["results"]):
    r = ch8["results"][key]; fcs["U"].append(r["U"])
    for q in ("c1", "fano", "c3c1"):
        fcs[q].append(r[q]["value"]); fcs[q+"_delta"].append(r[q]["delta"]); fcs[q+"_partial_sums"].append(r[q]["partial_sums"])
        if q != "c1": fcs[q+"_Nstar"].append(r[q]["N_star"])
for tag, d in (("maxdeg6_N4", ch8b), ("maxdeg7_N4", ch8c)):
    if d:
        fcs["maxdeg_check"][tag] = {k: dict(fano_partial_sums=v["fano"]["partial_sums"], c3c1_partial_sums=v["c3c1"]["partial_sums"], c1_partial_sums=v["c1"]["partial_sums"]) for k, v in d["results"].items()}
json.dump(fcs, open('fig_fcs_errorbars.json', 'w'), indent=1)

# ------------------------------------------------------------------ fig_kerr_convergence.pdf
from engine import destroy, steady
kp = kerr["parameters"]; kappa, Dc, eta = kp["kappa"], kp["Delta"], kp["eta"]
cn = np.array(kerr["coefficients"]); Nmax = kerr["Nmax"]
Us = np.linspace(0, 0.3, 61)
Nc = 32; A = destroy(Nc); Ad = A.conj().T
exact = np.array([np.trace(Ad@A@steady(Dc*Ad@A + 1j*(eta*Ad - np.conj(eta)*A) + U/2*Ad@Ad@A@A, [np.sqrt(kappa)*A])).real for U in Us])
mf = []
for U in Us:
    r = np.roots([U**2, 2*Dc*U, Dc**2 + kappa**2/4, -abs(eta)**2]) if U > 0 else [kp["alpha2"]]
    r = np.array(r); mf.append(np.min(r[np.abs(np.imag(r)) < 1e-9].real))
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(COLW, 4.3))
ax1.plot(Us, exact, "k-", lw=1.6, label="exact")
ax1.plot(Us, mf, color="gray", ls=(0, (4, 2)), lw=1.0, label="mean field ($L=0$)")
pw = Us[:, None]**np.arange(Nmax+1)[None, :]
for c, N in zip(ORDER_COLORS, [1, 2, 4, 8]):
    ax1.plot(Us, (pw[:, :N+1]*cn[None, :N+1]).sum(1), color=c, lw=0.9, label=f"$N\\leq{N}$")
# optimally truncated value on a coarser grid, with delta as error bar
Uo, So, Do = [], [], []
from truncation import optimal_truncation
for U in Us[2::4]:
    v, k, d, _ = optimal_truncation(np.cumsum(cn*U**np.arange(Nmax+1)))
    Uo.append(U); So.append(v); Do.append(d)
ax1.errorbar(Uo, So, yerr=Do, fmt="o", color="k", mfc="white", ms=3, capsize=1.5, lw=0.7, label="optimal truncation $\\pm\\delta$")
ax1.set_ylim(0.6, 3.2); ax1.set_xlim(0, Us[-1])
ax1.set_xlabel(r"$U/\kappa$"); ax1.set_ylabel(r"$\langle a^\dagger a\rangle$")
lg = ax1.legend(ncol=1, loc="upper left", bbox_to_anchor=(0.06, 1.0)); [l.set_linewidth(1.4) for l in lg.get_lines()]
ax1.text(0.97, 0.05, "(a)", transform=ax1.transAxes, ha="right", va="bottom")
for c, (key, rec) in zip(["#1b9e77", "#7570b3", "#d95f02", "#e7298a"], sorted(kerr["results"].items())):
    err = np.array(rec["true_error_all_orders"]); U = rec["U"]
    ax2.semilogy(np.arange(Nmax+1), err, "o-", color=c, ms=2.6, lw=0.8, label=f"$U/\\kappa={U}$")
    ks = rec["N_star"]
    ax2.semilogy([ks], [err[ks]], marker="*", color=c, ms=9, mec="k", mew=0.5, ls="none")
    ax2.semilogy(np.arange(1, Nmax+1), rec["differences"], ls=(0, (1, 1.5)), color=c, lw=0.7)
ax2.plot([], [], marker="*", color="0.5", mec="k", mew=0.5, ms=8, ls="none", label="$N^*$ (rule)")
ax2.plot([], [], ls=(0, (1, 1.5)), color="0.5", lw=0.7, label=r"$|S_N-S_{N-1}|$")
ax2.set_xlabel(r"truncation order $N$"); ax2.set_ylabel(r"$|\langle a^\dagger a\rangle_N-\langle a^\dagger a\rangle|$")
ax2.set_xlim(-0.3, Nmax+0.3); ax2.set_ylim(1e-12, 1e3); ax2.set_xticks(range(0, Nmax+1, 2))
ax2.legend(loc="lower left", ncol=2)
ax2.text(0.45, 0.95, "(b)", transform=ax2.transAxes, ha="left", va="top")
fig.tight_layout(pad=0.3, h_pad=0.6)
fig.savefig("fig_kerr_convergence.pdf")

# ------------------------------------------------------------------ timing_results.json + fig_ring.pdf
T = {}
for f in glob.glob('timing_*.json'):
    d = J(f); T[os.path.basename(f)[7:-5]] = d
def fitexp(Ks, ts):
    Ks, ts = np.array(Ks, float), np.array(ts, float)
    p = np.polyfit(np.log(Ks), np.log(ts), 1); return float(p[0]), float(np.exp(p[1]))
timing = dict(machine_note="all curves measured on one core (SLURM --cpus-per-task=1, OMP/OPENBLAS/MKL threads = 1) of the same machine",
              parameters=dict(kappa=1.0, Delta=-1.0, J=0.4, eta=1.0, U=0.15, observable="on-site g2, site 0 (n and n2 series)"), curves={}, fits={})
t2 = T.get("tensor2", {}).get("times", {}); t4 = T.get("tensor4", {}).get("times", {})
K2 = sorted(int(k) for k in t2); K4 = sorted(int(k) for k in t4)
timing["curves"]["tensor_N2"] = dict(K=K2, t=[t2[str(k)]["g2_total"] for k in K2], g2=[t2[str(k)]["g2"] for k in K2])
timing["curves"]["tensor_N4"] = dict(K=K4, t=[t4[str(k)]["g2_total"] for k in K4], g2=[t4[str(k)]["g2"] for k in K4])
mom = T.get("momentum", {}).get("times", {})
for key in ("n_unpruned_N2", "g2_unpruned_N2", "n_pruned_N2", "g2_pruned_N2", "g2_pruned_N4"):
    Ks = sorted(int(k) for k in mom if key in mom[k])
    if Ks: timing["curves"]["momentum_" + key] = dict(K=Ks, t=[mom[str(k)][key] for k in Ks])
ex = {}
for k in list(T):
    if k.startswith("exact_K"):
        for K, v in T[k]["times"].items(): ex[int(K)] = v
Ke = sorted(ex)
timing["curves"]["exact_evolve_nmax5"] = dict(K=Ke, t=[ex[k]["evolve"] for k in Ke], g2=[ex[k]["g2"] for k in Ke])
for name, cv in timing["curves"].items():
    if len(cv["K"]) >= 3 and name != "exact_evolve_nmax5":
        e, a = fitexp(cv["K"], cv["t"]); timing["fits"][name] = dict(exponent=e, prefactor=a, K_range=[min(cv["K"]), max(cv["K"])])
        # large-K exponent (upper half of the range) where the asymptotic power law is cleaner
        n = len(cv["K"]); 
        if n >= 5:
            e2, _ = fitexp(cv["K"][n//2-1:], cv["t"][n//2-1:]); timing["fits"][name]["exponent_upper_half"] = e2
if len(Ke) >= 2:
    e, a = fitexp(Ke, [ex[k]["evolve"] for k in Ke]); timing["fits"]["exact_evolve_nmax5_powerlaw_fit"] = dict(exponent=e, note="exponential in K, power-law fit only indicative")
    # exponential fit: log t = a + b K
    b = np.polyfit(Ke, np.log([ex[k]["evolve"] for k in Ke]), 1)[0]; timing["fits"]["exact_evolve_nmax5_exponential"] = dict(per_site_factor=float(np.exp(b)), nmax2=25.0)
timing["paper_values"] = dict(ring_timing_json=J(os.path.join(_ROOT, 'notebook', 'data', 'ring_timing.json')))
json.dump(timing, open('timing_results.json', 'w'), indent=1)

# fig_ring
Us13 = np.array(sorted(float(k) for k in ring["fig_ring_a"]))
fa = [ring["fig_ring_a"][f"{U:.3f}"]["g2"] for U in Us13]
PS = np.array([r["partial_sums"] for r in fa]); exG = np.array([r["exact"] for r in fa])
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(COLW, 4.3))
ax1.plot(Us13, exG, "k-", lw=1.6, label="exact")
for c, N, ls in zip(ORDER_COLORS, [1, 2, 4], [(0, (4, 2)), ":", "-"]):
    ax1.plot(Us13, PS[:, N], color=c, lw=1.0, ls=ls, label=f"$N\\leq{N}$")
ax1.errorbar(Us13, [r["value"] for r in fa], yerr=[r["delta"] for r in fa], fmt="o", color="k", mfc="white", ms=3, capsize=1.5, lw=0.7,
             label="optimal truncation $\\pm\\delta$")
for U, r in zip(Us13, fa):
    if U > 0: ax1.annotate(str(r["N_star"]), (U, r["value"]), textcoords="offset points", xytext=(0, -9), ha="center", fontsize=5.5, color="0.3")
ax1.set_xlabel(r"$U/\kappa$"); ax1.set_ylabel(r"on-site $g^{(2)}(0)$"); ax1.set_xlim(0, Us13[-1])
lg = ax1.legend(loc="upper left", fontsize=6.5); [l.set_linewidth(1.4) for l in lg.get_lines()]
ax1.text(0.97, 0.06, "(a) $K=3$", transform=ax1.transAxes, ha="right", va="bottom")
cv = timing["curves"]; ft = timing["fits"]
Kref = np.array([2, 3, 4, 5, 6])
if Ke:
    te = np.array(cv["exact_evolve_nmax5"]["t"])
    ax2.loglog(Kref, te[-1]*(5.0**(2*(Kref-Ke[-1]))), color="#d95f02", ls=":", lw=0.9, label=r"$\propto n_{\max}^{2K}$")
    ax2.loglog(Ke, te, "s", color="#d95f02", ms=3.5, label="exact steady state")
if "momentum_n_unpruned_N2" in cv:
    m = cv["momentum_n_unpruned_N2"]; e = ft["momentum_n_unpruned_N2"]["exponent"]
    ax2.loglog(m["K"], m["t"], "o-", color="#1b9e77", ms=3.5, lw=1.0, mfc="white", label=f"momentum engine, $N\\leq2$ ($K^{{{e:.1f}}}$)")
for key, col, mk, lab in (("tensor_N2", "#7570b3", "^", "tensor engine, $N\\leq2$"), ("tensor_N4", "#e7298a", "v", "tensor engine, $N\\leq4$")):
    if key in ft:
        e = ft[key].get("exponent_upper_half", ft[key]["exponent"])
        ax2.loglog(cv[key]["K"], cv[key]["t"], mk + "-", color=col, ms=3.5, lw=1.0, label=lab + f" ($K^{{{e:.1f}}}$)")
ax2.set_xlabel(r"ring size $K$"); ax2.set_ylabel("runtime (s)")
ax2.set_xlim(1.8, 36); ax2.set_ylim(1e-3, 1e6)
ax2.set_xticks([2, 3, 4, 6, 8, 12, 16, 24, 32]); ax2.set_xticklabels([str(k) for k in [2, 3, 4, 6, 8, 12, 16, 24, 32]])
ax2.tick_params(axis="x", which="minor", bottom=False, top=False)
lg = ax2.legend(loc="upper left", fontsize=6); [l.set_linewidth(1.4) for l in lg.get_lines()]
ax2.text(0.97, 0.06, "(b)", transform=ax2.transAxes, ha="right", va="bottom")
fig.tight_layout(pad=0.3, h_pad=0.5)
fig.savefig("fig_ring.pdf")
print(json.dumps(timing["fits"], indent=1))
print("figures written")
