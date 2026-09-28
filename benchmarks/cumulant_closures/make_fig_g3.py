"""Revised Fig. 11 (fig_g3.pdf): K=3 ring g3(0), closures of order 2,3,4 vs diagrams N<=2,4 vs exact.
Reads cumulant_closure_results.json; if ../truncation/fig_g3_errorbars.json exists, overlays the optimally
truncated diagrams with successive-difference error bars in panel (a)."""
import os, json, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 8,
    "axes.labelsize": 8, "legend.fontsize": 6.5, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "axes.linewidth": 0.6, "lines.linewidth": 1.0,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True, "legend.frameon": False,
})
COLW = 3.375
R = json.load(open(os.path.join(HERE, 'cumulant_closure_results.json')))
tab = R['table_fig11']; U = np.array([r['U'] for r in tab])
nv = R['variables_per_order']
def tmed(m):
    ts = [r[m]['t'] for r in tab if m in r and r[m].get('t') is not None]
    return np.median(ts) if ts else np.nan
def tfmt(t):
    return f"{t:.2f} s" if t < 1 else (f"{t:.1f} s" if t < 10 else f"{t:.0f} s")
# colour threading: orange/purple/green as in the original figure; new closures take the next Dark2 colours
S = {'closure2': dict(c="#d95f02", ls="--", mk="s", lab="2nd-order cumulant"),
     'closure3': dict(c="#e7298a", ls=(0, (5, 1.5, 1, 1.5)), mk="D", lab="3rd-order cumulant"),
     'closure4': dict(c="#a6761d", ls=(0, (1.5, 1)), mk="p", lab="4th-order cumulant"),
     'diagrams2': dict(c="#7570b3", ls=":", mk="o", lab=r"diagrams $N\leq2$"),
     'diagrams4': dict(c="#1b9e77", ls="-", mk="o", lab=r"diagrams $N\leq4$")}
for m in ('closure2', 'closure3', 'closure4'):
    S[m]['lab'] += f" ({nv[m[-1]]} var., {tfmt(tmed(m))})"
for m in ('diagrams2', 'diagrams4'):
    S[m]['lab'] += f" ({tfmt(tmed(m))})"
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(COLW, 4.3), sharex=True)
exg3 = np.array([r['exact_g3'] for r in tab])
ax1.plot(U, exg3, "k-", lw=1.8, label="exact", zorder=5)
for m in ('closure2', 'closure3', 'closure4', 'diagrams2', 'diagrams4'):
    s = S[m]; y = [r[m]['g3'] for r in tab]
    ax1.plot(U, y, color=s['c'], ls=s['ls'], marker=s['mk'], ms=3.2 if s['mk'] != 'p' else 3.8, lw=1.0 if m != 'diagrams4' else 1.3,
             label=s['lab'])
eb = os.path.join(HERE, '..', 'truncation', 'fig_g3_errorbars.json')
if os.path.exists(eb):
    E = json.load(open(eb)); Ue = np.array(E['U']); sel = [i for i, u in enumerate(Ue) if any(abs(u-x) < 1e-9 for x in U)]
    if sel:
        ax1.errorbar(Ue[sel], np.array(E['g3'])[sel], yerr=np.array(E['g3_delta'])[sel], fmt="o", color="k", mfc="white", ms=3.2,
                     capsize=1.5, lw=0.7, zorder=6, label=r"diagrams, optimal truncation $\pm\delta$")
ax1.set_ylabel(r"$g^{(3)}(0)$")
lg = ax1.legend(loc="upper left", fontsize=5.8, handlelength=2.6, borderaxespad=0.3, labelspacing=0.25)
[l.set_linewidth(1.3) for l in lg.get_lines()]
ax1.text(0.97, 0.05, "(a) $K=3$", transform=ax1.transAxes, ha="right", va="bottom")
for m in ('closure2', 'closure3', 'closure4', 'diagrams4'):
    s = S[m]
    ax2.semilogy(U, [r[m]['err_g3'] for r in tab], color=s['c'], ls=s['ls'], marker=s['mk'], ms=3.4, lw=1.1)
    ax2.semilogy(U, [r[m]['err_g2'] for r in tab], color=s['c'], ls=s['ls'], marker=s['mk'], ms=3.0, lw=0.7, mfc="white")
ax2.set_xlabel(r"$U/\kappa$"); ax2.set_ylabel("relative error")
h = [Line2D([], [], color="0.3", marker="o", ms=3.4, ls="-", lw=1.1, label=r"$g^{(3)}$"),
     Line2D([], [], color="0.3", marker="o", ms=3.0, ls="-", lw=0.7, mfc="white", label=r"$g^{(2)}$")]
ax2.legend(handles=h, loc="lower right", fontsize=6.5, ncol=2, handlelength=2.2, columnspacing=1.2)
lo = min(min(r[m]['err_g2'] for r in tab for m in ('closure2', 'closure3', 'closure4', 'diagrams4')), 1e-4)
ax2.set_ylim(10**np.floor(np.log10(lo)-0.4), 0.3)
ax2.text(0.03, 0.95, "(b)", transform=ax2.transAxes, ha="left", va="top")
ax2.set_xlim(U[0]-0.008, U[-1]+0.008); ax2.set_xticks(U)
fig.tight_layout(pad=0.3, h_pad=0.4)
fig.savefig(os.path.join(HERE, "fig_g3.pdf")); fig.savefig(os.path.join(HERE, "fig_g3_preview.png"), dpi=250)
print("wrote fig_g3.pdf")
