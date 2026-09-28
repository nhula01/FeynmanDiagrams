import os
"""
Generate the numerical-validation figures for the worked examples.
Requires engine.py (diagram engine) in the same directory.
    python make_figures.py      ->  fig_kerr_convergence.pdf,
                                    fig_cubic_convergence.pdf,
                                    fig_emitter.pdf
"""
import numpy as np, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from engine import *

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 8,
    "axes.labelsize": 8, "legend.fontsize": 6.5, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "axes.linewidth": 0.6, "lines.linewidth": 1.0,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True, "legend.frameon": False,
})
COLW = 3.375
ORDER_COLORS = ["#d95f02", "#1b9e77", "#7570b3", "#e7298a", "#66a61e", "#a6761d"]
numbers = {}

I1 = np.eye(1)
def mode_ops(alpha, d=1, I=None):
    I = np.eye(d) if I is None else I
    b = cav(0, 1, I); bd = cav(1, 0, I); one = cav(0, 0, I)
    a = add(b, scale(alpha, one)); ad = add(bd, scale(np.conj(alpha), one))
    return a, ad, b, bd

# ============================================================ Example I: Kerr
def kerr():
    kappa, Dc, eta = 1.0, -1.0, 1.0
    z = kappa/2 + 1j*Dc; alpha = eta/z
    a, ad, b, bd = mode_ops(alpha)
    HK = scale(0.5, mul(mul(ad, ad), mul(a, a)))           # U = 1
    mod = Model(1, z, lambda X: 0*X, I1, HK, vertex_drop=4)
    Nmax = 16                                             # orders shown in Fig. 5
    Ncoeff = 40                                           # Sec. V.A root-test value quoted at N=40
    sn_all = mod.series(mul(ad, a), Ncoeff).real          # <a^dag a>^(N) / U^N
    sn = sn_all[:Nmax+1]
    Us = np.linspace(0, 0.3, 61)
    Nc = 32; A = destroy(Nc); Ad = A.conj().T
    exact = []
    for U in Us:
        H = Dc*Ad@A + 1j*(eta*Ad - np.conj(eta)*A) + U/2*Ad@Ad@A@A
        exact.append(np.trace(Ad@A@steady(H, [np.sqrt(kappa)*A])).real)
    exact = np.array(exact)
    # mean-field (tree-level, all orders): z ab + i U |ab|^2 ab = eta
    mf = []
    for U in Us:   # |abar|^2 from  n[(Dc+U n)^2 + kappa^2/4] = |eta|^2  (monostable here)
        r = np.roots([U**2, 2*Dc*U, Dc**2 + kappa**2/4, -abs(eta)**2]) if U > 0 else [abs(alpha)**2]
        r = np.array(r); mf.append(np.min(r[np.abs(np.imag(r)) < 1e-9].real))
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(COLW, 4.3))
    ax1.plot(Us, exact, "k-", lw=1.6, label="exact")
    ax1.plot(Us, mf, color="gray", ls=(0, (4, 2)), lw=1.0, label="mean field ($L=0$)")
    for c, N in zip(ORDER_COLORS, [1, 2, 4, 8]):
        pw = Us[:, None]**np.arange(Nmax+1)[None, :]
        ps = (pw[:, :N+1]*sn[None, :N+1]).sum(1)
        ax1.plot(Us, ps, color=c, lw=0.9, label=f"$N\\leq{N}$")
    ax1.set_ylim(0.6, 3.2); ax1.set_xlim(0, Us[-1])
    ax1.set_xlabel(r"$U/\kappa$"); ax1.set_ylabel(r"$\langle a^\dagger a\rangle$")
    lg = ax1.legend(ncol=1, loc="upper left", bbox_to_anchor=(0.06, 1.0))
    [l.set_linewidth(1.4) for l in lg.get_lines()]
    ax1.text(0.97, 0.05, "(a)", transform=ax1.transAxes, ha="right", va="bottom")
    # (b) error versus truncation order
    for c, U in zip(["#1b9e77", "#7570b3", "#d95f02", "#e7298a"], [0.05, 0.1, 0.2, 0.3]):
        H = Dc*Ad@A + 1j*(eta*Ad - np.conj(eta)*A) + U/2*Ad@Ad@A@A
        ex = np.trace(Ad@A@steady(H, [np.sqrt(kappa)*A])).real
        ps = np.cumsum(sn*U**np.arange(Nmax+1))
        err = np.abs(ps - ex)
        ax2.semilogy(np.arange(Nmax+1), err, "o-", color=c, ms=2.6, lw=0.8, label=f"$U/\\kappa={U}$")
        numbers[f"kerr_err_U{U}"] = [float(e) for e in err]
    ax2.set_xlabel(r"truncation order $N$"); ax2.set_ylabel(r"$|\langle a^\dagger a\rangle_N-\langle a^\dagger a\rangle|$")
    ax2.set_xlim(-0.3, Nmax+0.3); ax2.set_ylim(1e-12, 1e3); ax2.set_xticks(range(0, Nmax+1, 2))
    ax2.legend(loc="lower left", ncol=2)
    ax2.text(0.45, 0.95, "(b)", transform=ax2.transAxes, ha="left", va="top")
    fig.tight_layout(pad=0.3, h_pad=0.6)
    fig.savefig("fig_kerr_convergence.pdf")
    numbers["kerr_coeffs"] = [float(x) for x in sn_all]
    numbers["kerr_exact_vs_mf"] = {f"{U:.2f}": [float(e), float(m)] for U, e, m in zip(Us[::10], exact[::10], np.array(mf)[::10])}
    numbers["kerr_alpha"] = float(abs(alpha))

# ================================================== Example II: cubic (APS)
def cubic():
    from sparse_ss import steady_sparse
    gam, om, f = 2.108, 2.153, 1.17+4.617j
    z = gam/2 + 1j*om; eta = f/2; alpha = eta/z
    a, ad, b, bd = mode_ops(alpha)
    Hc = scale(-0.5j, add(mul(mul(ad, ad), a), scale(-1, mul(ad, mul(a, a)))))   # g = 1
    mod = Model(1, z, lambda X: 0*X, I1, Hc, vertex_drop=3)
    Nmax = 60
    sa = mod.series(a, Nmax)
    gs = np.linspace(0, 0.8, 41)
    ex = {}
    for Nc in (40, 80):
        A = destroy(Nc); Ad = A.conj().T; vals = []
        for g in gs:
            H = om*Ad@A + 1j*(eta*Ad - np.conj(eta)*A) + g*(-0.5j)*(Ad@Ad@A - Ad@A@A)
            vals.append(np.trace(A@steady_sparse(H, [np.sqrt(gam)*A])))
        ex[Nc] = np.array(vals)
    diff = np.abs(ex[40] - ex[80])
    gbad = float(gs[np.argmax(diff > 1e-3)])
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(COLW, 4.3), sharex=True)
    for ax, fn, lab, lab2 in [(ax1, np.real, r"$\mathrm{Re}\,\langle a\rangle$", "(a)"),
                              (ax2, np.imag, r"$\mathrm{Im}\,\langle a\rangle$", "(b)")]:
        ax.axvspan(gbad, gs[-1], color="0.92", lw=0)
        ax.plot(gs, fn(ex[80]), "k-", lw=1.6, label="exact, $n_{\\max}=80$")
        ax.plot(gs, fn(ex[40]), color="0.45", ls=(0, (1.5, 1.5)), lw=1.2, label="exact, $n_{\\max}=40$")
        for c, N in zip(ORDER_COLORS, [1, 2, 3, 6, 12, 24]):
            ps = np.array([np.sum(sa[:N+1]*g**np.arange(N+1)) for g in gs])
            ax.plot(gs, fn(ps), color=c, lw=0.9, label=f"$N\\leq{N}$")
        ax.set_ylabel(lab)
        ax.text(0.03, 0.93, lab2, transform=ax.transAxes, ha="left", va="top")
    ax1.set_ylim(0.6, 1.1); ax2.set_ylim(0.15, 1.0)
    ax2.set_xlim(0, gs[-1]); ax2.set_xlabel(r"coupling $g$")
    lg = ax1.legend(ncol=2, loc="lower left")
    [l.set_linewidth(1.4) for l in lg.get_lines()]
    ax2.text(gbad+0.01, 0.2, "cutoff\ndependent", fontsize=6.5, color="0.35")
    fig.tight_layout(pad=0.3, h_pad=0.4)
    fig.savefig("fig_cubic_convergence.pdf")
    numbers["cubic_gbad"] = gbad
    numbers["cubic_root_test"] = [float(abs(sa[N])**(-1/N)) for N in (10, 20, 30, 40, 50, 60)]
    numbers["cubic_first"] = [complex(sa[1]).real, complex(sa[1]).imag]
    numbers["cubic_alpha"] = [alpha.real, alpha.imag]
    for gg in [0.3, 0.5, 0.6]:
        k = int(round(gg/0.02))
        e = ex[80][k]
        numbers[f"cubic_err_g{gg}"] = [float(abs(np.sum(sa[:N+1]*gg**np.arange(N+1)) - e)) for N in (1, 3, 6, 12, 24)]
        numbers[f"cubic_cutoffdiff_g{gg}"] = float(diff[k])

# ============================================ Example III: emitter in cavity
sm = np.array([[0, 0], [1, 0]], complex); sp_ = sm.conj().T
sz = np.diag([1., -1.]).astype(complex); sx = sm + sp_; I2 = np.eye(2)
def atom(Da, Om, gam):
    Ha = Da/2*sz + Om/2*sx
    La = lambda X: 1j*(Ha@X - X@Ha) + gam*(sp_@X@sm - 0.5*(sp_@sm@X + X@sp_@sm))
    return La, steady(Ha, [np.sqrt(gam)*sm])
def em_series(kap, Dc, Da, Om, gam, g, Nmax, O, dressed=True):
    z = kap/2 + 1j*Dc
    A = cav(0, 1, I2); Ad = cav(1, 0, I2)
    Hx = scale(g, add(mul(Ad, cav(0, 0, sm)), mul(A, cav(0, 0, sp_))))
    if dressed:
        La, rho = atom(Da, Om, gam); Hint = Hx
    else:
        La, rho = atom(Da, 0.0, gam); Hint = add(Hx, cav(0, 0, Om/2*sx))
    return Model(2, z, La, rho, Hint, vertex_drop=2).series(O, Nmax)
def em_exact(kap, Dc, Da, Om, gam, g, Nc=8):
    a = np.kron(destroy(Nc), I2); ad = a.conj().T
    S = np.kron(np.eye(Nc), sm); Sd = S.conj().T
    H = Dc*ad@a + np.kron(np.eye(Nc), Da/2*sz + Om/2*sx) + g*(ad@S + a@Sd)
    rho = steady(H, [np.sqrt(kap)*a, np.sqrt(gam)*S])
    return np.trace(ad@a@rho).real
def dressed_modes(Da, Om, gam):
    """eigen-decomposition of the emitter Heisenberg generator; returns
    lambdas (3), eigenoperators X_k, coefficients of dsigma^- and <dsigma^+ X_k>"""
    La, rho = atom(Da, Om, gam)
    M = np.zeros((4, 4), complex)
    for j in range(4):
        E = np.zeros(4, complex); E[j] = 1
        M[:, j] = La(E.reshape(2, 2, order="F")).reshape(-1, order="F")
    w, V = np.linalg.eig(M)
    idx = [i for i in range(4) if abs(w[i]) > 1e-9]
    lam = w[idx]; X = [V[:, i].reshape(2, 2, order="F") for i in idx]
    s_m = np.trace(sm@rho); dsm = sm - s_m*I2
    coef = np.linalg.lstsq(np.array([x.reshape(-1) for x in X]).T, dsm.reshape(-1), rcond=None)[0]
    dsp = sp_ - np.conj(s_m)*I2
    corr = np.array([np.trace(dsp@x@rho) for x in X])
    return lam, coef, corr, s_m
def emitter():
    gam, kap, Da = 1.0, 0.5, 0.0
    # (a) cavity-filtered Mollow triplet
    Om, g = 5.0, 0.2
    Dcs = np.linspace(-8, 8, 161)
    ex = np.array([em_exact(kap, D, Da, Om, gam, g) for D in Dcs])
    n2 = np.array([np.sum(em_series(kap, D, Da, Om, gam, g, 2, cav(1, 1, I2)))for D in Dcs]).real
    n4 = np.array([np.sum(em_series(kap, D, Da, Om, gam, g, 4, cav(1, 1, I2)))for D in Dcs]).real
    lam, coef, corr, s_m = dressed_modes(Da, Om, gam)
    contrib = {}
    for D in Dcs:
        z = kap/2 + 1j*D
        for k in range(3):
            contrib.setdefault(k, []).append((2*g**2/kap*coef[k]*corr[k]/(np.conj(z) - lam[k])).real)
        contrib.setdefault("coh", []).append(g**2*abs(s_m)**2/abs(z)**2)
    contrib = {k: np.array(v) for k, v in contrib.items()}
    chk = contrib[0] + contrib[1] + contrib[2] + contrib["coh"]
    numbers["emitter_order2_check"] = float(np.max(np.abs(chk - n2)))
    order = np.argsort(lam.imag)   # mode 2: Im<0 ; mode 1: real ; mode 3: Im>0
    lab_of = {order[1]: 1, order[0]: 2, order[2]: 3}
    numbers["emitter_lambdas"] = {str(lab_of[k]): [lam[k].real, lam[k].imag] for k in range(3)}
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(COLW, 4.4))
    fills = {1: "#fdd49e", 2: "#c6dbef", 3: "#c7e9c0"}
    for k in range(3):
        ax1.fill_between(Dcs, 0, 1e3*contrib[k], color=fills[lab_of[k]], lw=0, alpha=0.9,
                         label=f"dressed line {lab_of[k]}")
    ax1.plot(Dcs, 1e3*ex, "k-", lw=1.6, label="exact")
    ax1.plot(Dcs, 1e3*n2, color="#d95f02", ls=(0, (4, 2)), lw=1.0, label="dressed, $N\\leq2$")
    ax1.plot(Dcs, 1e3*n4, color="#1b9e77", ls=":", lw=1.2, label="dressed, $N\\leq4$")
    ax1.set_xlim(Dcs[0], Dcs[-1]); ax1.set_ylim(0, None)
    ax1.set_xlabel(r"cavity detuning $\Delta_c/\gamma$")
    ax1.set_ylabel(r"$10^3\,\langle a^\dagger a\rangle$")
    ax1.legend(loc="upper right", fontsize=6)
    ax1.text(0.03, 0.93, "(a)", transform=ax1.transAxes, ha="left", va="top")
    numbers["emitter_a_maxrel_err_N2"] = float(np.max(np.abs(n2-ex)/ex))
    numbers["emitter_a_maxrel_err_N4"] = float(np.max(np.abs(n4-ex)/ex))
    # (b) bare versus dressed as the drive grows, normalized by a classical dipole
    g = 0.1; Oms = np.linspace(0.02, 2.0, 100)
    def ncl(O):
        z = kap/2; za = gam/2
        Mx = np.array([[z, 1j*g], [1j*g, za]]); rhs = np.array([0, -1j*O/2])
        av, cv = np.linalg.solve(Mx, rhs)
        return abs(av)**2
    nc = np.array([ncl(O) for O in Oms])
    exb = np.array([em_exact(kap, 0.0, Da, O, gam, g) for O in Oms])/nc
    nd = np.array([np.sum(em_series(kap, 0.0, Da, O, gam, g, 2, cav(1, 1, I2))) for O in Oms]).real/nc
    nd4 = np.array([np.sum(em_series(kap, 0.0, Da, O, gam, g, 4, cav(1, 1, I2))) for O in Oms]).real/nc
    Nb = 16
    sb = np.array([np.cumsum(em_series(kap, 0.0, Da, O, gam, g, Nb, cav(1, 1, I2), dressed=False)).real for O in Oms])/nc[:, None]
    Omc = np.sqrt(2*(Da**2 + gam**2/4))
    ax2.axvspan(Omc, Oms[-1], color="0.92", lw=0)
    ax2.axhline(1.0, color="0.5", ls=(0, (4, 2)), lw=1.3, label="classical dipole")
    ax2.plot(Oms, exb, "k-", lw=1.6, label="atom, exact")
    for c, N in zip(["#7570b3", "#e7298a", "#66a61e", "#a6761d"], [4, 8, 12, 16]):
        ax2.plot(Oms, sb[:, N], color=c, lw=0.9, label=f"bare, $N\\leq{N}$")
    ax2.plot(Oms, nd, color="#d95f02", ls=(0, (4, 2)), lw=1.1, label="dressed, $N\\leq2$")
    ax2.plot(Oms, nd4, color="#1b9e77", ls=":", lw=1.3, label="dressed, $N\\leq4$")
    ax2.axvline(Omc, color="0.5", lw=0.6, ls=":")
    ax2.set_xlim(0, Oms[-1]); ax2.set_ylim(0, 1.25)
    ax2.set_xlabel(r"Rabi frequency $\Omega/\gamma$"); ax2.set_ylabel(r"$\langle a^\dagger a\rangle/\langle a^\dagger a\rangle_{\mathrm{cl}}$")
    lg = ax2.legend(loc="center right", ncol=2, fontsize=6, bbox_to_anchor=(1.0, 0.62))
    [l.set_linewidth(1.4) for l in lg.get_lines()]
    ax2.text(0.03, 0.08, "(b)", transform=ax2.transAxes, ha="left", va="bottom")
    ax2.text(Omc+0.03, 0.05, r"$\Omega>\gamma/\sqrt{2}$", fontsize=7, color="0.35")
    numbers["emitter_b_dressed_maxrel_err"] = float(np.max(np.abs(nd-exb)/exb))
    numbers["emitter_b_dressed4_maxrel_err"] = float(np.max(np.abs(nd4-exb)/exb))
    numbers["emitter_b_ratio_at_Om2"] = float(exb[-1])
    fig.tight_layout(pad=0.3, h_pad=0.6)
    fig.savefig("fig_emitter.pdf")
    numbers["emitter_b_dressed_maxrel_err"] = float(np.max(np.abs(nd-exb)/exb))
    numbers["emitter_b_dressed4_maxrel_err"] = float(np.max(np.abs(nd4-exb)/exb))
    numbers["emitter_a_params"] = dict(gamma=gam, kappa=kap, Omega=5.0, g=0.2)

def finite():
    from finite_time import kerr_finite_time, kerr_exact_time
    kap, Dc, eta, U = 1.0, -1.0, 1.0, 0.1
    ts = np.linspace(0, 10, 101)
    ser, parts = kerr_finite_time(kap, Dc, eta, U, ts, 8, lambda a, ad: mul(ad, a), by_degree=True)
    ex = kerr_exact_time(kap, Dc, eta, U, ts)
    ps = np.cumsum(ser.real, axis=1)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(COLW, 4.3), sharex=True)
    ax1.plot(ts, ex, "k-", lw=1.6, label="exact")
    for c, N in zip(ORDER_COLORS, [0, 1, 2, 4, 8]):
        ax1.plot(ts, ps[:, N], color=c, lw=0.9, label=f"$N\\leq{N}$")
    ax1.set_ylabel(r"$\langle a^\dagger a\rangle(t)$"); ax1.set_ylim(0, 1.5)
    lg = ax1.legend(ncol=2, loc="lower right"); [l.set_linewidth(1.4) for l in lg.get_lines()]
    ax1.text(0.03, 0.93, "(a)", transform=ax1.transAxes, ha="left", va="top")
    cols = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"]
    for c, dgr in zip(cols, range(5)):
        if dgr in parts[1]:
            lab = "no outgoing line ($r_{00}$)" if dgr == 0 else (f"{dgr} outgoing line" + ("s" if dgr > 1 else ""))
            ax2.plot(ts, parts[1][dgr].real, color=c, lw=1.0, label=lab)
    ax2.plot(ts, ser[:, 1].real, "k--", lw=1.0, label="total, $N=1$")
    ax2.axhline(0, color="0.7", lw=0.5)
    ax2.set_xlabel(r"$\kappa t$"); ax2.set_ylabel(r"$\langle a^\dagger a\rangle^{(1)}(t)$")
    lg = ax2.legend(ncol=2, loc="lower left", fontsize=6); [l.set_linewidth(1.4) for l in lg.get_lines()]
    ax2.text(0.03, 0.93, "(b)", transform=ax2.transAxes, ha="left", va="top")
    ax2.set_xlim(0, ts[-1])
    fig.tight_layout(pad=0.3, h_pad=0.4)
    fig.savefig("fig_kerr_finite_time.pdf")
    k = np.argmin(abs(ts-5)); numbers["finite_t5"] = [float(ex[k])] + [float(ps[k, N]) for N in (0, 1, 2, 4, 8)]
    numbers["finite_maxerr_N8"] = float(np.max(np.abs(ps[:, 8]-ex)))
    numbers["finite_maxerr_N4"] = float(np.max(np.abs(ps[:, 4]-ex)))
    numbers["finite_t10"] = [float(ex[-1]), float(ps[-1, 8])]

def atomcavity():
    from atom_cavity import series_cavity, exact_cavity
    kap, gam, g, eta = 1.0, 1.0, 0.3, 0.02
    obs = {"n": lambda a, ad: mul(ad, a), "n2": lambda a, ad: mul(mul(ad, ad), mul(a, a))}
    Ds = np.linspace(-2.0, 2.0, 81)
    orders = [2, 4, 8, 16, 32]
    Nmax = max(orders)
    exn, exg2, pn, pg2 = [], [], [], []
    for D in Ds:
        sa, al = series_cavity("atom", kap, gam, D, D, eta, g, Nmax, obs)
        n, n2, _ = exact_cavity("atom", kap, gam, D, D, eta, g)
        cn = np.cumsum(sa["n"].real); cn2 = np.cumsum(sa["n2"].real)
        exn.append(n/abs(al)**2); exg2.append(n2/n**2)
        pn.append(cn/abs(al)**2); pg2.append(cn2/cn**2)
    exn, exg2, pn, pg2 = map(np.array, (exn, exg2, pn, pg2))
    z = kap/2 + 1j*Ds; za = gam/2 + 1j*Ds
    Tcl = np.abs(z*za/(z*za + g**2))**2
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(COLW, 4.3), sharex=True)
    ax1.plot(Ds, exn, "k-", lw=1.6, label="atom, exact")
    ax1.plot(Ds, Tcl, color="0.5", ls=(0, (4, 2)), lw=1.3, label="classical dipole")
    for c, N in zip(ORDER_COLORS, [2, 4, 8]):
        ax1.plot(Ds, pn[:, N], color=c, lw=0.9, label=f"atom, $N\\leq{N}$")
    ax1.set_ylim(0, 1.05); ax1.set_ylabel(r"$\langle a^\dagger a\rangle/|\alpha|^2$")
    lg = ax1.legend(loc="lower right", fontsize=6); [l.set_linewidth(1.4) for l in lg.get_lines()]
    ax1.text(0.03, 0.93, "(a)", transform=ax1.transAxes, ha="left", va="top")
    ax2.axhline(1.0, color="0.5", ls=(0, (4, 2)), lw=1.3, label="classical dipole")
    ax2.plot(Ds, exg2, "k-", lw=1.6, label="atom, exact")
    for c, N in zip(ORDER_COLORS[1:], [4, 8, 16, 32]):
        ax2.plot(Ds, pg2[:, N], color=c, lw=0.9, label=f"atom, $N\\leq{N}$")
    ax2.set_ylim(0.86, 1.06); ax2.set_xlim(Ds[0], Ds[-1])
    ax2.set_xlabel(r"laser detuning $\Delta/\gamma$"); ax2.set_ylabel(r"$g^{(2)}(0)$")
    lg = ax2.legend(loc="lower right", ncol=2, fontsize=6); [l.set_linewidth(1.4) for l in lg.get_lines()]
    ax2.text(0.03, 0.93, "(b)", transform=ax2.transAxes, ha="left", va="top")
    fig.tight_layout(pad=0.3, h_pad=0.4)
    fig.savefig("fig_atom_cavity.pdf")
    k0 = np.argmin(abs(Ds)); k3 = np.argmin(abs(Ds-0.3))
    numbers["atomcav"] = dict(g2_res=float(exg2[k0]), g2_res_N32=float(pg2[k0, 32]), g2_res_N16=float(pg2[k0, 16]),
                              g2_max=float(exg2.max()), D_g2max=float(Ds[np.argmax(exg2)]),
                              maxdev_n_atom_cl=float(np.max(np.abs(exn/Tcl-1))),
                              maxerr_g2_N32=float(np.max(np.abs(pg2[:, 32]-exg2))),
                              maxerr_g2_N16_offres=float(np.max(np.abs(pg2[np.abs(Ds) > 0.5, 16]-exg2[np.abs(Ds) > 0.5]))))

def ring():
    """Driven-dissipative Kerr ring: diagrams versus exact numerics."""
    import json, time, sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import multimode as MM
    P = dict(kappa=1.0, Delta=-1.0, J=0.4, eta=1.0)
    # (a) accuracy versus U at K = 3, against a converged exact steady state
    Us = np.linspace(0, 0.3, 13)
    ex_U = json.load(open("ring_exactU.json"))
    n_ser, g2_ser = [], []
    for U in Us:
        r = MM.Ring(K=3, U=U, **P)
        a0, ad0 = r.site_ops(0)
        n = np.cumsum(r.series(MM.mul(ad0, a0), 2)).real
        n2 = np.cumsum(r.series(MM.mul(MM.mul(ad0, ad0), MM.mul(a0, a0)), 2)).real
        n_ser.append(n); g2_ser.append(n2/n**2)
    n_ser = np.array(n_ser); g2_ser = np.array(g2_ser)
    exU = np.array([ex_U[f"{U:.3f}"]["n0"] for U in Us])
    exG = np.array([ex_U[f"{U:.3f}"]["n0n0"]/ex_U[f"{U:.3f}"]["n0"]**2 for U in Us])
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(COLW, 4.3))
    ax1.plot(Us, exG, "k-", lw=1.6, label="exact")
    for c, N, ls in zip(ORDER_COLORS, [0, 1, 2], ["-", (0, (4, 2)), ":"]):
        ax1.plot(Us, g2_ser[:, N], color=c, lw=1.0, ls=ls, label=f"$N\\leq{N}$")
    ax1.set_xlabel(r"$U/\kappa$"); ax1.set_ylabel(r"on-site $g^{(2)}(0)$")
    ax1.set_xlim(0, Us[-1])
    lg = ax1.legend(loc="upper left", fontsize=6.5); [l.set_linewidth(1.4) for l in lg.get_lines()]
    ax1.text(0.97, 0.06, "(a) $K=3$", transform=ax1.transAxes, ha="right", va="bottom")
    tim = json.load(open("ring_timing.json"))
    Kd = np.array(sorted(int(k) for k in tim["diagrams"]))
    td = np.array([tim["diagrams"][str(k)] for k in Kd])
    Ke = np.array(sorted(int(k) for k in tim["exact"]))
    te = np.array([tim["exact"][str(k)] for k in Ke])
    Kall = np.arange(2, 11)
    ax2.semilogy(Kall, te[-1]*(5.0**(2*(Kall-Ke[-1]))), color="#d95f02", ls=":", lw=0.9,
                 label=r"$\propto n_{\max}^{2K}$")
    ax2.semilogy(Ke, te, "s-", color="#d95f02", ms=3.5, lw=1.0, label="exact steady state")
    ax2.semilogy(Kd, td, "o-", color="#1b9e77", ms=3.5, lw=1.0, label="diagrams, $N\\leq2$")
    ax2.set_xlabel(r"ring size $K$"); ax2.set_ylabel("runtime (s)")
    ax2.set_ylim(1e-1, 1e9); ax2.set_xlim(1.7, 10.3); ax2.set_xticks(range(2, 11))
    lg = ax2.legend(loc="upper left", fontsize=6.5); [l.set_linewidth(1.4) for l in lg.get_lines()]
    ax2.text(0.97, 0.06, "(b)", transform=ax2.transAxes, ha="right", va="bottom")
    fig.tight_layout(pad=0.3, h_pad=0.5)
    fig.savefig("fig_ring.pdf")
    numbers["ring"] = dict(g2_err_N2=float(np.max(np.abs(g2_ser[:, 2]-exG))),
                           g2_err_N1=float(np.max(np.abs(g2_ser[:, 1]-exG))),
                           n_err_N2=float(np.max(np.abs(n_ser[:, 2]-exU))))

if __name__ == "__main__":
    import sys
    todo = sys.argv[1:] or ["kerr", "cubic", "emitter"]
    for name in todo:
        globals()[name]()
    old = json.load(open("numbers.json")) if __import__("os").path.exists("numbers.json") else {}
    old.update(numbers)
    json.dump(old, open("numbers.json", "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in numbers.items() if "coeffs" not in k}, indent=1, default=str))
