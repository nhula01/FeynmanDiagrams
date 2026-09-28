"""Table V and the eight-site chain: counting cumulants order by order, optimal truncation.
usage: fcs_trunc.py [k1|k3|chain8|chain8check]"""
import sys, json, os, time, numpy as np
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
sys.path.insert(0, os.path.join(_ROOT, 'engines'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from truncation import truncation_record
RP = json.load(open(os.path.join(_ROOT, 'notebook', 'data', 'revision_params.json')))
what = sys.argv[1]
H = 0.05                                   # finite-difference step in chi (same as the paper)
CHIS = [-2*H, -H, 0.0, H, 2*H]

def cumulant_partial_sums(theta_table):
    """theta_table[c] = array of theta partial sums (orders 0..Nmax) at chi = CHIS[c];
    returns c1, c2/c1, c3/c1 partial-sum sequences (order by order)"""
    th = np.array(theta_table).real          # (5, Nmax+1)
    c1 = (th[3]-th[1])/(2*H); c2 = (th[3]-2*th[2]+th[1])/H**2
    c3 = (th[4]-2*th[3]+2*th[1]-th[0])/(2*H**3)
    return c1, c2/c1, c3/c1

def records(c1, f, c3, exact=None, **extra):
    ex = exact or {}
    return dict(c1=truncation_record(c1, exact=ex.get("c1")), fano=truncation_record(f, exact=ex.get("fano")),
                c3c1=truncation_record(c3, exact=ex.get("c3c1")), **extra)

if what == "k1":
    from fcs import theta_series
    from tilted_exact import theta_lead, cumulants as cum_exact
    from engine import destroy, steady
    P = RP["table5_fcs"]["K1"]; kappa, Dc, eta = 1.0, P["Delta"], P["eta"]; Nmax = P["Nmax"]
    out = dict(parameters=dict(kappa=kappa, Delta=Dc, eta=eta, h=H, Nmax=Nmax, maxdeg=14), results={})
    for U in P["U"]:
        t = time.time()
        table = [np.cumsum(theta_series(kappa, Dc, eta, U, x, Nmax, maxdeg=14)) for x in CHIS]
        table18 = [np.cumsum(theta_series(kappa, Dc, eta, U, x, Nmax, maxdeg=18)) for x in CHIS]
        c1, f, c3 = cumulant_partial_sums(table); c1b, fb, c3b = cumulant_partial_sums(table18)
        exact = {}
        for Nc in (26, 32):
            A = destroy(Nc); Ad = A.conj().T
            Hx = Dc*Ad@A + 1j*(eta*Ad - np.conj(eta)*A) + U/2*Ad@Ad@A@A
            ce = cum_exact(lambda x: theta_lead(Hx, [A], [kappa], [x]), h=H)
            exact[Nc] = dict(c1=ce[0], fano=ce[1], c3c1=ce[2])
        rec = records(c1, f, c3, exact=exact[32], U=U, exact_cutoff_change={k: abs(exact[26][k]-exact[32][k]) for k in exact[32]},
                      maxdeg_change_at_Nmax=dict(c1=abs(c1[-1]-c1b[-1]), fano=abs(f[-1]-fb[-1]), c3c1=abs(c3[-1]-c3b[-1])))
        out["results"][f"{U:.2f}"] = rec
        for k in ("c1", "fano", "c3c1"):
            r = rec[k]; print(f"K=1 U={U} {k}: N*={r['N_star']} {r['value']:.5f}+-{r['delta']:.1e} true {r['true_error']:.1e} exact {r['exact']:.5f}  S_N={np.round(r['partial_sums'],5)}", flush=True)
        print(f"  [{time.time()-t:.0f}s]", flush=True)
    json.dump(out, open('fcs_k1_truncation.json', 'w'), indent=1)

elif what == "k3":
    from fcs_multi import TiltedChain
    P = RP["table5_fcs"]["K3"]; U = P["U"]; Nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    MAXDEG = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    det3 = [P["Delta"]]*3; eta3 = np.zeros(3); eta3[P["eta_site"]] = 1.0
    k3 = json.load(open(os.path.join(_ROOT, 'notebook', 'data', 'k3_theta.json')))
    th = np.array([k3[f"{x:.3f}"] for x in CHIS])
    c1e = (th[3]-th[1])/(2*H); c2e = (th[3]-2*th[2]+th[1])/H**2; c3e = (th[4]-2*th[3]+2*th[1]-th[0])/(2*H**3)
    exact = dict(c1=c1e, fano=c2e/c1e, c3c1=c3e/c1e)
    out = dict(parameters=dict(K=3, Delta=P["Delta"], J=P["J"], U=U, eta_site=P["eta_site"], h=H, exact_cutoff=P["Fock_cutoff"]), exact=exact, results={})
    def build(chi):
        ch = TiltedChain(3, 1.0, det3, P["J"], eta3, U, ring=True); ch.set_tilt(P["eta_site"], chi, 1.0); return ch
    # (i) common maxdeg for all orders
    t = time.time()
    table = [np.cumsum(build(x).theta_series(Nmax, MAXDEG)) for x in CHIS]
    c1, f, c3 = cumulant_partial_sums(table)
    out["results"][f"maxdeg{MAXDEG}"] = records(c1, f, c3, exact=exact, maxdeg=MAXDEG, Nmax=Nmax, time=time.time()-t)
    for k in ("c1", "fano", "c3c1"):
        r = out["results"][f"maxdeg{MAXDEG}"][k]
        print(f"K=3 maxdeg={MAXDEG} {k}: N*={r['N_star']} {r['value']:.5f}+-{r['delta']:.1e} true {r['true_error']:.1e} exact {r['exact']:.5f}  S_N={np.round(r['partial_sums'],5)}", flush=True)
    print(f"  [{time.time()-t:.0f}s]", flush=True)
    # (ii) the paper's convention: order N with maxdeg N+2  (values at the final order only)
    paper = {}
    for N in range(2, min(Nmax, 4)+1):
        t = time.time()
        table = [np.cumsum(build(x).theta_series(N, N+2)) for x in CHIS]
        c1p, fp, c3p = cumulant_partial_sums(table)
        paper[f"N{N}"] = dict(c1=float(c1p[-1]), fano=float(fp[-1]), c3c1=float(c3p[-1]), maxdeg=N+2, time=time.time()-t)
        print(f"K=3 paper convention N<={N} maxdeg={N+2}: c1 {c1p[-1]:.5f} Fano {fp[-1]:.5f} c3/c1 {c3p[-1]:.5f} [{time.time()-t:.0f}s]", flush=True)
    out["paper_convention"] = paper
    json.dump(out, open(f'fcs_k3_truncation_maxdeg{MAXDEG}.json', 'w'), indent=1)

elif what in ("chain8", "chain8check"):
    from fcs_multi import TiltedChain
    P = RP["chain8_fcs"]; K = P["K"]; det = P["detunings"]; eta = np.array(P["eta"]); J = P["J"]
    MAXDEG = 5 if what == "chain8" else int(sys.argv[2]); Nmax = 6 if what == "chain8" else int(sys.argv[3])
    Us = [u for u in P["U_values"] if u > 0] if what == "chain8" else [float(x) for x in sys.argv[4:]]
    out = dict(parameters=dict(K=K, detunings=det, J=J, eta=list(eta), h=H, maxdeg=MAXDEG, Nmax=Nmax, count_site=P["count_site"]), results={})
    for U in Us:
        t = time.time()
        def build(chi):
            ch = TiltedChain(K, 1.0, det, J, eta, U, ring=P["ring"]); ch.set_tilt(P["count_site"], chi, 1.0); return ch
        table = [np.cumsum(build(x).theta_series(Nmax, MAXDEG)) for x in CHIS]
        c1, f, c3 = cumulant_partial_sums(table)
        out["results"][f"{U:.2f}"] = records(c1, f, c3, U=U, time=time.time()-t)
        for k in ("c1", "fano", "c3c1"):
            r = out["results"][f"{U:.2f}"][k]
            print(f"K=8 U={U} maxdeg={MAXDEG} {k}: N*={r['N_star']} {r['value']:.5f}+-{r['delta']:.1e}  S_N={np.round(r['partial_sums'],5)}", flush=True)
        print(f"  [{time.time()-t:.0f}s]", flush=True)
        json.dump(out, open(f'fcs_chain8_maxdeg{MAXDEG}_N{Nmax}.json', 'w'), indent=1)
