"""assemble fcs_rs_raw.json into cumulant partial sums"""
import sys, json, os, numpy as np
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
sys.path.insert(0, os.path.join(_ROOT, 'engines'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from series import series_record
RP = json.load(open(os.path.join(_ROOT, 'notebook', 'data', 'revision_params.json')))
H = 0.05; CHIS = [-2*H, -H, 0.0, H, 2*H]
raw = json.load(open('fcs_rs_raw.json'))
R = {}
for k, v in raw.items():
    kind, U, chi, Nmax, md, rs = k.split("|"); R[(kind, round(float(U), 3), round(float(chi), 3), int(Nmax), int(md), rs == "True")] = np.array([a + 1j*b for a, b in v["theta"]])
def cps(kind, U, Nmax, md, rs=True):
    keys = [(kind, round(U, 3), round(x, 3), Nmax, md, rs) for x in CHIS]
    if not all(k in R for k in keys): return None
    th = np.array([np.cumsum(R[k]).real for k in keys])
    c1 = (th[3]-th[1])/(2*H); c2 = (th[3]-2*th[2]+th[1])/H**2; c3 = (th[4]-2*th[3]+2*th[1]-th[0])/(2*H**3)
    return c1, c2/c1, c3/c1
def recs(t, exact=None, **extra):
    ex = exact or {}
    return dict(c1=series_record(t[0], exact=ex.get("c1")), fano=series_record(t[1], exact=ex.get("fano")), c3c1=series_record(t[2], exact=ex.get("c3c1")), **extra)
out = dict(note="full Rayleigh-Schroedinger recursion (TiltedChainRS in fcs_rs.py); fd step h=0.05 in chi as in the paper; rs=False reproduces the notebook code")
# K=1
old = json.load(open('fcs_k1_truncation.json'))
out["K1"] = {}
for U in (0.05, 0.10):
    ex = {k: old["results"][f"{U:.2f}"][k]["exact"] for k in ("c1", "fano", "c3c1")}
    d = {}
    for md in (14, 18):
        for rs in (True, False):
            t = cps("k1", U, 8, md, rs)
            if t: d[f"maxdeg{md}_{'fullRS' if rs else 'notebook'}"] = recs(t, exact=ex)
    out["K1"][f"{U:.2f}"] = d
# K=3
k3old = json.load(open('fcs_k3_truncation_maxdeg8.json')); ex3 = k3old["exact"]
out["K3"] = dict(exact=ex3, exact_note="exact theta(chi) from notebook/data/k3_theta.json (Fock cutoff 9), same finite differences")
for md in (8, 9):
    t = cps("k3", 0.05, 6, md)
    if t: out["K3"][f"maxdeg{md}_N6"] = recs(t, exact=ex3)
pc = {}
for N in (2, 3, 4):
    t = cps("k3", 0.05, N, N+2)
    if t: pc[f"N{N}_maxdeg{N+2}"] = dict(c1=float(t[0][-1]), fano=float(t[1][-1]), c3c1=float(t[2][-1]),
                                         err_c1=float(abs(t[0][-1]-ex3["c1"])), err_fano=float(abs(t[1][-1]-ex3["fano"])), err_c3c1=float(abs(t[2][-1]-ex3["c3c1"])))
out["K3"]["paper_convention_fullRS"] = pc
out["K3"]["paper_convention_notebook"] = k3old.get("paper_convention")
# chain8
ch8U = [u for u in RP["chain8_fcs"]["U_values"] if u > 0]
out["chain8"] = dict(parameters=RP["chain8_fcs"])
for md, Nm in ((5, 6), (6, 6), (7, 4)):
    d = {}
    for U in ch8U:
        t = cps("chain8", U, Nm, md)
        if t: d[f"{U:.2f}"] = recs(t, U=U)
    if d: out["chain8"][f"maxdeg{md}_N{Nm}"] = d
json.dump(out, open('fcs_rs_truncation.json', 'w'), indent=1)
def show(tag, r):
    print(tag, "  ".join(f"{q}: N={r[q]['N']} {r[q]['value']:.5f}" + (f" err {r[q]['error']:.1e}" if 'error' in r[q] else "") for q in ("c1", "fano", "c3c1")))
for U, d in out["K1"].items():
    for k, r in d.items(): show(f"K1 U={U} {k}", r)
for k in ("maxdeg8_N6", "maxdeg9_N6"):
    if k in out["K3"]: show(f"K3 {k}", out["K3"][k])
print("K3 paper conv fullRS", json.dumps(pc))
for k, d in out["chain8"].items():
    if k == "parameters": continue
    for U, r in d.items(): show(f"chain8 {k} U={U}", r)
