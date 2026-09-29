"""Kerr cavity photon number: partial sums to N=16 and the exact reference; the error of every partial sum."""
import sys, json, os, numpy as np
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
sys.path.insert(0, os.path.join(_ROOT, 'engines'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from engine import Model, cav, add, scale, mul, destroy, steady
from series import series_record
P = json.load(open(os.path.join(_ROOT, 'notebook', 'data', 'revision_params.json')))["kerr1"]
kappa, Dc, eta = P["kappa"], P["Delta"], P["eta"]
z = kappa/2 + 1j*Dc; alpha = eta/z
I1 = np.eye(1)
b, bd, one = cav(0, 1, I1), cav(1, 0, I1), cav(0, 0, I1)
a = add(b, scale(alpha, one)); ad = add(bd, scale(np.conj(alpha), one))
HK = scale(0.5, mul(mul(ad, ad), mul(a, a)))          # U = 1
mod = Model(1, z, lambda X: 0*X, I1, HK, vertex_drop=4)
Nmax = 16
cn = mod.series(mul(ad, a), Nmax).real                # coefficients c_N of <a^dag a> = sum c_N U^N
def kerr_exact(U, Nc):
    A = destroy(Nc); Ad = A.conj().T
    H = Dc*Ad@A + 1j*(eta*Ad - np.conj(eta)*A) + U/2*Ad@Ad@A@A
    return np.trace(Ad@A@steady(H, [np.sqrt(kappa)*A])).real
out = dict(parameters=dict(kappa=kappa, Delta=Dc, eta=eta, alpha2=abs(alpha)**2), coefficients=[float(x) for x in cn], Nmax=Nmax, results={})
for U in P["U_values_fig6b"]:
    ex32, ex40 = kerr_exact(U, 32), kerr_exact(U, 40)
    S = np.cumsum(cn*U**np.arange(Nmax+1))
    rec = series_record(S, exact=ex40, U=U, exact_cutoff_check=abs(ex32-ex40))
    rec["order_of_smallest_error"] = int(np.argmin(rec["error_all_orders"]))
    out["results"][f"{U:.2f}"] = rec
    print(f"U={U}: S_16={rec['value']:.6f} error {rec['error']:.2e}; smallest error {min(rec['error_all_orders']):.2e} "
          f"at N={rec['order_of_smallest_error']}; exact={ex40:.6f}")
json.dump(out, open('kerr_truncation.json', 'w'), indent=1)
