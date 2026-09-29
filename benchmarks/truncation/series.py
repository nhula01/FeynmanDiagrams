"""Records of perturbative partial sums at a stated order.

Every series in the paper is reported at a fixed, stated order: the highest order computed unless a
lower one is named.  No order is selected from the series itself.

    series_record(partial_sums, exact=None, order=None, **extra)

returns a dict with
    partial_sums       S_0, ..., S_Nmax  (real, or [re, im] where the imaginary part is not negligible)
    N                  the order reported (default Nmax)
    value              S_N
and, when an exact reference is given,
    exact              the reference value
    error              |S_N - exact|
    error_all_orders   |S_n - exact| for n = 0, ..., Nmax
"""
import numpy as np


def _plain(x):
    z = complex(x)
    return z.real if abs(z.imag) < 1e-12 * max(1.0, abs(z)) else [z.real, z.imag]


def series_record(partial_sums, exact=None, order=None, **extra):
    S = np.asarray(partial_sums)
    if S.ndim != 1 or len(S) < 1:
        raise ValueError("need a one-dimensional sequence of partial sums")
    N = len(S) - 1 if order is None else int(order)
    rec = dict(partial_sums=[_plain(x) for x in S], N=N, value=float(np.real(S[N])), **extra)
    if exact is not None:
        rec["exact"] = float(exact)
        rec["error"] = float(abs(np.real(S[N]) - exact))
        rec["error_all_orders"] = [float(abs(np.real(x) - exact)) for x in S]
    return rec
