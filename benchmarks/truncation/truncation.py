"""Optimal truncation of an asymptotic diagrammatic series.

Rule (used everywhere in the revised manuscript):
    Given partial sums S_0, S_1, ..., S_Nmax, stop at the order N* >= 1 at which the
    successive difference d_N = |S_N - S_{N-1}| is smallest over the computed orders.
    If the differences are still decreasing at the highest order, N* = Nmax.
    The error estimate is delta = d_{N*} = |S_{N*} - S_{N*-1}|.
"""
import numpy as np

def optimal_truncation(partial_sums, min_order=1):
    """partial_sums: sequence S_0..S_Nmax (real or complex).
    Returns (value, order, error, info) with
        value = S_{N*}, order = N*, error = |S_{N*} - S_{N*-1}|,
        info  = dict(differences=[d_1..d_Nmax], interior=bool)  where interior=True means the
                minimum lies strictly before the last computed order (asymptotic turnover seen)."""
    S = np.asarray(partial_sums)
    if S.ndim != 1 or len(S) < 2:
        raise ValueError("need at least two partial sums")
    d = np.abs(np.diff(S))                  # d[N-1] = |S_N - S_{N-1}|,  N = 1..Nmax
    cand = np.arange(1, len(S))
    mask = cand >= min_order
    k = int(cand[mask][np.argmin(d[mask])])
    value = S[k]
    if np.iscomplexobj(S) and abs(np.imag(value)) < 1e-12 * max(1.0, abs(value)):
        value = float(np.real(value))
    return value, k, float(d[k-1]), dict(differences=[float(x) for x in d], interior=k < len(S)-1)

def truncation_record(partial_sums, exact=None, **extra):
    """dict with partial sums, N*, delta and (if exact given) the true error at N* and at every order"""
    val, k, err, info = optimal_truncation(partial_sums)
    S = np.asarray(partial_sums)
    d = info["differences"]
    rec = dict(partial_sums=[complex(x).real if abs(complex(x).imag) < 1e-12*max(1, abs(x)) else [complex(x).real, complex(x).imag] for x in S],
               N_star=k, value=float(np.real(val)), delta=err, differences=d,
               interior_minimum=info["interior"],
               relative_delta=float(err/abs(np.real(val))) if val != 0 else None,
               # plateau flag: the smallest difference is not even a factor two below the one before it,
               # i.e. the terms have stopped decreasing -- the series is not usable at this coupling
               plateau=bool(k >= 2 and err > 0.5*d[k-2]), **extra)
    if exact is not None:
        rec["exact"] = float(exact)
        rec["true_error"] = float(abs(np.real(val) - exact))
        rec["true_error_all_orders"] = [float(abs(np.real(x) - exact)) for x in S]
        rec["delta_over_true_error"] = float(err/abs(np.real(val) - exact)) if abs(np.real(val) - exact) > 0 else None
    return rec
