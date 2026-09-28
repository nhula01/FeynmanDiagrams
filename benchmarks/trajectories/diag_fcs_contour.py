"""Diagrammatic counting cumulants of the eight-site chain, order by order, with exact chi-derivatives.

theta_N(chi) (partial sums of Eq. (thetaN) through order N, fcs_multi.TiltedChain) is evaluated on the
circle chi = r e^{i phi}; its Taylor coefficients (FFT) give c1, c2, c3 without the O(h^2) error of the
five-point finite differences used in the paper (h = 0.05).  The finite-difference values are computed
alongside for comparison.  usage: diag_fcs_contour.py maxdeg Nmax U1,U2,...  -> diag_chain8_maxdeg{m}.json
"""
import sys, os, json, math, time, numpy as np
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # the code/ directory
sys.path.insert(0, os.path.join(_ROOT, 'engines'))
from fcs_multi import TiltedChain
RP = json.load(open(os.path.join(_ROOT, 'notebook', 'data', 'revision_params.json')))['chain8_fcs']
maxdeg = int(sys.argv[1]); Nmax = int(sys.argv[2]); Us = [float(x) for x in sys.argv[3].split(',')]
K, det, J, eta = RP['K'], RP['detunings'], RP['J'], np.array(RP['eta'])

def partial_sums(U, chi):
    ch = TiltedChain(K, 1.0, det, J, eta, U, ring=RP['ring']); ch.set_tilt(RP['count_site'], chi, 1.0)
    return np.cumsum(ch.theta_series(Nmax, maxdeg))

def ratios(c):
    return [float(c[0]), float(c[1]/c[0]), float(c[2]/c[0])]

out = dict(description=__doc__, maxdeg=maxdeg, Nmax=Nmax, params=RP, results={})
for U in Us:
    t0 = time.time(); rec = dict(U=U)
    for r, M in ((0.1, 16), (0.05, 12)):
        z = r*np.exp(2j*np.pi*np.arange(M)/M)
        th = np.array([partial_sums(U, x) for x in z])          # (M, Nmax+1)
        co = np.fft.fft(th, axis=0)/M
        cn = np.array([(co[n]/r**n)*math.factorial(n) for n in range(1, 4)]).real   # (3, Nmax+1)
        rec[f'contour_r{r}'] = [ratios(cn[:, N]) for N in range(Nmax+1)]
    h = 0.05; tv = np.array([partial_sums(U, x) for x in np.array([-2, -1, 0, 1, 2])*h]).real
    c1 = (tv[3]-tv[1])/(2*h); c2 = (tv[3]-2*tv[2]+tv[1])/h**2; c3 = (tv[4]-2*tv[3]+2*tv[1]-tv[0])/(2*h**3)
    rec['fd_h0.05'] = [[float(c1[N]), float(c2[N]/c1[N]), float(c3[N]/c1[N])] for N in range(Nmax+1)]
    rec['time'] = time.time()-t0
    out['results'][f'{U:.2f}'] = rec
    print(U, 'contour', np.round(rec['contour_r0.1'][3], 5), np.round(rec['contour_r0.1'][4], 5), 'fd', np.round(rec['fd_h0.05'][3], 5), f"{rec['time']:.0f}s", flush=True)
    json.dump(out, open(f'diag_chain8_maxdeg{maxdeg}.json', 'w'), indent=1)

