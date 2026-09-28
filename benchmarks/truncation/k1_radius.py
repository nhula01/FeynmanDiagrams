"""K=1 exact contour cumulants versus radius (branch-point check) and location of the nearest level crossing."""
import sys, os, json, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import table5_contour as T5
out = {}
for U in (0.05, 0.10):
    H, aa = T5.ops([26], [T5.DC], 0.0, [T5.ETA], U, False); L, Jc, D = T5.liouv(H, aa, 0); L = L.toarray(); Jd = Jc.toarray()
    for r, M in ((0.03, 32), (0.05, 32), (0.08, 32), (0.08, 64), (0.12, 32), (0.15, 32)):
        z = T5.circle(r, M); t = time.time()
        vals = T5.fill([T5.lead_dense(L + (np.exp(z[j])-1)*Jd) for j in T5.half(M)], M)
        cn = T5.taylor(vals, r, nmax=6)
        out[f"U{U:.2f}_r{r}_M{M}"] = dict(U=U, r=r, M=M, c=[[x.real, x.imag] for x in cn], c1=cn[1].real, fano=(cn[2]/cn[1]).real, c3c1=(cn[3]/cn[1]).real,
                                         c4c1=(cn[4]/cn[1]).real, c5c1=(cn[5]/cn[1]).real, c0_abs=abs(cn[0]), time=time.time()-t)
        print(U, r, M, f"c1 {cn[1].real:.8f} F {(cn[2]/cn[1]).real:.8f} c3/c1 {(cn[3]/cn[1]).real:.8f} c5/c1 {(cn[5]/cn[1]).real:.3f} |c0| {abs(cn[0]):.1e}", flush=True)
        json.dump(out, open('k1_radius.json', 'w'), indent=1)
    # gap between the two eigenvalues with largest real part along the negative/positive real and imaginary chi axes
    gaps = {}
    for ang in (0, np.pi/2, np.pi):
        g = []
        for rr in np.linspace(0.02, 0.4, 20):
            w = np.linalg.eigvals(L + (np.exp(rr*np.exp(1j*ang))-1)*Jd); w = w[np.argsort(-w.real)]
            g.append([float(rr), float(abs(w[0]-w[1])), float(w[0].real-w[1].real)])
        gaps[f"{ang:.3f}"] = g
    out[f"U{U:.2f}_gaps"] = gaps; json.dump(out, open('k1_radius.json', 'w'), indent=1)
print("done")
