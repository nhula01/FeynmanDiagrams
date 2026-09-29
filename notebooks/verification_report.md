# Verification of the numbers quoted in the paper

238 checks passed, 0 failed, 3 informational

Tables V and SI compared against the snapshots in notebooks/paper_reference/.

| section | quantity | paper | code | ok |
|---|---|---|---|---|
| Table V | exact: 0.8596 1.0749 1.2511 | 0.8596 1.0749 1.2511 | 0.8596 1.0749 1.2511 | PASS |
| Table V | mean field: 0.8576 1.0000 1.0000 | 0.8576 1.0000 1.0000 | 0.8576 1.0000 1.0000 | PASS |
| Table V | Gaussian FCS: 0.8596 1.0730 1.2296 | 0.8596 1.0730 1.2296 | 0.8596 1.0730 1.2296 | PASS |
| Table V | diagrams $N\leq8$: 0.8596 1.0749 1.2510 | 0.8596 1.0749 1.2510 | 0.8596 1.0749 1.2510 | PASS |
| Table V | exact: 0.9460 1.2702 2.3585 | 0.9460 1.2702 2.3585 | 0.9460 1.2702 2.3585 | PASS |
| Table V | mean field: 0.9327 1.0000 1.0000 | 0.9327 1.0000 1.0000 | 0.9327 1.0000 1.0000 | PASS |
| Table V | Gaussian FCS: 0.9453 1.2292 1.7759 | 0.9453 1.2292 1.7759 | 0.9453 1.2292 1.7759 | PASS |
| Table V | diagrams $N\leq8$: 0.9458 1.2629 2.1848 | 0.9458 1.2629 2.1848 | 0.9458 1.2629 2.1848 | PASS |
| Table V | exact: 1.1458 1.1653 1.6393 | 1.1458 1.1653 1.6393 | 1.1458 1.1653 1.6393 | PASS |
| Table V | mean field ($U=0$ displacement): 1.0381 1.0000 1.0000 | 1.0381 1.0000 1.0000 | 1.0381 1.0000 1.0000 | PASS |
| Table V | diagrams $N\leq2$: 1.1422 1.1385 1.4572 | 1.1422 1.1385 1.4572 | 1.1422 1.1385 1.4572 | PASS |
| Table V | diagrams $N\leq3$: 1.1450 1.1558 1.5520 | 1.1450 1.1558 1.5520 | 1.1450 1.1558 1.5520 | PASS |
| Table V | diagrams $N\leq4$: 1.1456 1.1620 1.5997 | 1.1456 1.1620 1.5997 | 1.1456 1.1620 1.5997 | PASS |
| Table V | diagrams $N\leq6$: 1.1458 1.1649 1.6321 | 1.1458 1.1649 1.6321 | 1.1458 1.1649 1.6321 | PASS |
| Table V | row count | 14 | 14 | PASS |
| SM Table SI | U=0.10 g3: 1.312582 1.2(-2) 8.5(-4) 1.8(-4) 1.7(-2) 1.0(-3) 1.8(-6) | 1.312582 1.2(-2) 8.5(-4) 1.8(-4) 1.7(-2) 1.0(-3) 1.8(-6) | 1.312582 1.2(-2) 8.5(-4) 1.8(-4) 1.7(-2) 1.0(-3) 1.8(-6) | PASS |
| SM Table SI | U=0.10 g2: 1.093641 1.0(-3) 9.6(-5) 1.3(-6) 2.9(-3) 1.1(-4) 1.5(-7) | 1.093641 1.0(-3) 9.6(-5) 1.3(-6) 2.9(-3) 1.1(-4) 1.5(-7) | 1.093641 1.0(-3) 9.6(-5) 1.3(-6) 2.9(-3) 1.1(-4) 1.5(-7) | PASS |
| SM Table SI | U=0.10 n: 0.296458 4.1(-5) 4.6(-7) 2.8(-8) 8.6(-4) 3.2(-5) 1.7(-8) | 0.296458 4.1(-5) 4.6(-7) 2.8(-8) 8.6(-4) 3.2(-5) 1.7(-8) | 0.296458 4.1(-5) 4.6(-7) 2.8(-8) 8.6(-4) 3.2(-5) 1.7(-8) | PASS |
| SM Table SI | U=0.15 g3: 1.538157 3.1(-2) 4.7(-3) 1.8(-3) 5.3(-2) 6.5(-3) 1.6(-4) | 1.538157 3.1(-2) 4.7(-3) 1.8(-3) 5.3(-2) 6.5(-3) 1.6(-4) | 1.538157 3.1(-2) 4.7(-3) 1.8(-3) 5.3(-2) 6.5(-3) 1.6(-4) | PASS |
| SM Table SI | U=0.15 g2: 1.150869 3.9(-3) 6.2(-4) 1.4(-5) 9.8(-3) 7.5(-4) 1.4(-5) | 1.150869 3.9(-3) 6.2(-4) 1.4(-5) 9.8(-3) 7.5(-4) 1.4(-5) | 1.150869 3.9(-3) 6.2(-4) 1.4(-5) 9.8(-3) 7.5(-4) 1.4(-5) | PASS |
| SM Table SI | U=0.15 n: 0.302914 2.5(-4) 4.7(-6) 1.3(-7) 3.2(-3) 2.6(-4) 2.0(-6) | 0.302914 2.5(-4) 4.7(-6) 1.3(-7) 3.2(-3) 2.6(-4) 2.0(-6) | 0.302914 2.5(-4) 4.7(-6) 1.3(-7) 3.2(-3) 2.6(-4) 2.0(-6) | PASS |
| SM Table SI | U=0.20 g3: 1.831380 6.1(-2) 1.3(-2) 8.3(-3) 1.1(-1) 2.0(-2) 3.0(-3) | 1.831380 6.1(-2) 1.3(-2) 8.3(-3) 1.1(-1) 2.0(-2) 3.0(-3) | 1.831380 6.1(-2) 1.3(-2) 8.3(-3) 1.1(-1) 2.0(-2) 3.0(-3) | PASS |
| SM Table SI | U=0.20 g2: 1.216912 9.6(-3) 2.2(-3) 1.3(-5) 2.2(-2) 2.4(-3) 2.8(-4) | 1.216912 9.6(-3) 2.2(-3) 1.3(-5) 2.2(-2) 2.4(-3) 2.8(-4) | 1.216912 9.6(-3) 2.2(-3) 1.3(-5) 2.2(-2) 2.4(-3) 2.8(-4) | PASS |
| SM Table SI | U=0.20 n: 0.310913 9.4(-4) 7.9(-6) 1.0(-5) 8.5(-3) 1.2(-3) 5.5(-5) | 0.310913 9.4(-4) 7.9(-6) 1.0(-5) 8.5(-3) 1.2(-3) 5.5(-5) | 0.310913 9.4(-4) 7.9(-6) 1.0(-5) 8.5(-3) 1.2(-3) 5.5(-5) | PASS |
| SM Table SI | U=0.25 g3: 2.181138 8.8(-2) 2.0(-2) 2.7(-2) 1.8(-1) 3.2(-2) 1.9(-2) | 2.181138 8.8(-2) 2.0(-2) 2.7(-2) 1.8(-1) 3.2(-2) 1.9(-2) | 2.181138 8.8(-2) 2.0(-2) 2.7(-2) 1.8(-1) 3.2(-2) 1.9(-2) | PASS |
| SM Table SI | U=0.25 g2: 1.290142 1.7(-2) 5.3(-3) 1.7(-4) 3.9(-2) 3.7(-3) 1.7(-3) | 1.290142 1.7(-2) 5.3(-3) 1.7(-4) 3.9(-2) 3.7(-3) 1.7(-3) | 1.290142 1.7(-2) 5.3(-3) 1.7(-4) 3.9(-2) 3.7(-3) 1.7(-3) | PASS |
| SM Table SI | U=0.25 n: 0.321042 2.3(-3) 3.1(-4) 7.2(-5) 1.8(-2) 3.3(-3) 5.8(-4) | 0.321042 2.3(-3) 3.1(-4) 7.2(-5) 1.8(-2) 3.3(-3) 5.8(-4) | 0.321042 2.3(-3) 3.1(-4) 7.2(-5) 1.8(-2) 3.3(-3) 5.8(-4) | PASS |
| SM Table SI | U=0.30 g3: 2.527280 8.9(-2) 1.8(-2) 6.8(-2) 2.2(-1) 2.3(-2) 5.1(-2) | 2.527280 8.9(-2) 1.8(-2) 6.8(-2) 2.2(-1) 2.3(-2) 5.1(-2) | 2.527280 8.9(-2) 1.8(-2) 6.8(-2) 2.2(-1) 2.3(-2) 5.1(-2) | PASS |
| SM Table SI | U=0.30 g2: 1.363850 2.3(-2) 1.2(-2) 2.4(-3) 5.4(-2) 9.1(-4) 2.6(-3) | 1.363850 2.3(-2) 1.2(-2) 2.4(-3) 5.4(-2) 9.1(-4) 2.6(-3) | 1.363850 2.3(-2) 1.2(-2) 2.4(-3) 5.4(-2) 9.1(-4) 2.6(-3) | PASS |
| SM Table SI | U=0.30 n: 0.333841 3.7(-3) 1.6(-3) 1.0(-7) 3.3(-2) 7.0(-3) 3.0(-3) | 0.333841 3.7(-3) 1.6(-3) 1.0(-7) 3.3(-2) 7.0(-3) 3.0(-3) | 0.333841 3.7(-3) 1.6(-3) 1.0(-7) 3.3(-2) 7.0(-3) 3.0(-3) | PASS |
| SM Table SI | rows checked | 15 | 15 | PASS |
| Sec. V.A / Fig. 5 | U=0.05: error at N=16, 4e-12 | 4e-12 | 4.25171e-12 | PASS |
| Sec. V.A / Fig. 5 | U=0.10: error at N=16, 2e-7 | 2e-07 | 1.81281e-07 | PASS |
| Sec. V.A / Fig. 5 | U=0.05, 0.10: error still decreasing at N=16 | 1 | 1 | PASS |
| Sec. V.A / Fig. 5 | U=0.20: error at N=8, 0.061 | 0.061 | 0.0613993 | PASS |
| Sec. V.A / Fig. 5 | U=0.20: error at N=16, 0.37 | 0.37 | 0.373967 | PASS |
| Sec. V.A / Fig. 5 | U=0.20: minimum error (at N=9) | 0.019 | 0.0193687 | PASS |
| Sec. V.A / Fig. 5 | U=0.20: order of minimum error | 9 | 9 | PASS |
| Sec. V.A / Fig. 5 | U=0.20: minimum relative error, 1.3% | 1.3 | 1.30449 | PASS |
| Sec. V.A / Fig. 5 | U=0.30: minimum error | 0.16 | 0.162664 | PASS |
| Sec. V.A / Fig. 5 | U=0.30: order of minimum error | 5 | 5 | PASS |
| Sec. V.A / Fig. 5 | reported order of the stored series (N=16) | 16 | 16 | PASS |
| Sec. V.A / Fig. 5 | U=0.20: exact <a^dag a> | 1.485 | 1.48477 | PASS |
| Sec. V.A / Fig. 5 | U=0.30: exact <a^dag a> | 2.43 | 2.43311 | PASS |
| Sec. V.A / Fig. 5 | U=0.20: mean-field <a^dag a> | 1.215 | 1.21499 | PASS |
| Sec. V.A / Fig. 5 | U=0.30: mean-field <a^dag a> | 3.76 | 3.7571 | PASS |
| Sec. V.A / Fig. 5 | U=0.30: mean-field error | 1.3 | 1.32399 | PASS |
| Sec. V.A / Fig. 5 | root test |c_N|^(-1/N) at N=10 | 0.27 | 0.273252 | PASS |
| Sec. V.A / Fig. 5 | root test |c_N|^(-1/N) at N=40 | 0.15 | 0.151482 | PASS |
| Sec. VI / Fig. 8 | n: second-order relative error at U=0.05, 0.01% | 0.01 | 0.00978073 | PASS |
| Sec. VI / Fig. 8 | n: second-order relative error at U=0.15, 0.3% | 0.3 | 0.321512 | PASS |
| Sec. VI / Fig. 8 | n: second-order relative error at U=0.30, 3% | 3 | 3.31171 | PASS |
| Sec. VI / Fig. 8 | first/second-order error ratio on n, U=0.05-0.3: min (two) | 2 | 1.97829 | PASS |
| Sec. VI / Fig. 8 | first/second-order error ratio on n, U=0.05-0.3: max (twelve) | 12 | 11.5513 | PASS |
| Sec. VI / Fig. 8 | U=0.15: exact n_0 | 0.3029 | 0.302914 | PASS |
| Sec. VI / Fig. 8 | U=0.15: exact g2 | 1.151 | 1.15087 | PASS |
| Sec. VI / Fig. 8 | U=0.15: n_0 through second order | 0.3019 | 0.30194 | PASS |
| Sec. VI / Fig. 8 | U=0.15: g2 through second order | 1.14 | 1.13964 | PASS |
| Sec. VI / Fig. 8 | U=0.15: g2 through fourth order | 1.15 | 1.15001 | PASS |
| Sec. VI / Fig. 8 | U=0.15: g2 through tenth order, 1.15089 | 1.15089 | 1.15089 | PASS |
| Sec. VI / Fig. 8 | U=0.15: exact g2, 1.15087 | 1.15087 | 1.15087 | PASS |
| Sec. VI / Fig. 8 | U=0.15: tenth-order g2 error, 1.7e-5 | 1.7e-05 | 1.67522e-05 | PASS |
| Sec. VI / Fig. 8 | reported order of the ring series (N=10) | 10 | 10 | PASS |
| Sec. VI / Fig. 8 | U=0.10: tenth-order g2 error, 2e-7 | 2e-07 | 1.69103e-07 | PASS |
| Sec. VI / Fig. 8 | U=0.30: tenth-order g2 error, 4e-3 | 0.004 | 0.00389287 | PASS |
| Sec. VI / Fig. 8 | U=0.30: fourth-order g2 error, 1e-3 (smaller than tenth) | 0.001 | 0.000950349 | PASS |
| Sec. VI / Fig. 8 | g2: a lower order beats the tenth at some U>=0.2 (e.g. 0.2, 0.225, 0.3) | 1 | 1 | PASS |
| Sec. VI / Fig. 8 | g3: a lower order beats the tenth from U=0.175 (first U>=0.1) | 0.175 | 0.175 | PASS |
| Sec. VI / Fig. 8 | g3: ... and at every U>=0.175 | 1 | 1 | PASS |
| Sec. VI / Fig. 8 | exact steady state, 5 levels: growth per site (45) | 45 | 44.629 | PASS |
| Sec. VI / Fig. 8 | exact steady state, 5 levels, K=4: 328 s | 328 | 328.243 | PASS |
| Sec. VI / Fig. 8 | momentum engine N<=2, K=6: 2.4 s | 2.4 | 2.43876 | PASS |
| Sec. VI / Fig. 8 | momentum engine N<=2: exponent 5.7 | 5.7 | 5.65156 | PASS |
| Sec. VI / Fig. 8 | tensor engine N<=2, K=32: 2.1 s | 2.1 | 2.13121 | PASS |
| Sec. VI / Fig. 8 | tensor engine N<=2: exponent above K=16, 5.9 | 5.9 | 5.93716 | PASS |
| Sec. VI / Fig. 8 | tensor engine N<=4, K=14: 30 s | 30 | 30.2494 | PASS |
| Sec. VI / Fig. 8 | tensor engine N<=4: exponent 6.0 | 6 | 5.99493 | PASS |
| Sec. VI / Fig. 8 | tensor engine N<=4, contracted commutator, K=14: 6.1 s | 6.1 | 6.10118 | PASS |
| Sec. VI / Fig. 8 | label bound C(2K+Dmax,Dmax) ~ K^Dmax: Dmax at N=2 for d0=4 (K^4) | 4 | 4 | PASS |
| Sec. VI / Fig. 8 | label bound C(2K+Dmax,Dmax) ~ K^Dmax: Dmax at N=3 for d0=4 (K^4) | 4 | 4 | PASS |
| Sec. VI / Fig. 8 | label bound C(2K+Dmax,Dmax) ~ K^Dmax: Dmax at N=4 for d0=4 (K^6) | 6 | 6 | PASS |
| Sec. VI / Fig. 8 | label bound C(2K+Dmax,Dmax) ~ K^Dmax: Dmax at N=5 for d0=4 (K^6) | 6 | 6 | PASS |
| Sec. VI / Fig. 8 | measured distinct-label count of the tensor engine, g2, N=2, K>=8: local exponent 3.60 (bound: K^4) | 4 | 3.6002 | info |
| Sec. VI / Fig. 8 | measured distinct-label count of the tensor engine, g2, N=3, K>=8: local exponent 3.55 (bound: K^4) | 4 | 3.55329 | info |
| Sec. VI / Fig. 8 | measured distinct-label count of the tensor engine, g2, N=4, K>=8: local exponent 5.18 (bound: K^6) | 6 | 5.17857 | info |
| Sec. VI.B / Fig. 9 / SM Table SI | closure variables at order 2 | 27 | 27 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | closure variables at order 3 | 83 | 83 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | closure variables at order 4 | 209 | 209 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | closure order 2: 0.13 s per evaluation | 0.13 | 0.134304 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | closure order 3: 0.6 s | 0.6 | 0.624568 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | closure order 4: 5 s | 5 | 5.13944 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | diagrams N<=2, pruned symbolic recursion: 0.21 s per coupling | 0.21 | 0.21406 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | diagrams N<=4, pruned symbolic recursion: 1.1 s per coupling | 1.1 | 1.09812 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | sparse-matrix form, orders 0-2, all couplings: 0.03 s | 0.03 | 0.0298066 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | sparse-matrix form, orders 0-4, all couplings: 0.11 s | 0.11 | 0.113036 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | sparse-matrix form, orders 0-10, all couplings: 4.2 s | 4.2 | 4.16458 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI |    ... below one closure-4 evaluation | 1 | 1 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | closures cost 0.6 to 5 times the pruned fourth-order diagrams: min | 0.6 | 0.568759 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | closures cost 0.6 to 5 times the pruned fourth-order diagrams: max | 5 | 4.6802 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.1: exact g3 | 1.3126 | 1.31258 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.1: fourth-order diagrams reduce the g3 error of HFB by a factor 12 | 12 | 11.5671 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.1: HFB g3 error 1.2% | 1.2 | 1.15746 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.1: fourth-order diagram g3 error 0.10% | 0.1 | 0.100064 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.15-0.3: HFB/diagrams-4 g3 error ratio, min 2.8 | 2.8 | 2.78969 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.15-0.3: HFB/diagrams-4 g3 error ratio, max 4.8 | 4.8 | 4.82439 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.1: diagrams-4 g3 error / best closure, 5.5 | 5.5 | 5.4863 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.3: diagrams-4 g3 error / best closure, 1.25 | 1.25 | 1.25435 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.100: tenth-order g3 relative error 1.8e-06 | 1.8e-06 | 1.7709e-06 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.150: tenth-order g3 relative error 1.6e-04 | 0.00016 | 0.000159916 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.200: tenth-order g3 relative error 3.0e-03 | 0.003 | 0.00300392 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.250: tenth-order g3 relative error 1.9e-02 | 0.019 | 0.0186955 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.300: tenth-order g3 relative error 5.1e-02 | 0.051 | 0.0508276 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.100: tenth order better than best closure on g3 by 103 | 103 | 102.993 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.150: tenth order better than best closure on g3 by 11 | 11 | 11.0547 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.200: tenth order better than best closure on g3 by 2.8 | 2.8 | 2.7548 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.225: tenth order better than best closure on g3 by 1.8 | 1.8 | 1.81268 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.250: tenth order better than best closure on g3 by 1.05 | 1.05 | 1.04951 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.275: tenth order worse than best closure on g3 by 1.75 | 1.75 | 1.74693 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.300: tenth order worse than best closure on g3 by 2.8 | 2.8 | 2.75089 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | tenth order ahead of the best closure on g3 up to U=0.25 (abstract, intro, Sec. VI.B) | 0.25 | 0.25 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | tenth order ahead of the fourth-order closure on g2 only up to U=0.125 | 0.125 | 0.125 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | tenth order ahead of the fourth-order closure on n only up to U=0.1 | 0.1 | 0.1 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.125: tenth order ahead on g2 by 2.8 | 2.8 | 2.80273 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.15: closure-4 ahead on g2 by 1.06 | 1.06 | 1.05503 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.2: closure-4 ahead on g2 by 22 | 22 | 22.222 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.15-0.25: closure ahead on n by 6 to 16 (min) | 6 | 5.74414 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | U=0.15-0.25: closure ahead on n by 6 to 16 (max) | 16 | 15.5457 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | orders 0-10 of n, n2, n3 at K=3, symbolic recursion: 54 s per coupling (median of the raw runs) | 54 | 53.7803 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | sparse-matrix construction: exponent at second order over K=6-12 (4.9) | 4.9 | 4.88853 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | sparse-matrix construction: exponent at fourth order over K=6-12 (6.4) | 6.4 | 6.35248 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | sparse-matrix construction at K=12: second order 7.5 s | 7.5 | 7.54407 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | sparse-matrix construction at K=12: fourth order 240 s | 240 | 235.938 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | SM Table (pruning), none: total labels | 183618 | 183618 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | SM Table (pruning), none: labels after vertex 7 | 81793 | 81793 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | SM Table (pruning), degree_bound: total labels | 13927 | 13927 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | SM Table (pruning), degree_bound: labels after vertex 7 | 28 | 28 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | SM Table (pruning), degree_bound_momentum: total labels | 4663 | 4663 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | SM Table (pruning), degree_bound_momentum: labels after vertex 7 | 10 | 10 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | SM Table (pruning): same coefficient c8 in all three columns | 0 | 8.88178e-15 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | SM Table (pruning): time without selection, 86 s | 86 | 86.0378 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | SM Table (pruning): time with both rules, 1.3 s | 1.3 | 1.26956 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | largest label set of the sparse form at Nmax=10 (3702) | 3702 | 3702 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | labels after vertices 0-6 for n3 at Nmax=6: 136, 408, 691, 312, 72, 10, 1 | 0 | 0 | PASS |
| Sec. VI.B / Fig. 9 / SM Table SI | no plotted correlator below one (min g2, g3) | 1 | 1.09044 | PASS |
| Sec. VII / Table V text | K=1, U=0.05: eighth-order Fano error 8e-6 | 8e-06 | 8.1996e-06 | PASS |
| Sec. VII / Table V text | K=1, U=0.05: eighth-order c3/c1 error 1.5e-4 | 0.00015 | 0.000146968 | PASS |
| Sec. VII / Table V text | K=1, U=0.10: Gaussian error on c3/c1, 25% | 25 | 24.7032 | PASS |
| Sec. VII / Table V text | K=1, U=0.10: eighth-order error on c3/c1, 7% | 7 | 7.36281 | PASS |
| Sec. VII / Table V text | K=3: c3/c1 error at N<=3, 0.087 | 0.087 | 0.0872337 | PASS |
| Sec. VII / Table V text | K=3: c3/c1 error at N<=4, 0.040 | 0.04 | 0.0395696 | PASS |
| Sec. VII / Table V text | K=3: c3/c1 error at N<=6, 7e-3 | 0.007 | 0.00712254 | PASS |
| Sec. VII / Table V text | K=3: Fano error at N<=6, 3.4e-4 | 0.00034 | 0.000340635 | PASS |
| Sec. VII / Table V text | K=3: c3/c1 error falls through every order computed | 1 | 1 | PASS |
| Sec. VII / Fig. 10 | U=0.1: deviation of Fano factor from Poisson, 24% | 24 | 24.0074 | PASS |
| Sec. VII / Fig. 10 | U=0.1: deviation of c3/c1 from Poisson, 101% | 101 | 101.34 | PASS |
| Sec. VII / Fig. 10 | U=0.1: Gaussian share of the Fano deviation, 93% | 93 | 93.3602 | PASS |
| Sec. VII / Fig. 10 | U=0.1: Gaussian share of the c3/c1 deviation, 74% | 74 | 73.6363 | PASS |
| Sec. VII / Fig. 10 | U=0.02: diagrams minus Gaussian on c3/c1, 0.002 | 0.002 | 0.002153 | PASS |
| Sec. VII / Fig. 10 | U=0.10: diagrams minus Gaussian on c3/c1, 0.27 | 0.27 | 0.267169 | PASS |
| Sec. VII / Fig. 10 | non-Gaussian part: c3/c1 vs Fano ratio, min 17 | 17 | 16.7604 | PASS |
| Sec. VII / Fig. 10 | non-Gaussian part: c3/c1 vs Fano ratio, max 'forty' (code gives 38) | 40 | 37.9151 | PASS |
| Sec. VII / Fig. 10 | chain: reported order N=6 (highest computed), all cumulants | 6 | 6 | PASS |
| Sec. VII / Fig. 10 | U=0.1: c3/c1 increments shrink by ~0.65 per order | 0.65 | 0.626828 | PASS |
| Sec. VII / Fig. 10 | U=0.1: geometric tail of c3/c1, about 0.1 | 0.1 | 0.088556 | PASS |
| Sec. VII / Fig. 10 | U=0.04: number of trajectories | 1216 | 1216 | PASS |
| Sec. VII / Fig. 10 | U=0.04: total counting time, 1.2e6/kappa | 1.2 | 1.216 | PASS |
| Sec. VII / Fig. 10 | U=0.04: Fock cutoffs [10, 6, 5, 4, 3, 3, 2, 2] | 0 | 0 | PASS |
| Sec. VII / Fig. 10 | U=0.08: number of trajectories | 3072 | 3072 | PASS |
| Sec. VII / Fig. 10 | U=0.08: total counting time, 3.1e6/kappa | 3.1 | 3.072 | PASS |
| Sec. VII / Fig. 10 | U=0.08: Fock cutoffs [12, 7, 5, 4, 3, 3, 2, 2] | 0 | 0 | PASS |
| Sec. VII / Fig. 10 | U=0.10: number of trajectories | 3552 | 3552 | PASS |
| Sec. VII / Fig. 10 | U=0.10: total counting time, 3.6e6/kappa | 3.6 | 3.552 | PASS |
| Sec. VII / Fig. 10 | U=0.10: Fock cutoffs [14, 7, 5, 4, 3, 3, 2, 2] | 0 | 0 | PASS |
| Sec. VII / Fig. 10 | U=0.1 trajectories: c1 = 0.9309(6) | 0.9309 | 0.930898 | PASS |
| Sec. VII / Fig. 10 |    se(c1) = 0.0006 | 0.0006 | 0.000584998 | PASS |
| Sec. VII / Fig. 10 | U=0.1 trajectories: c2/c1 = 1.2363(49) | 1.2363 | 1.23628 | PASS |
| Sec. VII / Fig. 10 |    se = 0.0049 | 0.0049 | 0.00492687 | PASS |
| Sec. VII / Fig. 10 | U=0.1 trajectories: c3/c1 = 2.054(49) | 2.054 | 2.05357 | PASS |
| Sec. VII / Fig. 10 |    se = 0.049 | 0.049 | 0.0490818 | PASS |
| Sec. VII / Fig. 10 | U=0.1 diagrams N<=6 = chain8 sixth-order partial sum (c3/c1) | 0 | 0 | PASS |
| Sec. VII / Fig. 10 | U=0.1 diagrams N<=6: c1 = 0.9306 | 0.9306 | 0.930573 | PASS |
| Sec. VII / Fig. 10 | U=0.1 diagrams N<=6: c2/c1 = 1.2401 | 1.2401 | 1.24007 | PASS |
| Sec. VII / Fig. 10 | U=0.1 diagrams N<=6: c3/c1 = 2.013 | 2.013 | 2.0134 | PASS |
| Sec. VII / Fig. 10 | U=0.1 diagrams within one standard error (max |z|) | 1 | 0.82 | PASS |
| Sec. VII / Fig. 10 | U=0.1 third order: c1 = 0.9282 | 0.9282 | 0.928191 | PASS |
| Sec. VII / Fig. 10 | U=0.1 third order: c2/c1 = 1.2075 | 1.2075 | 1.20749 | PASS |
| Sec. VII / Fig. 10 | U=0.1 third order: c3/c1 = 1.749 | 1.749 | 1.74881 | PASS |
| Sec. VII / Fig. 10 | U=0.1 third order low by 4.6 s.e. on c1 | 4.6 | 4.63 | PASS |
| Sec. VII / Fig. 10 | U=0.1 third order low by 5.8 s.e. on c2/c1 | 5.8 | 5.84 | PASS |
| Sec. VII / Fig. 10 | U=0.1 third order low by 6.2 s.e. on c3/c1 | 6.2 | 6.21 | PASS |
| Sec. VII / Fig. 10 | U=0.1 Gaussian: c1 = 0.9305 | 0.9305 | 0.930473 | PASS |
| Sec. VII / Fig. 10 | U=0.1 Gaussian: c2/c1 = 1.2241 | 1.2241 | 1.22413 | PASS |
| Sec. VII / Fig. 10 | U=0.1 Gaussian: c3/c1 = 1.746 | 1.746 | 1.74623 | PASS |
| Sec. VII / Fig. 10 | U=0.1 Gaussian low by 2.5 s.e. on Fano | 2.5 | 2.47 | PASS |
| Sec. VII / Fig. 10 | U=0.1 Gaussian low by 6.3 s.e. on c3/c1 | 6.3 | 6.26 | PASS |
| Sec. VII / Fig. 10 | U=0.1 non-Gaussian part of c3/c1, trajectories 0.31 | 0.31 | 0.307339 | PASS |
| Sec. VII / Fig. 10 | U=0.1 non-Gaussian part of c3/c1, diagrams 0.27 | 0.27 | 0.267169 | PASS |
| Sec. VII / Fig. 10 | U=0.08 trajectories: c1 = 0.9009(6) | 0.9009 | 0.900887 | PASS |
| Sec. VII / Fig. 10 | U=0.08 trajectories: c2/c1 = 1.1541(49) | 1.1541 | 1.15414 | PASS |
| Sec. VII / Fig. 10 | U=0.08 trajectories: c3/c1 = 1.585(45) | 1.585 | 1.58499 | PASS |
| Sec. VII / Fig. 10 | U=0.08 diagrams: c1 = 0.9002 | 0.9002 | 0.90024 | PASS |
| Sec. VII / Fig. 10 | U=0.08 diagrams: c2/c1 = 1.1613 | 1.1613 | 1.16133 | PASS |
| Sec. VII / Fig. 10 | U=0.08 diagrams: c3/c1 = 1.617 | 1.617 | 1.61745 | PASS |
| Sec. VII / Fig. 10 | U=0.08 diagrams within 1.5 s.e. (max |z|) | 1.5 | 1.47 | PASS |
| Sec. VII / Fig. 10 | U=0.08 Gaussian: c3/c1 = 1.499 | 1.499 | 1.49925 | PASS |
| Sec. VII / Fig. 10 | U=0.08 Gaussian low by 1.9 s.e. on c3/c1 | 1.9 | 1.92 | PASS |
| Sec. VII / Fig. 10 | U=0.04 trajectories: c1 = 0.8507(9) | 0.8507 | 0.850677 | PASS |
| Sec. VII / Fig. 10 | U=0.04 trajectories: c2/c1 = 1.0635(71) | 1.0635 | 1.0635 | PASS |
| Sec. VII / Fig. 10 | U=0.04 trajectories: c3/c1 = 1.286(58) | 1.286 | 1.28553 | PASS |
| Sec. VII / Fig. 10 | U=0.04 trajectory error on c3/c1, 0.058 | 0.058 | 0.0576781 | PASS |
| Sec. VII / Fig. 10 | U=0.04 diagrams: c1 = 0.8497 | 0.8497 | 0.849669 | PASS |
| Sec. VII / Fig. 10 | U=0.04 diagrams: c2/c1 = 1.0597 | 1.0597 | 1.05974 | PASS |
| Sec. VII / Fig. 10 | U=0.04 diagrams: c3/c1 = 1.196 | 1.196 | 1.19622 | PASS |
| Sec. VII / Fig. 10 | U=0.04 non-Gaussian part of c3/c1, 0.013 | 0.013 | 0.0131321 | PASS |
| Sec. VII / Fig. 10 | U=0.04 diagrams within 1.6 s.e. | 1.6 | 1.55 | PASS |
| Sec. VII / Fig. 10 | U=0.04 Gaussian within 1.8 s.e. | 1.8 | 1.78 | PASS |
| Sec. VII / Fig. 10 | U=0.04: counted c1 vs kappa<n_0>, |z| (SM) | 1.1 | 1.09918 | PASS |
| Sec. VII / Fig. 10 | U=0.08: counted c1 vs kappa<n_0>, |z| (SM) | 1 | 0.970716 | PASS |
| Sec. VII / Fig. 10 | U=0.10: counted c1 vs kappa<n_0>, |z| (SM) | 0.3 | 0.296816 | PASS |
| Sec. VII / Fig. 10 | U=0.1 cutoff-8 run: c1 = 0.9283(10) | 0.9283 | 0.928288 | PASS |
| Sec. VII / Fig. 10 | U=0.1 cutoff-8 run: c2/c1 = 1.2098(84) | 1.2098 | 1.20978 | PASS |
| Sec. VII / Fig. 10 | U=0.1 cutoff-8 run: c3/c1 = 1.742(85) | 1.742 | 1.74233 | PASS |
| Sec. VII / Fig. 10 | cutoff-8 run low by 2.6 s.e. (c1) [in the cutoff-8 run's own s.e.] | 2.6 | 2.58975 | PASS |
| Sec. VII / Fig. 10 | cutoff-8 run low by 3.2 s.e. (c2/c1) | 3.2 | 3.17177 | PASS |
| Sec. VII / Fig. 10 | cutoff-8 run low by 3.7 s.e. (c3/c1) | 3.7 | 3.68055 | PASS |
| Table V | K=1, U=0.10: exact c3/c1 (finite differences would give 2.4310) | 2.3585 | 2.35849 | PASS |
| Table V | K=3, U=0.05: exact c3/c1 (lab-basis finite differences would give 1.6116) | 1.6393 | 1.63925 | PASS |
| Table V | K=1, U=0.05: mean-field c1 (finite differences would give 0.8587) | 0.8576 | 0.857573 | PASS |
| Table V | K=1, U=0.10: mean-field c1 (finite differences would give 0.9399) | 0.9327 | 0.932697 | PASS |
| SM, validation | K=3 ring (Table V): trajectories c1 = 1.1456(2) | 1.1456 | 1.14555 | PASS |
| SM Figs. S1, S5, S7 / Fig. 7 | emitter (a): max relative error of fourth order, 0.37% | 0.37 | 0.3665 | PASS |
| SM Figs. S1, S5, S7 / Fig. 7 | emitter (b): dressed fourth-order max relative error, 2.0% | 2 | 2.01511 | PASS |
| SM Figs. S1, S5, S7 / Fig. 7 | Kerr switch-on: max error of N<=8 partial sum, 1.46e-3 | 0.00146 | 0.00146283 | PASS |
| SM Figs. S1, S5, S7 / Fig. 7 | atom in cavity: g2(0) on resonance exact, 0.8936 | 0.8936 | 0.893561 | PASS |
| SM Figs. S1, S5, S7 / Fig. 7 | atom in cavity: g2(0) on resonance N<=32, 0.8936 | 0.8936 | 0.893572 | PASS |
| SM Figs. S1, S5, S7 / Fig. 7 | cubic root test |c_N|^(-1/N) at N=10 | 1.074 | 1.0739 | PASS |
| SM Figs. S1, S5, S7 / Fig. 7 | cubic root test |c_N|^(-1/N) at N=20 | 0.864 | 0.864296 | PASS |
| SM Figs. S1, S5, S7 / Fig. 7 | cubic root test |c_N|^(-1/N) at N=30 | 0.789 | 0.788702 | PASS |
| SM Figs. S1, S5, S7 / Fig. 7 | cubic root test |c_N|^(-1/N) at N=40 | 0.758 | 0.758249 | PASS |
| SM Figs. S1, S5, S7 / Fig. 7 | cubic root test |c_N|^(-1/N) at N=50 | 0.707 | 0.707382 | PASS |
| SM Figs. S1, S5, S7 / Fig. 7 | cubic root test |c_N|^(-1/N) at N=60 | 0.664 | 0.663508 | PASS |
