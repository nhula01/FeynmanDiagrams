# Higher-order cumulant closures for the K = 3 Kerr ring (Supplemental Material note)

## Model and observables

We consider the uniform three-site ring of Sec. [ring] with Delta = -1, J = 0.4, eta = 1, kappa = 1 and U between 0 and 0.3
(parameter set `ring3` of revision_params.json). The observables are the on-site correlators
g2(0) = <a0^dag^2 a0^2>/<a0^dag a0>^2 and g3(0) = <a0^dag^3 a0^3>/<a0^dag a0>^3.

## Construction of the closures

The closure code works for a general bosonic lattice with Lindblad dissipation, and it does not use translation symmetry.
A normal-ordered monomial prod_i a_i^dag^{m_i} prod_i a_i^{n_i} on K sites is stored as the multi-index pair (m, n), and
operators are sparse dictionaries of such keys. Products are brought back to normal order site by site with
a^n a^dag^m = sum_k k! C(n,k) C(m,k) a^dag^{m-k} a^{n-k}. For every monomial O of total order 1 <= m+n <= N we evaluate
d<O>/dt = <i[H,O] + sum_j kappa (a_j^dag O a_j - {a_j^dag a_j, O}/2)>. The result is a linear combination of moments of order up to N+2.
For the Kerr ring these come from the quartic interaction and the hopping, which conserves order.
Each moment of order p > N is then replaced by its cumulant expansion with every joint cumulant of order larger than N set to zero:
<X_1...X_p> = -sum_{pi, |pi| >= 2} (-1)^{|pi|-1} (|pi|-1)! prod_{B in pi} <X_B>, where the sum runs over set partitions of the p factors.
Blocks that are still of order larger than N are expanded recursively, so the closed moments end up as polynomials in the retained
variables. At order N = 2 this is the Gaussian (time-dependent Hartree-Fock-Bogoliubov) theory of the main text. At N = 3 and N = 4
it is the standard higher-order truncation of Refs. [kira2008cluster, plankensteiner2022quantumcumulants]. For K = 3 the closures
contain 27, 83 and 209 complex moment variables at N = 2, 3, 4. The observables <a0^dag^2 a0^2> and <a0^dag^3 a0^3> are either
variables or, when their order exceeds N, evaluated from the same cumulant polynomials. At N = 2 this reduces g3 to Wick's theorem.

The steady state is obtained by integrating the closed equations from the product coherent state of the linear (U = 0) problem to
t = 300/kappa (DOP853, rtol 1e-10, atol 1e-12). The end point is then polished with a Newton-type root search (MINPACK hybrd). We
accept a solution when the residual of the closed equations is below 1e-8. At every U the residual was below 3e-12. We check
linear stability through the finite-difference Jacobian of the closed equations at the fixed point, and confirm that the solution
is translation invariant although this was never imposed.

## Exact reference

The converged reference is the steady state of the full Liouvillian. We write a_j = alpha_j + b_j, with alpha the classical
(mean-field) steady state, and truncate the fluctuation Fock space at total fluctuation number sum_j n_j(b) <= Nt. The displacement
is a unitary change of frame, so the Liouvillian is unchanged and only the truncation differs. Because the state is close to
coherent (<sum_j b_j^dag b_j> = 0.038 at U = 0.3), the truncation error falls much faster with Nt than it does in the undisplaced
basis. With P the projector on the truncated space, every normal-ordered product obeys P a^dag^k a^l P = (P a^dag P)^k (P a P)^l,
so the Hamiltonian, the jump terms and the observables are the exact projections. The density matrix is evolved to its fixed point
with a fixed-step fourth-order Runge-Kutta scheme applied matrix-free (sparse operators acting on a dense rho). The step is set from a
power-iteration estimate of the spectral radius. The Runge-Kutta map is a polynomial p(hL) with p(0) = 1, so its fixed point is the
exact null vector of L for any stable step. We iterate until max|L rho|/max|rho| < 1e-11 and raise Nt, warm-starting from the
previous cutoff, until successive g3 values agree to 1e-8 (1e-10 for U <= 0.125) or a time limit is reached. The final cutoffs are
Nt = 11 (U = 0.1), 13 (0.15), 16 (0.2), 19 (0.25) and 21 (0.3). The larger of the last two changes in g3 is at most 1e-8 up to U = 0.25
and 8.5e-7 at U = 0.3. That is at least three orders of magnitude below any method error discussed here. At U = 0.02 and 0.025 the
result agrees with an independent calculation in the undisplaced per-site cube n_j < 10 to 2.4e-7 in g3, which is within that
calculation's own cutoff error (4e-6). A reference with a lower cutoff gives g3 of about 2.52 at U = 0.3, against the converged
2.52728.

## Validation of the closure code

At N = 2 our general code reproduces cumulant.py (ring_cumulant_ss) to 4e-13 in n and 5e-12 in g2 over all 14 couplings.
It also reproduces hfb_chain.py in n and g2. It does not reproduce the g3 printed by hfb_chain.py, because that routine has a
combinatorial error in the Gaussian three-photon moment: it uses 18 n|m|^2 in <b^dag^3 b^3>, where Wick's theorem gives
<b^dag^3 b^3> = 6 n^3 + 9 n |m|^2 (n = <b^dag b>, m = <b b>). We checked the coefficient against a displaced squeezed thermal state
in QuTiP. The error makes the Gaussian g3 too large by 6e-5 at U = 0.1, 1.5e-3 at U = 0.2 and 1.2e-2 at U = 0.3. The second-order
closure curve in Fig. 9 uses the corrected value. At U = 0 all three orders reproduce the coherent state
(|g2 - 1|, |g3 - 1| < 2e-14, and n = |alpha|^2 to 1e-16). At U = 0.02 the relative error in g3 is 3.3e-4, 1.1e-6 and 3.6e-8 at
N = 2, 3, 4, so each order gains between one and three decades in the weakly nonlinear limit.

## Stability of the higher-order closures

The known failure mode of cumulant truncations beyond second order is a loss of positivity that produces unstable or unphysical
fixed points. We did not find it in this parameter range. The N = 3 and N = 4 closures converge from the coherent-state start at every
U up to 0.3, the imaginary parts of <a0^dag a0>, <a0^dag^2 a0^2> and <a0^dag^3 a0^3> vanish to 1e-13, and the slowest Jacobian eigenvalue stays at
Re lambda = -0.50 to -0.53 kappa. The breakdown that does appear is a loss of convergence in the order: from U = 0.25 on
the fourth-order g3 is worse than the third-order one. The g3 errors at U = 0.25 are 2.0e-2 at N = 3 and 2.7e-2 at N = 4; at U = 0.3 they
are 1.8e-2 and 6.8e-2. The g2 of the fourth-order closure remains the most accurate of all methods up to U = 0.25.

## Accuracy and cost at the couplings of Fig. 9

Relative errors against the exact steady state:

| U/kappa | exact g3 | closure N=2 | closure N=3 | closure N=4 | diagrams N<=2 | diagrams N<=4 | diagrams N<=10 |
|---|---|---|---|---|---|---|---|
| 0.10 | 1.312582 | 1.2e-2 / 1.0e-3 | 8.5e-4 / 9.6e-5 | 1.8e-4 / 1.4e-6 | 1.7e-2 / 2.9e-3 | 1.0e-3 / 1.1e-4 | 1.8e-6 / 1.5e-7 |
| 0.15 | 1.538157 | 3.1e-2 / 3.9e-3 | 4.7e-3 / 6.2e-4 | 1.8e-3 / 1.4e-5 | 5.3e-2 / 9.8e-3 | 6.5e-3 / 7.5e-4 | 1.6e-4 / 1.4e-5 |
| 0.20 | 1.831380 | 6.1e-2 / 9.6e-3 | 1.3e-2 / 2.2e-3 | 8.3e-3 / 1.3e-5 | 1.1e-1 / 2.2e-2 | 2.0e-2 / 2.4e-3 | 3.0e-3 / 2.8e-4 |
| 0.25 | 2.181138 | 8.8e-2 / 1.7e-2 | 2.0e-2 / 5.3e-3 | 2.7e-2 / 1.7e-4 | 1.8e-1 / 3.9e-2 | 3.2e-2 / 3.7e-3 | 1.9e-2 / 1.7e-3 |
| 0.30 | 2.527280 | 8.9e-2 / 2.3e-2 | 1.8e-2 / 1.2e-2 | 6.8e-2 / 2.4e-3 | 2.2e-1 / 5.4e-2 | 2.3e-2 / 9.1e-4 | 5.1e-2 / 2.6e-3 |

Each entry gives the error on g3, then the error on g2. The wall-clock time per evaluation on one core is 0.13 s (N = 2),
0.6 s (N = 3) and 5 s (N = 4) for the closures, counting generation of the equations and the steady-state solve. It is 1.7 s
(through second order) and 34 to 38 s (through fourth order) for the diagrams, counting the n, g2 and g3 series. Both sets of
timings were taken in the same process, one after the other. A single-thread repetition later in the day, on an oversubscribed
node, scattered by up to a factor of ten between identical calls but did not change the ordering. The exact steady state takes
90 to 600 s on one core in the displaced basis and grows exponentially with K.

Against the second-order closure, the fourth-order diagrams reduce the g3 error by a factor of 12 at U = 0.1 and by factors of
2.8 to 4.8 for U = 0.15 to 0.3, as stated in the main text. Against the higher-order closures the comparison reverses. On g3 the
better of the third- and fourth-order closures is more accurate than the fourth-order diagrams at every coupling of Fig. 9. The
diagram error is larger by a factor of 5.5 at U = 0.1, 3.7 at U = 0.15, 2.4 at U = 0.2, 1.6 at U = 0.25 and 1.25 at U = 0.3. On g2 the
fourth-order closure is better by factors of 20 to 190 up to U = 0.25. Only at U = 0.3 do the diagrams win on g2, by a factor of 2.6.
The third-order closure costs about 60 times less than the fourth-order diagrams and the fourth-order closure about 7 times less.

Carried to tenth order (the partial sums through N = 10 of n, n2 and n3 from ../truncation/, scored against the converged
reference above), the diagrams have a relative g3 error of 1.8e-6 at U = 0.1, 1.6e-4 at U = 0.15, 3.0e-3 at U = 0.2, 1.9e-2 at
U = 0.25 and 5.1e-2 at U = 0.3. This beats the best closure on g3 by factors of 103, 11, 2.8, 1.8 and 1.05 at U = 0.1, 0.15, 0.2,
0.225 and 0.25, and is worse by factors of 1.75 at U = 0.275 and 2.8 at U = 0.3. On g2 the tenth-order diagrams are better than the
fourth-order closure only up to U = 0.125, and on the photon number only up to U = 0.1. The tenth-order values and errors are stored
under diagrams_order10 in cumulant_closure_results.json.

Verdict: for the K = 3 ring the fourth-order diagrams do not beat the best cumulant closure. The best closure has a smaller g3 error
by a factor of 1.25 (U = 0.3) to 5.5 (U = 0.1) at 7 to 60 times lower cost, and the diagrams are ahead only on g2 at U = 0.3,
by a factor of 2.6. The tenth-order diagrams beat the closures on g3 for U <= 0.25, at a cost (4.2 s for all couplings in the
sparse-matrix form) below one fourth-order closure evaluation.

## Files

cumulant_closure.py (moment-hierarchy generator and closure), run_cumulants.py (closure and diagram driver),
exact_disp.py (displaced-basis exact reference), diag_run.py and closure_timing.py (diagram values and timing repeats),
assemble.py (writes cumulant_closure_results.json), make_fig_g3.py (fig_g3.pdf). All are in this directory.
