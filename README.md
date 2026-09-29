# Simulation code

All code for "Operator-language Feynman rules for driven-dissipative quantum systems:
from mean field to non-Gaussian photon correlations" (P. Ehlers, P. H. Nguyen, D. Soh).

## Layout

    engines/                   the diagram engines and reference solvers (importable modules)
    notebooks/                 one notebook per section of the paper (00-10): the calculation, the regenerated figure or
                               table, and the numbers the text quotes for it; verify_numbers.py checks every number
    notebook/                  the original companion notebook, its cached exact data (data/) and its output figures
    benchmarks/
      cumulant_closures/       Fig. 9: cumulant closures of order 2-4 against the displaced-frame exact reference
      truncation/              Figs. 5 and 8, Table V: fixed-order partial sums, contour cumulants, runtimes, label counts
      trajectories/            Fig. 10: quantum-jump counting statistics of the eight-site chain
        configs/               production and validation run configurations
        est/                   counting-cumulant estimates computed from the jump records
        diag/, exact/          diagram (contour) and exact tilted-Liouvillian references for the chain and for K = 1, 3

Every path in the scripts is resolved relative to this directory, so the tree can be moved as a whole.

## engines/

Library modules

    engine.py                  single-mode diagram engine: vertex table, shift relations, propagators, partial sums
    multimode.py               multimode (momentum-space) engine for rings and chains
    tensor_engine.py           tensor implementation for large K; label enumeration and contraction
    sparse_ring.py             the ring recursion X_N = V G0 X_{N-1} stored once as a sparse matrix over the reachable
                               zero-momentum labels; one run gives every order, observable and coupling (the vertex is
                               proportional to U, so the N-th order term is c_N U^N). check_sparse_ring.py validates it
                               against multimode.Ring with the same pruning (agreement 1e-16)
    lattice.py, chain_diagrams.py   lattice geometry, disorder realizations, chain diagram sets
    fcs.py, fcs_multi.py       counting-field (tilted) engines, single mode and multimode; full
                               Rayleigh-Schroedinger recursion for theta_N; contour_cumulants, cumulants_contour
    gaussian_fcs.py, gaussian_fcs_riccati.py   Gaussian (HFB) counting statistics; tilted Riccati equations
    tilted_exact.py            exact tilted Liouvillian for K = 1, 2, 3; theta_lead(..., complex_out=True)
    cumulant.py, hfb_chain.py, hfb_lattice.py   cumulant closures (order 2 = HFB) for cavity, chain, lattice
    ring_exact.py, ring_ss_evolve.py, sparse_ss.py   exact steady states by diagonalization or sparse time evolution
    finite_time.py, twotime.py   finite-time rules in Laplace space; two-time correlations
    atom_cavity.py             two-level atom in a cavity: dressed lines and emitter moves

Figure and table scripts

    make_figures.py, gen_chain_figs.py

Checks

    check_wick.py              Wick inversion on a squeezed-thermal steady state
    check_inversion.py         eigenoperator inversion
    k3_validate.py             K = 3 ring against the exact steady state

## notebooks/ — reproducing the paper section by section

    00_engine_checks                 engine validation (SM Table SII)
    01_kerr_cavity_fig5              Sec. V.A, Fig. 5
    02_kerr_switch_on_figS1          finite-time rules, SM Fig. S1
    03_cubic_nonlinearity_figS5      SM Fig. S5
    04_atom_in_cavity_figS7          SM Fig. S7
    05_driven_atom_fig7              Sec. V.B, Fig. 7
    06_kerr_ring_fig8                Sec. VI, Fig. 8 (engine cross-checks, K = 3 ring, disorder, timings)
    07_cumulant_closures_fig9        Sec. VI.B, Fig. 9, SM Table SI, and the sparse-matrix engine
    08_counting_statistics_tableV    Sec. VI.B, Table V
    09_eight_site_chain_fig10        Sec. VI.B, Fig. 10 and the quantum-jump benchmark
    10_verification                  runs verify_numbers.py

Each notebook has a Part A that runs the calculation of that section (FAST = True: reduced orders and cutoffs, seconds
to a minute each; FAST = False: the orders and cutoffs of the paper) and a Part B that regenerates the paper's figure or
table from the result files in benchmarks/ and prints the numbers the text quotes next to the values in the files; each
Part B names the script that recomputes the raw inputs. Run them in order from notebooks/ (about 15 minutes in total).

    python3 notebooks/verify_numbers.py --report notebooks/verification_report.md

checks every number quoted in Secs. V-VII of the paper, Tables V and SI, and the quantum-jump section of the
Supplemental Material against the result files (232 checks). Tables V and SI are compared with the LaTeX sources if
PAPER_DIR points to them, otherwise with the snapshots in notebooks/paper_reference/. engines/numbers.json holds the
quoted numbers for Figs. 5 and 7 and SM Figs. S1, S5, S7; make_figures.py rewrites it.

## notebook/

    paper_simulations_run.ipynb    executed, with outputs
    paper_simulations_clean.ipynb  the same without outputs
    paper_simulations.html         read without Jupyter
    data/parameters.json           parameters of every figure and table
    data/revision_params.json      parameters of the benchmark runs
    data/*.json                    cached exact results

Run the notebook from this directory; it imports the modules from ../engines. FAST = True (default)
reproduces the results at reduced order and cutoff in a few minutes; FAST = False uses the orders and
cutoffs of the paper. Its fig_kerr_convergence.pdf is superseded by the one from benchmarks/truncation/.

## benchmarks/

Each directory holds a methods note (the text of the corresponding Supplemental Material section), the
scripts, the raw outputs, and the result files from which the figures and tables were drawn.

    cumulant_closures/   run_cumulants.py, exact_disp.py, diag_run.py  ->  assemble.py  ->  make_fig_g3.py
                         raw inputs diag_U*.json, exactD_U*.json, closure_partial.json, closure_timing_1core.json;
                         result cumulant_closure_results.json
    truncation/          series.py, kerr_trunc.py, exact_pool.py, exact_disp_pool.py, fcs_nz.py, fcs_nz_md.py, k1_radius.py,
                         table5_contour.py, label_counts.py, timing_v2.py, timing_v3.py, sparse_ring_timing.py
                         (raw outputs *_raw*.json, timing*_*.json, pass3/)
                         ->  assemble_truncation.py  (writes truncation_results.json, table5_corrected.json,
                             table5_fragment.tex, timing_results.json, fig_kerr_convergence.pdf, fig_ring.pdf)
    trajectories/        traj_qjump.py (waiting-time MCWF solver, numpy or cupy backend),
                         prod_gpu.py / prod_cpu_pinned.py (production; batch_*.npz and cpubatch_*.npz),
                         traj_fcs_estimator.py (two-window estimator with bootstrap), traj_validate.py,
                         traj_fcs_final.py (comparison and Fig. 10)
                         ->  regenerate_traj_fcs.sh  (writes traj_fcs_results.json, chain8_contour.json, fig_fcs.pdf)
                         The jump records themselves (605 MB) are in the data deposit, not in this repository;
                         regenerate_traj_fcs.sh expects them under trajectories/runs/<run>/; without them it re-assembles Fig. 10
                         from the estimates in est/.

Other scripts in these directories are checks and intermediate versions kept for completeness.

## Environment

Python 3.12 is recommended. The production benchmark environment recorded in
`notebook/data/revision_params.json` used NumPy 2.5 and SciPy 1.18; the code is
compatible with NumPy 2.x / SciPy 1.x. Matplotlib is used for figures and SymPy
for the symbolic Wick-rule check. CuPy 14.2 is optional and is needed only for
the GPU backend of `traj_qjump.py`. Install the CPU/notebook dependencies with
`python -m pip install -r requirements.txt`.
Set OMP_NUM_THREADS=1 for the pinned CPU production runs.

## Sign convention

    Schrodinger picture :  L(rho)   = -i[H,rho] + D(rho)
    Adjoint (Heisenberg):  L^dag(O) = +i[H,O]   + D^dag(O)

Both appear in engines/; the plus sign is correct only for the adjoint. Every tilted (counting-field)
calculation is checked against c1 = kappa <a^dag a>, which fails if the Schrodinger sign is wrong.
