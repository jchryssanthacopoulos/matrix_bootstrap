# To do

*Updated 2026-09-18. Current plan of record: `docs/trace_bootstrap_plan.md` (exact-finite-$N$ trace bootstrap) — see the section "Plan of record" below.*

## Setup
- [ ] Create a Python environment (numpy, scipy, cvxpy or an SDP solver) so that `src/` can be run; re-verify the ED numbers recorded in `research/notes/` (single-matrix $E_0(N_\Psi)$ for $N\le4$; 3-matrix $N=2$ BPS counts 243:486:243).
- [ ] Decide on an SDP backend (SCS via cvxpy for prototyping; SDPB/SDPA-GMP for production — Lin–Zheng 2025 App. A warns that these problems are not strictly feasible).

## Theory (docs/derivations.md)
- [ ] Write out the refined definition of R-charge concentration for the matrix models (irreducible complexes labelled by $N_\Psi$ mod 3, $SU(N)$ irrep, flavor) and derive what the single-matrix window (Chen 2.22) implies per complex.
- [ ] Derive the refined index $\mathrm{Tr}[(-1)^F q^{N_\Psi}]$ per $SU(N)$ irrep for the 3-matrix model as a unitary matrix integral with character insertion (Chen §3, Anninos–Silva §4); evaluate at $N=2$ and compare with ED.
- [x] Derive the exact single-trace form of $H=\{Q,\bar Q\}$ for the 3-matrix model and check which parts are Casimirs — done 2026-09-16, `docs/derivations.md` D2 (general cyclic $C$; verified numerically for $(N,p)$ up to $(3,2)$ and $(2,3)$; $H$ is not a Casimir function; flavor symmetry is $\mathbb Z_3$, not $S_3$). Verification: `scripts/verify_trace_hamiltonian.py`.
- [ ] Large-$N$ moments $\mathrm{tr}H,\mathrm{tr}H^2$ per charge sector via planar diagrams (Gaitan et al. §4 method).
- [ ] State precisely when ground-state positivity may be imposed sector-by-sector (only for sector-preserving operators).

## Numerics (after theory)
- [x] Reproduce the single-matrix covariant bootstrap at $N=2,3$; lifted sectors at $N=3$ tight (2026-09-16, sector bootstrap D3: exact in all sectors at $N=3$).
- [x] Port to the 3-matrix model at $N=2$ (flavor is $\mathbb Z_3$, not $S_3$); validated against ED — $k\le3$ exact (2026-09-17, D4).
- [~] $N=3$ lifted sectors for the 3-matrix model: **not via the Hilbert-space bootstrap** (closed by design, $d_k$ up to $2\times10^7$; the "$U(3)$ highest weights" / "never form $H_k$ densely" extensions are retired). Low sectors via the trace engine (plan §4, V3–V4); the window via sparse Lanczos (plan §8).

## Findings from re-running the old analyses (2026-09-15)
- [x] ED numbers in `research/notes/` reproduced (single-matrix spectra $N=2,3$; 3-matrix $N=2$: 972 BPS, 243:486:243, 168 levels).
- [x] Section-4 bootstrap validations reproduced (`mm_bootstrap.py`, `sm_charge_profile.py`).
- [x] Single-trace $H$: constant corrected to $-\tfrac32N(N^2-1)$ (docs/derivations.md D1); the fitted $-9(N-1)^2$ fails at $N=4$.
- [x] Covariant single-trace bootstrap at $N=3$: ground-state positivity runs (36 s) and level 4 runs (3 min) — both leave lifted sectors at the trivial bound 0. Diagnosis in the addendum to `research/notes/matrix_syk_covariant_bootstrap.md`.
- [x] (done, see above)

## From the second batch of papers (2026-09-15)
- [ ] Implement the Turiaci–Witten chaos test: singular values of $Q_k:\mathcal H_k\to\mathcal H_{k+3}$ per $(k,\text{irrep})$, $r$-ratio vs $\beta=2$ surmise; try $p=2,N=3$ and edge sectors of $p=3,N=3$ by sparse Lanczos on $Q_k^\dagger Q_k$.
- [ ] Compute BPS counts per charge and per $SU(N)$ irrep for the 3-matrix model at $N=2$ and compare with $\cos(\pi k/3)$; state the $N\ge3$ prediction (Q1d) precisely.
- [ ] Refined index $\mathrm{Tr}[(-1)^Fe^{2\pi irN_\Psi/3}\chi_R]$ per irrep (FGMS (5.5) + character insertion) for both models; check saturation at $N=2$.
- [x] Rebuild the covariant bootstrap with adjoint-valued operators (2026-09-16: `src/sector_bootstrap.py`, formulation D3). Single matrix: exact in all sectors at $N=3$, 15/17 at $N=4$; edge sector converging slowly with level. See `research/notes/sector_bootstrap_results.md`.
- [x] Symmetry-reduced version (2026-09-17: `src/symmetry_reduction.py`, `symmetry=True` in `SectorSDP`, derivation D4). Invariant $\rho_k=\oplus_R\sigma_R\otimes\mathbf 1$ via the isotypic decomposition of the $k$-block under $SU(2)_{\rm gauge}\times\mathbb Z_3$ (Ward identities are then automatic), $\mathbb Z_3$ Fourier letters grade the cones; streaming Gram accumulation. Three-matrix $N=2$, $k=3$ level 3 now builds in $<1$ GB (was 13–37 GB). $N\ge3$ gauge reduction not implemented (reducer needs $U(N)$ highest weights).
- [x] Prune linearly dependent EOM rows before the SDP (done 2026-09-16; $23\,970\to226$). Remaining solver-memory driver is the PSD cone size (interior point $\sim O(m^2)$); use SCS or shrink cones by symmetry.
- [x] Island/archipelago scans implemented (`scripts/scan_sector.py`, `scripts/plot_scans.py`); results in `research/notes/sector_bootstrap_results.md` §7.
- [x] (duplicate of the item above; confirmed done 2026-09-17 — in the symmetry-reduced path the EOM span is accumulated as a $2D\times2D$ Gram matrix and an orthonormal basis of it is used directly, so no redundant rows ever reach the solver.)
- [ ] Memory: the cvxpy→SCS pipeline triples the coefficient data (peak 18.3 GB for $k=3$ level 4 vs 6.4 GB estimated). To go further (level 4 with $L_{\rm sing}=4$, $k=4$ level 4, or $N=3$) feed SCS/SDPA directly from the coefficient arrays, or store coefficients in float32 for the solver stage; recalibrate `plan_sector` (currently: peak ≈ 3× estimate at level 4).
- [ ] $k=4$ (almost-BPS, $E_0=0.08796$): level 3 gives 0; level 4 needs the memory item above. Also try level $(4,3,4)$ with only the $q=0$ cones at length 4 as a cheaper intermediate.
- [~] Faster convergence in window-edge sectors (assessed 2026-09-17, D4.4): gauge Ward identities are implied by the invariant functional (verified: zero extra rows); $Q$-descendant words at a fixed length are linear combinations of words already present and cannot enlarge a Gram cone, so they help only as a sparser *subset* at a higher length — low priority. Still open: SDPB/SDPA-GMP for the 'inaccurate' SCS solves.
- [ ] (Later) KMS/thermal bootstrap of the 3-matrix model for $E(\beta)$ and the near-BPS density of states (Cho et al. (1.6), (2.11)).

## From derivation D2 (2026-09-16)
- [ ] Two-particle sector of the 3-matrix model analytically ($H_4$ acts non-trivially for the first time); compare with $N=2$ ED value $E_0(2)=5.1654$.
- [ ] Traceless ($su(N)$) variant of the 3-matrix model: does it concentrate at $N=2$ ($2^9$ states)? Decide which variant the large-$N$ engine should target (trace modes couple for $p=3$, unlike $p=1$).
- [ ] Use $E_0(0)=16N^3-15N$ and the one-particle levels $16N^3-53N,\,16N^3-20N,\,16N^3-16N$ as exact anchors when validating any bootstrap of the 3-matrix model.


## Plan of record (2026-09-18): exact-finite-$N$ trace bootstrap — `docs/trace_bootstrap_plan.md`
Same constraint set as the (now exact at $N=2$) sector bootstrap, re-expressed as identities among trace-word expectation values with $N$ as a numerical coefficient; cost set by level, not by $N$. Reach: $k\lesssim4$–$5$ at any $N$ (fixed-$k$ large-$N$ regime), not the half-filling window.
- [ ] **M1** `src/trace_algebra.py`: letters, open/trace words, monomials with fermionic parity; canonical cyclic form via (P.1) with double-trace terms; trace-product reordering with merging terms; $[H,\cdot]$ from D2.11/D2.13; sector rows $\phi((N_\Psi-k)X)=0$; finite-$N$ relations (antisymmetriser over $N+1$ indices). **Harness:** every emitted identity verified as an explicit operator equation at $N=2$ ($p=3$) and $N=2,3$ ($p=1$) with `word_matrix`. Acceptance: no failures on $10^3$ random words up to length 6.
- [ ] **M2** `src/trace_bootstrap.py`: closure-generated variables with $N$-polynomial coefficients; evaluation at numerical $N$; independent rows / free-variable parametrisation (sparse, numerical; exact elimination fallback); singlet + adjoint-projected cones; **direct** SCS/Clarabel arrays or SDPA-format export (no cvxpy in the hot path — this also unblocks the Hilbert-space $k=4$ level-4/5 run); $\phi_{\rm GS}$ feasibility diagnostic. Acceptance: single-matrix Casimir ladder tight at $N=3$, level $\le3$.
- [ ] **M3** three-matrix validation: $N=2$: $98,22$; $5.16536$ (level 3); $1.22706$ (level 4) with `finiteN_len` reported, and the deficit without finite-$N$ relations; $N=3$: $387,273$.
- [ ] **M4** first new numbers: $E_0(k;N)$ for $k=2,3,4$, $N=3..10$, level 4 (5 if affordable); ED cross-checks in the small sectors ($d_k=\binom{3N^2}{k}$); test the "level $\approx k+1$" hypothesis.
- [ ] **M5** large $N$ ($10^2,10^3$; 't Hooft normalisation): fixed-$k$ scaling; $k=0,1$ must match D2 exactly; fit $k\ge2$.
- [ ] **M6** write-up: D5 (fermionic trace algebra, harness-fixed signs) in `docs/derivations.md`; results note; research questions.
- [ ] Every run: memory estimate first, watchdog, budget 10 GB (20 GB ceiling).

## Parallel track: sparse Lanczos at $N=3$ (plan §8)
- [ ] Matrix-free $H$ on $N_\Psi$ sectors of the $2^{27}$ space with $\mathbb Z_3$ and $U(3)$-Cartan reductions; $E_0(k)$ for all $k$ (window width 3 vs 4 at $N=3$ — the next concentration data point); per-irrep BPS counts (Q1b–d); Turiaci–Witten singular-value statistics of $Q_k$ (Q2). Memory estimate before running ($2\times10^7$-dimensional vectors, $0.3$ GB each).
