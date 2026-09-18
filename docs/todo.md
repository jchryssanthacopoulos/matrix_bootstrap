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
- [x] **M1** `src/trace_algebra.py` + `tests/test_trace_algebra.py` (2026-09-18). Single primitive = adjacent swap of a routed term (letters in operator order + routing permutation; anticommutator branches rewire the routing by union-find on index edges, closed loops give $N$). Built on it: canonical cyclic rotation, trace-factor ordering, products, daggers, $[H,\cdot]$, finite-$N$ antisymmetriser relations with letter or block entries, and self-symmetry reductions ($\mathrm{Tr}[\ell\ell]=0$, $TT=\tfrac12\{T,T\}$ for odd $T$, periodic words). Coefficients are polynomials in $N$. **Harness:** 1000 random routed terms up to length 6 (worst error $4\times10^{-15}$), traces/products/daggers/commutators, $Q$, $N_\Psi$, $H$(D2.11) and $\{Q,\bar Q\}=H$ *computed by the algebra* — all at machine precision at $N=2$ ($p=3$, both letter bases) and $N=2,3$ ($p=1$). Timings: $[H,\mathrm{Tr}\,w]$ with $|w|=4$ in 2 ms; Gram entries of length 8–10 in $<1$ ms.
- [x] **M2** `src/trace_bootstrap.py` + `scripts/run_trace_bound.py` (2026-09-18): closure-generated monomial variables, reality/EOM/sector/finite-$N$ rows, singlet + adjoint-projected cones, direct Clarabel/SCS interface (no cvxpy), row normalisation, 't Hooft column scaling, pivoted-QR row pruning for Clarabel, `check_exact` diagnostic. V1: Casimir ladder exact at $N=2,3$ (and $N=4$, $k=1$) at level 2.
- [x] **M3** (2026-09-18): $N=2$, $k=2$, level 3: 3.263 / 3.438 / **5.16536** with `finiteN_len` 0/4/6 (Hilbert-space value reproduced); anchors $16N^3-15N$, $16N^3-53N$ reproduced to $10^{-8}$ for $N=2..30$ (Clarabel) and to $10^{-10}$ at $N=100$ (SCS); $N=1000$ only to $10^{-7}$ (SCS inaccurate, not certified). Not done: $k=3$ at $N=2$ (needs level 4).
- [~] **M4** first new numbers (2026-09-18, `research/notes/trace_bootstrap_results.md` §5–6): $E_0(k;N)=16N^3-(15+38k)N$ for $k\le N-1$, tight at level 2 for $k=2$ ($N=3,4,5,10,100$) and $k=3$ ($N=4,5,10$), confirmed by few-body ED (`src/fewbody.py`: $159,660,1545,508$); first interacting sector $k=N$: $E_0(3;3)\ge70.35$ at level 3 (exact 75.79), unchanged by finite-$N$ relations to length 6 → a level-4 deficit. Open: level-4 assembly (restrict length-4 words to the informative cones; eliminate monomials through sector rows before assembly; measure memory first); $k=4$ at $N=4$ and $k=N$ at $N=4$ ($E_0(4;4)$ by few-body ED, dim 194 580, Lanczos).
- [ ] **D5 candidate**: prove $E_0(k;N)=16N^3-(15+38k)N$ for $k\le N-1$ (k antisymmetrised copies of the lowest adjoint one-particle mode annihilated by the quartic interaction; obstruction at $k=N$) and characterise the $k=N$ sector.
- [ ] **M5** large $N$ ($10^2,10^3$; 't Hooft normalisation): fixed-$k$ scaling; $k=0,1$ must match D2 exactly; fit $k\ge2$.
- [ ] **M6** write-up: D5 (fermionic trace algebra, harness-fixed signs) in `docs/derivations.md`; results note; research questions.
- [ ] Every run: memory estimate first, watchdog, budget 10 GB (20 GB ceiling).

## Parallel track: concentration at $N=3$ by cohomology ranks (plan §8)
- [ ] Exact BPS counts per $(k,\text{irrep},z)$ at $N=3$ from ranks of $Q_k$ on the $N_\Psi$ sectors of the $2^{27}$ space (mod-$p$ elimination or zero modes of $Q_k^\dagger Q_k$), with $\mathbb Z_3$ / $U(3)$-Cartan reductions; window width 3 vs 4; compare with the refined index (Q1b), test $1{:}2{:}1$ (Q1d); Turiaci–Witten singular-value statistics (Q2) from the same $Q_k$. Memory estimate before running.
- Division of labour recorded in plan §8: concentration ⇒ cohomology + index; near-BPS energies and $N$-scaling ⇒ bootstrap.
