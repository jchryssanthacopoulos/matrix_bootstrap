# To do

*Updated 2026-09-15.*

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
- [ ] Reproduce the single-matrix covariant bootstrap at $N=2,3$; push lifted sectors at $N=3$ to tightness (higher level + ground-state positivity).
- [ ] Port to the 3-matrix model at $N=2$ with $S_3$-flavor blocking; validate against ED.
- [ ] Attempt $N=3$ lifted sectors for the 3-matrix model.

## Findings from re-running the old analyses (2026-09-15)
- [x] ED numbers in `research/notes/` reproduced (single-matrix spectra $N=2,3$; 3-matrix $N=2$: 972 BPS, 243:486:243, 168 levels).
- [x] Section-4 bootstrap validations reproduced (`mm_bootstrap.py`, `sm_charge_profile.py`).
- [x] Single-trace $H$: constant corrected to $-\tfrac32N(N^2-1)$ (docs/derivations.md D1); the fitted $-9(N-1)^2$ fails at $N=4$.
- [x] Covariant single-trace bootstrap at $N=3$: ground-state positivity runs (36 s) and level 4 runs (3 min) — both leave lifted sectors at the trivial bound 0. Diagnosis in the addendum to `research/notes/matrix_syk_covariant_bootstrap.md`.
- [ ] Redesign the covariant bootstrap with adjoint-valued (open-index) operators organised into $SU(N)$ irreps; test on $N=3$ lifted sectors.

## From the second batch of papers (2026-09-15)
- [ ] Implement the Turiaci–Witten chaos test: singular values of $Q_k:\mathcal H_k\to\mathcal H_{k+3}$ per $(k,\text{irrep})$, $r$-ratio vs $\beta=2$ surmise; try $p=2,N=3$ and edge sectors of $p=3,N=3$ by sparse Lanczos on $Q_k^\dagger Q_k$.
- [ ] Compute BPS counts per charge and per $SU(N)$ irrep for the 3-matrix model at $N=2$ and compare with $\cos(\pi k/3)$; state the $N\ge3$ prediction (Q1d) precisely.
- [ ] Refined index $\mathrm{Tr}[(-1)^Fe^{2\pi irN_\Psi/3}\chi_R]$ per irrep (FGMS (5.5) + character insertion) for both models; check saturation at $N=2$.
- [ ] Rebuild the covariant bootstrap with adjoint-valued operators and $U(N)$ tensor-structure decomposition (Cho–Gabai–Sandor–Yin (3.5)); test on the $N=3$ single-matrix lifted sectors.
- [ ] (Later) KMS/thermal bootstrap of the 3-matrix model for $E(\beta)$ and the near-BPS density of states (Cho et al. (1.6), (2.11)).

## From derivation D2 (2026-09-16)
- [ ] Two-particle sector of the 3-matrix model analytically ($H_4$ acts non-trivially for the first time); compare with $N=2$ ED value $E_0(2)=5.1654$.
- [ ] Traceless ($su(N)$) variant of the 3-matrix model: does it concentrate at $N=2$ ($2^9$ states)? Decide which variant the large-$N$ engine should target (trace modes couple for $p=3$, unlike $p=1$).
- [ ] Use $E_0(0)=16N^3-15N$ and the one-particle levels $16N^3-53N,\,16N^3-20N,\,16N^3-16N$ as exact anchors when validating any bootstrap of the 3-matrix model.
