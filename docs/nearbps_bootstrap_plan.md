# Plan: bootstrapping the near-BPS spectrum from both sides (2026-10-06)

*Status: proposed. Builds on `research/notes/sector_ed_results.md` (exact energies up to the window at $N=2,3$),
`research/notes/trace_bootstrap_results.md` §10 (calibration and its failure near the window), and two papers read
on 2026-10-06: Cho–Gabai–Lin–Yeh–Zheng 2025 (`literature/notes/cho_gabai_lin_yeh_zheng_2025.md`) and Adams 2025
(`literature/notes/adams_2025.md`).*

## Goal

Rigorous bounds on the lowest non-BPS energies $E_0(k,\lambda)$ next to the BPS window, in the low-Casimir
non-singlet irreps where exact diagonalisation finds them (at $N=3$: $\mathbf 8$ at $k=10$,
$\mathbf{10}\oplus\overline{\mathbf{10}}$ at $k=11$). The bounds must hold at $N$ where exact diagonalisation cannot
reach half filling ($N\ge4$), and should eventually give the near-BPS spectral weight seen by simple probes, which
is what the super-Schwarzian predicts.

## Where we stand

| | lower bounds on $E_0(k)$ | upper bounds on $E_0(k)$ |
|---|---|---|
| free sectors $k\le\lfloor N^2/4\rfloor$ | exact (D16) | exact (D16) |
| intermediate sectors | level 3 captures 93%→13% at $N=3$, $k=3..6$ | D16 band bound, weak |
| **near the window** ($N=3$, $k\ge7$) | **0**: fake "BPS-like" functionals ($\phi(H)=0$) survive every constraint, including irrep rows and finite-$N$ relations | **none** |

## The central idea

**Use the BPS states as the reference states.** They are the exact *global* ground states of $H=\{Q,\bar Q\}\ge0$
($E=0$, $QB=\bar QB=0$), exactly as the singlet ground state is in CGLYZ's ungauged matrix QM. Their adjoint-gap
construction then carries over.
- **The bound.** For trial operators $O_j$ of charge $k-k_B$ that carry a BPS state $B$ (in window sector $k_B$)
  into sector $k$,
  $$E_0(k)\ \le\ \Delta^*\equiv\max\{\Delta:\ \mathcal N-\Delta\mathcal M\succeq0\},\qquad \mathcal M_{ij}=\langle B|O_i^\dagger O_j|B\rangle ,$$
  $$\mathcal N_{ij}=\langle B|O_i^\dagger HO_j|B\rangle=\langle B|[Q,O_i\}^\dagger[Q,O_j\}|B\rangle+\langle B|[\bar Q,O_i\}^\dagger[\bar Q,O_j\}|B\rangle .$$
- **The physics reading.** A near-BPS state is a BPS state dressed by an operator that is *almost* $Q$- and
  $\bar Q$-closed, and the gap is bounded by the smallest relative size of that operator's supersymmetry variation.
- **The non-singlet target comes for free.** A low-Casimir $B$ dressed by adjoint-valued words lands in
  low-Casimir irreps ($\mathbf 8\otimes\mathbf 8=\mathbf 1+\mathbf 8+\mathbf 8+\mathbf{10}+\overline{\mathbf{10}}+\mathbf{27}$),
  exactly where the near-window ground states sit.
- **Positivity with charged operators is allowed here.** Ground-state positivity holds for charged, non-singlet
  operators when the reference is a BPS state, but not for sector ground states.

So the BPS-anchored bootstrap supplies **upper** bounds. The sector bootstrap, given additional rows, supplies
**lower** bounds. Together they bracket the near-BPS gap.

## Phases

### Phase 1: calibrate at $N=3$ with exact BPS states (upper bounds; no SDP)

- **1a. Reference states.**
  - Extract explicit BPS states in window sectors $k_B=12,13,14$ with `src/sector_ed.py`.
  - Zero-weight blocks have 0.98–1.1M states. Zero modes come from Lanczos on $H$, or on $H+\mu\hat C_2$ to select
    irreps (the singlet and adjoint BPS states that matter here).
  - Estimate: at most about 3 GB, minutes to hours.
- **1b. Trial spaces.** $V_L=\mathrm{span}\{O\,B\}$: words $O$ up to level $L=1..4$ of charge $k-k_B$
  (for example $k=11$ from $k_B=12$ (charge $-1$) or from $k_B=14$ (charge $-3$; $\bar Q$ itself kills $B$)),
  restricted to weight-zero components.
- **1c. Rayleigh–Ritz.** Solve the generalised eigenproblem of $(\mathcal N,\mathcal M)$ on $V_L$, which bounds
  $E_0(k)$ from above, per irrep after Casimir resolution. Compare with the exact values 0.0614 ($k=11$),
  0.345 ($k=10$), 1.146 ($k=9$).
- **1d. Ground truth for Phase 4.** Exact correlators $G_O(\tau)=\langle B|O^\dagger e^{-\tau H}O|B\rangle$ and
  their spectral weights $|\langle n|O|B\rangle|^2$ on near-BPS states.
- **Decision.**
  - If the upper bounds come within roughly 10–30% of the exact values at $L\le3$–4, near-BPS states are
    well-approximated by low-level dressings of BPS states, which is itself a physics result. Go to Phase 3.
  - If not, the CGLYZ route does not transfer, which is also worth recording.

### Phase 2: lower bounds from the supermultiplet eigen-bootstrap (go/no-go first)

- **The rows.** For a sector ground state $\psi$ below the window, using the multiplet-averaged eigenprojector, add
  three families:
  - lower-member rows, $\phi(QX)=\phi(X\bar Q)=0$;
  - eigen rows, $\phi(X(H-E))=\phi((H-E)X)=0$ at fixed $E>0$, with $E$ scanned;
  - the **cohomology input** that sector $k$ has no BPS state, which excludes $E=0$.
- **Rigour.**
  - The lower-member assumption is removed by $E_0(k)=\min(L(k),L(k-3))$, where $L(j)$ is the lowest
    lower-member energy in sector $j$.
  - Exclusion between grid points of the $E$ scan is made rigorous through dual certificates with margins.
- **2a. Go/no-go at $N=2$, $k=4$** (exact 0.08796, singlet; Hilbert-space engine `src/sector_bootstrap.py`). Does
  the feasible set separate from $E=0$ at level $\le3$? It is about an hour of coding and well under 10 GB.
- **2b. If yes: $N=3$, $k=9..11$** in the trace engine.
  - Build the rows for $\phi(XH)$ and $\phi(X)$ once and combine them at assembly time for each $E$, so the
    program is not rebuilt for every grid point.
  - Add Casimir and cubic-Casimir rows for the ground-state irreps, and highest-weight rows if needed.
  - Estimate: at most about 6 GB per SDP, 20–40 grid points per sector.
- **2c. $N=4$.** Lower bounds at $k=21$ (window $22..26$) need the input "no BPS state below the window" in the
  relevant irreps. That is the band statement, verified only at high Casimir. **At $N\ge4$, lower bounds are
  conditional on it** unless BPS-exclusion certificates can be obtained.

### Phase 3: upper bounds at $N\ge4$ from a bootstrapped BPS functional

- **The reference functional.** Replace the exact $B$ by a bootstrapped BPS functional $\phi_B$ on window
  sector $k_B$. Impose:
  - BPS rows, $\phi_B(QX)=\phi_B(XQ)=\phi_B(\bar QX)=\phi_B(X\bar Q)=0$, as in Lin–Zheng 2024's supersymmetric
    ground-state bootstrap;
  - sector and Casimir rows;
  - positivity of Gram blocks of every charge.
- **The bound.** The largest $\Delta$ for which some feasible $\phi_B$ has $\mathcal N(\phi_B)-\Delta\mathcal M(\phi_B)\succeq0$
  on the trial block is a rigorous upper bound on $E_0(k)$, because the true $B$ is feasible. One feasibility SDP
  per $\Delta$. $\mathcal N$ is a Gram matrix of the words $[Q,O\}$ and $[\bar Q,O\}$, one level up.
- **Order of runs.** Calibrate at $N=3$ against Phase 1, then run $N=4,5$.
- **Deliverable.** Rigorous upper bounds on the near-window gaps as a function of $N$. If they fall with $N$ in the
  natural units, that is the first rigorous evidence in this model for a closing near-BPS gap, the direction the
  super-Schwarzian needs.
- **Main risk.** Fake BPS functionals in the window sector make $\Delta^*$ weak. Mitigations: irrep resolution,
  $\mathbb Z_3$ charge, higher level.

### Phase 4: BPS-anchored two-point correlators (near-BPS spectral weight)

- **The objects.** $G_O(\tau)=\langle B|O^\dagger e^{-\tau H}O|B\rangle=\sum_nw_ne^{-\tau E_n}$ (plus a constant
  if $O$ stays in the window), bounded through CGLYZ's dual "inequalities of motion":
  - reflection positivity;
  - $\partial_\tau\mathcal M=-\mathcal MD$, with $D$ from the trace engine's commutators;
  - ground-state positivity, valid because $B$ is a global ground state;
  - zero-time data from Phase 3, i.e. the "partial knowledge" variant, their App. A.6.
- **Calibration** at $N=3$ against Phase 1d, then $N=4$.
- **Read-out.** Fit the late-$\tau$ bounds to the super-Schwarzian density to estimate its scale against $N$. These
  are estimates, as in CGLYZ §4.5 and Adams.
- **Tooling.** Spline dual with Clarabel or MOSEK; polynomial matrix program with SDPB if needed.

### Phase 5 (optional): thermal route

- **Setup.** The grand-canonical Gibbs state $e^{-\beta(H-\mu N_\Psi)}$, for which KMS holds for every operator.
  Impose it exactly with QICS (Adams), then fit low-$T$ bounds to super-Schwarzian thermodynamics.
- **Dependency.** Only worth attempting once Phase 2 shows that fake BPS-like functionals can be removed, since
  thermal relaxations are likely to inherit them.

## Order and decision points

1. **Phases 1 and 2a in parallel.** Both are cheap and decisive, about 1–2 days.
2. **Phase 1 positive:** Phase 3 at $N=3$, then $N=4,5$, then Phase 4.
3. **Phase 2a positive:** Phase 2b, then 2c, giving two-sided brackets.
4. **Both negative:** stop pushing the bootstrap toward the window. Record why, and put the near-BPS effort into
   exact diagonalisation at $N=3$: singlet towers, and Turiaci–Witten level statistics of $Q$'s singular values.

## Software and memory

- **Phases 1–2:** existing code (`sector_ed.py`, `sector_bootstrap.py`, `trace_bootstrap.py`) plus new rows and a
  trial-space Rayleigh–Ritz.
- **Phase 4:** a new dual (inequalities of motion) module.
- **New dependencies:** MOSEK needs a licence; SDPB is optional; `qics` (pip) only for Phase 5.
- **Memory.** Every run is estimated first and guarded; the estimates above are all within 10 GB.

## What would count as success

1. Rigorous two-sided bounds on near-window gaps at $N=4$, beyond exact diagonalisation.
2. Rigorous upper bounds showing the near-BPS gap shrinks with $N$.
3. Estimates of the near-BPS spectral weight seen by simple probes, to compare with the super-Schwarzian density.

None is guaranteed. Phases 1 and 2a exist to find out quickly which, if any, is within reach.
