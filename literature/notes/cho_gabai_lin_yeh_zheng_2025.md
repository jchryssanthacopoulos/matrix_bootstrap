# Cho, Gabai, Lin, Yeh, Zheng (2025) — Bootstrapping Euclidean Two-point Correlators

**File:** `papers/cho_gabai_lin_yeh_zheng_2025.pdf` (56 pp., read in full 2026-10-06)

## Full citation

Minjae Cho, Barak Gabai, Henry W. Lin, Jessica Yeh, Zechuan Zheng, *Bootstrapping Euclidean Two-point
Correlators*, arXiv:2511.08560v3 [hep-th] (30 Jun 2026). Cited as ref. [28] by Klebanov–Lin–Meshcheriakov 2026 for
the adjoint gap of matrix QM. Probably the source of Lin's talk on "bootstrapping non-singlet physics" (not
verified).

## Main question

Can the quantum-mechanical bootstrap constrain *dynamical* observables, namely Euclidean two-point correlators
$\langle\bar{\mathcal O}_i(\tau)\mathcal O_j(0)\rangle$ in the ground state or a thermal state, with rigorous bounds
at continuous $\tau$? Can the bounds then be used to extract the spectrum and matrix elements of non-singlet
(adjoint) states?

## Physical system

- Warm-up: the anharmonic oscillator $H=\frac12p^2+\frac12x^2+gx^4$ (64).
- Main application: the **ungauged** one-matrix QM $H=N^2\,\mathrm{tr}(\frac12P^2+V(X))$ (2), with
  $V=\frac12X^2+gX^4$ (65), in the 't Hooft normalisation of footnote 3. All $U(N)$ representations are in the
  Hilbert space. The singlet sector is free fermions; the adjoint sector is governed at large $N$ by the
  Marchesini–Onofri equation (175).
- Couplings: $g=1$, and $g\to g_c=-\sqrt2/(6\pi)\approx-0.075026$ (the $c=1$ double-scaling point).

## Important definitions

- **Variables** (3)–(4). $\mathcal M_{ij}(\tau)=\langle\Omega|\bar{\mathcal O}_i(\tau)\mathcal O_j(0)|\Omega\rangle$,
  with $\bar{\mathcal O}(\tau)=e^{\tau H}\bar{\mathcal O}e^{-\tau H}$. Time-translation invariance gives
  $\mathcal M_{ij}(\tau)=\langle\Omega|\bar{\mathcal O}_i(\tau/2)\mathcal O_j(-\tau/2)|\Omega\rangle$.
- **Reflection positivity** (5): $\mathcal M(\tau)\succeq0$ for all $\tau\ge0$.
- **Heisenberg equations** (6)–(8). Write $[H,\mathcal O_i]=\mathcal O_kD_{ki}$, with the basis closed under
  commutation up to truncation. Then $\partial_\tau\mathcal M=-\mathcal MD$ and $D^\dagger\mathcal M-\mathcal MD=0$.
- **Ground-state positivity** (9)–(11). $\mathcal N_{ij}(\tau)=\langle\Omega|\bar{\mathcal O}_i(\tau/2)[H,\mathcal O_j(-\tau/2)]|\Omega\rangle\succeq0$,
  with $\mathcal N=-\partial_\tau\mathcal M=\mathcal MD$. It "uses the fact that the ground state has lower energy
  than the state $\mathcal O_j(-\tau/2)|\Omega\rangle$", so it needs $\Omega$ to be the *global* ground state of the
  space the operators reach.
- **Symmetry constraints** (13)–(15); time reversal makes $\mathcal M$ real symmetric (12).
- **Primal problem** (16): minimise $\mathrm{Tr}\,\hat O\mathcal M(T)$ subject to the constraints above. Convex,
  but infinite-dimensional because the variables are functions of $\tau$.
- **Dual problem** (17)–(18), with Lagrange multipliers $\lambda_A,\lambda_D,\lambda_G$ as functions of $\tau$. The
  Heisenberg equations become **"inequalities of motion"**:
  $$\lambda_AD^\dagger-D\lambda_A+\tfrac12(D\lambda_G+\lambda_GD^\dagger+D\lambda_D+\lambda_DD^\dagger)+\partial_\tau\lambda_D\succeq0,\qquad -\lambda_G\succeq0 .$$
  Any feasible dual gives a rigorous bound (weak duality).
- **Finite dual.**
  - B-splines: expand $\lambda$ in a positive basis; coefficient matrices PSD then suffice. Equation (21),
    solved with MOSEK.
  - Polynomials: $\tau=Ty/(1+y)$ gives a polynomial matrix program (22), (125)–(134), solved with SDPB.
- **Level.** $\ell(X)=1$, $\ell(P)=2$. "Level $L$" means $\mathcal M$ built from operators up to level $L/2$, and
  $\mathcal N$ to level $L+1$.
- **Thermal version** (23)–(40). KMS gives $\mathcal M_{ij}(\beta)=\mathcal M_{\bar j\bar i}(0)$ (24), with
  multipliers on the two intervals $[0,T]$ and $[T,\beta]$.
- **Zero-time input** (Table 1): full knowledge of $\mathcal M(0)$ (App. A.2/A.4), bounds from a one-point
  bootstrap (A.6, eq. (99)), or none (A.3/A.5).

## Main assumptions

- A stationary state: the ground state, or a thermal state.
- For the ground-state positivity used in the adjoint application, that the singlet ground state is the global
  ground state of the ungauged Hilbert space.
- Large-$N$ factorisation, only in the zero-time data (App. G) and in some analytic bounds (footnote 17).

## Main analytical results

- **Log-convexity** (41)–(44). The connected correlator satisfies $(\log G_c)''\ge0$, hence
  $G_c(\tau)\ge G_c(0)\,e^{-\mu\tau}$ with $\mu=-G'_c(0^+)/G_c(0)$. For MQM (46),
  $\langle\mathrm{tr}X(\tau)X(0)\rangle\ge\langle\mathrm{tr}X^2\rangle e^{-\mu\tau}$ with
  $\mu=1/(2\langle\mathrm{tr}X^2\rangle)$, which gives an **upper bound on the adjoint gap**,
  $\Delta_{\rm adj}\le\mu$. Saturated by the harmonic oscillator (48)–(50).
- **New derivation of energy–entropy balance** (51)–(54): log-convexity plus KMS gives
  $\log\frac{\langle\bar OO\rangle_\beta}{\langle O\bar O\rangle_\beta}\le\beta\frac{\langle\bar O[H,O]\rangle_\beta}{\langle\bar OO\rangle_\beta}$.
- **High-temperature limit** (55)–(63): the KMS two-point bootstrap reduces to the loop / Schwinger–Dyson
  bootstrap of the classical (matrix) integral.
- **Adjoint gap as a generalised eigenvalue problem** (71)–(80).
  - Spectral decomposition: $G_c(\tau)=\sum_n|\langle\Omega|\bar{\mathcal O}_{ba}|n,ab\rangle|^2e^{-\Delta_n\tau}$.
  - Single operators give $\Delta_{\rm gap}\le-G'_c(0)/G_c(0)$ (72), e.g. $\Delta_1\le1/(2\langle\mathrm{tr}X^2\rangle)$ (73)
    and $\Delta_2$ (76).
  - Superposing operators, $\Delta_{\rm gap}\le\bar\alpha\tilde{\mathcal N}\alpha/\bar\alpha\tilde{\mathcal M}\alpha$
    for every $\alpha$. So the best bound is the largest $\Delta$ with $\tilde{\mathcal N}(0)\succeq\Delta\tilde{\mathcal M}(0)$ (80),
    a rigorous **upper** bound (attributed to Lanzette–Zheng, in prep.).
  - Closed forms at level 6 (83) and level 8 (86).
- **Extremal functional** (91)–(92): the generalised eigenproblem $\mathcal M^{-1/2}\tilde{\mathcal N}\mathcal M^{-1/2}v=\Delta v$
  estimates the excited levels and matrix elements $|\langle\Omega|\mathcal O|n\rangle|^2$. These are estimates,
  not bounds.
- **Physics.**
  - The adjoint gap stays finite as $\mu\to0$ (§5). Rare low-energy "short" long strings sit below the typical
    $\log\mu$ scale.
  - Conjecture: $\Delta$ is non-analytic, $\propto C\mu^k\log\mu$ with $1\le k\le2$.

## Important equations

(3)–(11), (16), (18), (21)–(22), (24), (41)–(44), (46), (51)–(54), (71)–(80), (83), (86), (91)–(92), (99), (109),
(113), (125)–(134), (175)–(177), (178)–(184).

## Numerical methods

- Dual SDPs: MOSEK (B-splines; default $N=20$ knots, degree 10) and SDPB (polynomial matrix program, degree
  $d=16$).
- Ground-state zero-time data come from the exact large-$N$ free-fermion solution (App. F/G). Thermal zero-time
  data are bounded by the level-14 EEB bootstrap of Cho–Gabai–Sandor–Yin.
- App. A.1: operator-basis truncation, with constraint spaces built from $\mathcal D_\pm[M]=DM\pm MD^\dagger$.
- App. E: positivity imposed beyond $\tau=T$ gives no improvement (Table 4, to $10^{-10}$).
- App. D: a rigorous (spline) version of the Lawrence–McPeak–Neill time-dependent bootstrap.
- Arbitrary precision was needed for the level-16 gap.

## Relevant figures/results

- Table 2: $\langle x(5)x(0)\rangle$ for the anharmonic oscillator. The level-10 window is
  $[0.000315160730,\,0.000315201914]$, against the exact $0.000315160732$.
- Fig. 6: $\langle\mathrm{tr}X(\tau)X(0)\rangle$ bounds at $g=1$ and $g\approx g_c$.
- Figs. 8–11: adjoint-gap upper bounds against Marchesini–Onofri. The level-$L$ estimates converge to $10^{-14}$ by
  $L=16$.
- Table 3 at $g=-0.075026$:
  - $\Delta_1$: rigorous upper bound $0.741573662448591$ at level 16; Marchesini–Onofri gives $0.74158$; fits give
    $0.74436$ and $0.74199$.
  - $\Delta_2$: upper bound $1.26122$; Marchesini–Onofri gives $1.26130$.
  - $\Delta_3$: extremal-functional estimate $1.68192$ (not rigorous); Marchesini–Onofri gives $1.68218$.
- Figs. 13–14: thermal correlators, tighter than Monte Carlo (App. H).

## Limitations

- **One direction only.** The rigorous gap statements are upper bounds (variational in content). Levels and matrix
  elements beyond the gap come from fits or extremal functionals and are not rigorous. A *lower* bound on the MQM
  adjoint gap is cited from elsewhere (Gross–Klebanov 1990, their (4.16); footnote 10), not obtained by the method.
- **Easy test case.** One solvable matrix model (exact singlet data), bosonic, large $N$. Multi-matrix and chaotic
  models are future work.
- **A requirement on the reference state.** Ground-state positivity with non-singlet operators needs the
  reference to be the global ground state, true here because the model is ungauged.
- Convergence and strong duality are unproven in general.

## Relationship to our project

1. **The adjoint gap is the template for our near-BPS gap.**
   - In our model the global ground states are the **BPS states** $B$: $E=0$ exactly, and $QB=\bar QB=0$.
   - The lowest non-BPS energy in sector $k$ (irrep $\lambda$) is the analogue of their adjoint gap, reached by
     charged operators acting on $B$.
   - Their (80) becomes $E_0(k,\lambda)\le\max\{\Delta:\ \mathcal N-\Delta\mathcal M\succeq0\}$, with
     $\mathcal M_{ij}=\langle B|O_i^\dagger O_j|B\rangle$ and $\mathcal N_{ij}=\langle B|O_i^\dagger HO_j|B\rangle$.
   - Because $QB=\bar QB=0$,
     $$\mathcal N_{ij}=\langle B|[Q,O_i\}^\dagger[Q,O_j\}|B\rangle+\langle B|[\bar Q,O_i\}^\dagger[\bar Q,O_j\}|B\rangle .$$
     So the bound is the smallest ratio "norm of the supersymmetry variation / norm" over the trial operators.
     Near-BPS states are made by *almost $Q$-closed* operators acting on BPS states.
   - This gives rigorous **upper** bounds, the direction our sector bootstrap cannot supply.
2. **Ground-state positivity is valid here with charged operators.** It holds for every operator when the
   reference is a BPS state. That removes the restriction to charge-neutral operators that applies to sector ground
   states (noted 2026-10-06). For BPS references it follows from positivity of $Q$-exact words (Lin–Zheng 2024,
   App. A).
3. **A route to the near-BPS density.** Correlators $\langle B|O^\dagger e^{-\tau H}O|B\rangle$ decompose over the
   near-BPS states, with weights $|\langle n|O|B\rangle|^2$. Their dual "inequalities of motion" give rigorous
   bounds at continuous $\tau$, which is the natural way to reach the near-BPS density seen by simple probes. Within
   the window the correlator has a constant BPS piece plus decaying near-BPS terms.
4. **Tools.** Polynomial matrix programs with SDPB; B-splines with MOSEK; arbitrary precision for high-level gaps.
   Our Clarabel/SCS stack would need extending.
5. **Not a replacement for lower bounds.** Lower bounds still need the sector bootstrap with new constraints
   (eigen rows, lower-member rows, cohomology input).
