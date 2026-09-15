# Cho, Gabai, Sandor, Yin (2024) — Thermal Bootstrap of Matrix Quantum Mechanics

**File:** `papers/cho_gabai_sandor_yin_2024.pdf`

## Full citation

Minjae Cho, Barak Gabai, Joshua Sandor, Xi Yin, *Thermal Bootstrap of Matrix Quantum Mechanics*, arXiv:2410.04262v3 [hep-th] (27 Mar 2025). Cited as [34] in Lin–Zheng 2025 and [111] in Lin's TASI notes.

## Main question

Can the bootstrap determine *finite-temperature* observables of (ungauged) matrix quantum mechanics — at finite $N$ and in the planar limit — by combining stationarity, positivity and an SDP-representable form of the KMS condition?

## Physical system

Ungauged one-matrix QM $H=\mathrm{Tr}[\tfrac12P^2+V(X)]$, $V=\tfrac12X^2+\tfrac gNX^4$ (1.1), (3.1), traceless Hermitian $X,P$ with $[P_{ab},X_{cd}]=-i(\delta_{ad}\delta_{bc}-\frac1N\delta_{ab}\delta_{cd})$ (3.2); $N\le10$ and $N=\infty$; also $g<0$ metastable case and a preliminary two-matrix model (1.8), (6.1).

## Important definitions

- Thermal expectation $\langle O\rangle_\beta$; stationarity $\langle[H,O]\rangle_\beta=0$ (1.2); positivity $A_{ij}=\langle O_i^\dagger O_j\rangle_\beta\succeq0$ (1.3)–(1.4); $B_{ij}=\langle O_jO_i^\dagger\rangle_\beta$, $C_{ij}=\langle O_i^\dagger[H,O_j]\rangle_\beta$ (1.5).
- **KMS as a matrix inequality** (1.6): $\beta C-A^{1/2}\log(A^{1/2}B^{-1}A^{1/2})A^{1/2}\succeq0$, equivalent to the KMS condition (2.3) given stationarity; $\beta\to\infty$ limit is ground-state positivity $C\succeq0$.
- **Semidefinite relaxation of the matrix logarithm** (§2.1): Gauss–Radau rational approximants $r_{m,k}$ (2.4)–(2.5) bounding $\log$ from above; relative-entropy cone (2.8)–(2.10); SDP representation (2.11) with auxiliary matrices $T_j$, $Z_i$ (Theorem 3 of ref. [36]). Rigorous as long as the quadrature bounds $\log$ from above.
- **Adjoint-valued basis** in the ungauged theory: open-string words $O_i\in\{1,X,P,X^2,XP,\dots\}$ as $N\times N$ matrices (3.4); two-point functions decompose into $U(N)$ tensor structures (3.5), reducing to $\langle\mathrm{Tr}(O^\dagger O')\rangle$ and $\langle\mathrm{Tr}O^\dagger\,\mathrm{Tr}O'\rangle$ blocks (3.7)–(3.10). No gauge constraint $\langle\mathrm{tr}GO\rangle=0$ (footnote 11: the thermal state is $SU(N)$-invariant but $G\rho\ne0$).
- Planar limit: factorisation (3.12), with multi-trace stationarity (3.13) used to reduce variables; remaining nonlinear constraints simply omitted (few).

## Main assumptions

Thermal (KMS) state; ungauged model so that adjoint operators act on the Hilbert space; at large $N$, factorisation.

## Main analytical results

- SDP formulation (2.17), (3.11); two-sided bounds on $E(\beta)$ converging in word length $L$: finite $N\le10$ at $L=8$ (Fig. 2, agreement with Hamiltonian truncation at $N=2$ to $10^{-4}$–$10^{-2}$), planar $L=10$ with bounds differing by $\sim10^{-3}$ (Figs. 3–5).
- Interpolation between the high-$T$ perturbative expansion and the low-$T$ **long-string effective theory** (§4, App. B): the leading correction is set by the adjoint gap $\Delta_1$ (4.11), extracted from the bounds to 0.5–3 %.
- Metastable $g<0$: SDP infeasibility above $T_c^{\rm SDP}$ gives a rigorous upper bound on the critical temperature (Figs. 6–7).
- Gauged MQM: with only traced operators the KMS inequality trivialises in the planar limit ($A=B+O(1/N)$) — an unresolved obstacle (§6).

## Important equations

(1.2)–(1.6), (2.4)–(2.11), (2.17), (3.2), (3.4)–(3.13), (4.11), App. A tables.

## Numerical methods

MOSEK (double precision, fast) and SDPA-DD (double-double) when MOSEK is unstable or for infeasibility certificates; $(m,k)=(3,3)$ quadrature suffices. Costs (App. A): planar $L=10$: 175 variables, $126\times126$ positivity block, $120\times120$ $A,B,C$; symbolic setup 41 min and 46.5 GB; solver seconds–minutes; infeasibility scans hours.

## Relevant figures/results

Fig. 1 (anharmonic oscillator $E(T)$, gap to 0.3 %), Fig. 2 ($E/(N^2-1)$ vs $T$ for $N=2..10,\infty$), Figs. 3–5 (planar bounds vs analytic limits), Figs. 6–7 ($g<0$, $T_c$ bounds), Fig. 8 (two-matrix $L=4$).

## Limitations

Symbolic setup dominates cost and memory; gauged theories not treatable at large $N$ by this route; multi-matrix bases explode with word length.

## Relationship to our project

Three lessons. (1) It is the source and cleanest formulation of ground-state positivity as the $\beta\to\infty$ limit of KMS (the constraint that gives us *upper* bounds), and the only paper that handles the thermal problem rigorously — relevant if we want $E(\beta)$, the entropy, or the near-BPS density of states of the 3-matrix model (a bootstrap test of the super-Schwarzian). (2) Methodologically it is the closest precedent for what our covariant bootstrap must become: an **ungauged** model bootstrapped with **adjoint-valued (open-index) operators**, whose Gram matrix is decomposed into $U(N)$ tensor structures (3.5) — exactly the ingredient found missing in the project's gauge-invariant-only covariant bootstrap. (3) Practical: solver choice (MOSEK vs SDPA-DD), and the warning that the symbolic constraint generation, not the SDP, is the bottleneck.
