# Research questions

*Last updated 2026-09-15 (after the literature review; see `literature/synthesis.md`). Keep this file current.*

## Primary question

**Q1. Does Chen's 3-matrix "matrix SYK" model exhibit R-charge concentration for $N\ge3$?**

Model: $p=3$ $U(N)$-adjoint complex fermion matrices, $Q=\sum_{1\le i\le j\le k\le3}\mathrm{Tr}[\Psi_i\Psi_j\Psi_k]$, $H=\{Q,\bar Q\}$ (Chen 2025, eq. 4.3). Chen reports only that an $N=2$ simulation "suggests" concentration.

Sharp formulation (Chang–Chen–Sia–Yang 2024, Conjecture 1 and §2): in every irreducible cochain complex — labelled by $N_\Psi$ mod 3, the $SU(N)$ irrep, and flavor ($S_3$) quantum numbers — with macroscopic index, BPS states occupy a *single* degree. Observable proxy: the sector-resolved ground energy $E_0(N_\Psi,\mathbf r)$; concentration ⇔ $E_0>0$ outside a window of at most three consecutive $N_\Psi$ (one per $\mathbb Z_3$ class) per irrep. Baseline: the single-matrix model has $E_0=0$ on $N+1$ consecutive sectors, all in $r_*$ (Chen (2.22)) — non-concentrated for $N\ge3$.

Sub-questions:
- Q1a. Can the bootstrap give rigorous lower bounds $E_0(N_\Psi)>0$ in lifted sectors at $N$ beyond exact diagonalisation ($N=3$: $2^{27}$ states)?
- Q1b. Are refined indices $\mathrm{Tr}[(-1)^Fq^{N_\Psi}]$ per $SU(N)$ irrep saturated by BPS counts (index saturation ⇔ concentration, Chang et al. §5.1)?
- Q1c. How does the number of BPS states scale with $N$ ($e^{cN^2}$ multiplicity, as for a genuine black-hole-like sector — Chang–Lin 2024 Conjecture 3 — versus a large irrep with small multiplicity as in the single-matrix model)?
- Q1d. **Do the BPS counts in the three charges of a concentrated complex stay in ratio $1{:}2{:}1$ as $N$ grows?** This is the $\hat q=3$ super-Schwarzian prediction $N_{\rm BPS}(k)\propto\cos(\pi k/3)$ (Turiaci–Witten (3.10)), realised exactly in $\mathcal N=2$ SYK (FGMS (5.7)), and matched by the project's $243{:}486{:}243$ at $N=2$; the single-matrix binomial weights depart from it for $N\ge3$.

## Secondary questions

**Q2. Is the 3-matrix model chaotic?** Primary diagnostic: the $r$-ratio statistics of the *singular values of $Q_k$* between adjacent charge sectors $(k,k+3)$ at fixed $SU(N)$ irrep, compared with the Altland–Zirnbauer $\beta=2$ surmise (Turiaci–Witten §2.5; $\beta=1$ for the $CT$-self-conjugate multiplet at odd fermion number). Secondary: level statistics of $H$ in symmetry-resolved sectors at the largest accessible $N$; LMRS statistics of projected simple operators in the BPS subspace; information/entanglement entropy of BPS states (Chang et al. §3.4); whether gaps close as $N\to\infty$ (Casimir-type fermionic matrix models keep $O(1)$ gaps, Klebanov et al. 2018 §6.1).

**Q3. What is the minimal structure that breaks the Casimir shortcut?** Which cubic supercharges (number of flavors, which trace structures) have $H$ *not* a function of Casimirs? Is $p=3$ really the "magic number"? Do $p=2$ models with mixed traces already concentrate?

**Q4. Large-$N$ organisation of the fortuitous sector.** Fortuitous states are maximal-Casimir, non-singlet (single matrix) or partly singlet (3-matrix, per project ED at $N=2$). Does large-$N$ factorisation of single-trace correlators hold in the relevant states? If not, what replaces it (covariant/Wigner–Eckart blocks, Marchesini–Onofri-type non-singlet equations)?

**Q5. Near-BPS spectrum and the super-Schwarzian.** Quantitative targets from Turiaci–Witten (3.11) and FGMS (5.39)–(5.42): non-BPS multiplet gaps $E_0(q)=q^2/(4\hat q^2)$ in Schwarzian units (growing quadratically with distance of the average charge $q$ from the window), edge density $\rho_q\propto\sinh(2\pi\sqrt{E-E_0(q)})/E$, and $1/\sqrt\varepsilon$ threshold only for the $q=0$ multiplet. Unknown for the matrix model: the overall Schwarzian scale as a function of $N$ (is it $\propto1/N^2$, the natural guess with $N^2$ fermions?). (Chen's conjecture; untested.)

**Q6. Methodological.** For a purely fermionic model the operator algebra truncates; what is the smallest operator level at which the bootstrap is tight in (a) BPS sectors, (b) lifted sectors, as a function of $N$? Which constraints are essential (supercharge EOM, ground-state positivity restricted to sector-preserving operators, $SU(N)$ Ward identities)?

## Status of what is known

*Re-verified 2026-09-15 in the project `.venv` unless stated; see `docs/todo.md` for the two corrections found.*

- Single matrix: $E_0(N_\Psi)$ exact for $N\le4$; BPS window $[N(N-1)/2,N(N+1)/2]$ confirmed; covariant single-trace bootstrap tight at $N=2$ in all sectors and at $N=3$ in BPS sectors; lifted sectors at $N=3$ not yet tight.
- 3-matrix, $N=2$: 972 BPS states at $N_\Psi\in\{5,6,7\}$ with counts $243:486:243$; 168 distinct energies; gauge-singlet BPS multiplets present.
- 2026-09-16: sector-resolved bootstrap with adjoint-valued operators (D3) is exact for the single-matrix model in all sectors at $N=3$ and in 15/17 sectors at $N=4$; the window-edge sectors converge slowly in level ($0\to11.4\to14.7$ of $18$ at levels $\le3,4,5$). Three-matrix $N=2$: $k=0,1$ exact (analytic anchors), $k=2,3$ at the SUSY floor at level 2; level 3 feasible only for $k=2$ within memory.
- 2026-09-17: symmetry-reduced bootstrap (D4: invariant $\rho_k$ over $SU(2)_{\rm gauge}\times\mathbb Z_3$ isotypic blocks, $\mathbb Z_3$-graded cones, Cho et al. adjoint-projected channel). Three-matrix $N=2$: $E_0(2)=5.16536$ **reproduced exactly at level 3** (the plain adjoint channel stalls at 88 %); $k=3$: $E_0(3)\ge0.9622$ at level 3 (78 %, 9 GB) and **exact ($1.22706$) at level 4** (18 GB, at the memory ceiling); $k=4$ (almost BPS, $E_0=0.088$) not lifted at level 3. Partial answer to Q6: at $N=2$ the essential ingredients in lifted sectors are the sector restriction, the adjoint-valued words, the EOM and the singlet/adjoint split of the Gram matrix; Ward identities are automatic and GS positivity is not needed. Q1a status and plan (2026-09-18): the Hilbert-space bootstrap stops at $N=2$ by design ($d_k$ up to $2\times10^7$ at $N=3$). The plan of record (`docs/trace_bootstrap_plan.md`) re-expresses the same, now-validated constraint set at the trace level with $N$ as a numerical coefficient (exact at finite $N$, no factorisation); realistic reach is the fixed-$k$ ($k\lesssim4$–$5$), any-$N$ regime, where $k=0,1$ are the D2 anchors and $k\ge2$ is new. The half-filling window — the sector that decides concentration, i.e. the first lifted sector next to the BPS window, with many particles and a tiny gap ($k=4$ at $N=2$: $E_0=0.088$, unlifted at level 3) — is where every bootstrap is weakest; at $N=3$ it is a sparse-Lanczos problem (plan §8), at $N\ge4$ neither method reaches it and index/cohomology counting remains the tool. Concentration is a $1/N$-subleading effect on $O(N^3)$ energies, so the planar ('t Hooft) bootstrap is blind to it.
- Corrections: the single-trace $H$ constant is $-\tfrac32N(N^2-1)$ (docs/derivations.md D1); the covariant gauge-invariant-only bootstrap gives trivial bounds in lifted sectors at $N=3$ (levels 3 and 4, with or without ground-state positivity) and must be rebuilt with adjoint-valued operators.
