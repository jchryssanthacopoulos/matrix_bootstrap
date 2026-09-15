# Research questions

*Last updated 2026-09-15 (after the literature review; see `literature/synthesis.md`). Keep this file current.*

## Primary question

**Q1. Does Chen's 3-matrix "matrix SYK" model exhibit R-charge concentration for $N\ge3$?**

Model: $p=3$ $U(N)$-adjoint complex fermion matrices, $Q=\sum_{1\le i\le j\le k\le3}\mathrm{Tr}[\Psi_i\Psi_j\Psi_k]$, $H=\{Q,\bar Q\}$ (Chen 2025, eq. 4.3). Chen reports only that an $N=2$ simulation "suggests" concentration.

Sharp formulation (Chang–Chen–Sia–Yang 2024, Conjecture 1 and §2): in every irreducible cochain complex — labelled by $N_\Psi$ mod 3, the $SU(N)$ irrep, and flavor ($S_3$) quantum numbers — with macroscopic index, BPS states occupy a *single* degree. Observable proxy: the sector-resolved ground energy $E_0(N_\Psi,\mathbf r)$; concentration ⇔ $E_0>0$ outside a window of at most three consecutive $N_\Psi$ (one per $\mathbb Z_3$ class) per irrep. Baseline: the single-matrix model has $E_0=0$ on $N+1$ consecutive sectors, all in $r_*$ (Chen (2.22)) — non-concentrated for $N\ge3$.

Sub-questions:
- Q1a. Can the bootstrap give rigorous lower bounds $E_0(N_\Psi)>0$ in lifted sectors at $N$ beyond exact diagonalisation ($N=3$: $2^{27}$ states)?
- Q1b. Are refined indices $\mathrm{Tr}[(-1)^Fq^{N_\Psi}]$ per $SU(N)$ irrep saturated by BPS counts (index saturation ⇔ concentration, Chang et al. §5.1)?
- Q1c. How does the number of BPS states scale with $N$ ($e^{cN^2}$ multiplicity, as for a genuine black-hole-like sector, versus a large irrep with small multiplicity as in the single-matrix model)?

## Secondary questions

**Q2. Is the 3-matrix model chaotic?** Diagnostics: level statistics in symmetry-resolved sectors at the largest accessible $N$; LMRS statistics of projected simple operators in the BPS subspace; information/entanglement entropy of BPS states (Chang et al. §3.4); whether gaps close as $N\to\infty$ (Casimir-type fermionic matrix models keep $O(1)$ gaps, Klebanov et al. 2018 §6.1).

**Q3. What is the minimal structure that breaks the Casimir shortcut?** Which cubic supercharges (number of flavors, which trace structures) have $H$ *not* a function of Casimirs? Is $p=3$ really the "magic number"? Do $p=2$ models with mixed traces already concentrate?

**Q4. Large-$N$ organisation of the fortuitous sector.** Fortuitous states are maximal-Casimir, non-singlet (single matrix) or partly singlet (3-matrix, per project ED at $N=2$). Does large-$N$ factorisation of single-trace correlators hold in the relevant states? If not, what replaces it (covariant/Wigner–Eckart blocks, Marchesini–Onofri-type non-singlet equations)?

**Q5. Near-BPS spectrum and the super-Schwarzian.** In a concentrating sector, does the gap above the BPS states scale as a power of $1/N$ and does the density of states show a $\sinh\sqrt{E}$-type edge, as the $\mathcal N=2$ super-Schwarzian would predict? (Chen's conjecture; untested.)

**Q6. Methodological.** For a purely fermionic model the operator algebra truncates; what is the smallest operator level at which the bootstrap is tight in (a) BPS sectors, (b) lifted sectors, as a function of $N$? Which constraints are essential (supercharge EOM, ground-state positivity restricted to sector-preserving operators, $SU(N)$ Ward identities)?

## Status of what is known (from `research/notes/`, not re-verified this session)

- Single matrix: $E_0(N_\Psi)$ exact for $N\le4$; BPS window $[N(N-1)/2,N(N+1)/2]$ confirmed; covariant single-trace bootstrap tight at $N=2$ in all sectors and at $N=3$ in BPS sectors; lifted sectors at $N=3$ not yet tight.
- 3-matrix, $N=2$: 972 BPS states at $N_\Psi\in\{5,6,7\}$ with counts $243:486:243$; 168 distinct energies; gauge-singlet BPS multiplets present.
- No Python environment with numpy/scipy is currently installed in this checkout; results above cannot be re-run until one is set up.
