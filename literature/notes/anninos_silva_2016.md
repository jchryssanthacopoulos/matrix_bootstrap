# Anninos, Silva (2016) — Solvable Quantum Grassmann Matrices

**File:** `papers/anninos_silva_2016.pdf`

## Full citation

Dionysios Anninos, Guillermo A. Silva, *Solvable Quantum Grassmann Matrices*, J. Stat. Mech. 1704 (2017) 043102, arXiv:1612.03795v1 [hep-th] (12 Dec 2016). Cited by Chen 2025 as ref. [51].

## Main question

Solve exactly (at finite $N$ and at large $N$) quantum mechanical models of many non-locally interacting fermions with vector-like and matrix-like index structure and quartic interactions, by mapping the thermal partition function to an ordinary bosonic (matrix) integral; determine the thermal phase structure and real-time correlators.

## Physical system

- Vector model: $N$ complex fermions $\psi^I$, $Z[\beta]=\int D\psi D\bar\psi\,e^{-\oint[\bar\psi\dot\psi-\frac{1}{4N\gamma}(\bar\psi\psi)^2]}$ (2.1); $\hat H=-\frac1{4N\gamma}(\bar\psi\psi)^2+\frac1{4\gamma}\bar\psi\psi-\frac N{16\gamma}$ (2.15).
- Matrix model: $NL$ complex fermions $\psi^{iA}$, $U(N)\times U(L)$ bifundamental, $\alpha\equiv L/N$; $H=-(4L\gamma)^{-1}\sum\bar\psi^{Ai}\psi^{iB}\bar\psi^{Bj}\psi^{jA}$ + quadratic normal-ordering terms (4.1); Hilbert space $2^{NL}$.

## Important definitions

- Hubbard–Stratonovich field $\lambda(\tau)$ (2.2) / Hermitian matrix $M_{ij}(\tau)$ (4.2), which behaves as an emergent $(0+1)$-d gauge field: reparametrisation invariance (2.4), non-compact $U(1)$ gauge symmetry (3.25), $U(N)$ gauge symmetry (4.4); determinant depends only on the Polyakov loop $W=\mathrm{tr}\,Pe^{i\oint M}$ (4.5).
- Fermion two-point functions as Wilson lines (3.9)–(3.10), (5.6)–(5.7).

## Main assumptions

- Quartic interactions; $\gamma>0$ (vector) and $\tilde\gamma=-\gamma>0$ (matrix, obtained by $\mu_i\to i\mu_i$); at finite $N$ expressions are analytic in $\gamma$.
- Large-$N$ saddle analysis with $\alpha$ fixed.

## Main analytical results

1. Vector model exact: $Z[\beta]=\sum_{n=0}^N\binom Nn e^{\beta(N-2n)^2/(16N\gamma)}$ (2.13), $E_n=-\frac{(N-2n)^2}{16N\gamma}$, $d_n=\binom Nn$ (2.14). Doubly-degenerate ground state $E_g=-N/(16\gamma)$; $2^N$ states peaked at $n=N/2$, $E=0$. Large-$N$ thermal transition at $\beta\approx8\gamma$ from $O(1)$ to $O(N)$ entropy (Fig. 1); saddle equation $\tanh\frac{\beta\lambda_0}2=4\gamma\lambda_0$ (3.14).
2. Real-time correlator (3.2); low $T$: single oscillation with $\Delta E=1/4\gamma$ (3.16); high $T$: Gaussian decay $G\sim e^{-t^2/32N\gamma(\gamma-\beta/8)}$ (3.20) with recurrences at $t\sim4\pi N\gamma$ — "a flavor of integrability"; at $\beta=8\gamma$ decay after $t\sim N^{1/4}\beta$ (3.23). Low-energy effective theory: a hydrodynamic $U(1)$ phase $\phi(\tau)$ with $S=N\gamma\int(\partial_\tau\phi)^2$ (3.27).
3. Matrix model exact partition function as a unitary-type eigenvalue integral
   $$Z=\mathcal N\int[DU]e^{L\gamma\,\mathrm{tr}\oint(U\dot U^\dagger)^2}\int\prod_id\mu_i\prod_{i<j}\sin^2\frac{\mu_i-\mu_j}2\prod_i\cos^L\frac{\mu_i}2e^{-L\gamma\mu_i^2}\quad(4.10)$$
   or, after $\mu\to i\mu$, with $\sinh^2$ and $\cosh^L$ (4.12)–(4.13) — the modified Vandermonde characteristic of finite-temperature gauge theory / Chern–Simons. Checks: $L=1,N=2$: $Z=e^{5\beta/8\tilde\gamma}+3e^{\beta/8\tilde\gamma}$ (4.14); $L=N=2$ (4.15); $L=2,N=3$ (4.16); state counts $2^{NL}$ verified.
4. Phase structure (§5.1): eigenvalue potential $V(\mu)=L\sum_i[\tilde\gamma\mu_i^2-\log\cosh\frac{\mu_i}2]$ (5.1); high $T$: single-cut distribution (5.3) of Chern–Simons type; low $T$: double well at $\mu=\pm1/4\tilde\gamma$, connected for small $\alpha$, disconnected for large $\alpha$ ($\alpha>2/(1-8\tilde\gamma\log2)$).
5. Matrix correlators: $G(\delta\tau)=\int dy\,\rho(y)\frac{e^{-\delta\tau y}}{1+e^{-y}}$ (5.10); high $T$: $G\approx I_1(2\sqrt t\,\delta\tau)/(2\sqrt t\,\delta\tau)$ (5.11), decaying as $\delta t^{-3/2}$ — much faster than the vector model.

## Important equations

(2.1), (2.13)–(2.15), (3.2), (3.20), (4.1), (4.10)–(4.16), (5.1), (5.3), (5.10)–(5.11).

## Numerical methods

Exact finite-$N$ sums/plots; numerical evaluation of eigenvalue-density integrals; comparison of exact vs large-$N$ correlators (Figs. 2–4).

## Relevant figures/results

Fig. 1 (energy/entropy vs $\beta$, transition at $\beta=8\gamma$), Fig. 2 (recurrences), Fig. 3 (critical-temperature correlator), Fig. 4 (matrix correlator decay).

## Limitations

- Detailed large-$N$ phase diagram of the matrix model deferred; correlators only at high $T$.
- Non-supersymmetric; Hamiltonian is a (deformed) Casimir-type quartic — integrable flavour; no chaos.

## Relationship to our project

Shows how a purely fermionic matrix model's thermodynamics reduces *exactly* to a unitary/eigenvalue integral — the same mathematical object Chen uses to count $r_*$ multiplicities (Chen's (3.2)/(3.9)). It is a useful reference for (a) computing exact finite-$N$ partition functions / state counts of our fermionic models, (b) what "non-chaotic" looks like in real-time correlators (oscillations, early recurrences), which provides a baseline against which SYK-like decay in the 3-matrix model would be judged, and (c) the fact that fermion number sectors dominate thermodynamics near $N_\Psi\approx$ half filling — the same region where the fortuitous states of the matrix SYK model are expected to concentrate.
