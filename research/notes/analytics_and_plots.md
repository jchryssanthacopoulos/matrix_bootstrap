# Analytics and plots

Two questions: what can be done analytically / enumeratively, and what plots are informative. Results below are verified against ED in-session.

## Part 1 — analytic and quasi-analytic results

### Single matrix (exactly solvable) — closed forms

**1. The BPS distribution is exactly binomial.** From Chen's $Z_{\rm BPS}(q)=3^{N(N-1)/2}(1+q)^N q^{N(N-1)/2}$, the number of BPS states at fixed $R$-charge is
$$
\#\mathrm{BPS}(N_\Psi)=3^{\frac{N(N-1)}2}\binom{N}{\,N_\Psi-\tfrac{N(N-1)}2\,},\qquad
\tfrac{N(N-1)}2\le N_\Psi\le \tfrac{N(N+1)}2 .
$$
Verified exactly ($N=2$: $3,6,3$; $N=3$: $27,81,81,27$). This makes "non-concentration" precise: the distribution is a **binomial of width $N$** (support $N+1$ sectors, standard deviation $\sqrt N/2$). Concentration in the SYK sense means replacing $\binom{N}{\cdot}$ by a binomial of $O(1)$ width.

**2. Edge energy law.** For small charge the maximal Casimir grows linearly, $C_2^{\max}(N_\Psi)=N\,N_\Psi$ (fill $N_\Psi$ "longest" roots, each contributing $C_2(\mathrm{adj})=N$), giving the closed form
$$
E_0(N_\Psi)=3N(N^2-1)-9N\,N_\Psi = 3N\,(N^2-1-3N_\Psi),\qquad 0\le N_\Psi\le N
$$
(and mirrored near the top). Verified: $N=4$ gives $144,108,72,36$ at $N_\Psi=1,2,3,4$; $N=3$ gives $45,18$. So near the edge $E_0$ falls linearly with slope $9N$ per unit charge.

**3. The whole profile is representation theory.** Because $E_0(N_\Psi)=3N(N^2-1)-9\,C_2^{\max}(N_\Psi)$, the entire $E_0(N_\Psi)$ curve at **any $N$** is a Casimir-maximization over Young diagrams in $\wedge^{N_\Psi}\mathbf{adj}$ — an enumerative computation, no ED needed. This is an analytic route to the single-matrix baseline at $N=5,6,\dots$ (a cross-check for the ED scripts) and to the large-$N$ intensive curve $E_0(\nu)/N^3$.

### Three-matrix SYK (the target) — enumerative structure

**4. A strategically important fact: the SYK model has gauge-singlet BPS states.** Decomposing the $N=2$ BPS states by gauge $SU(2)$ Casimir:

| $N_\Psi$ | singlet ($j{=}0$) | triplet ($j{=}1$) | $j{=}2$ | $j{=}3$ | total |
|---|---|---|---|---|---|
| 5 | 9 | 27 | 18 | 9 | 243 |
| 6 | 18 | 54 | 36 | 18 | 486 |
| 7 | 9 | 27 | 18 | 9 | 243 |

(counts are *multiplicities* of each gauge irrep). Unlike the single matrix — whose fortuitous states are a **single** non-singlet rep $\mathbf r_*$ — the three-matrix model has genuine **gauge-invariant** ($j=0$) fortuitous states (36 of them). **Consequence:** the standard singlet-sector matrix bootstrap (Lin–Zheng / Laliberte–McPeak, which assumes a near-singlet master field) *applies directly to the three-matrix target*. The master-field obstruction that forced the covariant formulation for the single matrix is **much milder for the actual model of interest.** This meaningfully de-risks the SDPB plan: bootstrap the singlet-sector $E_0(N_\Psi)$ of the three-matrix model with off-the-shelf technology.

**5. Edge energies are exact at any $N$ (small sectors).** The low-charge sectors are tiny, so $E_0$ near the edges is computable even where full ED is impossible. $E_0(N_\Psi=0)=\lVert Q|0\rangle\rVert^2$ exactly; new numbers at $N=3$: $E_0(0)=387$, $E_0(1)=273$ (and mirror at $N_\Psi=27,26$). These pin the *ends* of the $N=3$ profile analytically; the scripts fill in the center.

**6. Structure of the $243:486:243=3^5(1{:}2{:}1)$ pattern.** The $1{:}2{:}1$ in charge is the same binomial-in-$N$ ($=\binom{2}{\cdot}$) as the single matrix; the model-dependence is entirely in the multiplicity ($3^1\to3^5$) and the recentering ($N_\Psi\!\sim\!1\to\!\sim\!6$). Table (4) exhibits the gauge content behind the $3^5$.

### Further analytic handles (not yet done)
- **Refined Witten index** $\mathrm{Tr}[(-1)^F q^{N_\Psi}]$ / $Q$-cohomology Euler characteristics: give rigorous *lower bounds* on BPS counts per charge for the three-matrix model, purely combinatorially.
- **Large-$N$ intensive curve** $E_0(\nu)/N^3$ for the single matrix from $C_2^{\max}$, as the analytic backdrop against which the bootstrap's SYK curve is read.

## Part 2 — plots

Four example figures are generated in `matrix_syk_figures.pdf/.png` (`make_plots.py`), from verified data:

- **(A) $R$-charge-resolved ground energy** $E_0(N_\Psi)$, single matrix $N=2,3,4$ vs SYK $N=2$: the flat-zero stretch is the BPS window. The headline observable.
- **(B) BPS window width vs $N$**: single matrix is the exact line $N+1$; the concentration claim is that the SYK points stay flat as the line climbs. This is the whole program in one panel.
- **(C) BPS distribution**: single matrix binomials visibly widen with $N$; SYK is pinned. Concentration made visual.
- **(D) Spectra**: single matrix (4 integer levels, hugely degenerate) vs SYK (168 irrational levels) — solvable vs interacting.

### The "islands" plot (the analog of the paper figures you mentioned)
The bootstrap papers plot the **allowed region** in a plane of two observables (e.g. $E$ vs $\Tr\langle X_1X_2\rangle$), the physical theory living inside the boundary. Our direct analog, once SDPB is running, is a **2-D allowed region at fixed $R$-charge**: e.g.
$$
\big(E,\ \langle \Tr[\Psi\bar\Psi\Psi\bar\Psi]\rangle\big)\quad\text{or}\quad
\big(\langle \hat C_2\rangle,\ \langle \Tr[\Psi^2\bar\Psi^2]\rangle\big)\ \text{at fixed }N_\Psi,
$$
swept over operator level to show the region shrinking onto the true point (which we know exactly for the single matrix — a perfect calibration of the island against ground truth). For the three-matrix singlet sector these islands are genuine predictions.

### Other plots worth having
- **Bootstrap convergence**: lower bound on $E_0(N_\Psi)$ vs operator level $L$, one curve per charge, approaching the exact value — the standard "it converges" figure, and our test of whether the lifted sectors tighten.
- **$C_2^{\max}(N_\Psi)$ profile** for the single matrix: shows the plateau at $N(N^2-1)/3$ = the BPS window, with the linear edges of law (2).
- **Level-spacing ratio / spectral form factor** for the SYK model (chaos), once a large-enough symmetry-resolved sector is reachable — the chaos analog of the concentration plots.
