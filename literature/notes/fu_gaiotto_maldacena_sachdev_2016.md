# Fu, Gaiotto, Maldacena, Sachdev (2016) — Supersymmetric SYK models

**File:** `papers/fu_gaiotto_maldacena_sachdev_2016.pdf`

## Full citation

Wenbo Fu, Davide Gaiotto, Juan Maldacena, Subir Sachdev, *Supersymmetric Sachdev–Ye–Kitaev models*, Phys. Rev. D 95 (2017) 026009, arXiv:1610.08917v3 [hep-th] (4 Mar 2017).

## Main question

Define supersymmetric ($\mathcal N=1,2$) generalisations of SYK with $H=Q^2$ or $\{Q,\bar Q\}$ for a random polynomial supercharge; solve them at large $N$; determine whether SUSY is broken; derive the super-Schwarzian low-energy theory; count ground states via a refined Witten index.

## Physical system

- $\mathcal N=1$: $N$ Majoranas, $Q=i\sum_{i<j<k}C_{ijk}\psi^i\psi^j\psi^k$ (1.1), $\overline{C^2}=2J/N^2$ (1.3), $H=Q^2=E_0+\sum J_{ijkl}\psi^4$ (1.4)–(1.5); generalisation $Q\sim\psi^{\hat q}$.
- $\mathcal N=2$ (§V): $N$ complex fermions, $Q=i\sum_{i<j<k}C_{ijk}\psi^i\psi^j\psi^k$, $\bar Q=i\sum\bar C^{ijk}\bar\psi_i\bar\psi_j\bar\psi_k$ (5.1), $Q^2=\bar Q^2=0$ automatically; $H=\{Q,\bar Q\}=|C|^2+\sum J_{ijkl}\psi\psi\bar\psi\bar\psi$ (5.3); $U(1)_R$ with $\psi$ of charge $1/\hat q$ ($=1/3$), $Q$ of charge 1; $\mathbb Z_{\hat q}$ subgroup commutes with $Q$.

## Important definitions

- Auxiliary boson $b_i$ linearising SUSY (2.2); $G$–$\Sigma$ superspace formulation; IR dimensions $\Delta_\psi=1/(2\hat q)$, $\Delta_b=\Delta_\psi+\tfrac12$.
- **Refined Witten index** $W_r=\mathrm{Tr}[(-1)^Fe^{2\pi irQ_R}]=(1-e^{2\pi ir/\hat q})^N$ (5.5), maximal at $r=(\hat q\pm1)/2$: $\log|W|=N\log[2\cos\frac{\pi}{2\hat q}]$ (5.6) = the large-$N$ ground-state entropy ⇒ index saturated, SUSY unbroken even non-perturbatively.
- **$\mathcal N=2$ super-Schwarzian** (5.27)–(5.29): $\mathrm{SU}(1,1|1)$ symmetry; low-energy thermodynamics $E/N=\alpha_s(2\pi^2/\beta^2J+\mu^2/J)$, $Q_R/N=4\alpha_s\mu/J$, $S-S_0=\pi\sqrt{8\alpha_sE/JN-(Q_R/N)^2}$ (5.39)–(5.40); one-loop $Z\sim e^{S_0}e^{\frac{N\alpha_s}{J}(2\pi^2/\beta+2\beta\mu^2)}$ with no $\beta$-dependent prefactor (equal boson/fermion zero modes) (5.41); density of states $D(E,Q_R)\propto e^{S_0+\Delta S}/(\Delta S)^2$ (5.42); sum over windings of the $U(1)$ phase (5.43).

## Main assumptions

Gaussian random couplings, large $N$, melonic dominance; conformal ansatz in the IR.

## Main analytical results

- $\mathcal N=1$: unbroken SUSY at large $N$ but broken non-perturbatively, $E_0\sim e^{-\alpha N}$ (ED §III); super-Schwarzian with bosonic $h=2$ and fermionic $h=3/2$ modes; ladder kernel spectrum with an extra $(h=1,\tfrac32)$ pair.
- $\mathcal N=2$: exact zero-energy ground states with $e^{N\log 2\cos(\pi/2\hat q)}$ degeneracy; **ED counts for $\hat q=3$** (5.7): $D(N,0)=2\cdot3^{N/2-1}$, $D(N,\pm\tfrac13)=3^{N/2-1}$ for even $N$ (zero elsewhere) — i.e. BPS states in exactly three adjacent $U(1)_R$ charges in ratio **1:2:1**; for odd $N$, $D(N,\pm\tfrac16)=3^{(N-1)/2}$ ($N\equiv3$ mod 4) or $D(N,\pm\tfrac16)=3^{(N-1)/2}$ and $D(N,\pm\tfrac12)=1$ or $3$ ($N\equiv1$ mod 4). Footnote: ground states ↔ $Q$-cohomology.

## Important equations

(1.1)–(1.5), (5.1)–(5.7), (5.27)–(5.29), (5.39)–(5.43).

## Numerical methods

Exact diagonalisation of $\mathcal N=1$ ($N\le$ moderate) and $\mathcal N=2$ models; numerical solution of the large-$N$ Schwinger–Dyson equations; ladder-kernel eigenvalues.

## Relevant figures/results

ED ground-state energies vs $N$ ($\mathcal N=1$); the degeneracy table (5.7).

## Limitations

Disordered couplings; only $\hat q=3$ counts reported; no statement about non-generic (structured) supercharges.

## Relationship to our project

The parent model: Chen's $Q=\mathrm{Tr}\Psi^3$ and the 3-matrix supercharge are *non-random, sparse* special cases of (5.1) with $N\to pN^2$ fermions. Three direct uses: (1) the refined index (5.5) is exactly computable for the matrix models in each $SU(N)$ irrep and gives the rigorous BPS lower bounds we want; (2) the ED pattern (5.7) — three adjacent charges with counts $1{:}2{:}1$ — is what a concentrating, super-Schwarzian-governed model produces, and it is precisely the $243{:}486{:}243$ pattern found for the 3-matrix model at $N=2$ (project ED, reproduced 2026-09-15), while the single-matrix model gives binomial $\binom N{\cdot}$ weights (Chen (2.20)) that agree with $1{:}2{:}1$ only at $N=2$; (3) the $\mathcal N=2$ super-Schwarzian formulas (5.39)–(5.42) are the concrete low-energy predictions (gap scaling with charge away from the BPS window, $\sqrt{E}$ edge) that a large-$N$ bootstrap of the 3-matrix model could test.
