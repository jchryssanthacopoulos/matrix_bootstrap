# Gaitan, Klebanov, Pakrouski, Pallegar, Popov (2020) — Hagedorn Temperature in Large N Majorana Quantum Mechanics

**File:** `papers/gaitan_klebanov_pakrouski_pallegar_popov_2020.pdf`

## Full citation

G. Gaitan, I. R. Klebanov, K. Pakrouski, P. N. Pallegar, F. K. Popov, *Hagedorn Temperature in Large N Majorana Quantum Mechanics*, Phys. Rev. D 101 (2020) 126002, arXiv:2002.02066v3 [hep-th] (8 Apr 2020). Cited by Chen 2025 as ref. [54].

## Main question

What are the large-$N$ densities of states and thermodynamics of Majorana fermion models with orthogonal symmetry groups, of vector type (one large rank) and of matrix type (two large ranks)? Do they show Hagedorn behaviour?

## Physical system

$O(N_1)\times O(N_2)\times O(N_3)$ Majorana model (1.1) specialised to: vector-type $O(N)\times O(2)^2$ ($N_2=N_3=2$) and $O(N)\times SO(4)$ with $H=\frac g4\epsilon_{IJKL}\psi_{aI}\psi_{aJ}\psi_{a'K}\psi_{a'L}$ (1.2); matrix-type $O(N_1)\times O(N_2)\times O(2)$ (4.1) = complex fermionic matrix model. Large $N$ with $gN=\lambda$ fixed.

## Important definitions

- Hilbert space decomposition $\mathcal H=\bigoplus_{\mu\subset\mu_{\max}}[\mu]_{O(N)}\otimes[(\mu_{\max}/\mu)^T]_{O(4)}$ (2.1) (Cauchy identities, App. A).
- $SO(4)\cong SU(2)_+\times SU(2)_-$ generators $K^\pm_i$ (2.4); energies $E(Q_+,Q_-,q_+,q_-)=g[Q_+(Q_++2)-Q_-(Q_-+2)+2q_-^2-2q_+^2]$ (2.7); $SO(N)$ irrep dimension (2.8).
- Hagedorn temperature $T_H=\lambda$ where $\log\rho\approx-|E|/\lambda$.

## Main assumptions

- Large $N$; for the matrix case, planar dominance of $\mathrm{tr}H^n$ by disconnected pieces.

## Main analytical results

1. Vector models: exact spectrum via Casimirs (2.5)–(2.7); ground state $E_0=-gN(N+2)$; density of states near $E_0$ $\rho\sim(8gE)^{-1/2}$ (2.10); for $|E|/\lambda=O(1)$, $\log\rho\approx-|E|/\lambda$ + slowly varying ⇒ specific heat diverges as $(T_H-T)^{-2}$ in the strict large-$N$ limit; finite $N$ smooths the transition (arbitrarily high $T$ reachable).
2. **Fermionic matrix model** $O(N)^2\times O(2)$ (§4): Hamiltonian in terms of charges $H=-2g\big(4C_2^{SU(N_1)}-C_2^{SO(N_1)}+C_2^{SO(N_2)}+\frac2{N_1}Q^2+(N_2-N_1)Q-\frac14N_1N_2(N_1+N_2)\big)$ (4.2); density of states near $E\approx0$ is Gaussian,
   $$\log\rho(E)=N^2\log2-\tfrac12\Big(\frac E{\lambda N}\Big)^2\quad(4.3),\qquad\log\rho=N^2\log2-\tfrac12x^2-\tfrac1{12N^2}x^4+\dots,\ E=gN^2x\quad(4.12)$$
   derived from moments $\mathrm{tr}H^n/\mathrm{tr}1$ computed by planar Feynman diagrams with propagator $\langle\psi_{ab}\psi_{a'b'}\rangle=\frac12\delta\delta$ and $H$ as a single-trace vertex: $\sigma_E^2=g^2(N^4-N^3)\approx(\lambda N)^2$ (4.6)–(4.7); connected $\mathrm{tr}H^4$ (4.11). Spectrum splits into even/odd branches (Fig. 5–6). Lanczos for $N_1=N_2=6$ ($7\times10^{10}$ states).

## Important equations

(1.1)–(1.2), (2.1), (2.5)–(2.10), (4.1)–(4.3), (4.6)–(4.7), (4.11)–(4.12).

## Numerical methods

Exact spectra from group theory for $N\le10$; Lanczos for the $6\times6$ matrix model; Gaussian/quartic fits to $\log\rho$; specific-heat curves.

## Relevant figures/results

Fig. 1 (cusp and linear $\log\rho$ in the vector model), Fig. 5–7 (Gaussian $\log\rho$ for $N=8,9,10$ matrix model), specific-heat figure.

## Limitations

- Non-supersymmetric, Casimir-solvable models; no chaos; the matrix model's smooth Gaussian density is generic for any model with a single-trace $H$ and finite Hilbert space, so it does not discriminate integrable from chaotic.

## Relationship to our project

The moment method of §4 — $\langle E^n\rangle=\mathrm{tr}H^n/\mathrm{tr}1$ via planar diagrams with a free-fermion propagator — applies verbatim to the matrix SYK models: for $H=\{Q,\bar Q\}$ with $Q=\sum\mathrm{Tr}\Psi_i\Psi_j\Psi_k$ the infinite-temperature moments are computable at large $N$ and give the bulk density of states (the $N$-scaling of the mean and width of $\rho(E)$ for $H=\{Q,\bar Q\}$ remains to be derived). This is a cheap analytic cross-check for any bootstrap or ED result on the *bulk* spectrum, and $\mathrm{tr}[H]$, $\mathrm{tr}[H^2]$ restricted to charge sectors are natural normalisation constraints for the bootstrap. The paper also shows how the density of states in fermionic matrix models is organised by Casimirs — useful when decomposing our spectra into $(N_\Psi,\text{irrep})$ blocks and when estimating whether a given sector can host a Schwarzian-like $\rho(E)\sim\sinh\sqrt E$ edge.
