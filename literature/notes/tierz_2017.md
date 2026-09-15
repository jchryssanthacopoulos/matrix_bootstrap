# Tierz (2017) — Polynomial solution of quantum Grassmann matrices

**File:** `papers/tierz_2017.pdf`

## Full citation

Miguel Tierz, *Polynomial solution of quantum Grassmann matrices*, J. Stat. Mech. 1705 (2017) 053203, arXiv:1703.02454v2 [hep-th] (8 May 2017). Cited by Chen 2025 as ref. [52].

## Main question

Compute exactly, for arbitrary $N$ and given $L$, the thermal partition function of the Anninos–Silva fermionic matrix model, using $q$-deformed orthogonal polynomials, and read off the full spectrum and degeneracies.

## Physical system

The Anninos–Silva matrix model: $NL$ complex fermions $\psi^{iA}$, $U(N)\times U(L)$, quartic $H=-(4L\gamma)^{-1}\sum\bar\psi^{Ai}\psi^{iB}\bar\psi^{Bj}\psi^{jA}$ (+ quadratic normal-ordering terms) (2.1); $\dim\mathcal H=2^{NL}$. Starting point is the Anninos–Silva eigenvalue integral (2.2)–(2.3).

## Important definitions

- Rewriting (3.1): $\cosh^L(\mu_i/2)$ is a characteristic-polynomial insertion $\prod_i(\lambda+x_i)^L$ in a Stieltjes–Wigert ensemble with weight $e^{-L\tilde\gamma\log^2x}$; spectral parameter $\lambda=-\exp(1/4\tilde\gamma+N/2L\tilde\gamma)=-q^{-L/2-N}$, $q=\exp(-1/(2\tilde\gamma L))$.
- Stieltjes–Wigert polynomials (3.4), monic form (3.6); $q$-Pochhammer, $q$-binomials.

## Main assumptions

- Follows Anninos–Silva in working with the hyperbolic ($\sinh$) version after $\mu\to i\mu$; contour-rotation subtleties argued to be absent for $L>0$.

## Main analytical results

- General $L\times N$: $Z$ is an $L\times L$ Wronskian of Stieltjes–Wigert polynomials $S_N,\dots,S_{N+L-1}$ and derivatives (3.3).
- $L=1$: $\hat Z_{1\times N}=q^{-N^2-N/2}\sum_{k=0}^N\binom Nk_q q^{k^2-kN}$ (3.8) — a single SW polynomial evaluated at $\lambda$; symmetric reduction for odd/even $N$.
- For fixed $L$ the number of distinct energies is quadratic in $N$ ⇒ exponential degeneracies; mid-spectrum levels are exponentially degenerate (unimodality of $q$-binomials), ground/first excited states have polynomial degeneracy, top levels nearly non-degenerate (Hagedorn-like structure in the sense of dual resonance models).
- High-temperature limit → Gaussian matrix model; SW → Hermite polynomials; duality swapping $N\leftrightarrow L$ with Wick rotation of the spectral parameter; generic $L\times N$ in terms of a generalised Hermite polynomial (§6); large-$N$ at fixed $q$ related to Rogers–Ramanujan identities (Appendix).

## Important equations

(2.1)–(2.3), (3.1)–(3.4), (3.6), (3.8).

## Numerical methods

Symbolic evaluation of determinants/polynomials in Mathematica (accompanying notebook); e.g. $L=5,N=10$ ($2^{50}$ states) is trivial.

## Relevant figures/results

Explicit spectra for small $L,N$ (in the text and notebook); confirmation of Anninos–Silva's $1\times2$, $2\times2$, $2\times3$ cases.

## Limitations

- Specific to the quartic Anninos–Silva model, whose solvability rests on the bosonic-integral representation; no supersymmetry; no dynamics/correlators.

## Relationship to our project

Illustrates that thermal partition functions of fermionic matrix models can be exact and compact even when the Hilbert space is astronomically large — a possible strategy for *index/partition-function-level* data on our models (e.g. refined Witten indices $\mathrm{Tr}(-1)^Fq^{N_\Psi}$ per $SU(N)$ irrep as unitary-matrix integrals with character insertions, à la Chen §3), which provide rigorous BPS-count bounds to compare with the bootstrap. Also a reminder that exponential degeneracies (here a symptom of integrability) must be distinguished from the exponential BPS degeneracy that fortuity produces.
