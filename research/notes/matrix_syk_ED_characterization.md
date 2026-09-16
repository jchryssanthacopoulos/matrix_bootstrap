# Matrix SYK model — III. What $N=2$ exact diagonalization actually says

*Target model (Chen arXiv:2511.00790, §4): three $U(N)$ fermion matrices with the cubic supercharge $Q=\sum_{1\le i\le j\le k\le 3}\mathrm{Tr}[\Psi_i\Psi_j\Psi_k]$, $H=\{Q,\bar Q\}$. At $N=2$ this is $pN^2=12$ complex fermions, $\dim=2^{12}=4096$ — the largest ED-accessible case (the next, $N=3$, is $2^{27}$). Files: `msyk_model.py`, `msyk_ed.py`, `msyk_ed_run.txt`.*

The model is verified as a genuine $\mathcal N=2$ theory: $Q^2=0$, $H=H^\dagger$, $[H,N_\Psi]=0$, $[H,J^a]=0$ (gauge $SU(N)$), $[N_\Psi,Q]=3Q$.

## 1. R-charge concentration — confirmed and quantified

Chen states that a numerical simulation of this model "suggests R-charge concentration." Resolving the $E=0$ (BPS) states by fermion number $N_\Psi$ makes it sharp. Out of **972 BPS states** (23.7% of the Hilbert space), *every one* sits in $N_\Psi\in\{5,6,7\}$:

| $N_\Psi$ | 0–4 | 5 | 6 | 7 | 8–12 |
|---|---|---|---|---|---|
| #BPS | 0 | 243 | 486 | 243 | 0 |

- Centered exactly at $\langle N_\Psi\rangle = 6 = n/2$, with standard deviation $0.707$.
- BPS states occupy only **3 of the 13** available charge sectors.
- The counts are $243:486:243 = 3^5\times(1:2:1)$ — a clean, suggestive pattern (recall the single-matrix BPS multiplicity was $3^{N(N-1)/2}=3$ for $N=2$; here $3^5$ appears).

**Contrast with the solvable single-matrix model** ($N=2$): there the BPS states spread across their whole allowed window $\tfrac{N(N-1)}2\le N_\Psi\le\tfrac{N(N+1)}2$ with binomial weights $Z_{\mathrm{BPS}}\propto(1+q)^N q^{N(N-1)/2}$ (for $N=2$: counts $3,6,3$ over $N_\Psi=1,2,3$ — spread over the entire allowed range). The matrix SYK model instead **pins** its BPS states to the center — the defining feature of the black-hole-like ("meat", not "skeleton") fortuity that the single-matrix model lacks.

This is the concrete, quantitative version of Chen's remark: concentration is not approximate here, it is exact (width 3 out of 13, tightly peaked, symmetric $1{:}2{:}1$).

## 2. The spectrum is genuinely interacting

The full spectrum has only **168 distinct energies** among 4096 states, and the nonzero levels are **irrational** ($E=0.0880,\,0.2260,\,0.3572,\dots$), spanning $[0,104.6]$. This is qualitatively unlike the single-matrix model, whose Hamiltonian is the affine-Casimir $H=3N(N^2-1)-9\hat C_2$ with only the integer levels $\{0,18,45,72\}$. The matrix SYK Hamiltonian is **not** a function of a single Casimir — it is a genuinely interacting, non-solvable model. That non-algebraic spectrum is the prerequisite for chaos.

## 3. Chaos at $N=2$ is not diagnosable — and that is the point

Attempting random-matrix level statistics in the central $N_\Psi=6$ sector (dim 924) runs into a wall: that sector contains only **~50 distinct energies**. The huge degeneracy is symmetry — gauge $SO(3)$ $\times$ flavor $S_3$ $\times$ SUSY multiplets — and after resolving gauge $SU(2)$ the pooled level-spacing ratio is $\langle r\rangle\approx0.32$, below even the Poisson value $0.386$: a sign that too few independent levels remain for a meaningful Wigner–Dyson test, not evidence of integrability.

The honest conclusion: **$N=2$ is too small to see spectral chaos.** Genuine SYK-like chaos (and the super-Schwarzian low-energy behavior Chen anticipates) is a **large-$N$** statement — precisely the regime where ED is impossible and where the bootstrap is the only available tool.

## 4. Consequence for the project — the beyond-ED target

$N=2$ ED gives everything at $N=2$; it cannot say anything at large $N$. So the genuinely new, "can't-get-it-any-other-way" statement must come from a **large-$N$ bootstrap**. The cleanest such target, made concrete by §1, is:

> **Rigorous large-$N$ lower bounds on the ground-state energy in each R-charge sector, $E_0(N_\Psi)$.** A strictly positive lower bound in sectors away from the center *proves* the absence of BPS states there — i.e. proves R-charge concentration at large $N$, where no ED or simulation can reach.

This requires upgrading the bootstrap from the explicit-matrix (ED-scale) implementation validated in Note II to the **planar / single-trace** formulation: variables are large-$N$-factorized single-trace correlators of words in $\Psi_i,\bar\Psi_i$; constraints are moment positivity, the fermionic algebra with large-$N$ factorization, gauge (Gauss law), the $S_3$/$O(3)$ flavor decomposition (as in Lin–Zheng's $O(D)$ blocking), the SUSY relations, and ground-state positivity for two-sided control. The $N=2$ ED data above is the validation anchor.

### Sources (project library)
- Y. Chen, *Fortuity with a Single Matrix*, arXiv:2511.00790 (§4 proposal; R-charge concentration claim).
- C.-M. Chang, Y. Chen, B. S. Sia, Z. Yang, *Fortuity in SYK Models*, arXiv:2412.06902 (R-charge concentration).
- H. W. Lin, Z. Zheng, *High-Precision Bootstrap of Multimatrix Quantum Mechanics*, arXiv:2507.21007 (planar single-trace bootstrap, symmetry blocking).
- S. Laliberte, B. McPeak, *Bootstrapping supersymmetric (matrix) quantum mechanics*, arXiv:2510.01356 (SUSY constraints, ground-state positivity).

---

## Addendum (2026-09-16): flavor symmetry correction

§3 refers to a "flavor $S_3$" symmetry. This is incorrect: Chen's supercharge $\sum_{i\le j\le k}\mathrm{Tr}[\Psi_i\Psi_j\Psi_k]$ contains the cyclic class of $\mathrm{Tr}[\Psi_1\Psi_2\Psi_3]$ but not that of $\mathrm{Tr}[\Psi_1\Psi_3\Psi_2]$, so only **cyclic** flavor permutations ($\mathbb Z_3$) are symmetries. Verified at $N=2$: $\|[H,U_{\rm cyc}]\|=0$, $\|[H,U_{\rm swap}]\|_\infty=6$ for all three transpositions; the full one-body commutant of $H$ is exactly gauge $su(2)\oplus u(1)_R$ (`docs/derivations.md` D2.9). Symmetry-sector bookkeeping (e.g. for level statistics) should use $SU(N)\times U(1)_R\times\mathbb Z_3$.
