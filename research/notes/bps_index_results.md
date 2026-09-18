# Refined Witten index of the three-matrix model — results (2026-09-19)

*Code: `src/bps_index.py`, `scripts/compute_bps_index.py`. Data: `results/data/refined_index_N{2,3,4}_p3.json` (map "c|w|lambda" → $I$). Derivation and conventions: `docs/derivations.md` D6. Everything here is exact (integer arithmetic); nothing depends on the dynamics beyond the symmetries of $Q$.*

## 1. What is computed

For the complex $(c,\lambda,\omega)$ — $c=k\bmod3$ (the only part of $N_\Psi$ commuting with $Q$), $\lambda$ a $U(N)$ irrep with $\sum\lambda_i=0$, $\omega$ the $\mathbb Z_3$ flavor charge —
$$I_{c,\lambda,\omega}=\sum_{k\equiv c\,(3)}(-1)^k\,n(k,\lambda,\omega),$$
with $n(k,\lambda,\omega)$ the multiplicity of $\lambda$ in the $(k,\omega)$ subspace of $\wedge^\bullet(\mathbb C^3\otimes\mathbf{adj})$. $|I|$ is a lower bound on the number of BPS multiplets in the complex, with equality iff its cohomology sits in one degree (index saturation, CCSY §5.1) — which is what concentration asserts. $\lambda=(j,-j)$ is spin $j$ at $N=2$.

## 2. Gate: $N=2$ against exact diagonalisation

All 36 complexes with BPS states are **saturated**: index = BPS multiplet count, complex by complex (`compute_bps_index.py --N 2 --check_ed`). The counts are flavor-independent and, per flavor charge and spin $j=0,1,2,3$: $3,9,6,3$ at $k=5$ and $k=7$; $6,18,12,6$ at $k=6$ — totals $243,486,243$. Sector dimensions per $(k,j,\omega)$ also match the isotypic-block dimensions of the ED reducer, so the $U(N)$ (with coupled trace modes) and Fourier-flavor conventions are verified.

## 3. Results at $N=2,3,4$

**Class totals** $I_c=\sum_\lambda I_{c,\lambda,\omega}\dim\lambda$ summed over flavor, and the closed form (D6): $I_c(N)=\tfrac23\,3^{3N^2/2}\cos\!\big(\tfrac{2\pi c}3+\tfrac{\pi N^2}2\big)$:

| $N$ | $c=0$ | $c=1$ | $c=2$ | pattern |
|---|---|---|---|---|
| 2 | $486$ | $-243$ | $-243$ | $2{:}1{:}1$ |
| 3 | $0$ | $-3^{13}$ | $+3^{13}$ | $0{:}1{:}1$ |
| 4 | $2\cdot3^{23}$ | $-3^{23}$ | $-3^{23}$ | $2{:}1{:}1$ |

This is the Turiaci–Witten profile $\cos\big(\pi(k-k_*)/3\big)$ about half filling $k_*=\tfrac32N^2$: a three-sector window $\{k_*-1,k_*,k_*+1\}$ with weights $1{:}2{:}1$ for even $N$, a two-sector window $\{k_*\pm\tfrac12\}$ with weights $1{:}1$ for odd $N$ — realised kinematically here exactly as in $\mathcal N=2$ SYK (FGMS 5.7), *provided* the BPS states sit at those degrees. The $N=2$ ED window $\{5,6,7\}$ is the even case.

**Structure of the refined index (all three $N$):**
* Within a class every complex has the **same sign** (12/12, 75/75, 639/639): no cancellations between complexes.
* The three flavor charges carry exactly equal weight ($3^{12}$ each per class at $N=3$; $2\cdot3^{22}$ / $3^{22}$ at $N=4$).
* At $N=3$ the index vanishes **in every complex** of class $c=0$ (150 non-zero complexes, all in $c=1,2$). Under concentration this says the class-$0$ sectors $k=12,15$ (and all other $k\equiv0$) carry no BPS states at $N=3$ and the window is $\{13,14\}$.
* Multiplicities per irrep (flavor $\omega=0$): largest $|I|=18$ ($N=2$, adjoint), $891$ ($N=3$), $347\,328$ ($N=4$); number of irreps per class $4,\,25,\,213$; the largest complexes sit in big self-conjugate irreps ($(4,0,-4)$ dim 125 and $(3,0,-3)$ dim 64 at $N=3$; $(6,2,-2,-6)$ dim 15 625 and $(5,2,-2,-5)$ dim 7 020 at $N=4$). The **gauge-singlet** index is $6,\,81,\,1944$ for $N=2,3,4$ (class 0 / 1 / 0), i.e. singlet BPS states exist and their number grows.

**Contrast with the single-matrix model** (same code, $p=1$): per class exactly **one** irrep — the maximal-Casimir irrep of Chen's solution: adjoint $(1,-1)$ at $N=2$ with multiplicities $1,1,2$; $(2,0,-2)$ (dim 27) at $N=3$ with $9,9$; $(3,1,-1,-3)$ (dim 729) at $N=4$ with $3,3,6$ — and singlet index $0$. Total BPS numbers are exponential in $N^2$ in both models ($3^{3N^2/2}$ vs $3^{N^2/2}$), but in the single-matrix model they form one huge irrep with $O(1)$–$O(10)$ multiplicity, whereas in the three-matrix model they are spread over hundreds of irreps with multiplicities themselves growing like $e^{cN^2}$ ($c\approx0.75$ from $18\to891\to3.5\times10^5$), including singlets. That is the "many microstates per irrep" structure Chang–Lin associate with black-hole-like BPS sectors (Q1c), now established at the level of the index for $N\le4$.

## 4. What the index cannot say, and the hand-over

The index fixes the number of BPS multiplets per complex (assuming saturation) but not the degree $k$ at which they sit; concentration is precisely the statement that each complex has a single degree. The $N=2$ gate shows saturation there, where it is forced. For $N\ge3$ the location of the window needs either cohomology ranks (plan §8) or the bootstrap **exclusion** test (positivity + sector + $\phi(XQ)=\phi(\bar QX)=0$ infeasible ⇒ no BPS states in sector $k$): if exclusion empties every complex down to one degree, the cohomology equals $|I|$ there and concentration is proved with the counts above. The index has thus turned the concentration question at $N=3$ into a finite list: show that $k\equiv0\ (3)$ sectors and all but one degree of each $c=1,2$ complex are empty.

Also available now for free: the index per complex at $N=5$ is within reach with a memory-lean implementation ($\sim10$ GB with int64 arrays; $N\le4$ takes seconds), and the large-$N$ asymptotics of the per-irrep distribution (saddle point of the character integral) is the natural analytic continuation (Q1c growth exponent, Q1d profile).
