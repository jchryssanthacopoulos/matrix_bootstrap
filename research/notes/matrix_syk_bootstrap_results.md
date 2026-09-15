# Bootstrapping matrix SYK models — II. Python validation on $Q=\mathrm{Tr}[\Psi^3]$

*Companion to the formulation note. Implements formulation 2.3(a) — the eigenstate bootstrap evaluated on/against the fixed reference state $|\lambda\rangle$ — in Python, and validates it on Chen's exactly-solvable single-matrix model. Everything here is $N=2$ (Fock dimension $2^{N^2}=16$); the exact answers from the formulation note are the ground truth.*

Files: `mm_model.py` (model + reference state + ED), `mm_bootstrap.py` (bootstrap SDP), `bootstrap_run_N2.txt` (raw run log). Requires `numpy`, `scipy`, `cvxpy`.

## Method as implemented

The bootstrap variable is a **linear functional** $\phi(O)=\mathrm{Tr}[\rho\,O]$ with $\rho$ a Hermitian $16\times16$ matrix that is **not** assumed positive — so $\rho$ is just the most general linear functional on operators, and any bound returned is a genuine truncated bootstrap bound, not exact diagonalization. The truncation is entirely in *which* constraints we impose:

- **Positivity** of the moment matrix $M_{IJ}=\phi(B_I^\dagger B_J)\succeq0$, where $B$ is a set of fermion words (deduplicated to a linearly independent subset). Optionally the supercharge $Q,\bar Q$ are appended to $B$.
- **Equations of motion** $\phi([H,C])=0$ for all words $C$ up to a level (192 independent EOM constraints at $N=2$).
- **Normalization** $\phi(1)=1$.
- (Formulation 2.3(a) targeting) optional *sharp* constraints $\phi(N_\Psi)=q_*$, $\phi(N_\Psi^2)=q_*^2$ and $\phi(\hat C_2)=C_2^{\max}$, $\phi(\hat C_2^2)=(C_2^{\max})^2$, which force the functional onto the highest-weight BPS state $|\lambda\rangle$.

Objective: minimize $\phi(H)$ (rigorous lower bound on $E_0$), or maximize $\phi(H)$ / $\phi(\hat C_2)$.

The moment matrix is built vectorized ($M=\mathrm{reshape}(A\,\mathrm{vec}\,\rho)$) to keep memory small; the SDP is solved with SCS.

## Results ($N=2$: exact $E_0=0$, $E_\text{top}=18$, $C_2^{\max}=2$)

**1. The truncation is real.** With *words only* as the operator basis:

| operator basis | moment size $m$ | $\min\phi(H)$ |
|---|---|---|
| words, level 1 | 9 | $-\infty$ (unbounded) |
| words, level 2 | 37 | $-2.250$ |

A level-2, words-only bootstrap returns $-2.25$, strictly below the true $E_0=0$ — confirming this is a genuine relaxation, not diagonalization in disguise.

**2. SUSY positivity pins $E_0=0$.** Appending the supercharge $Q,\bar Q$ to the operator basis makes $\phi(H)=\phi(\bar QQ)+\phi(Q\bar Q)$ a sum of moment-matrix diagonals, so positivity forces $\phi(H)\ge0$; combined with the existence of the true BPS state this pins the bound:

| operator basis | $m$ | $\min\phi(H)$ |
|---|---|---|
| level 1 $+\{Q,\bar Q\}$ | 11 | $-7\times10^{-12}$ |
| level 2 $+\{Q,\bar Q\}$ | 39 | $+7\times10^{-7}$ |

This is exactly the sum-of-squares mechanism anticipated in §2.5 of the formulation note.

**3. Top of the spectrum.** Maximizing instead of minimizing, $\max\phi(H)=18.0000=3N(N^2-1)$ (level 2 $+\{Q,\bar Q\}$) — the singlet ($C_2=0$) states at the top of the spectrum.

**4. Maximal Casimir.** $\max\phi(\hat C_2)=2.00000=N(N^2-1)/3$ already at level 1 $+\{Q,\bar Q\}$. Since $H=3N(N^2-1)-9\hat C_2$, bounding the energy from below is identically bounding the Casimir from above; the bootstrap recovers the largest representation $\mathbf r_*$.

**5. The constraints are correct.** Feeding the *exact* $|\lambda\rangle$ functional $\rho_\lambda=|\lambda\rangle\langle\lambda|$ into the constraints: min eigenvalue of the moment matrix $=-2\times10^{-16}$ ($\ge0$ ✓) and $\max|\phi([H,C])|=0$ ✓. The true ground state is feasible, so no constraint is spuriously over-tightening.

**6. Formulation 2.3(a) targeting.** Imposing sharp $R$-charge $q_*=N(N-1)/2=1$ and sharp Casimir $C_2^{\max}=2$ selects the highest-weight BPS multiplet and returns $\phi(H)=3\times10^{-8}\approx0$ — the bootstrap lands on $|\lambda\rangle$'s energy.

## Reading of the result

For this solvable model the bootstrap validation is essentially analytic: **positivity + the SUSY algebra $\{Q,\bar Q\}=H$ ⟹ $E\ge0$ ⟹ $\langle\hat C_2\rangle\le N(N^2-1)/3$**, saturated by the fortuitous states. The numerics confirm this and, importantly, show that *without* the supercharge in the basis the bound is loose ($-2.25$), so the machinery — operator enumeration, Pauli/normal-ordering via explicit matrices, moment-matrix positivity, EOM, and 2.3(a) state targeting — is doing real work and is implemented correctly.

## R-charge-resolved ground energy $E_0(N_\Psi)$ — the concentration observable

The physically pointed observable (chosen for the matrix SYK target) is the ground-state energy in each fixed R-charge sector, $E_0(N_\Psi)$: the set of $N_\Psi$ with $E_0=0$ is exactly the BPS window, and its width vs. $N$ is what distinguishes concentrated (SYK) from spread (single-matrix) fortuity. For the single-matrix model this is $E_0(N_\Psi)=3N(N^2-1)-9\,C_2^{\max}(N_\Psi)$, so bootstrapping it = bounding the maximal Casimir at fixed fermion number. Imposing sharp charge $\phi(N_\Psi)=k$, $\phi(N_\Psi^2)=k^2$ and maximizing $\phi(\hat C_2)$ gives a rigorous lower bound on $E_0(k)$ (`sm_charge_profile.py`):

| $N_\Psi$ | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| $E_0$ exact (ED) | 18 | 0 | 0 | 0 | 18 |
| bootstrap lower bound, level 1 $+Q$ | 0 | 0 | 0 | 0 | 0 |
| bootstrap lower bound, level 2 $+Q$ | **18** | **0** | **0** | **0** | **18** |

At level 2 $+\{Q,\bar Q\}$ the bound is already tight: it reproduces the BPS window $N_\Psi\in[N(N-1)/2,N(N+1)/2]=[1,3]$ (where $E_0=0$) and the lifted sectors $\{0,4\}$ (where $E_0=18$) exactly. This validates the exact observable we will carry to the matrix SYK model, on the solvable case.

## Scaling note and next steps

- $N=2$ ($\dim=16$) runs comfortably. $N=3$ ($\dim=512$) is memory-heavy in this flat implementation because operator coordinate vectors have length $\dim^2$; the fix is to (i) work within fixed $R$-charge blocks and (ii) use an SVD-reduced operator basis rather than the full $\dim^2$ space. $N=3$ exact targets ($E\in\{0,18,45,72\}$, $C_2^{\max}=8$) are already in hand from ED (`mm_model.py`) for cross-checking once the blocked SDP is in place.
- With the machinery validated, the next step is the **matrix SYK target**: swap the single $\Psi$ for $p=3$ matrices with $Q=\sum_{i\le j\le k}\mathrm{Tr}[\Psi_i\Psi_j\Psi_k]$, add the $O(p)$ flavor blocking, and (crucially) drop the Casimir shortcut — there $\min\phi(H)$ is a real unknown and the bound is the physics output. At $N=2$ ($2^{12}=4096$) it remains ED-checkable.

### Sources (project library)
- Y. Chen, *Fortuity with a Single Matrix*, arXiv:2511.00790.
- S. Laliberte, B. McPeak, *Bootstrapping supersymmetric (matrix) quantum mechanics*, arXiv:2510.01356.
- H. W. Lin, Z. Zheng, *High-Precision Bootstrap of Multimatrix Quantum Mechanics*, arXiv:2507.21007.
