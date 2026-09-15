# Matrix SYK — IV. Toward the large-N bootstrap: baseline, trace Hamiltonian, and the real obstruction

*Building the single-trace ("planar") bootstrap engine, validating on the single-matrix model, targeting $E_0(N_\Psi)$ per R-charge sector. Files: `sm_charge_ED.py`, `sm_trace_hamiltonian.py`.*

## 1. The concentration question, made precise (single-matrix baseline)

Sparse ED of the single-matrix model $Q=\mathrm{Tr}[\Psi^3]$ at $N=2,3,4$ (up to $\dim=2^{16}$) gives the R-charge-resolved ground energy $E_0(N_\Psi)$ exactly. The BPS window (where $E_0=0$) is:

| $N$ | BPS window $[\tfrac{N(N-1)}2,\tfrac{N(N+1)}2]$ | # sectors | $\dim$ |
|---|---|---|---|
| 2 | $[1,3]$ | 3 | 16 |
| 3 | $[3,6]$ | 4 | 512 |
| 4 | $[6,10]$ | 5 | 65 536 |

The window has **$N+1$ sectors — it grows linearly with $N$.** Outside it, $E_0(N_\Psi)$ climbs in steps of 18 (the Casimir ladder). So the single-matrix model provably does **not** R-charge-concentrate in the SYK sense (an $O(1)$-width window).

This sharpens the target. The $p=3$ matrix SYK model at $N=2$ has its BPS states in $N_\Psi\in\{5,6,7\}$ — also 3 sectors, *identical* to the single-matrix $N=2$ window. The two models are only distinguished by how the window scales:

> **The concentration question is: does the SYK BPS window stay $O(1)$ as $N$ grows, while the single-matrix window grows as $N$?** This cannot be answered by ED (the next SYK point, $N=3$, is $2^{27}$) — it is exactly the large-$N$ statement the bootstrap must supply.

## 2. Engine foundation: $H$ in single-trace form

A single-trace bootstrap never builds the $2^{N^2}$-dimensional operators; it works with traces of words in $\Psi,\bar\Psi$. Matching to the explicit Hamiltonian at $N=2$ and $N=3$ (residual $<10^{-11}$) gives the exact identity

$$
\boxed{\,H = \tfrac{9}{2}\Big(\mathrm{Tr}[\Psi\Psi\bar\Psi\bar\Psi] + \mathrm{Tr}[\bar\Psi\bar\Psi\Psi\Psi]\Big) - 9(N-1)^2\,}
$$

Two features matter for the bootstrap:

- The single-trace quartic coefficient $9/2$ is **$N$-independent**, and the remainder is a **pure constant** (the naive double-trace pieces from $\hat C_2$ collapse via $\mathrm{Tr}[\Psi\bar\Psi]+\mathrm{Tr}[\bar\Psi\Psi]=N^2\mathbb 1$).
- Hence the dynamics is **purely single-trace** — the most favorable structure for large-$N$: the equations of motion $\langle[H,O]\rangle=0$ close on single-trace correlators without genuine double-trace couplings.

This is the reusable input the engine needs, and it will be re-derived the same way for the 3-matrix model (where flavor structure enters but the method is identical).

## 3. The real obstruction — and why it is the physics, not a bug

There is a genuine reason the standard planar bootstrap does not simply drop in here, and it is worth stating plainly because it *is* the content of fortuity.

The fortuitous ($E=0$) states carry the **maximal color Casimir**, $C_2\sim N^3$ — they are as far from a gauge singlet as a state can be. Standard large-$N$ technology (Lin–Zheng, Laliberte–McPeak, and the collective-field / master-field picture generally) assumes the relevant state is a gauge singlet described by a classical master field, for which single-trace expectation values factorize. That assumption **fails** for the fortuitous sector:

- the ground state is non-singlet (the singlet sector of this model sits at the *top* of the spectrum, $E=3N(N^2-1)$);
- a maximal-Casimir state has no classical master-field description, so naive single-trace factorization is not guaranteed.

This is precisely Chen's own point that fortuitous states are *not* captured by a classical large-$N$ saddle. So the bootstrap of the fortuitous sector is genuinely less-charted than the published singlet-sector matrix bootstraps — the difficulty is physical, not incidental.

**Consequences for the engine design** (the two viable routes):

1. **Covariant / R-charge-resolved moment matrix.** Work at fixed $N_\Psi$ with gauge-*covariant* operators (open color indices, adjoint-valued), imposing positivity of the multiplet-averaged moment matrix $\langle O_i^\dagger O_j\rangle$ (formulation 2.3(b)) together with $\langle[H,O]\rangle=0$ and the gauge Ward identities. This is the natural home for $E_0(N_\Psi)=3N(N^2-1)-9\,C_2^{\max}(N_\Psi)$ and it does not presuppose a master field. Since the dynamical $H$ is single-trace (§2), the EOM stay tractable.
2. **Direct Casimir-maximization** as an independent check: $E_0(N_\Psi)$ for the single matrix is *exactly* $3N(N^2-1)-9\,C_2^{\max}(N_\Psi)$, a representation-theory optimization (largest Casimir irrep in $\wedge^{N_\Psi}\mathbf{adj}$) solvable at any $N$ — the rigorous answer any bootstrap must reproduce, and the analytic baseline for the concentration scaling.

## 4. Status and next step

Delivered and verified: the exact $E_0(N_\Psi)$ baseline ($N\le4$); the single-trace Hamiltonian; and a precise, well-posed large-$N$ target (window scaling). The engine's next concrete piece is the **covariant, fixed-$N_\Psi$ moment-matrix bootstrap** (route 1), validated against the single-matrix $E_0(N_\Psi)$ above, then transported to the 3-matrix model to test whether its window stays $O(1)$ — the statement obtainable no other way.

### Sources (project library)
- Y. Chen, *Fortuity with a Single Matrix*, arXiv:2511.00790.
- C.-M. Chang, Y. Chen, B. S. Sia, Z. Yang, *Fortuity in SYK Models*, arXiv:2412.06902.
- H. W. Lin, Z. Zheng, *High-Precision Bootstrap of Multimatrix Quantum Mechanics*, arXiv:2507.21007.
- S. Laliberte, B. McPeak, *Bootstrapping supersymmetric (matrix) quantum mechanics*, arXiv:2510.01356.
