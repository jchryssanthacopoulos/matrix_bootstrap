# Sector-resolved bootstrap with adjoint-valued operators — results (2026-09-16)

*Formulation: `docs/derivations.md` D3. Code: `src/fermion_matrix_model.py`, `src/sector_bootstrap.py`. Data: `results/data/sector_bootstrap_2026-09-16.json`. All numbers below were produced in-session on 2026-09-16; the scratch logs were lost when the machine had to be restarted after a memory blow-up (see §5), so they are transcribed from the session record. Every run is reproducible from the script with the stated parameters.*

## 1. What changed relative to `matrix_syk_covariant_bootstrap.md`

Two ingredients, both diagnosed on 2026-09-15:
1. **Sector restriction.** The functional lives on the fixed-$N_\Psi$ block ($\rho_k$ Hermitian on $\mathcal H_k$), instead of a functional on the whole space with $\phi(N_\Psi)=k$, $\phi(N_\Psi^2)=k^2$.
2. **Adjoint channel.** Positivity of $\phi(\mathrm{Tr}[w^\dagger w'])$ over *open-index* words, in addition to the singlet channel $\phi((\mathrm{Tr}w)^\dagger\mathrm{Tr}w')$.

Everything else (EOM, optional sector-preserving ground-state positivity, $Q,\bar Q$ in the basis, reduction to the span of touched operators) is as before.

## 2. Single-matrix model (exact solution known)

| $N$ | level | result |
|---|---|---|
| 2 | $L_{\rm adj}=L_{\rm sing}=2$ | exact in all 5 sectors (trivially: $H_k=0$ in BPS sectors) |
| 3 | $L_{\rm adj}=L_{\rm sing}=2$ | **exact in all 10 sectors**: $72,45,18,0,0,0,0,18,45,72$; $r\le5$ SDP variables; runtime $<1$ s per sector |
| 4 | $L_{\rm adj}=L_{\rm sing}=2$ | exact in 15 of 17 sectors ($180,144,108,72,36,\cdot,0,0,0,0,0,\cdot,36,\dots$); the two failures are $k=5,11$ (exact $18$, bound $0$) — the sectors **adjacent to the BPS window** $[6,10]$ |

**Ablation at $N=3$** (sectors $k=1,2,3$; exact $45,18,0$):

| constraints | $k=1$ | $k=2$ | $k=3$ |
|---|---|---|---|
| full (adjoint + singlet + GS + $Q$) | 45 | 18 | 0 |
| no adjoint channel | 45 | 0 | 0 |
| no singlet channel (hence no $Q$) | 24 | −12 | −36 |
| no ground-state positivity | 45 | 18 | 0 |
| no $Q,\bar Q$ | 45 | 18 | −9 |

Reading: the sector restriction alone already fixes $k=1$ (the old formulation gave 0 there); the adjoint channel is what fixes $k=2$; $Q,\bar Q$ supply the SUSY floor in BPS sectors; ground-state positivity is not needed at $N=3$.

**Level scan in the $N=4$ edge sector $k=5$** (exact 18): levels 2 and 3 give the trivial bound 0 (with or without GS positivity); level 4 lifts it: $1.91$ ($L_{\rm adj}=4$, $L_{\rm sing}=2$), $5.72$ ($4,4$), $11.42$ ($4,4$ + GS positivity, $r=46$, 142 s with SCS). Level 5 did not complete (§5). So at $N=4$ the edge of the window is resolved only slowly in level, and GS positivity does help there.

**Important structural fact discovered:** in the single-matrix model *all* EOM $\phi([H,X])=0$ with gauge-invariant $X$ are identically satisfied, because every gauge-invariant operator commutes with $\hat C_2$ and $H=3N(N^2-1)-9\hat C_2$. The exactness at $N=3$ therefore comes from positivity + normalisation + sector restriction alone, and the single-matrix model does **not** test the EOM machinery. That machinery is exercised for the first time by the three-matrix model.

## 3. Three-matrix model at $N=2$ (Chen's $C$; ED available)

Level $L_{\rm adj}=L_{\rm sing}=2$, $L_{\rm eom}=3$, SCS:

| $k$ | $d_k$ | bound | exact | $r$ | EOM |
|---|---|---|---|---|---|
| 0 | 1 | 98 | 98 | 1 | 0 |
| 1 | 12 | 22 | 22 | 18 | 489 |
| 2 | 66 | 0 | 5.165 | 267 | 1008 |
| 3 | 220 | 0 | 1.227 | 368 (604 with GS) | 1008 |

$k=0,1$ reproduce the exact analytic anchors of D2.8 ($16N^3-15N$, $16N^3-53N$). The first genuinely interacting sectors ($k=2,3$, where $H_4$ acts) sit on the SUSY floor at level 2. Level 3 was attempted and is what exhausted the machine's memory (§5); it has not been completed. $k=4$ (lifted, $d=495$) and the BPS sectors $5,6$ were not run.

## 4. Assessment

- The formulation is correct and, where it converges, exact: a purely gauge-invariant set of expectation values, organised by open-index Gram matrices, bounds the energy of *non-singlet* sector ground states tightly. This settles the "master-field obstruction" worry of the earlier notes for the single-matrix model at $N\le4$ (away from the window edge).
- The physically decisive sectors — those just outside the BPS window — are the hardest: at $N=4$ the bound moves off zero only at level 4. For the three-matrix model at $N=2$ the analogous sectors ($k=2,3,4$) are unresolved at level 2. Whether level 3–4 suffices there is the open empirical question, and it is now a memory/engineering question rather than a formulation question.

## 5. Memory post-mortem and required fixes

Two independent mistakes in the first implementation caused a machine-wide out-of-memory event on 2026-09-16 while running (a) the $N=4$, level-5 single-matrix edge sector in the background and (b) the three-matrix $N=2$ level-3 sectors in the foreground:

1. **Open words were built on the full $2^n$-dimensional space** and only then restricted to the sector. A length-$L$ word maps each basis state to $\lesssim N^{L-1}$ states, so at $N=4$, $L=5$ a single matrix entry has $\sim10^7$ non-zeros, and there are $N^2$ entries per word and dozens of words. **Fix:** letters map $\mathcal H_k\to\mathcal H_{k\pm1}$; build every word as a product of *rectangular sector-to-sector blocks*, never touching the full space.
2. **The reduction formed a dense complex $n_h\times n_h$ Gram matrix** of all touched operators. For the three-matrix model at level 3 ($n_h\approx3\times10^4$) this is $\sim12$ GB before intermediates (the sparse product `F.conj() @ F.T` materialises a second copy). **Fix:** deduplicate operators, accumulate the real Gram matrix in row blocks into a preallocated `float64` array, and **estimate the footprint before building** ($8n_h^2$ bytes for $\mathcal G$ plus $\sim3\times$ for `eigh`, plus $24\,\mathrm{nnz}(F)$ bytes), refusing to run above a user-set budget.

Both fixes were implemented later on 2026-09-16 (`word_blocks`, `_reduced_basis` with two routes — Gram over the $n_h$ operators or over the $2d_k^2$ real components, whichever is smaller — a `budget_gb` guard and a `dry_run` mode). The estimator was found to be ~20 % low when nnz$(F)$ dominates; keep a margin.

**A false shortcut, recorded so it is not repeated.** Restricting to *normal-ordered* open words (all $\Psi$ left of all $\bar\Psi$) is **not** lossless: reordering fermion operators keeps the index-contraction pattern of the original word, so e.g. $(\Psi\bar\Psi\Psi)_{ij}=N\Psi_{ij}-\Psi_{im}\Psi_{nj}\bar\Psi_{mn}$, and the last term is not a matrix product of reordered letters. With that restriction the $N=4$, $k=5$ bound collapsed from $11.4$ (level 4) to $0$ (levels 5–6). All results quoted in this note use all orderings.

## 6. Runs after the memory fixes (2026-09-16, all orderings, single process each, budget 6 GB)

| model | sector | level ($L_{\rm adj},L_{\rm sing}$, GS) | bound | exact | $r$ | peak RSS |
|---|---|---|---|---|---|---|
| single matrix $N=4$ | $k=5$ | (4,4,on) | 11.417 | 18 | 46 | 3.9 GB |
| single matrix $N=4$ | $k=5$ | (5,4,on) | 14.734 (SCS "inaccurate") | 18 | 68 | 6.9 GB |
| three-matrix $N=2$ | $k=2$ | (3,3,off) | 4.548 (SCS "inaccurate", 88 min, 23 970 EOM) | 5.165 | 270 | 5.7 GB |

Three-matrix $k=2$: level 2 gave the SUSY floor 0; level 3 recovers 88 % of the exact value — the EOM (non-trivial here, unlike the single-matrix model) are doing real work. The same SDP with CLARABEL (interior point) gave $4.559$ in 37 min but at a **9.9 GB peak** (the KKT system for $\sim2.4\times10^4$ equality rows); the two solvers agree to $0.2\,\%$, i.e. within their tolerances (both report "inaccurate"), so the level-3 value is $\approx4.55$. Lesson: the `budget_gb` guard covers only the reduction step; interior-point solver memory scales with the number of equality rows and must be budgeted separately (or the redundant EOM rows must be pruned — many of the 23 970 rows are linearly dependent, since $r=270$). For now SCS is the memory-safe choice.

Edge-sector convergence at $N=4$: $0\ (L\le3)\to11.4\ (L=4)\to14.7\ (L=5)$ — monotone but slow; level 6 is estimated at $>8$ GB with the current $d^2=1.9\times10^7$-column sparse representation and was not run.

Not feasible within 6 GB with the present reduction: three-matrix $N=2$, $k=3$ at level 3 (dry-run estimate 13–22 GB). This requires the symmetry-reduced formulation (gauge-invariant $\rho_k$ via Ward identities, $\mathbb Z_3$-Fourier flavor basis), which shrinks both the operator space ($d_k^2\to\sum_Rm_R^2$) and the word count.

