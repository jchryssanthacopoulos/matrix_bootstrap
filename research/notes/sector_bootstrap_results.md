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


## 7. Allowed-region ("island / archipelago") scans (2026-09-16, evening)

Tooling: `scripts/scan_sector.py` (fix $\phi(H)=E$, test feasibility — optionally with the eigenstate constraints $\phi(OH)=E\phi(O)$ — and min/max chosen observables; SCS), `scripts/plot_scans.py`. Data: `results/data/scan_*.json`, `results/data/island_widths_p3_N2_k2_L2.txt`. Figures: `results/figures/sector_scan_single_matrix_N3.png`, `results/figures/sector_scan_three_matrix_N2_k2.png`.

**Single matrix, $N=3$ (calibration).**
- *Archipelago.* With the eigenstate constraints, the feasible energies are: $N_\Psi=2$, level 2 → exactly $\{18,45\}$; $N_\Psi=3$, level 2 → the whole interval $[0,72]$ (four irreps, too few independent Casimir moments at this level to separate four weights); $N_\Psi=3$, level 3 → exactly $\{0,18,45,72\}$. Islands have zero width on a grid of step 0.25 (they are points), i.e. the eigenstate constraints pin the levels completely once the level suffices.
- *$(E,\langle\hat C_2\rangle)$ region.* Already at level 2, without eigenstate constraints, the region has zero width and coincides with the exact line $C_2=(72-E)/9$ — because the operator identity $H=72-9\hat C_2$ lies in the touched span. (The plain scan's feasible $E$-range in $N_\Psi=2$ is $[18,72]$ although the sector's spectrum is $\{18,45\}$: without eigenstate constraints, mixtures and the unphysical top are not excluded.)

**Three-matrix model, $N=2$, $N_\Psi=2$ ($\dim 66$, exact spectrum $5.165,\ 12,\ 13.49,\ 18.33,\ 20,\ 22.34,\ 78,\dots$; 15 distinct levels).**
- *Regions at level 2* (figure): $(E,\langle\hat C_2^{\rm gauge}\rangle)$ and $(E,\langle N_s\rangle)$ regions contain every exact eigenstate (a non-trivial correctness check of the whole pipeline on an interacting model), extend down to the SUSY floor $E\approx0$ (consistent with the level-2 bound $E_0\ge0$) and up to a loose top ($120.6$ vs exact $\max E\approx98$). The Casimir region is loose (allows $C_2$ up to 12 where the sector has at most 6). $\langle N_s\rangle$ is not conserved, so its exact "points" are ranges within each level; the low-lying levels have $\langle N_s\rangle\gtrsim1$ — mild support for the guess (D2.8) that low-energy states are rich in the flavor-symmetric mode, but at $N=2$ only. $\langle\bar QQ\rangle$ was dropped from the panels: for $N_\Psi<3$, $Q\bar Q\equiv0$ in the sector, so $\langle\bar QQ\rangle=E$ identically.
- *Archipelago at level 2.* Remarkably, with the eigenstate constraints the level-2 problem is feasible **only at the exact eigenvalues**: feasible at $E_0=5.16536436$, $12$, $13.4907418$ (each to machine precision), infeasible at $E=4$, at the midpoint between the first two levels, and everywhere beyond $\pm1.5\times10^{-5}$ of each level (bisection). So the eigenstate constraints already determine the spectrum of this sector at level 2, whereas the convex ground-state *bound* at the same level is 0. The price is that $E$ enters bilinearly: one can only scan, the islands are points, and a blind grid would miss them (the single-matrix grid found them only because the eigenvalues were integers). This is the same phenomenon as Lin–Zheng 2024 App. B / Han–Hartnoll–Kruthoff Fig. 1, in its extreme finite-dimensional form (cf. HHK footnote [20]: with enough constraints the feasible set collapses onto the exact eigenstates).

**Implications.** (i) The convex, scannable object of interest remains $E_0(k)$ (and its irrep refinements) — that is what large-$N$ needs. (ii) The eigenstate-constrained scan is a powerful *finite-$N$* tool: one could locate levels by a root-finding strategy on an infeasibility measure rather than a grid, or use it to certify that a candidate $E$ is *not* an eigenvalue. (iii) For the physics question the relevant sectors are the window-edge ones; there the level needed grows with $N$ (single matrix: level 4–5 at $N=4$), and for the three-matrix model at level $\ge3$ the PSD cones ($84\times84$) make interior-point solvers memory-hungry (CLARABEL exceeded 17 GB on the $k=2$ level-3 problem after EOM pruning; SCS solves it in $\sim1.5$ h at $<6$ GB). Symmetry reduction of the Gram blocks ($\mathbb Z_3$ flavor, $SU(N)$ irreps) is the way to shrink the cones.

**Solver notes.** CLARABEL fails numerically on some of these (not strictly feasible) SDPs where SCS succeeds; `SectorSDP.solve` now retries with a bounded SCS run. EOM rows are pruned to an independent set before the SDP (e.g. $23\,970\to226$ at level 3), which removes the equality-row memory problem; the remaining cost is the PSD cone size.

## 8. Symmetry-reduced bootstrap (2026-09-17): three-matrix $k=2$ exact at level 3, first $k=3$ bound

*Formulation: `docs/derivations.md` D4. Code: `src/symmetry_reduction.py`, `symmetry=True` / `adjoint_projected=True` in `src/sector_bootstrap.py`, runner `scripts/run_sector_bound.py` (one JSON line per run in `results/data/sector_bounds_symmetric_2026-09-17.jsonl`), memory guard `scripts/watchdog.sh`. All runs: Chen's $C$, $N=2$, $p=3$, Fourier flavor letters, budget 10 GB, watchdog kill at 12 GB, single process.*

### 8.1 What was implemented

1. **Invariant functional.** $\rho_k=\oplus_R\sigma_R\otimes\mathbf 1$ over the isotypic components of the $k$-block under $SU(2)_{\rm gauge}\times\mathbb Z_3$; every operator is reduced to its $D=\sum_Rm_R^2$ invariant components the moment it is formed (streaming), and the touched span and the EOM span are accumulated as $2D\times2D$ Gram matrices. Nothing of size $d^2$ or $n_h\times D$ is stored. $D$: 6, 90, 666, 2727, 6156, 8112 for $k=1..6$ (versus $d^2$ = 144, 4356, 48 400, 245 025, 627 264, 853 776).
2. **$\mathbb Z_3$ Fourier letters** $\Psi^{(m)}$ carry a $\mathbb Z_3$ charge $z$; Gram cones are graded by $(q,z)$, so each level-3 cone has $\le28$ words instead of $84$.
3. **Adjoint-projected channel** (Cho et al. 2024 eq. 3.5): for an invariant $\rho$ the singlet and adjoint parts of $\langle(w_{ij})^\dagger w'_{kl}\rangle$ are separately PSD, so $\phi(\mathrm{Tr}[w^\dagger w'])-\tfrac1N\phi((\mathrm{Tr}w)^\dagger\mathrm{Tr}w')\succeq0$ is imposed in place of the plain adjoint channel. This was not available before because the old functional was not gauge invariant.
4. **Conditioning.** Each Gram cone is rescaled by a diagonal congruence (unit-norm diagonal coefficient rows). CLARABEL converges only with chordal decomposition *and* equilibration disabled (`solver_options`); SCS is always run with a bounded iteration count.
5. **Diagnostics.** `SectorSDP.exact_coordinates()` evaluates every constraint on the multiplet-averaged exact ground state. This caught a real bug (a round-off "Hermitian part" of an anti-Hermitian commutator was normalised into a spurious EOM row, making the SDP infeasible); the cut is now relative to $\|X\|$ in all code paths (D4.3).

Regression: unreduced and reduced SDPs agree on every earlier number (single matrix $N=2,3$ all sectors; three-matrix $k=1$: 22, $k=2$ level 2: 0). With `ward=True` no extra rows appear (Ward identities are implied by invariance).

### 8.2 Results (three-matrix Chen model, $N=2$; exact $E_0(k)$ from ED: $98,\,22,\,5.16536,\,1.22706,\,0.08796,\,0,0,0,\dots$)

| $k$ | level $(L_{\rm adj},L_{\rm sing},L_{\rm eom})$ | channel | bound | exact | status | $r$ / $D$ | EOM rows raw→indep. | build / solve | peak RSS |
|---|---|---|---|---|---|---|---|---|---|
| 1 | (2,2,3) | adj-proj | 22.00000 | 22.00000 | optimal | 6 / 6 | 0 | 1 s / 0 s | 0.15 GB |
| 2 | (2,2,3) | adj-proj | 0.00000 | 5.16536 | optimal | 87 / 90 | 468→46 | 1 s / 0 s | 0.19 GB |
| 2 | (3,3,4) | plain adjoint | 4.55415 | 5.16536 | SCS inaccurate | 90 / 90 | 14 862→66 | 14 s / 255 s | 0.94 GB |
| 2 | (3,3,4) | **adj-proj** | **5.16536** | 5.16536 | **optimal** | 90 / 90 | 15 102→66 | 15 s / 14 s | 1.0 GB |
| 3 | (2,2,3) | adj-proj | 0.00000 | 1.22706 | optimal | 122 / 666 | 468→60 | 3 s / 0 s | 0.39 GB |
| 3 | (3,3,4) | adj-proj | **0.95350** | 1.22706 | SCS inaccurate | 666 / 666 | 20 718→583 | 70 s / 1910 s | 6.3 GB |
| 3 | (3,3,4) | adj-proj, CLARABEL | 0.96218 | 1.22706 | CLARABEL inaccurate | 666 / 666 | 20 718→583 | 69 s / 68 s | 9.0 GB |
| 3 | (3,3,4) + GS | adj-proj, CLARABEL | **0.96218** | 1.22706 | **optimal** | 666 / 666 | 20 718→583 | 69 s / 71 s | 9.2 GB |
| 3 | **(4,3,4)** | adj-proj | **1.22706** | 1.22706 | **SCS optimal** | 666 / 666 | 304 146→592 | 993 s / 5286 s | **18.3 GB** (ru_maxrss; 17.4 GB sampled) |
| 4 | (3,3,4) | adj-proj | 0.00000 | 0.08796 | SCS optimal (CLARABEL failed) | 1218 / 2727 | 20 718→649 | 502 s / 167 s | 10.1 GB |

Cross-checks on the $k=2$, level-3 problem: plain adjoint channel with CLARABEL 5.1737 (inaccurate, chordal off), 5.16536 (optimal, chordal + equilibration off, 22 s), SCS 5.16536 (optimal, 13 s); the plain-channel value 4.554 agrees with the unreduced runs of §6 (4.548 SCS / 4.559 CLARABEL, 88 min / 37 min, 5.7 / 9.9 GB) — the same bound, now in 4 min and under 1 GB.

### 8.3 Reading of the results

* **$k=2$ is solved exactly at level 3** once the adjoint-projected channel is used: $E_0(2)=5.16536$ to solver precision. The plain adjoint channel saturates at 88 % at the same level. Since $r=D=90$ already at level 3, the functional is fully general within the invariant class; the remaining gap of the plain channel was a *constraint* deficit, and the singlet/adjoint split supplies exactly the missing information. This is the first three-matrix sector outside the trivially exact ones ($k=0,1$: SUSY floor / one-particle levels of D2) obtained exactly by the bootstrap.
* **$k=3$ at level 3 for the first time:** $E_0(3)\ge0.9622$ (exact 1.22706, 78 %). Level 2 gives 0 (SUSY floor) as for $k=2$. Three solves agree: SCS 0.9535 (inaccurate, 32 min), CLARABEL 0.96218 (inaccurate, 68 s), CLARABEL with GS blocks 0.96218 (*optimal*, 71 s) — so 0.9622 is the accurate level-3 optimum and GS positivity adds nothing here, as in the single-matrix model. $r=D=666$ is saturated at level 3, so the remaining gap is a constraint deficit and the next lever is level 4 words.
* **$k=3$ is solved exactly at level 4** (adjoint words of length $\le4$, singlets $\le3$): $E_0(3)=1.22706$ to solver precision, SCS status *optimal*. The 45 cones (largest $169,168,168,111,\dots$) add $3\times10^5$ raw EOM rows but only 9 new independent ones (583→592) — the new information is almost entirely in the positivity cones, not the EOM. Together with $k=2$ this gives the three-matrix pattern: **each lifted sector becomes exact one level after it leaves the SUSY floor** ($k=2$: $0\to$ exact at levels 2→3; $k=3$: $0\to78\%\to$ exact at levels 2→3→4), provided the singlet/adjoint split is used.
* **Cost warning.** The level-4 run peaked at 18.3 GB against a plan estimate of 6.4 GB: the estimate covers the stored Gram entries and coefficient arrays with a factor 3 for cvxpy, but SCS's own copies of the $3.4\times10^5\times666$ constraint matrix (real embedding, plus its normalised and factorised forms) roughly triple that again. Rule of thumb from these runs: **peak ≈ 3 × plan estimate** for SCS at level 4; CLARABEL is worse (9 GB already at level 3 for $r=666$). Level 4 with $L_{\rm sing}=4$, or level 4 at $k=4$ ($D=2727$), is therefore *not* feasible under 20 GB with the present cvxpy pipeline; it would need a direct SCS/SDPA interface without cvxpy's intermediate copies, or a coefficient-side compression (the $r=666$ coordinates of the $1.7\times10^5$ entries are not sparse).
* **$k=4$ at level 3:** bound 0 against $E_0(4)=0.08796$ — the level-3 constraints cannot lift the almost-BPS sector off the SUSY floor; here $r=1218<D=2727$, so both the functional and the constraints are still unsaturated. (The build took 8 min because `reduce` costs $O(d\cdot D)$ per operator at $d=495$; the SDP itself took 3 min.)
* **Costs.** The build is now negligible (70 s for $k=3$ level 3, which was estimated at 13–37 GB before the reduction and never ran). The cost has moved entirely to the SDP solve: the number of variables is $\le D$, but the total cone size $\sum m^2\approx10^4$ (level 3) $\to3\times10^5$ (level 4) with dense $r$-vectors per entry. $k=3$ level 4 with $L_{\rm sing}=3$ is estimated at $\sim7$ GB of coefficient data (with cvxpy overhead), $L_{\rm sing}=4$ at $\sim12$ GB — inside the 20 GB ceiling but not the 10 GB one. $k=4$ level 3 ($D=2727$) is estimated at 2 GB and is the natural next run; $k=5,6$ are BPS ($E_0=0$) and the SUSY floor already makes the bound exact there, as at $k=0$.
