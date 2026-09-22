# BPS-exclusion bootstrap — results log

*Idea: R-charge concentration is an exclusion statement. A BPS state satisfies $Q\psi=\bar Q\psi=0$, hence the linear rows $\phi(XQ)=\phi(QX)=\phi(X\bar Q)=\phi(\bar QX)=0$ for every operator $X$; if the sector-$k$ functional with these rows (plus positivity, sector rows, reality) is infeasible, sector $k$ has no BPS state. Combined with the refined index (D6, `bps_index_results.md`), emptying every degree but one in a complex determines its BPS count exactly. Code: `bps=True` in `TraceSDP` (`src/trace_bootstrap.py`), margin solve `bps_margin`, CLI `scripts/run_bps_exclusion.py`. The margin $t^*=\max t$ s.t. all cones $\succeq t\mathbf 1$: $t^*\approx0$ = not excludable at this level (singular cones make $t^*>0$ impossible), $t^*<0$ = excluded.*

## 1. $N=2$ (window $\{5,6,7\}$; $E_0(k)=98,22,5.165,1.227,0.088,0,0$ for $k=0..6$)

| $k$ | level 2, no BPS rows | level 2, BPS rows | level 3, no BPS rows | level 3, BPS rows | energy route (for comparison) |
|---|---|---|---|---|---|
| 0 | 0 | $-0.922$ **excl.** | 0 | $-3.428$ **excl.** | exact at level 2 |
| 1 | 0 | $-0.453$ **excl.** | 0 | $-1.904$ **excl.** | exact at level 2 |
| 2 | 0 | $-0.120$ **excl.** | 0 | $-0.741$ **excl.** | floor at level 2, exact at level 3 |
| 3 | 0 | 0 | 0 | $-0.0167$ **excl.** | floor at level 3, exact at level 4 (+ finite-$N$ relations) |
| 4 | 0 | 0 | 0 | 0 | not lifted through level 4 |
| 5,6 | 0 | 0 | 0 | 0 | BPS (correct) |

The exact multiplet-averaged BPS functional satisfies every BPS row at $k=5,6,7$ (residuals $10^{-12}$), so the rows are valid; without the BPS rows the margin is $0$ everywhere (the energy ground state is always feasible), so all exclusions come from the BPS rows. Observations: (i) exclusion is strictly stronger than the energy bound — $k=2$ is excluded at level 2 and $k=3$ at level 3, each one level earlier than the energy route and without finite-$N$ relations; (ii) the required level still grows: **exclusion at level $k$** in this data, vs level $k+1$ for energies; (iii) $k=4$ ($E_0=0.088$, first sector outside the window) is not excluded at level 3 with the basic row family (32 / 1284 BPS rows at levels 2 / 3). Row-family enlargements (longer charge-$\mp3$ words, split-energy rows $\phi(X\bar QQ)=0$) tested next.

**Row-family enlargements tested on the non-excluded cases** ($k=3$ at level 2, $k=4$ at level 3): longer charge-$\mp3$ traced words (`L_bps=5`), split-energy rows $\phi(X\bar QQ)=\phi(XQ\bar Q)=0$ for all cone operators $X$ (`bps_split`), longer EOM words, and the finite-$N$ relations to length 6 — **none changes the margin** (all $0$). At fixed level the limit is the operator content of the cones, not the number of BPS rows. At $N=2$ the rule is therefore "exclusion at level $k$".

## 2. $N=3$ and $N=4$ at level 2 (2026-09-19; `results/data/bps_exclusion_2026-09-19.jsonl`)

Same constraint set as at $N=2$, evaluated at $N=3,4$ (seconds per sector):

| $N$ | excluded at level 2 | first non-excluded | margins |
|---|---|---|---|
| 2 | $k\le2$ | 3 | $-0.92,-0.45,-0.12$ |
| 3 | $k\le6$ | 7 | $-3.70,-2.89,-2.14,-1.49,-1.00,-0.52,-0.04$ (linear in $k$) |
| 4 | $k\le11$ | 12 | $-9.42\to-0.20$, again linear in $k$ |

So the **reach of a fixed level grows with $N$**, roughly like $0.7N^2$ (the margin decreases linearly in $k$ with an $N$-dependent slope) — the $N=2$ extrapolation "level $\approx k$" was misleading because at $N=2$ every sector is close to the window. By particle–hole symmetry the same sectors are excluded from the top ($k\ge3N^2-6$ at $N=3$, $k\ge37$ at $N=4$). At $N=3$ this already excludes BPS states from $k\le6$ and $k\ge21$; the index says the $k\equiv0$ classes carry zero net index, so the remaining question for the window $\{13,14\}$ is the exclusion of $k=7..12$ and $15..20$ (equivalently, by symmetry, $k=7..12$). Level 3 sweeps at $N=3$ ($k=7..13$) and $N=4$ ($k=12..24$) are running.

## 3. Level 3 at $N=3,4$ (2026-09-19/20; same data file)

| $N$ | level 2 reach | level 3 reach | level-3 margins at the edge | window predicted by index + TW profile |
|---|---|---|---|---|
| 2 | $k\le2$ | $k\le3$ | $-0.0167$ ($k=3$) | $\{5,6,7\}$ (ED: correct) |
| 3 | $k\le6$ | $k\le6$ | $-1.32$ ($k=6$), $0$ ($k=7..13$) | $\{13,14\}$ |
| 4 | $k\le11$ | $k\le12$ | $-4.18$ ($k=11$), $-1.25$ ($k=12$), $0$ ($k=13..24$) | $\{23,24,25\}$ |

Level 3 sharply amplifies the margins of the already-excluded sectors but extends the reach by at most one sector. The non-excluded range at $N=3$ is $7\le k\le20$: the class-1 BPS states ($3^{13}$ of them) therefore sit at $k\in\{7,10,13\}$ and, by particle–hole symmetry, the class-2 ones at $\{14,17,20\}$; the index alone cannot distinguish these, all of which give the same $0{:}1{:}1$ profile. (The level-3 solves at $N=4$ took 1–3 h each with Clarabel at 5 GB — the AlmostSolved statuses are interior-point tolerance flags, the margins are far from 0 where exclusion is claimed.)

**Assessment.** Exclusion is real and cheap (level 2, seconds per sector at any $N$) and reaches roughly half-way to the window ($6/13.5$ at $N=3$, $12/24$ at $N=4$), well beyond the free sectors $k\le N-1$; but the remaining half needs either higher level (not available) or a different idea. Candidates: (i) few-body / sparse ED at $N=3$ for $k=5..7$ to see whether the gap is closing where exclusion stops; (ii) irrep-resolved exclusion (Casimir rows) — the index says which irreps carry the BPS states, so excluding *those* irreps in a sector suffices and is a much smaller problem; (iii) the cohomology-rank computation for the window itself.

## 4. Irrep-resolved exclusion (2026-09-21/22): the frontier in the $(k,\hat C_2)$ plane

*The quadratic Casimir is a trace polynomial (D7), so $\phi((\hat C_2-c_\lambda)X)=0$ are linear rows of the same kind as the sector rows; with the BPS rows this tests each $(k,\hat C_2)$ cell separately. Data: `results/data/bps_exclusion_irrep_2026-09-21.jsonl`.*

**$N=2$ (ED truth: BPS only at $k=5,6,7$).** Level 2: the Casimir rows sharpen the margins of the already-excluded sectors $k\le2$ (e.g. $k=0$: $-0.96,-1.00,-1.08,-1.23$ for $j=0..3$) but exclude nothing new. Level 3: $k=3$ is excluded in every spin with margins growing with $j$ ($-0.018,-0.091,-0.357$ for $j=1,2,3$), and — new — **$(k=4,\ j=3)$ is excluded** ($-0.033$) although the whole sector $k=4$ is not ($j=0,1,2$: $0$). The reach extends from the high-Casimir end.

**$N=3$, level 2** ($k=7$ is the first non-excludable sector). Per Casimir value at $k=7$: $C_2=48,42,38,36$ excluded ($-0.20,-0.10,-0.036,-0.004$), $C_2\le35$ not; at $k\ge8$ nothing at level 2. Finite-$N$ relations (length $\le6$) do not change any margin.

**$N=3$, level 3** (7 min and $\le7$ GB per cell):

| cell $(k,C_2)$ | class | margin | verdict |
|---|---|---|---|
| $(7,35)$ | 1 | $-0.962$ | excluded |
| $(7,30)$ | 1 | $-0.681$ | excluded |
| $(7,24)$ | 1 | $-0.446$ | excluded (this is $(4,0,-4)$, the largest complex, $540$ multiplets/flavour) |
| $(8,48)$ | 2 | $-0.876$ | excluded |
| $(8,42)$ | 2 | $-0.281$ | excluded |
| $(8,38)$ | 2 | $-0.114$ | excluded |
| $(9,48)$ | 0 | $-0.340$ | excluded (index already zero in class 0) |
| $(10,48)$ | 1 | $-0.050$ | **excluded** |
| $(10,42)$ | 1 | $0$ | not excluded |

So the frontier moves down in $C_2$ at fixed $k$ with level ($k=7$: from $C_2\ge36$ at level 2 to $\le24$ at level 3) and moves out in $k$ at fixed high $C_2$ ($C_2=48$: $k=7$ at level 2, $k=10$ at level 3). Margins are roughly linear in $C_2$ at fixed $k$ ($k=7$: $-0.96,-0.68,-0.45$ at $35,30,24$; $k=8$: $-0.88,-0.28,-0.11$ at $48,42,38$).

**Localisation of the $C_2=48$ complex $(6,0,-6)$** (index $\mp27$ per flavour in classes $1,2$; $81$ multiplets per class in all). Class 1 ($k\equiv1$): $k=1,4$ excluded as whole sectors, $k=7,10$ by the cells above; $k=19,22,25$ are particle–hole images of the excluded class-2 cells $k=8,5,2$. Remaining: $k\in\{13,16\}$. The cell $(11,48)$ (class 2, image $k=16$) decides: if excluded, the class-1 complex sits at $k=13$ and the class-2 complex at $k=14$ — the first BPS multiplets at $N=3$ localised to a single R-charge, with count fixed by the index. Running, together with $(13,48)$ as the consistency check (must come back *not* excluded).

**Confirmation of the deciding $N=3$ cell $(11,48)$.** Clarabel with pruned rows at tolerance $10^{-9}$: $t^*=-0.002556$; Clarabel with all rows and stronger regularisation: $-0.002556$ (six digits agree); SCS (independent first-order solver, $\epsilon=10^{-7}$, stopped at its 40 000-iteration cap, 2.4 h): $-0.002329$ — same sign and magnitude, the 10 % difference being SCS's residual inaccuracy. The controls $(13,48)$ and $(12,48)$ are feasible ($t^*=0$ at tolerance $10^{-9}$), $(12,48)$ being in class $0$ where the index vanishes (empty or cancelling; exclusion cannot tell). With three solves on two solvers agreeing: **the complexes $(c{=}1,(6,0,-6),\omega)$ sit at $k=13$ and $(c{=}2,(6,0,-6),\omega)$ at $k=14$, each with exactly $27$ BPS multiplets (dimension $125$) per flavour charge** — the first BPS multiplets of the three-matrix model at $N=3$ localised to a single R-charge, at the half-filling position the Turiaci–Witten profile predicts. (A fully rigorous version would extract and verify the dual certificate in high precision; the margin is $\sim10^{3}$ times the solver tolerance and reproduced across solvers, so we regard it as established numerically.)

## 5. $N=4$: the top-Casimir complex (2026-09-22)

The maximal-Casimir irrep at $N=4$ is $\lambda=(9,3,-3,-9)$, $C_2=120$, $\dim=117\,649$, self-conjugate, with index $162:{-}81:{-}81$ per flavour in classes $0:1:2$ (the $2{:}1{:}1$ pattern). It occurs only in sectors $k=18..30$ (multiplicities $1,4,22,76,165,264,312,264,\dots$), so the candidate degrees are $\{18,21,24,27,30\}$, $\{19,22,25,28\}$, $\{20,23,26,29\}$ for the three classes, and particle–hole symmetry ($k\to48-k$) maps class 0 to itself and classes 1↔2. Level 3, one constraint set, $\sim5$ min and 6 GB per cell ($N=4$ is *not* more expensive than $N=3$):

| $k$ | class | $t^*$ | verdict |
|---|---|---|---|
| 18 | 0 | $-0.459$ | excluded (hence $30$) |
| 19 | 1 | $-0.101$ | excluded (hence $29$) |
| 20 | 2 | $-0.012$ | excluded (hence $28$) |
| 21 | 0 | $0$ | not excluded |
| 22 | 1 | $0$ | not excluded |
| 23, 24, 25 | 2, 0, 1 | $0$ | feasible (controls, as required) |

Result: the top complex is narrowed to $k\in\{22,25\}$ (class 1), $\{23,26\}$ (class 2) and $\{21,24,27\}$ (class 0); the frontier at $C_2=120$ lies between $k=20$ and $21$, i.e. $k_*-4$ at $N=4$ versus $k_*-2.5$ at $N=3$. Level 3 excludes $k\le20$ in this irrep against $k\le12$ for the whole sector. The last step (one more unit in $k$) is beyond level 3.

**Assessment.** Irrep resolution is the principled extension it promised to be: the exclusion frontier lives in the $(k,C_2)$ plane, is monotone in both variables, and penetrates deepest at maximal Casimir — where the single-matrix intuition (BPS = maximal Casimir) says the "least fortuitous" states live. It completes the localisation of the top complex at $N=3$ and comes within one sector of doing so at $N=4$. The physically weightier complexes (Casimirs $15$–$40$, multiplicities $10^2$–$10^5$) remain partially localised; reaching them needs either level 4 (solver-limited) or a new source of constraints — the cubic Casimir (which separates irreps sharing $C_2$ and is also a trace polynomial) is cheap to add but does not obviously deepen the frontier; the finite-$N$ relations do not help at $N\ge3$.
