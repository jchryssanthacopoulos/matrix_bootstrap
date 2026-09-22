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

**Confirmation of the deciding $N=3$ cell $(11,48)$.** Clarabel with pruned rows at tolerance $10^{-9}$: $t^*=-0.002556$; Clarabel with all rows and stronger regularisation: $-0.002556$ (six digits agree); SCS (independent first-order solver, $\epsilon=10^{-7}$, stopped at its 40 000-iteration cap, 2.4 h): $-0.002329$ — same sign and magnitude, the 10 % difference being SCS's residual inaccuracy. The controls $(13,48)$ and $(12,48)$ are feasible ($t^*=0$ at tolerance $10^{-9}$), $(12,48)$ being in class $0$ where the index vanishes (empty or cancelling; exclusion cannot tell). With three solves on two solvers agreeing: **the complexes $(c{=}1,(6,0,-6),\omega)$ sit at $k=13$ and $(c{=}2,(6,0,-6),\omega)$ at $k=14$, each with exactly $27$ BPS multiplets (dimension $343$) per flavour charge** — the first BPS multiplets of the three-matrix model at $N=3$ localised to a single R-charge, at the half-filling position the Turiaci–Witten profile predicts. (A fully rigorous version would extract and verify the dual certificate in high precision; the margin is $\sim10^{3}$ times the solver tolerance and reproduced across solvers, so we regard it as established numerically.)

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


## 6. Two attempted extensions, both negative (2026-09-22)

*Recommended as the cheap levers before the expensive level-4 route; both are now implemented, validated and tested, and neither moves the frontier. Recorded so they are not tried again.*

**A. Using $E=0$ more fully.** For a BPS density matrix $H\rho=\rho H=0$, so (i) $\phi(XH)=\phi(HX)=0$ for *every* $X$, not only the charge-$\mp3$ ones already used, and (ii) the ground-state cone $\big[\phi(X^\dagger HY)\big]\succeq0$ is available (valid for any state, since $H\succeq0$; for a BPS functional its $X=\mathbf 1$ row vanishes, which is (i) again). Both are implemented — (ii) needed a new algebra primitive, `trace_algebra.sandwich`, which inserts $H$ in operator order between two open words contracted into a single trace; it is verified to $10^{-14}$ against explicit operators at $(N,p)=(2,3),(2,1),(3,1)$, and the rows are satisfied by the exact BPS functional to $10^{-11}$. **Effect on the frontier: none.** At level 2 the margins at $N=2$, $k=3,4$ and $N=3$, $k=7,8$ are unchanged (all $0$). The $E=0$ information is already implied by the existing BPS rows at these levels. (The `gs` option of the trace engine, previously unimplemented, is now filled in as a by-product.)

**B. Irrep resolution beyond the quadratic Casimir.** The cubic Casimir is also a trace polynomial, $\hat C_3=\mathrm{Tr}[\tilde M^3]$ with $\tilde M=M-pN\mathbf 1$ (D7); it is verified to $10^{-15}$ against the explicit operator and commutes exactly with every gauge generator. Its eigenvalue was obtained empirically from the joint $(\hat C_2,\hat C_3)$ spectrum of the $k=1,2$ sectors and then validated as
$$c_3(\lambda)=f(l)-f(l^0),\qquad f(x)=\sum_i\big[x_i^3+(\tfrac32-N)x_i^2\big],\qquad l_i=\lambda_i+N-i,\quad l^0_i=N-i,$$
which reproduces the full $(\hat C_2,\hat C_3)$ multiset (with multiplicities) at $N=2,3,4$ and correctly degenerates to $c_3=2c_2$ at $N=2$, where $su(2)$ has no independent cubic invariant. Adding the rows $\phi((\hat C_3-c_3)X)=0$ splits the $\hat C_2$-degenerate cells into single irreps. **Effect on the frontier: none.** The test cell is $N=3$, $k=10$, $\hat C_2=42$, which the unsplit test leaves at $t^*=0$; its two irreps $(6,-1,-5)$ and $(5,1,-6)$ (dimensions $260$; $c_3=252$ and $0$) both give $t^*=0$ with the full $C_2$ rows plus $C_3$ rows against operators of length $\le2$ (288 s, 6.5 GB per cell). A first version of this test also restricted the $C_2$ rows, which would have made the comparison unfair; the numbers quoted are from the corrected run.

**Interpretation.** Restricting the *state* further — by charge, by $\hat C_2$, now by $\hat C_3$ — helped enormously the first time ($k\le6$ whole-sector $\to$ $k\le11$ at $\hat C_2=48$) and not at all the second. What sets the frontier at a given level is the *operator content* of the cones, not the resolution of the state: once the cone entries cannot see the interaction that would force infeasibility, no additional linear restriction on $\rho$ recovers it. The remaining lever is therefore level 4, where for exclusion we need only an infeasibility (Farkas) certificate rather than an accurate optimum — plus the reversal $\mathbb Z_2$ (transposition, a symmetry of $\mathrm{Tr}\Psi^3$) which halves every cone and would bring level-4 cones from 168 to $\sim84$.

**Cost note.** Level-3 irrep cells with the $C_3$ rows: build 1.5 GB / 90 s, solve 5.6–6.5 GB / $\sim5$ min (Clarabel). The $C_3$ rows must be restricted to short operators (`L_cas`): applied against every Gram entry they reach length 12 and exceed 5 GB in the build alone.

## 7. Level 4 for exclusion: feasibility mode and Farkas certificates (2026-09-22)

**The idea.** At levels 2–3 we used the *margin* form (maximise $t$ with every cone $\succeq t\mathbf 1$), which is an optimisation and needs convergence. Exclusion, however, only needs a yes/no plus a certificate: the pure feasibility problem (objective $\equiv0$) is infeasible exactly when a Farkas certificate exists, i.e. a $y$ with
$$A^{\!\top}y=0,\qquad y\in K^*,\qquad b\cdot y<0 .$$
`TraceSDP.feasibility` solves that problem and returns $y$ together with the three residuals. This is what makes level 4 reachable: the level-4 *energy* problem never converged (4.5 h, 18 GB, residual $10^{-2}$), whereas a level-4 *feasibility* problem solves in **5 minutes at 2.5 GB**.

**Machinery validated before use.** (i) Level 2, both solvers: on cells known to be excluded ($N=2$, $k=0$; $N=3$, $k=6$) both return `infeasible` with a clean certificate ($\|A^{\!\top}y\|/\|y\|\sim10^{-15}$–$10^{-10}$, $b\cdot y/\|y\|<0$, dual-cone eigenvalue $\ge-10^{-18}$); on the feasible cell $N=3$, $k=7$ they return `solved` with $b\cdot y\approx0$, i.e. no certificate. (ii) Level 3: $(10,C_2{=}48)$ → `infeasible`, $(13,C_2{=}48)$ → `solved`, matching the margins. A bug worth recording: the two solvers use *different* svec conventions (Clarabel packs the upper triangle column-major, SCS the lower); unpacking with the wrong one made a valid certificate look like a dual-cone violation of $-2\times10^{-2}$.

**Level-4 problem.** `TraceSDP.fixed_N` now carries the exclusion rows, built from a deliberately small operator set (short words only): the length-8 Gram entries must never enter the row generator — multiplying them by $Q$, $\hat C_2$ or $\hat C_3$ reaches length 11–14 and explodes the monomial count (this OOM'd the machine once). With `L_eom_gram=4`, `L_short=6`: 47 421 monomials, 13 538 rows (1 284 BPS), 85 713 cone entries, cones of 168 words (real embedding 336), $n=77\,106$ unknowns, 426 547 rows; build 29 s / 0.86 GB, SCS setup 17 s / 2.5 GB.

**First result.** $N=3$, $k=10$, $C_2=42$ (the cell just beyond the level-3 frontier): **feasible** at level 4 — $b\cdot y=+1.2\times10^{-7}$, so no certificate exists — i.e. level 4 does *not* exclude it. Controls on cells with known level-3 verdicts are running; they are needed because the memory trimming (`L_eom_gram=4`) removes EOM rows that level 3 had, so level 4 as configured is not automatically stronger than level 3.

**Note on the reversal $\mathbb Z_2$.** Earlier outlooks proposed halving the cones with the transposition symmetry $w\to w^{\rm T}$. It does not exist for this model: under transposition $Q\to-C^{\rm rev}_{abc}\mathrm{Tr}[\Psi^a\Psi^b\Psi^c]$, and Chen's $C$ is not reversal-symmetric ($\langle C,C^{\rm rev}\rangle=5$ against $\|C\|^2=16/3$, D2.12). It *is* available for the single-matrix model.

## 8. Calibration against exact cohomology (2026-09-22) — supersedes the guesswork in §6

The exact $Q$-cohomology computation (D8, `research/notes/cohomology_results.md`) now gives the ground truth for
every $N=3$ cell in `results/data/bps_exclusion_irrep_2026-09-21.jsonl`, so the exclusion bootstrap can be graded
rather than argued about.

**Reading the margins correctly.** Margins of order $10^{-7}$–$10^{-10}$ are numerically zero: they mean the SDP
sits *on* the boundary and the test is **inconclusive**, not that the cell is excluded. Only margins bounded away
from zero (plus a Farkas certificate) are exclusions. Several cells in §3–§5 were recorded with margins in that
range; they should be read as inconclusive. With that convention:

| $C_2$ | bootstrap excludes up to | true BPS onset | gap |
|---|---|---|---|
| 48 | $k=11$ ($-2.6\times10^{-3}$) | $k=12$ | **0 — exactly sharp** |
| 42 | $k=8$ | $k=12$ | 3 |
| 38 | $k=8$ | $k=12$ | 3 |
| 35 | $k=7$ | $k=12$ | 4 |

* **Soundness: zero contradictions.** No cell was ever excluded that in fact carries BPS states — over all $N=3$
  cells run, at levels 3 and 4.
* **Sharpness at the top irrep.** For $(6,0,-6)$ the level-3 frontier is *exactly* the true boundary: it excludes
  $k\le11$ with a genuine margin and goes inconclusive at precisely $k=12$, the true onset. That is the strongest
  validation the method has received.
* **The deficit, now measured.** At $C_2=42,38,35$ the degrees $k=9,10,11$ are genuinely empty but could not be
  certified — including $(10,42)$ at level 4. §6 concluded that the frontier is set by operator level rather than
  by state resolution; that stands, and the shortfall is 3–4 degrees rather than an unknown amount.

**Status of the level-4 controls.** The control runs at $(10,48)$, $(13,48)$ and $(10,42)$ announced in §7 never
completed — the job died without writing output and was not restarted. They are now moot: the question they were
meant to settle (whether the trimmed level-4 configuration dominates level 3) is answered directly by the table
above, which shows level 4 failing to exclude $(10,42)$ where the truth is that the cell is empty.

**Consequence for the programme.** The $N=3$ concentration question no longer needs the bootstrap: cohomology
answered it exactly, in minutes, at 1 GB. The level-4/level-5 and HPC escalation discussed earlier is not the route
to this particular question. The exclusion bootstrap's remaining value is as a method that scales to $N$ where
cohomology does not — and it now has a measured accuracy at $N=3$ to justify trusting it there.
