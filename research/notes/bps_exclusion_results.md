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
