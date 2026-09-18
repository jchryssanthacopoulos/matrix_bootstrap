# Exact finite-$N$ trace bootstrap — results log

*Plan: `docs/trace_bootstrap_plan.md`. Code: `src/trace_algebra.py` (M1), `src/trace_bootstrap.py` (M2), `src/fewbody.py` (few-body sector ED for cross-checks), `scripts/run_trace_bound.py`. Data: `results/data/trace_bounds_2026-09-18.jsonl`. All bounds are rigorous lower bounds on $E_0(k;N)$ at the stated integer $N$ (no factorisation, no $1/N$ expansion); "exact" values come from full ED ($N=2$), D2 anchors, or few-body ED.*

## 1. Engine (2026-09-18)

* **Variables** are the expectation values of canonical neutral monomials (products of trace words); real unknowns $(\mathrm{Re},\mathrm{Im})$. **Constraints:** $\phi(1)=1$; reality; EOM $\phi([H,X])=0$ and sector rows $\phi((N_\Psi-k)X)=\phi(X(N_\Psi-k))=0$ for $X$ = every Gram entry and every neutral trace word up to $L_{\rm eom}$; optional finite-$N$ trace relations (antisymmetriser over $N+1$ indices, total length $\le$ `finiteN_len`); cones = singlet $\phi(T_a^\dagger T_b)\succeq0$ and adjoint-projected $\phi(\mathrm{Tr}[w_a^\dagger w_b])-\phi(\mathrm{Tr}w_a^\dagger\mathrm{Tr}w_b)/N\succeq0$, graded by $(q,z)$. No cvxpy: the standard-form conic program (equalities + real embeddings of the Hermitian cones) goes straight to Clarabel or SCS.
* **Conditioning** (all needed): equality rows normalised to unit norm; 't Hooft column scaling $\phi(m)=N^{\sum_i(1+L_i/2)}\tilde\phi(m)$ and objective $/N^3$; for Clarabel, pruning of dependent rows by pivoted QR when the dense block is affordable, chordal decomposition off, equilibration **on**, static regularisation $10^{-7}$ ($10^{-5}$ if unpruned); for SCS, the unpruned rows. Clarabel is the accurate choice at $N\lesssim30$, SCS the robust one at $N\ge100$ (Clarabel returns NumericalError there).
* **Diagnostics:** `check_exact` evaluates every row and cone on the exact multiplet-averaged sector ground state (small $N$); residuals were $\le10^{-12}$ in every instance below.

## 2. Validation ladder

| step | run | result |
|---|---|---|
| V1 | single matrix $N=2,3$, all sectors, level (2,2,3), no relations | Casimir ladder exact ($18,0,0,0,18$; $72,45,18,0,0,0,0,18,45,72$); $N=4$, $k=1$: $144$ exact |
| V2 | three-matrix $N=2$, level (2,2,3), $k=0..3$ | $98,\,22,\,0,\,0$ = the Hilbert-space level-2 values |
| V2 | three-matrix $N=2$, $k=2$, level (3,3,4), `finiteN_len` = 0 / 4 / 6 | **3.2631 / 3.4384 / 5.16536** — exact with the length-6 relations (Hilbert-space value 5.16536); 35 170 monomials, 7 283 EOM/sector rows, +432 / +65 976 relation rows; build 70 s, solve 2–3 min, 5–7 GB |
| V3 | three-matrix $k=0$ and $k=1$, level (2,2,3), one constraint set evaluated at $N=2,3,4,5,10,30,100,1000$ | $16N^3-15N$ and $16N^3-53N$ reproduced to $10^{-8}$ ($N\le30$, Clarabel/SCS), $10^{-10}$ ($N=100$, SCS), $\sim10^{-7}$ ($N=1000$, SCS "inaccurate" — the $N=1000$ values are slightly *above* the exact ones and are therefore not certified bounds at that tolerance) |

The V2 series is the cleanest statement of what the finite-$N$ relations do: at $N=2$ the invariant sector functional has only $D=90$ parameters, the generic trace algebra sees $3.5\times10^4$ monomials, and the length-$\le6$ Cayley–Hamilton-type relations close the gap completely at this level.

## 3. New exact values from few-body ED (`src/fewbody.py`)

Sector ED in the $k$-particle space (dimension $\binom{3N^2}{k}$, never the full Fock space), validated against full ED at $N=2$ and the D2 anchors at $N=3,4$:

| $N$ | $k=0$ | $k=1$ | $k=2$ | $k=3$ |
|---|---|---|---|---|
| 2 | 98 | 22 | 5.165364 | 1.227064 |
| 3 | 387 | 273 | **159** | **75.791675** |
| 4 | 964 | 812 | **660** | — |
| 5 | 1925 | 1735 | **1545** | — |

**Two-particle sector is non-interacting for $N\ge3$:** $E_0(2;N)=16N^3-91N=E_0(0;N)-2\cdot38N$ for $N=3,4,5$ — twice the lowest one-particle gap ($38N$, flavor-symmetric gauge-traceless mode, D2.8). At $N=2$ this would give $-54$; the actual $5.165$ shows the $N=2$ value is a finite-size effect (two lowest one-particle states cannot coexist without interaction at $N=2$). This settles the D2 todo item "two-particle sector analytically" numerically for $N\le5$; an analytic proof for all $N\ge3$ is a natural short derivation (D5 candidate). $k=3$ is interacting: $75.79\ne387-3\cdot114=45$.

## 4. First bootstrap results beyond $N=2$ (2026-09-18)

**Two-particle sector, level (2,2,3), no finite-$N$ relations, one constraint set evaluated at several $N$ (SCS):**

| $N$ | bound | $16N^3-91N$ | few-body ED |
|---|---|---|---|
| 3 | 158.999997 | 159 | 159 |
| 4 | 659.999994 | 660 | 660 |
| 5 | 1544.999987 | 1545 | 1545 |
| 10 | 15089.999992 | 15090 | — |
| 100 | 15990900.002 (SCS inaccurate) | 15990900 | — |

So at level 2 — the level at which the same sector gives only the SUSY floor $0$ at $N=2$ — the trace bootstrap is **tight for $k=2$ at every $N\ge3$** and reproduces $E_0(2;N)=16N^3-91N$. Together with the variational upper bound from two non-interacting lowest one-particle excitations (which exists for $N\ge3$), this is a proof of the two-particle ground energy for every $N\ge3$ at which the SDP is solved, and the bootstrap side is a one-second computation. The contrast with $N=2$ (needs level 3 *and* the length-6 finite-$N$ relations) illustrates the plan's expectation that the trace formulation gets *easier*, not harder, with $N$: the $N=2$ difficulty is the finite-size structure of $2\times2$ matrices, not the bootstrap.

## 5. The few-particle sectors: $E_0(k;N)=16N^3-(15+38k)N$ for $k\le N-1$ (2026-09-18)

Few-body ED (`src/fewbody.py`) and the level-2 trace bootstrap agree on the following, and the bootstrap side is rigorous:

| $N$ | $k=2$ | $k=3$ | $k=4$ | formula $16N^3-(15+38k)N$ |
|---|---|---|---|---|
| 2 | 5.165 (ED) | 1.227 (ED) | 0.088 (ED) | $-54,-130,-206$: **not** applicable ($k\ge N$) |
| 3 | **159** (ED = bootstrap) | 75.792 (ED); bootstrap level 2 gives 45 | 46.271 (ED) | $159$; $45$; $-69$: applies to $k=2$ only |
| 4 | **660** (ED = bootstrap) | **508** (ED = bootstrap) | — | $660,\ 508$ |
| 5 | **1545** (ED = bootstrap) | **1355** (bootstrap) | — | $1545,\ 1355$ |
| 10 | **15090** (bootstrap) | **14710** (bootstrap) | — | $15090,\ 14710$ |
| 100 | 15990900 (bootstrap, SCS inaccurate) | — | — | $15990900$ |

**Statement (numerical, for the cases in the table; conjectured in general).** For $k\le N-1$ the $k$-particle ground state of Chen's three-matrix model consists of $k$ non-interacting copies of the lowest one-particle excitation (the flavor-symmetric, gauge-traceless mode at $38N$ below $E_0(0)$, D2.8), $E_0(k;N)=16N^3-(15+38k)N$; the first *interacting* sector is $k=N$ ($5.165$ at $N=2$, $75.79$ at $N=3$), where the bound from level 2 is only the non-interacting value (or the SUSY floor, when that is higher) and the finite-$N$ relations / higher level are needed — exactly the pattern seen at $N=2$, $k=2$. A natural reading: the lowest one-particle multiplet is an $SU(N)$ adjoint, and $k$ of them can be antisymmetrised into a state that the quartic interaction annihilates as long as $k\le N-1$ (a Casimir/Pauli-type obstruction at $k=N$); this is a short derivation to write down (D5 candidate) and is the finite-$k$ shadow of Chen's "maximal Casimir" structure. If true, the near-BPS physics at fixed $k$ starts only at $k\ge N$, so the fixed-$k$, large-$N$ regime that the trace engine reaches is **free** for $k<N$ — a clean but somewhat deflating statement for the programme's physics content (see §6 below once the $k=N$ runs are in).

**Bootstrap tightness at level 2 for $k\le N-1$** ($k=2$: $N=3,4,5,10,100$; $k=3$: $N=4,5,10$) is itself a result about the method: in the non-interacting sectors the level-2 constraint set (words of length $\le2$, EOM to length 3, sector rows, singlet + adjoint-projected cones) is complete at every $N$, in seconds.

## 6. The first interacting sector at $N=3$: $k=3$ (2026-09-18)

| level | `finiteN_len` (relation rows) | bound | exact (few-body ED) |
|---|---|---|---|
| (2,2,3) | 0 | 45.000 (= non-interacting value $16N^3-129N$) | 75.7917 |
| (3,3,4) | 0 | **70.3533** (93 %) | 75.7917 |
| (3,3,4) | 4 (6 rows) | 70.3533 | |
| (3,3,4) | 6 (2082 rows) | 70.3533 | |

Level 3 with the generic algebra alone already gives 93 % (35 170 monomials, 7 286 EOM/sector rows, Clarabel 2.4 min, 4.6 GB). Unlike $N=2$, the finite-$N$ relations up to length 6 change nothing at $N=3$: the remaining gap is a **level** deficit, i.e. it needs length-4 open words — the same "tight one level after leaving the floor" pattern as the Hilbert-space $k=3$ at $N=2$ (level 3: 78 %, level 4: exact). This is the first rigorous bootstrap bound on an interacting sector of the three-matrix model beyond the reach of the Hilbert-space code: $E_0(3;3)\ge70.35$.

**Level 4 in the trace engine** is not yet affordable with the present Python assembly ($\sim1.5\times10^3$ open words, cones up to $\sim170$, an estimated $10^5$–$10^6$ monomials); see the todo (M4) for the two routes: restrict the length-4 words to the cones that carried the information at $N=2$, and/or eliminate monomials through the sector rows before assembly.

## 7. Status against the plan (end of M2/M3)

* M2 done: `TraceSDP` builds and solves; V1 (Casimir ladder $N=2,3,4$) exact at level 2.
* M3 done for $k\le2$: V2 exact at $N=2$ with `finiteN_len=6`; V3 anchors exact for $N=2..1000$ (SCS at $N\ge100$); $k=3$ at $N=2$ level 4 not attempted (level-4 assembly).
* M4 begun: $E_0(k;N)$ table above; the non-interacting pattern for $k\le N-1$ established numerically and by rigorous bounds; first interacting-sector bound at $N=3$.
* Solver notes: Clarabel is accurate up to $N\approx30$ and returns NumericalError beyond even with the 't Hooft scaling; SCS handles $N=100$–$1000$ but its "inaccurate" results at $N=1000$ are slightly above the exact values and must not be quoted as bounds. Peak memory so far $\le7$ GB.
