# Exact $Q$-cohomology of the three-matrix model at $N=3$ (2026-09-22)

Status: **computed exactly**, not bootstrapped and not conjectured. Every number below is an exact integer
(rank over $\mathbb{F}_p$, cross-checked against a second prime and floating point), and every one is validated
against the independently computed refined Witten index. Code: `src/cohomology.py`; data:
`results/data/cohomology_N3_top.json`, `results/data/cohomology_N3_flavour.json`; derivation: D8 in `docs/derivations.md`.

## 1. What made this tractable

The obstruction was always size: at $N=3$, $p=3$ the Fock space is $2^{27}$ and the middle sector has
$\dim\mathcal H_{13}=2\times10^7$, far beyond a rank computation. The observation that removes the problem is that
$Q=\sum_{abc}C_{abc}\,\mathrm{Tr}[\Psi^a\Psi^b\Psi^c]$ is a $U(N)$ **singlet**, so it commutes with the Cartan and
preserves the $U(N)$ weight. The complex therefore splits over weights,

$$\cdots \xrightarrow{\;Q\;} W_\lambda(k) \xrightarrow{\;Q\;} W_\lambda(k+3) \xrightarrow{\;Q\;}\cdots ,$$

and the weight spaces are tiny: for $\lambda=(6,0,-6)$ the largest is $\dim W=126$ against $2\times10^7$ for the
full sector. The whole computation for the eight highest weights runs in **under two minutes at 1.0 GB peak**.

Turning weight-space cohomology into per-irrep BPS counts is a triangular inversion,
$\dim H^k(W_\lambda)=\sum_{\mu\ge\lambda}K(\mu,\lambda)\,h^k_\mu$, with $K$ the Kostant weight multiplicity. Higher
weights have *smaller* weight spaces, so the set of weights that is affordable is automatically closed upward and
the peel is self-contained — no truncation is involved.

## 2. Result: BPS states of the eight highest irreps at $N=3$

These are exactly the eight irreps of the $N=3$ Fock space with $C_2\ge35$ (verified exhaustively against
`bps_index.irrep_multiplicities`), so the peel is complete for all of them. $h^k_\lambda$ = number of BPS
multiplets of irrep $\lambda$ at degree $k=N_\Psi$, summed over the three $\mathbb Z_3$ flavour charges.

| irrep $\lambda$ | $C_2$ | $\dim\lambda$ | $h^{12}$ | $h^{13}$ | $h^{14}$ | $h^{15}$ | total |
|---|---|---|---|---|---|---|---|
| $(6,0,-6)$   | 48 | 343 |  27 |  81 |  81 |  27 |  216 |
| $(6,-1,-5)$  | 42 | 260 |  54 | 162 | 162 |  54 |  432 |
| $(5,1,-6)$   | 42 | 260 |  54 | 162 | 162 |  54 |  432 |
| $(6,-2,-4)$  | 38 | 162 |  81 | 243 | 243 |  81 |  648 |
| $(4,2,-6)$   | 38 | 162 |  81 | 243 | 243 |  81 |  648 |
| $(6,-3,-3)$  | 36 |  55 |  27 |  81 |  81 |  27 |  216 |
| $(3,3,-6)$   | 36 |  55 |  27 |  81 |  81 |  27 |  216 |
| $(5,0,-5)$   | 35 | 216 |  38 | 486 | 486 |  38 | 1048 |

**The cohomology is supported on exactly four consecutive degrees, $k\in\{12,13,14,15\}$, and vanishes identically
elsewhere** — in particular $h^k=0$ for all $k\le11$ and all $k\ge16$, for every one of these irreps.

Two structural facts:

* $h^{12}=h^{15}$ and $h^{13}=h^{14}$ in every row. This is the particle–hole involution $k\mapsto pN^2-k=27-k$,
  which maps $12\leftrightarrow15$ and $13\leftrightarrow14$; it is a consistency check, not an input.
* The ratio $1:3:3:1$ holds for the maximal-Casimir family ($C_2=48,42,38,36$) but **fails** for $(5,0,-5)$, where
  it is $38:486$. The $38$ was re-verified with a second prime and in floating point; it is not a rounding artefact.
  So the $1:3:3:1$ shape is a property of the top of the Casimir spectrum, not a general law.

## 3. The concentration question, stated per complex

This is the point that the raw table hides. The complexes are labelled by $c=k\bmod 3$ (since $[N_\Psi,Q]=3Q$),
together with the irrep $\lambda$ and the flavour charge $w$. The support $\{12,13,14,15\}$ distributes as

| complex | degrees carrying cohomology | index | verdict |
|---|---|---|---|
| $c=1$ | $k=13$ only | $-h^{13}$ | **index-saturated: concentrated in a single degree** |
| $c=2$ | $k=14$ only | $+h^{14}$ | **index-saturated: concentrated in a single degree** |
| $c=0$ | $k=12$ **and** $k=15$ | $h^{12}-h^{15}=0$ | not saturated; index-invisible |

Flavour-refined, this was checked for all $3\times3\times3=27$ complexes $(\lambda,w,c)$ of the top three irreps:
every $\chi$ equals the refined index exactly (`cohomology_N3_flavour.json`). For $(6,0,-6)$ each flavour charge
carries $h^{12},h^{13},h^{14},h^{15}=9,27,27,9$.

So at $N=3$:

* **Where the index is macroscopic ($c=1,2$), concentration holds exactly** — the cohomology sits in one degree and
  the refined index counts the BPS states on the nose. This is index saturation in the sense of
  Chang–Chen–Sia–Yang §5.1.
* **Where the index vanishes ($c=0$), there are nonetheless BPS states** — $149{,}526$ of them across these eight
  irreps, sitting in two degrees ($k=12$ and $k=15$) in perfectly cancelling pairs. They are invisible to the index.

This is **not** a counterexample to the CCSY concentration conjecture, which is conditioned on a macroscopic index;
it is a quantitative statement about how much the index misses. Across the eight irreps, $758{,}808$ BPS states
in total, of which $\mathbf{19.7\%}$ are index-invisible. Any argument that uses the index as a proxy for the BPS
count at finite $N$ is wrong by that much here.

## 4. Calibration of the exclusion bootstrap against exact truth

Every $N=3$ cell we ran in `results/data/bps_exclusion_irrep_2026-09-21.jsonl` can now be graded. Margins of
order $10^{-7}$–$10^{-10}$ are numerically zero and are counted as **inconclusive**, not as exclusions (a genuine
exclusion needs a margin bounded away from zero plus a Farkas certificate).

| $C_2$ | bootstrap excludes up to | true BPS onset | gap |
|---|---|---|---|
| 48 | $k=11$ (margin $-2.6\times10^{-3}$) | $k=12$ | **0 — exactly sharp** |
| 42 | $k=8$  | $k=12$ | 3 |
| 38 | $k=8$  | $k=12$ | 3 |
| 35 | $k=7$  | $k=12$ | 4 |

**Zero contradictions**: no cell was ever excluded that actually carries BPS states. At the top irrep the level-3
exclusion frontier is *exactly* the true boundary — it excludes $k\le11$ and goes inconclusive at precisely $k=12$,
where the BPS states really begin. That is a strong independent validation of the exclusion machinery.

It also quantifies the machinery's weakness: at $C_2=42,38,35$ the cells $k=9,10,11$ are genuinely empty but the
level-3 (and, at $(10,42)$, level-4) bootstrap could not certify it. The earlier reading that the frontier is set
by operator level rather than state resolution (`bps_exclusion_results.md` §6) is confirmed, and we now know the
exact size of the deficit rather than guessing at it.

## 4b. $N=2$, and a conjecture on the window

The same machinery at $N=2$ ($p=3$, Fock space $2^{12}$) is exact and instantaneous, and reproduces the ED
multiplet counts:

| irrep (spin $j$) | $\dim$ | $h^5$ | $h^6$ | $h^7$ |
|---|---|---|---|---|
| $(3,-3)$, $j=3$ | 7 |  9 | 18 |  9 |
| $(2,-2)$, $j=2$ | 5 | 18 | 36 | 18 |
| $(1,-1)$, $j=1$ | 3 | 27 | 54 | 27 |
| $(0,0)$,  $j=0$ | 1 |  9 | 18 |  9 |

BPS states per degree: $243:486:243$. This $1{:}2{:}1$ was already known at $N=2$ from ED and the index
(`bps_index_results.md`, Q1d); what is new is that the cohomology **reproduces it independently**, and shows it
holds irrep by irrep, not merely in the total.

Setting the two $N$ side by side:

| $N$ | window | width | ratios (top irreps) | centre |
|---|---|---|---|---|
| 2 | $k=5,6,7$ | $3$ | $1:2:1$ | $6=pN^2/2$ |
| 3 | $k=12,13,14,15$ | $4$ | $1:3:3:1$ | $13.5=pN^2/2$ |

**Conjecture (stated as such, not proven):** the BPS window is centred on half filling $k=pN^2/2$, has width
$N+1$, and the degree profile goes as $\binom{N}{j}$. *Outcome (§4c, §4d): the **width** half has held for every
irrep computed at $N=2,3,4$ — 25 irreps in total, no exceptions. The **profile** half is false in general: it
breaks at $(5,0,-5)$ for $N=3$ and at the two $C_2=107$ irreps for $N=4$, in both cases by suppressing the window
edges while leaving the centre intact.*

The caveat is sharp and already visible: at $N=3$ the $\binom{3}{j}=1{:}3{:}3{:}1$ profile holds for the seven
highest-Casimir irreps but **fails** for $(5,0,-5)$, where the profile is $38:486:486:38$. At $N=2$ no such failure
occurs — all four irreps obey $1{:}2{:}1$. So either the profile is a property of the top of the Casimir spectrum
that happens to be exhaustive at $N=2$, or the total over *all* irreps obeys it while individual irreps need not.
Distinguishing these requires the low-Casimir $N=3$ weights, which are not yet computed (§6). The width-$(N+1)$
part of the conjecture is on firmer ground than the profile part; $N=4$ (§4c) confirmed both for the maximal
weight.

## 4c. $N=4$: the method reaches the maximal weight (2026-09-23)

**It works.** The maximal weight $\lambda=(9,3,-3,-9)$ ($C_2=120$, $\dim\lambda=117{,}649$) runs in **6 s at
0.33 GB** — against $\dim\mathcal H_{24}=3.2\times10^{13}$ for the full sector. And because nothing sits above a
maximal weight, **no peeling is needed**: $h^k=\dim H^k(W_\lambda)$ directly.

One change was required. The original `weight_basis` enumerated occupation matrices with
`itertools.product(range(p+1), repeat=N*N)` $=4^{16}\approx4.3\times10^9$ at $N=4$ — unusable regardless of how
small the weight space is. It is now a pruned row-by-row recursion (`cohomology.occupation_matrices`), validated
exhaustively: the $N=3$ weight spaces still partition every sector ($\sum_\lambda\dim W_\lambda(k)=\binom{27}{k}$
for all $k$), and at $N=4$ the counts match the D6 generating function for every $k$ and the four highest weights.
All previously reported $N=2$ and $N=3$ numbers are unchanged under the swap (test suite re-run).

| $k$ | 22 | 23 | 24 | 25 | 26 |
|---|---|---|---|---|---|
| $h^k$ (all charges) | 81 | 324 | 486 | 324 | 81 |
| per flavour charge   | 27 | 108 | 162 | 108 | 27 |

$81:324:486:324:81 = 81\times\binom{4}{j}$: at the maximal weight, both halves of the §4b conjecture hold at
$N=4$ — width $5=N+1$, centred on $k=24=pN^2/2$, profile $\binom{N}{j}$. *Extending to twelve further irreps
(§4d) then kept the width but broke the profile at $C_2=107$, so read the profile claim as holding at the top of
the spectrum only.*

### The important part: saturation fails differently at $N=4$

| $c=k\bmod3$ | degrees | multiplets (per charge) | $|$index$|$ | verdict |
|---|---|---|---|---|
| 0 | $\{24\}$      | 162 | 162 | **saturated** |
| 1 | $\{22,25\}$   | 135 |  81 | **not saturated, and the index does not vanish** |
| 2 | $\{23,26\}$   | 135 |  81 | **not saturated, and the index does not vanish** |

Checked for all nine $(c,w)$ complexes; every $\chi$ equals the refined index.

At $N=3$ the failures of saturation were confined to the $c=0$ complexes, which have *zero* index — so one could
say the index was reliable wherever it was non-zero. **That escape closes at $N=4$.** Here the $c=1$ and $c=2$
complexes carry a non-vanishing index ($-81$ per charge) and still have cohomology in two degrees, $135$ multiplets
against $|I|=81$. Across the whole irrep: $1296$ BPS multiplets ($152{,}473{,}104$ states) against $972$ counted by
the index — a $25\%$ excess.

So index saturation is **not** a reliable proxy at finite $N$ even for complexes with macroscopic index. This is
not a refutation of CCSY, whose conjecture is asymptotic in $N$ and whose $O(1)$ corrections are invisible at
$N=4$; but any *finite-$N$* argument that reads BPS counts off the index is wrong here by a quarter, and the error
does not shrink between $N=3$ ($19.7\%$, top eight irreps) and $N=4$ ($25\%$, top irrep).

### Reach at $N=4$

Weight-space sizes for the highest weights: $(9,3,-3,-9)$ → 924; $(9,3,-4,-8)$, $(8,4,-3,-9)$, $(9,2,-2,-9)$ →
4752; then 11286, 14696, 25542, 30294. Only the first is comfortable for dense elimination. Going past the single
maximal irrep at $N=4$ therefore needs a sparse or blocked rank routine — the ranks are the wall now, not the
enumeration.

## 4d. $N=4$ pushed to thirteen irreps with a sparse/blocked rank routine (2026-09-23)

The ranks were the wall; two changes moved it by a factor of ~36 000 in raw $n^3$ work (largest weight space
$924\to30\,294$), and the whole run below took **~30 minutes at 6.1 GB peak**.

**1. Blocked exact rank (`rank_mod_p_blocked`).** Right-looking LU with partial pivoting, skipping rank-deficient
columns, where the trailing-submatrix update is one `L21 @ U12` matmul that numpy hands to BLAS. It stays *exact*
by working modulo a prime small enough that float64 products are exact integers:
$\text{block}\times(P-1)^2<2^{53}$ with $P\approx2^{20}$. (Two bugs on the way, both caught by the reference
routine: the panel initially normalised the pivot row across the panel only, leaving $U$ inconsistent with the unit
diagonal of $L_{11}$; and normalising before rather than after the elimination, which do not commute. The final
version does standard LU with no rescaling at all.)

**2. $\mathbb Z_3$ flavour blocking (`complex_cohomology_blocked`).** $Q$ commutes with the flavour rotation, so
each weight space splits into three blocks of $\approx n/3$ — a $9\times$ saving in time and memory, and it
resolves the flavour charge for free. Done exactly over $\mathbb F_P$ by choosing $P\equiv1\pmod 3$ so that a
primitive cube root of unity exists; the eigenbasis is built from the orbits of the rotation (size 1 or 3, with
fermionic signs), giving
$(Q_w)_{O',O}=\sum_r a_{u_r}d_r\,\omega^{wr}$ for $Q|s_0\rangle=\sum_u a_u|u\rangle$.
Validated at $N=2,3,4$: totals and per-charge counts match the dense computation, and the block ranks sum to the
dense rank for every $k$.

### Result: all thirteen $N=4$ irreps with $C_2\ge107$

| irrep | $C_2$ | $\dim$ | $h^{22},\dots,h^{26}$ | profile |
|---|---|---|---|---|
| $(9,3,-3,-9)$ | 120 | 117649 | 81, 324, 486, 324, 81 | $81\binom4j$ |
| $(9,3,-4,-8)$ | 114 | 91000 | 162, 648, 972, 648, 162 | $2\cdot81\binom4j$ |
| $(8,4,-3,-9)$ | 114 | 91000 | 162, 648, 972, 648, 162 | $2\cdot81\binom4j$ |
| $(9,2,-2,-9)$ | 114 | 94640 | 162, 648, 972, 648, 162 | $2\cdot81\binom4j$ |
| $(9,3,-5,-7)$ | 110 | 57456 | 243, 972, 1458, 972, 243 | $3\cdot81\binom4j$ |
| $(9,1,-1,-9)$ | 110 | 61236 | 243, 972, 1458, 972, 243 | $3\cdot81\binom4j$ |
| $(7,5,-3,-9)$ | 110 | 57456 | 243, 972, 1458, 972, 243 | $3\cdot81\binom4j$ |
| $(8,4,-4,-8)$ | 108 | 69825 | 324, 1296, 1944, 1296, 324 | $4\cdot81\binom4j$ |
| $(9,3,-6,-6)$ | 108 | 19635 | 81, 324, 486, 324, 81 | $81\binom4j$ |
| $(9,0,0,-9)$ | 108 | 21175 | 81, 324, 486, 324, 81 | $81\binom4j$ |
| $(6,6,-3,-9)$ | 108 | 19635 | 81, 324, 486, 324, 81 | $81\binom4j$ |
| $(8,3,-2,-9)$ | 107 | 80640 | **114, 1572, 2916, 1572, 114** | **not binomial** |
| $(9,2,-3,-8)$ | 107 | 80640 | **114, 1572, 2916, 1572, 114** | **not binomial** |

Every entry passes four checks: weight dimensions against the D6 multiplicities, Euler characteristic against the
refined index *per flavour charge* as well as per class, non-negativity of the peeled multiplicities, and the
particle–hole symmetry $h^k=h^{48-k}$ (unimposed). Exact total: **3 008 825 568 BPS states**, of which
**19.0% are invisible to the index**.

Three things are uniform across all thirteen:

* **The window is exactly $k=22\ldots26$ in every case** — width $5=N+1$, centred on $k=24=pN^2/2$. Not one irrep
  deviates.
* **$c=0$ is saturated** (cohomology at $k=24$ alone) in all thirteen.
* **$c=1$ and $c=2$ are never saturated** — two degrees each, $\{22,25\}$ and $\{23,26\}$, with non-vanishing
  index. The $N=4$ failure of saturation reported in §4c is not special to the maximal weight; it is the rule here.

### The binomial profile breaks, exactly as at $N=3$

Eleven of thirteen irreps follow $g\binom4j$. The two at $C_2=107$ do not: $114, 1572, 2916, 1572, 114$ against
the $486,1944,2916,1944,486$ that the centre value would imply. This is the precise analogue of $(5,0,-5)$ at
$N=3$ ($38,486,486,38$ against $162,486,486,162$), and the deviation has the same shape in both cases: **the centre
is unaffected and the window edges are suppressed**. So the $\binom Nj$ profile is a property of the high-Casimir
top of the spectrum, not a law — while the *window* $k\in[pN^2/2\pm N/2]$ has survived every irrep tested at
$N=2,3,4$. The width conjecture is on much firmer ground than the profile conjecture, and it is the width that
matters for concentration.

## 4e. Why the window is wider than the index-based prediction, and why saturation must fail for $N\ge3$

Before the cohomology, the expected windows were $\{13,14\}$ at $N=3$ and $\{23,24,25\}$ at $N=4$
(`docs/todo.md`, BPS-exclusion entry). Those came from the refined index **plus the assumption of index
saturation**: the index has three classes $c=k\bmod3$ (because $[N_\Psi,Q]=3Q$), and if each class concentrates
in one degree, the window is the set of those degrees — two at odd $N$ (one class has vanishing index), three at
even $N$. The measured windows are $\{12,\dots,15\}$ and $\{22,\dots,26\}$, one degree wider on each side.

The reason is pure pigeonhole, and it reproduces every measurement including *which* class fails:

| $N$ | window (width $N+1$) | $c=0$ | $c=1$ | $c=2$ | consequence |
|---|---|---|---|---|---|
| 2 | $5,6,7$ | $\{6\}$ | $\{7\}$ | $\{5\}$ | every class saturates |
| 3 | $12,13,14,15$ | $\{12,15\}$ | $\{13\}$ | $\{14\}$ | one class cannot saturate |
| 4 | $22,\dots,26$ | $\{24\}$ | $\{22,25\}$ | $\{23,26\}$ | two classes cannot saturate |
| 5 | $35,\dots,40$ | $\{36,39\}$ | $\{37,40\}$ | $\{35,38\}$ | *(predicted)* all three fail |

A window of width $N+1$ distributed over a grading of period 3 can put one degree in each class only when
$N+1\le3$. So:

* **$N=2$ saturates for kinematic reasons**, not because anything special holds — width 3 against period 3 is the
  one case where each class gets exactly one degree. The clean $N=2$ saturation recorded earlier is therefore not
  evidence that saturation persists.
* **At $N=3$ the surplus lands in $c=0$**, whose index vanishes identically, so the extra states are wholly
  index-invisible and the non-zero-index classes still saturate. That is exactly what §3 measured.
* **At $N=4$ the surplus lands in $c=1,2$**, whose index does *not* vanish — so saturation fails where the index
  is macroscopic. That is the stronger failure reported in §4c/§4d.
* **If the width conjecture holds, strict single-degree concentration is impossible for every $N\ge3$**, and at
  $N=5$ all three classes would fail.

The logic runs one way only: the width is the empirical input (25 irreps at $N=2,3,4$, no exceptions, and it
requires the cohomology to occupy *all* $N+1$ degrees — which it does in every irrep computed), and the class
structure is then forced. Nothing here is a proof of the width itself.

Note also that the index is not missing a runaway fraction: $19.7\%$ of BPS states invisible at $N=3$ (top eight
irreps), $19.0\%$ at $N=4$ (top thirteen) — roughly flat over the range computed.

**One clarification about the bootstrap.** The exclusion bootstrap never asserted that $k=12,15$ were empty at
$N=3$; its margins there were $\sim10^{-8}$, i.e. numerically zero and inconclusive (§4). The narrowing to
$k=13/14$ came from the index-saturation assumption layered on top, not from a certificate. The bootstrap's actual
exclusions — $k\le11$ at $C_2=48$ — are all confirmed.

## 4f. $N=5$, and what the result does and does not say about concentration (2026-09-23)

**$N=5$ maximal weight $(12,6,0,-6,-12)$, 49 s at 0.41 GB** (max weight space 6435). Every §4e prediction holds:

| $k$ | 35 | 36 | 37 | 38 | 39 | 40 |
|---|---|---|---|---|---|---|
| $h^k$ | 243 | 1215 | 2430 | 2430 | 1215 | 243 |

Window width $6=N+1$, centred on $k=37.5=pN^2/2$; profile $3^5\binom5j$; particle–hole symmetry unimposed and
satisfied; and — the part that had not been tested — **all three classes carry two degrees each**
($c=0:\{36,39\}$, $c=1:\{37,40\}$, $c=2:\{35,38\}$), so none can saturate, exactly as the pigeonhole argument
of §4e predicts for $N\ge5$. The $c=0$ index cancels to zero; $c=1,2$ have $|I|=2187$ against 2673 states.

### The maximal irrep reproduces Chen's non-concentrated single-matrix baseline

Verified at $N=2,3,4,5$, the maximal irrep's BPS partition function is

$$Z_{\rm BPS}(q)=3^N(1+q)^N q^{k_{\min}},\qquad\text{total } 6^N \text{ multiplets.}$$

Chen 2025 (2.20) for the **single-matrix** model is $3^{N(N-1)/2}(1+q)^Nq^{N(N-1)/2}$ — the *same* $(1+q)^N$
factor, hence the same $N+1$ sectors and the same binomial weights. Chen (2.22) and note item 4 cite exactly that
$N+1$ spread as the reason the single-matrix model has **no R-charge concentration**. By the project's own
criterion (`docs/research_questions.md`, Q1: "concentration ⇔ $E_0>0$ outside a window of at most three
consecutive $N_\Psi$"), the three-matrix model's maximal irrep is therefore **not concentrated**, at every $N$ we
can compute.

### Why this is not yet a refutation of the conjecture

CCSY condition concentration on complexes of **macroscopic index**. The maximal irrep has $6^N$ BPS multiplets —
$e^{1.79N}$, exponential in $N$ but not in $N^2$ — and every irrep computed here at $N\ge3$ sits in that same thin
high-Casimir top. The macroscopic-index complexes are elsewhere entirely:

| $N$ | maximal irrep multiplicity | largest multiplicity in the Fock space | located at |
|---|---|---|---|
| 3 | 512 | 511 488 | $(2,0,-2)$, $C_2=8$ |
| 4 | 4 096 | 6 588 297 216 | $(3,1,-1,-3)$, $C_2=20$ |

So increasing $N$ at the top of the spectrum — $N=6$ and beyond — would keep extending the same pattern without
ever reaching a complex the conjecture actually speaks about. **The informative direction is down in Casimir, not
up in $N$.** Reachability at $N=3$ with the blocked routine (full map computed 2026-09-23):

* $C_2=24$, $(4,0,-4)$: multiplicity $110\,080$, weight space 38 763 → blocked 12 921, ~1.3 GB — comfortably
  affordable, and $215\times$ the multiplicity of the top irrep.
* $C_2=20$, $(3,1,-4)$ and $(4,-1,-3)$: multiplicity $161\,280$, blocked 25 528, ~10–13 GB and hours — needs the
  float64-block memory reduction first.
* $C_2\le15$ (including the peak at $(2,0,-2)$): blocked dimension $\ge61\,000$, out of reach.

A width-4 window at $C_2=24$ would make the negative result much harder to attribute to the special kinematics of
the maximal irrep; a narrower window there would overturn the present picture.

## 4g. Descending in Casimir at $N=3$: the window never narrows (2026-09-23)

Thirteen irreps, complete Casimir tiers from 48 down to 24, all checks passing (weight dimensions, Euler
characteristic per class *and* per flavour charge, non-negativity, particle–hole). Data:
`results/data/cohomology_N3_full.json`.

| $C_2$ | Fock multiplicity | support | width | $h^{12}$ | $h^{13}$ | centre/edge | index-invisible |
|---|---|---|---|---|---|---|---|
| 48 | 512 | $\{12..15\}$ | 4 | 27 | 81 | 3.0 | 25.0% |
| 42 (×2) | 2 560 | $\{12..15\}$ | 4 | 54 | 162 | 3.0 | 25.0% |
| 38 (×2) | 4 608 | $\{12..15\}$ | 4 | 81 | 243 | 3.0 | 25.0% |
| 36 (×2) | 2 560 | $\{12..15\}$ | 4 | 27 | 81 | 3.0 | 25.0% |
| 35 | 15 360 | $\{12..15\}$ | 4 | 38 | 486 | 12.8 | 7.3% |
| 30 (×2) | 35 840 | $\{12..15\}$ | 4 | 36 | 810 | 22.5 | 4.3% |
| 27 (×2) | 35 840 | $\{12..15\}$ | 4 | 18 | 648 | 36.0 | 2.7% |
| **24** | **110 080** | $\{12..15\}$ | 4 | 64 | 1620 | 25.3 | 3.8% |

**The headline: the support is exactly four consecutive degrees in every one of the thirteen**, including
$(4,0,-4)$ at $C_2=24$, whose Fock multiplicity of $110\,080$ is $215\times$ that of the maximal irrep. Strict
single-degree concentration therefore fails not only at the thin top of the spectrum but at the most macroscopic
complex currently reachable. This is the robust negative result, and it is exactly what the pigeonhole of §4e
forces once $N\ge3$.

**Correction to the first version of this section.** On the twelve irreps down to $C_2=27$ the edge multiplicity
$h^{12}$ was falling ($38,36,18$) and the index-invisible fraction was falling monotonically
($25\%\to7.3\%\to4.3\%\to2.7\%$); I flagged the possibility that $h^{12}$ would reach zero, collapsing the
support to $\{13,14\}$ and restoring exact concentration. **It does not.** At $C_2=24$ the edge rises to $64$ and
the invisible fraction back to $3.8\%$. Neither quantity is monotone in $C_2$, and there is no sign of the window
closing.

What survives is weaker but still real: the *profile* is far more peaked away from the top ($3.0$ at
$C_2\ge36$ versus $25$–$36$ below it), so the index is a good proxy for the BPS **count** (96–97%) even though it
is wrong about the **support** everywhere. Over all thirteen irreps: $1\,914\,064$ BPS states, $10.1\%$ invisible
to the index — and that aggregate is again dominated by the flat-ratio top, so it should not be quoted as a
property of the spectrum.

**Still untested:** $C_2\le20$, where the multiplicity keeps growing to the peak $511\,488$ at $(2,0,-2)$,
$C_2=8$. Reachability (§4f): $C_2=20$ needs blocked dimension 25 528 (hours, ~5 GB with the float64 fix);
$C_2\le15$ needs $\ge61\,000$ and is out of reach. So "the window is 4 at every $N=3$ complex" is established
over multiplicities $512$–$110\,080$, a factor of 215, but not over the final factor of ~5 to the peak.

## 4h. Is the negative result an artefact of under-refinement?  No — the flavour symmetry is $\mathbb Z_3$, not $S_3$ (2026-09-23)

A finer grading can only narrow the support of each complex, so before concluding anything negative one must check
that every available symmetry has been used. `docs/research_questions.md` Q1 describes the complexes as labelled by
"$N_\Psi$ mod 3, the $SU(N)$ irrep, and flavor ($S_3$) quantum numbers" — but for Chen's $C$ the flavour symmetry
group is only $\mathbb Z_3$:

    perm (0,1,2) even : C invariant
    perm (1,2,0) even : C invariant
    perm (2,0,1) even : C invariant
    perm (0,2,1) odd  : NOT invariant (max |dC| = 1/3)
    perm (1,0,2) odd  : NOT invariant
    perm (2,1,0) odd  : NOT invariant

The cyclic symmetrisation of $C_{abc}=\mathbb 1_{a\le b\le c}$ is invariant under the even permutations $A_3$ only;
the transpositions are not symmetries, so there is no $S_3$ refinement available and no finer grading than the
$(k\bmod3,\ \lambda,\ w)$ we already use. Q1's "$S_3$" should read $\mathbb Z_3$. **The negative result is
therefore not an artefact of insufficient refinement.**

## 4i. Answer to Q1, with scope

**For every complex computed, the model does not exhibit R-charge concentration.** The BPS cohomology occupies
$N+1$ consecutive degrees in all 27 irreps computed across $N=2,3,4,5$, with no exception; since the complexes are
graded by $k\bmod3$, at least one class carries $\ge2$ degrees whenever $N\ge3$, violating CCSY's single-degree
criterion.

**The load-bearing observation concerns Chen's evidence.** At $N=2$ the window has width $3$ against a period-3
grading, so each class receives exactly one degree and concentration holds *for kinematic reasons alone*. That is
the regime of the unpublished $N=2$ simulation which Chen 2025 §4 cites as "suggesting" concentration. From $N=3$
onwards it cannot hold. The sole evidence for the matrix-SYK proposal therefore comes from the one value of $N$
where the conclusion is automatic and uninformative.

Reinforcing this, the maximal irrep satisfies $Z_{\rm BPS}=3^N(1+q)^Nq^{k_{\min}}$ (§4f), carrying the same
$(1+q)^N$ as Chen's own single-matrix result (2.20) — the model Chen (2.22) explicitly identifies as *not*
concentrated.

**Scope, and what would change the answer.**

1. *Coverage.* At $N=3$ the descent reaches Fock multiplicity $110\,080$ against a peak of $511\,488$ — a factor
   $\sim5$ short. At $N=4$, 13 of 213 irreps, all high-Casimir. The low-Casimir bulk is untested and needs a
   sparse (not merely blocked) rank routine.
2. *Asymptotics.* CCSY's conjecture is a large-$N$ statement conditioned on macroscopic index; finite-$N$ data
   cannot formally refute it. What can be said is that there is no sign of concentration switching on, and that
   the mechanism defeating it ($N+1>3$) strengthens with $N$.
3. *The approximate statement survives.* Away from the high-Casimir top the profile is strongly peaked and the
   index captures 96–97% of the BPS states (§4g). The index is nearly right about the **count**; it is the
   **support** claim that fails.

**Assessment.** Chen's conjecture that $p=3$ yields a concentrating "matrix SYK" model is not supported by this
computation: at the top of the spectrum the model reproduces the solvable single-matrix structure, and the $N=2$
evidence behind the proposal is the unique case where concentration is kinematically forced. This is stated as a
finding over the computed range, not as a theorem; $C_2\le20$ at $N=3$ is the next thing that would firm it up.

## 4j. The window is independent of $C$, and the no-go is structural (2026-09-23)

$Q^2=0$ holds for **any** $C$ here: $Q=C_{abc}\mathrm{Tr}[\Psi^a\Psi^b\Psi^c]$ is cubic in creation operators, so
$Q^2$ contracts the symmetric $C\otimes C$ against a totally antisymmetric six-fermion product and vanishes. Any
cubic tensor is therefore a legitimate model, and the machinery can screen them in seconds. At $N=3$, maximal
weight:

| $p$ | $C$ | window | width | profile |
|---|---|---|---|---|
| 1 | single matrix (Chen) | $\{3,4,5,6\}$ | 4 $=N+1$ | $\binom3j$ |
| 2 | Chen-type | $\{6,\dots,12\}$ | **7 $=2N+1$** | $\binom6j$ |
| 3 | Chen's | $\{12,\dots,15\}$ | 4 | $\binom3j$ |
| 4 | Chen-type | $\{15,\dots,21\}$ | **7** | $27\binom6j$ |
| 3 | random, seed 1 | $\{12,\dots,15\}$ | 4 | $\binom3j$ |
| 3 | random, seed 2 | $\{12,\dots,15\}$ | 4 | $\binom3j$ |

Three consequences.

1. **$p=1$ reproduces Chen (2.22)/(2.21) exactly** — window $N(N-1)/2\le N_\Psi\le N(N+1)/2$ with weights
   $\binom Nj$. An independent end-to-end validation of the pipeline against a published closed form.
2. **The cohomology at the maximal weight is *generically* independent of $C$** — but not universally.
   **Corrected 2026-09-23** (the first version of this item claimed outright $C$-independence on the strength of
   two random seeds; that was an over-generalisation). Sampling 30 random cyclic tensors at $N=3$, $p=3$:

   | outcome | count | width |
   |---|---|---|
   | $h=(27,81,81,27)$ — Chen's exact profile | 26/30 | 4 |
   | $h=(1,9,36,84,126,126,84,36,9,1)$ — the whole weight space | 4/30 | 10 |

   The second case is not $Q=0$ as an operator (those tensors act non-trivially on 20/20 generic states); $Q$
   simply annihilates the *maximal-weight* complex, leaving all of it BPS. So the generic answer is Chen's, and
   the exceptional set gives a **wider** window, never a narrower one. Either way no choice of $C$ improves
   matters, which is what the no-go needs — but the correct statement is "generic", not "always".
3. **Even $p$ is worse**: width $2N+1$ rather than $N+1$. Whatever the parity dependence is, no case gets close to
   the width $\le3$ that concentration demands.

**The no-go, stated structurally.** Concentration requires each complex — graded by $k\bmod q$ with $q=3$ here —
to carry a single degree, so it requires the BPS window to be no wider than the grading period. In SYK the window
is $O(q)$ and $N$-independent, so this is satisfiable. In these matrix models the window grows with $N$ ($N+1$ for
odd $p$, $2N+1$ for even $p$) while the period stays at 3. Concentration is therefore impossible for *any* cubic
$C$ at *any* $p$ once $N$ is large enough — not a property of Chen's particular choice.

**Caveat.** Items 2 and 3 are measured at the maximal weight at $N=3$ only. The genericity in item 2 rests on a
30-tensor sample, not a proof; the width at lower Casimir was checked only for Chen's $C$ (§4g, all 13 irreps
width 4). Whether the parity pattern $N+1$ / $2N+1$ persists in $N$ and down the spectrum is open.

## 4k. The traceless $su(N)$ variant: concentration $\iff$ $\mathrm{rank}\,G\le2$ (2026-09-23)

Chen (2.17)/(2.21) identifies the $(1+q)^N$ — i.e. the window — with the $N$ diagonal modes $\Psi_{ii}$, each
freely occupied. $su(N)$ has $N-1$ Cartan directions instead of $N$, so the window should lose exactly one.
Implemented in `src/su_model.py` (basis $E_{ij}$, $i\ne j$, plus $H_m=E_{mm}-E_{m+1,m+1}$; couplings
$g=C_{abc}\mathrm{Tr}[T_\alpha T_\beta T_\gamma]$; the $N$-slot diagonal pool of the $u(N)$ enumeration becomes a
single pool of $p(N-1)$ zero-weight modes). Mode enumeration validated by
$\sum_\lambda\dim W_\lambda(k)=\binom{p(N^2-1)}{k}$; all ranks exact mod two primes.

| model | $\mathrm{rank}\,G$ | $\lambda^*$ | window | width | $h$ | concentrates? |
|---|---|---|---|---|---|---|
| $u(2)$ | 2 | $(3,-3)$ | $\{5,6,7\}$ | 3 | $9,18,9$ | **yes** |
| $su(2)$ | 1 | $(3,-3)$ | $\{3,4,5,6\}$ | 4 | $1,3,3,1$ | no (degenerate, $Q\equiv0$) |
| $u(3)$ | 3 | $(6,0,-6)$ | $\{12..15\}$ | 4 | $27,81,81,27$ | no |
| **$su(3)$** | **2** | $(6,0,-6)$ | $\{11,12,13\}$ | **3** | $9,18,9$ | **YES** |
| $u(4)$ | 4 | $(9,3,-3,-9)$ | $\{22..26\}$ | 5 | $81,\dots$ | no |
| $su(4)$ | 3 | $(9,3,-3,-9)$ | $\{21..24\}$ | 4 | $13,81,81,13$ | no |
| $u(5)$ | 5 | $(12,6,0,-6,-12)$ | $\{35..40\}$ | 6 | $243,\dots$ | no |

**Width $=\mathrm{rank}(G)+1$** in six of the seven; the exception is $su(2)$, where $Q$ annihilates the maximal
weight space entirely so the "window" is just that space, a degenerate small-model artefact. Since concentration
requires the window to fit inside the grading period $q=3$:

$$\text{concentration}\iff\mathrm{rank}(G)\le2 .$$

$u(2)$ — the case behind Chen's $N=2$ evidence — and $su(3)$ are exactly the rank-2 models, and they give
*identical* BPS content $9,18,9$, the $1{:}2{:}1$ Turiaci–Witten profile of Q1d.

**Consequence.** Going traceless buys exactly one unit of rank and therefore rescues $N=3$ and nothing else. The
obstruction is not the trace mode as such — it is one Cartan direction among $\mathrm{rank}\,G$ of them. The window
tracks the rank while the grading period is pinned at 3 by the cubic supercharge, and the large-$N$ limit *is* the
limit of growing rank. An escape would need the supercharge degree $q$ to grow with the rank, which is outside the
class of cubic "matrix SYK" models Chen proposes.

**Caveat.** Measured at the maximal weight only, and $su(N)$ has not been pushed down in Casimir the way $u(3)$ was
in §4g. The rank rule rests on seven models at one weight each.

## 4l. Screening candidate supercharges: the escape route closes (2026-09-23)

With `cohomology.apply_Q_degree` / `maximal_weight_window` any odd-degree supercharge
$Q=\sum C_{a_1\ldots a_q}\mathrm{Tr}[\Psi^{a_1}\cdots\Psi^{a_q}]$ can be screened in seconds. Only **odd** $q$ is
admissible: $\mathrm{Tr}[\Psi^q]$ picks up $(-1)^{q-1}$ under a cyclic shift so it vanishes identically for even
$q$, and $Q^2=0$ is automatic exactly when $q$ is odd (exchanging the two $q$-tuples is a permutation of $2q$
fermions with $q^2$ transpositions, odd iff $q$ is odd, against a totally antisymmetric product). Validated: at
$q=3$, $p=1$ it reproduces Chen (2.22) exactly for $N=2,3,4$.

| axis | range screened | result |
|---|---|---|
| gauge algebra | $u(2..5)$, $su(2..4)$ | width $=\mathrm{rank}(G)+1$ (6 of 7; $su(2)$ degenerate) |
| flavours $p$ | 1, 2, 3, 4 at $N=3$ | width 4 for odd $p$, **7** for even $p$ |
| coupling $C$ | 30 random cyclic tensors | 26/30 give Chen's $27,81,81,27$; 4/30 give width **10** |
| degree $q$ | 3, 5, 7 | see below |

**The $q$ axis is the only lever that could work, and it closes.** Concentration needs the window to fit in the
grading period, width $\le q$; since the window is $\mathrm{rank}+1$ this means $q\ge\mathrm{rank}+1$. But:

* For $q\le N$ the supercharge acts non-trivially and the window is $N+1>q$. Fails.
* For $q>N$ a closed index loop of length $q$ must repeat an index, Pauli kills those terms, and $Q$ **degenerates
  on the maximal-weight complex**: at $N=3$, $p=3$, $q=5$, all 8 sampled tensors give $Q\equiv0$ there, so the
  entire weight space is BPS and the window *widens* from 4 to 10. Fails worse.

The condition $q\ge\mathrm{rank}+1$ therefore lands in precisely the regime where $Q$ stops acting. Raising the
supercharge degree buys a longer grading period but destroys the differential that was supposed to cut the window
down — the two effects work against each other.

**Caveats.** The degeneracy at $q>N$ is established at $N=3$; whether the squeeze is general needs a cell with
$q\le N$ and $q>3$ (running: $p=2$, $N=5$, $q=5$). The $p=1$ sweep across $q=3,5,7$ showed width $N+1$
independent of $q$, but there $Q$ vanishes on the maximal weight for *every* $q$ (this is Chen's (2.22), where the
whole maximal-weight space is BPS), so it only exhibits the pigeonhole and is not evidence of genuine
$q$-independence.

## 5. What this changes

* The $N=3$ concentration question is **answered exactly** for the high-Casimir end of the spectrum. It no longer
  needs the bootstrap, and it did not need level 4, level 5, or HPC.
* The width of the BPS window at $N=3$ is **4 consecutive degrees** ($12$–$15$), centred on half filling
  $k=13.5$. Earlier notes describing the top complex as "localised to $k=13/14$" were reading the index, which
  sees only $c=1,2$; the full cohomology adds the cancelling pair at $k=12,15$.
* **Q1b needs amending.** Index saturation was recorded as holding in every complex at $N=2$. At $N=3$ it holds
  in the $c=1,2$ complexes but **fails in every $c=0$ complex**, which has vanishing index and non-zero cohomology.
  Saturation is therefore not automatic once $N>2$.
* **Q1d needs amending.** The Turiaci-Witten $\cos(\pi(k-k_*)/3)$ profile describes the *index*, and at $N=3$ it
  gives $0{:}1{:}1$ over the classes $c=0,1,2$ — correctly. But the *BPS counts* over the same classes are
  $149{,}526 : 304{,}641 : 304{,}641 \approx 1 : 2 : 2$, not $0{:}1{:}1$. The $c=0$ class is not empty; it is merely
  invisible to the index. (Caveat: summed over the eight $C_2\ge35$ irreps only, not the full spectrum.)
* The bootstrap's role narrows to what cohomology cannot reach: near-BPS *energies*, $N$-scaling, and the
  large-$N$ limit. Exact cohomology is limited by weight-space size, which grows fast away from the top.

## 6. Limitations

* Only the eight irreps with $C_2\ge35$ are done. The remaining $N=3$ weights run from $\dim W\approx4000$ up to
  $1.1\times10^6$ at $\lambda=(0,0,0)$; the latter needs a sparse rank routine, not dense elimination.
* At $N=4$, thirteen irreps down to $C_2=107$ are done (§4d). The next tier ($C_2=104$, weight spaces 60 984) is
  ~2.6 h each and ~12 GB, past the standing memory budget without a further idea; the bulk of the $N=4$ spectrum
  (weight spaces up to $\sim10^9$) remains far out of reach, so all $N\ge3$ statements here are about the
  high-Casimir top.
* The flavour refinement uses a dense projector built from $U^t$; for the larger weight spaces it should be
  replaced by an orbit construction (orbits of the $\mathbb Z_p$ action have size 1 or $p$, giving the eigenvectors
  directly) before it is pushed further.
