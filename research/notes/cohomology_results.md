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

**Conjecture (stated as such, not proven, and supported by exactly two data points):** the BPS window is centred
on half filling $k=pN^2/2$, has width $N+1$, and the degree profile goes as $\binom{N}{j}$.

The caveat is sharp and already visible: at $N=3$ the $\binom{3}{j}=1{:}3{:}3{:}1$ profile holds for the seven
highest-Casimir irreps but **fails** for $(5,0,-5)$, where the profile is $38:486:486:38$. At $N=2$ no such failure
occurs — all four irreps obey $1{:}2{:}1$. So either the profile is a property of the top of the Casimir spectrum
that happens to be exhaustive at $N=2$, or the total over *all* irreps obeys it while individual irreps need not.
Distinguishing these requires the low-Casimir $N=3$ weights, which are not yet computed (§6). The width-$(N+1)$
part of the conjecture is on firmer ground than the profile part, and $N=4$ would be the first real test of both.

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
* $N=4$ is untouched by this method so far. The maximal weight $(9,3,-3,-9)$ should be comparably cheap and is the
  obvious next target; the interior weights are not.
* The flavour refinement uses a dense projector built from $U^t$; for the larger weight spaces it should be
  replaced by an orbit construction (orbits of the $\mathbb Z_p$ action have size 1 or $p$, giving the eigenvectors
  directly) before it is pushed further.
