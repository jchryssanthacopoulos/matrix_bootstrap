# The $U(n)^3$ fermionic quiver as a project: literature check and fortuity test (2026-10-07)

**Progress report (2026-10-08):** `research/tex/quiver_progress_report.tex` → `research/pdfs/quiver_progress_report.pdf` (12 pp.; model, literature, all exact and numerical results of §§1–8, prospects and next steps). New figure for it: `results/figures/quiver_concentration.{pdf,png}` (`scripts/plot_quiver_concentration.py`): singlet dimensions and BPS states per degree, and the lowest non-BPS singlet energy per degree, at $(2,2)$ and $(2,3)$.

## Status summary (updated 2026-10-08; details in §§1–8)

| property | status | evidence |
|---|---|---|
| Escapes the adjoint-model obstruction | **derived** | bifundamentals have no zero weights ($V_0=0$), so the window-widening mechanism of Chen-type models is absent; this does not by itself prove concentration |
| Concentration (singlet sector) | **verified** at $(2,2)$ (exact ranks); **numerical** at $(2,3)$ | 90 at $k=12$; 1680 at $k=18$ (§5) |
| Concentration at $n=3$ | **open** | $\vert I_0\vert(3,2)=1680$ exactly, with a sign consistent with half filling; zero-weight sectors of about $5\times10^{10}$ block direct tests |
| $p=1$ | trivial | every per-irrep index is $\pm1$; singlet index 0 |
| Macroscopic BPS entropy | **$p=2$: no, proved** ($I_0=(-1)^n(3n)!/(n!)^3$ for all $n$, D18); **$p\ge3$: yes, at the level of the index** ($\ln|I_0|=F^*(p)n^2+O(\ln n)$, $F^*(3)=1.2373$, $F^*(4)=2.5550$; 41% of the singlet log-count at $p=3$; rigorous up to standard log-gas large deviations for even $p$, plus a sign proviso for odd $p$) | §§7–8, D18 |
| Fortuity | **verified** at leading order at $(2,2)$, ordinary sense | 90/90 obstructed (§3); expected for this class; refined (Choi–Choi–Kim) criterion not checked |
| Chaos (non-BPS or BPS) | **untested** | — |
| Near-BPS gaps / super-Schwarzian | **untested**; one hint | the lowest penalised singlet eigenvalues fall to 0.054 next to the window at $(2,3)$ |
| Large-$n$ control | planar (not melonic); singlet sector suits the standard bootstrap | Witten eq. 3.10; synthesis §10 |
| Novelty | no prior study found | §1 |

Model (cohomology notes §4x): $A\in(n,\bar n,1)$, $B\in(1,n,\bar n)$, $C\in(\bar n,1,n)$ of $U(n)^3$, $p$ flavours per
edge, all fermionic creation operators, $Q=\sum_{abc}C_{abc}\mathrm{Tr}(A^aB^bC^c)$, $H=\{Q,Q^\dagger\}$. Known before
this note: at $(n,p)=(2,2)$ (random integer couplings, seed 3) there are exactly 90 BPS singlets, all at $k=12$
(half filling), saturating the singlet index.

## 1. Literature check (web search, 2026-10-07)

Searched for purely fermionic quiver / bifundamental supersymmetric quantum mechanics with a cubic supercharge,
its BPS cohomology, fortuity, and bootstraps. **No paper found that studies this model.** Absence from searches is
not proof; the closest relatives are:

| relative | what it shares | what differs |
|---|---|---|
| Gurau–Witten coloured models, $q=D+1$ fields of rank $D$ (Witten arXiv:1610.09758) | at $q=3$, $D=2$ the coloured structure *is* three bifundamentals with interaction $\psi^1\psi^2\psi^3$ | real fermions, no supersymmetry, studied for $D\ge3$ (melonic); $D=2$ is planar |
| Peng–Spradlin–Volovich 2016 (in `papers/`) | supersymmetric Gurau–Witten; their "purely mesonic model" has $Q\propto\chi_{ab}\chi_{bc}\chi_{ca}$, a cyclic triangle | $\mathcal N=1$, rank-4 tensors (melonic); they note the $q=4$ quarks admit **no** cubic invariant, which is exactly what $D=2$ bifundamentals do admit |
| Chang–Colin-Ellerin–Rangamani 2018, *On melonic supertensor models* (arXiv:1806.09903) | $\mathcal N=2$ tensor QM | dynamical bosons; SUSY broken in the IR |
| Denef quiver QM; Bena–Berkooz–de Boer–El-Showk–Van den Bleeken, *Scaling BPS solutions and pure-Higgs states* (arXiv:1205.5023) | three-node quiver with closed loop and cyclic superpotential; exponentially many "pure-Higgs" BPS states, read as single-centre black-hole microstates | bosonic chiral multiplets, abelian/low-rank nodes, Higgs-branch cohomology rather than a fermionic Fock space |
| De Marco–Mullahasanoglu–Raeymaekers–Rossi–Sengör 2026 (arXiv:2607.15348); Sanli 2025 (arXiv:2509.07838) | large-$N$ cyclic quivers, emergent conformal symmetry, superconformal quiver index | same bosonic Denef setting |
| Biggs–Lin–Maldacena 2026 (in `papers/`) | disorder-free $\mathcal N=2$ fermionic model with SYK IR | single $SU(2)$ multiplet, melonic, global symmetry |
| Choi–Choi–Kim 2026, *Fortuity on the conifold* (arXiv:2609.19297) | fortuity in a bifundamental (Klebanov–Witten) gauge theory | 4d $\mathcal N=1$ field theory |
| Chang–Zhang 2025, D1–D5 (arXiv:2511.23294); Belin et al. 2025 ABJM (arXiv:2512.04146); Choi–Kim 2025 (arXiv:2512.12674) | fortuity and R-charge concentration beyond SYK/$\mathcal N=4$ | CFTs |

**Assessment.** The fermionic $\mathcal N=2$ "$D=2$ coloured" quiver appears to be unstudied. Its natural framing
is the matrix ($D=2$, planar) member of the supersymmetric coloured family, and the possible link to pure-Higgs
states is suggestive but unestablished (different fields, different branch).

**Update 2026-10-07 (papers read; notes in `literature/notes/`, synthesis §10).** All four papers suggested above
have now been read. What they establish for this model:
- **Witten.** The melonic argument needs $D\ge3$ (his eq. 3.10), so the quiver, the $D=2$ member, is planar. His
  "gauging is harmless" holds only for $D\ge3$.
- **Chang–Colin-Ellerin–Rangamani.** Purely fermionic $\mathcal N=2$ models with an odd-degree supercharge had no
  melonic tensor version (their eq. 6.4). Their bosonic supertensor models break supersymmetry at large $N$.
- **Bena et al.** For the bosonic abelian quiver: closed-form pure-Higgs counts and exponential growth, and a
  suggestion that the combinatorics has a fermionic origin. A cheap $n=1$ comparison is proposed in their note; the
  naive match fails.
- **Choi–Choi–Kim.** The refined "generalised monotone" criterion for bifundamental theories is a caveat on §3
  below.

## 2. Fortuity test: definition and leading-order obstruction

**Definition used** (Tierz 2026, as adopted in `literature/notes/tierz_2026.md`). The projection
$\pi_{3\to2}$ deleting every monomial with a third index on any node is a cochain map,
$\pi Q^{(3)}=Q^{(2)}\pi$, because $Q$ contains only creation operators. A rank-2 class is monotone only if it lies in
$\mathrm{Im}\,\Pi_{M\to2}$ for every $M>2$. The inclusion $\iota$ is a cochain map for $Q^\dagger$ (whose new terms
annihilate states without index-3 modes), and with harmonic representatives $\iota_*[z]\ne0$ in
$Q^\dagger$-cohomology iff $z\not\perp\mathrm{Im}\,\Pi$. So the $Q$- and $Q^\dagger$-versions of the
Chang–Chen–Sia–Yang definition are decided by the same subspace.

**Obstruction (derived here).**
- **Setup.** Write $Q^{(3)}=Q^{(2)}+Q_{\rm new}$ and $V=V_0\oplus V_+$, with $V_+$ the states carrying index-3 modes.
  $Q^{(3)}$ preserves $V_+$.
- **Lifting criterion.** A closed $z_0\in V_0$ is $\pi$ of a closed rank-3 state iff $Q_{\rm new}z_0$ is $Q^{(3)}$-exact
  in $V_+$. Shifting $z_0$ by $Q^{(2)}y_0$ does not change this, since $Q_{\rm new}Q^{(2)}y_0=-Q^{(3)}Q_{\rm new}y_0$
  using $(Q^{(3)})^2=0$. This is the connecting map of $0\to V_+\to V\to V_0\to0$, i.e. Chang–Lin's $I_N$ picture
  with $I_N=V_+$.
- **Grading.** Grade $V_+$ by the number $m$ of index-3 modes. $Q^{(2)}$ preserves $m$. A term of $Q_{\rm new}$ with
  one index equal to 3 adds two index-3 modes; one with two or three such indices adds three. Nothing raises $m$
  by one.
- **The $m=2$ condition.** The lift equation at $m=2$ is $Q^{(2)}z_2=-Q^{(1)}_{\rm new}z_0$. Its components are
  $\sum_bC_{abc}B^b_{jk}z_0$ (from $i=3$), $\sum_cC_{abc}C^c_{ki}z_0$ (from $j=3$) and $\sum_aC_{abc}A^a_{ij}z_0$
  (from $k=3$). For couplings of full rank in the summed index this reduces to:
  $$\text{every single letter }X\text{ applied to }z_0\text{ is }Q^{(2)}\text{-exact}\iff\langle\beta|X|z_0\rangle=0\ \ \forall\,\beta\in\mathrm{BPS}_{k+1}.$$
- **Conclusion.** Any class with a non-zero "BPS three-point function" with a single letter is not in
  $\mathrm{Im}\,\Pi_{3\to2}$, hence fortuitous. Vanishing is necessary for liftability, not sufficient: orders
  $m\ge3$ are untested.
- **Reduction by gauge covariance.** Testing $X_{11}$ for each edge and flavour suffices. The map
  $X\mapsto P_{\rm BPS}Xz_0$ is equivariant and the letters of an edge form one irrep.

Code: `scripts/run_quiver_fortuity.py`. Results: §3.

## 3. Result

`results/data/quiver_fortuity_n2_p2_k12.json`.
- **Run.** $(n,p)=(2,2)$, seed 3, $C_{abc}=[[[5,1],[1,2]],[[1,5],[5,3]]]$, 5238 s, peak 5.3 GB (guarded at 9 GB).

**Checks.**
- Exactly 90 singlet BPS states at $k=12$, reproducing §4x. They are the kernel of $[Q_{12};Q_9^T;\text{all lowering
  operators}]$ in the 16 436-state zero-weight block; the last zero eigenvalue of $S^TS$ is $1.4\times10^{-13}$
  and the first non-zero one is $4.0$.
- Every letter block at $k=13$ (10 980 states) has 644 BPS states, with gaps of $6\times10^{-13}$ against 7.5–7.8.
- $X z_0$ is $Q$-closed to $4\times10^{-11}$, as required.
- The coupling maps in the summed flavour index have rank 2 for all three edges, so the per-letter test is
  equivalent to the actual $m=2$ condition.

**Result.**

| letter | BPS states at $k=13$ in its block | rank of $\langle\beta|X|z_0\rangle$ | top singular values |
|---|---|---|---|
| $A^0_{11}$ | 644 | 54/90 | 0.93, 0.93, 0.91, 0.83 |
| $A^1_{11}$ | 644 | 54/90 | 0.98, 0.98, 0.97, 0.91 |
| $B^0_{11}$ | 644 | 54/90 | 1.00, 0.99, 0.99, 0.94 |
| $B^1_{11}$ | 644 | 54/90 | 0.95, 0.93, 0.92, 0.83 |
| $C^0_{11}$ | 644 | 54/90 | 1.00, 0.99, 0.99, 0.94 |
| $C^1_{11}$ | 644 | 54/90 | 0.95, 0.93, 0.92, 0.83 |
| **all six combined** | | **90/90** | |

- **All 90 singlet BPS classes are obstructed at leading order.** None lies in $\mathrm{Im}\,\Pi_{3\to2}$, so all
  are fortuitous in Tierz's sense; dually, their inclusions are trivial in rank-3 $Q^\dagger$-cohomology.
- Each letter alone obstructs a 54-dimensional subspace, and different letters obstruct different subspaces.
- The $B\leftrightarrow C$ coincidences of singular values reflect a symmetry of these couplings; not analysed.

**Caveats.**
1. **Floating-point ranks.** The per-letter singular values are $O(1)$, but the smallest singular value of the
   combined $540\times90$ matrix was not recorded. The 90/90 verdict uses a relative cut of $10^{-8}$, and an exact
   check mod $P$, or a rerun printing the spectrum, would make it airtight.
2. **Scope.** One $n$, one random coupling, and the singlet sector only.
3. **The result is expected rather than discriminating.**
   - The BPS window sits at half filling, which moves with $n$: $k=12$ at $n=2$, $k=27$ at $n=3$. A class at fixed
     $k=12$ is therefore very unlikely to survive to rank 3.
   - This mirrors Chang–Chen–Sia–Yang's argument that generic $\mathcal N=2$ SYK has only fortuitous BPS states.
   - So fortuity confirms the black-hole-like picture but does not single the quiver out. What would
     discriminate is concentration at larger $n$, chaos of the BPS and near-BPS sectors, and the near-BPS gaps.
4. **Definition-dependence** (Choi–Choi–Kim 2026).
   - The verdict is fortuity in the ordinary (Chang–Lin/Tierz) sense.
   - $U(n)^3$ singlets are all multi-trace, but determinant-type singlets such as $\det(ABC)$ have $n$-dependent
     multi-trace expansions. Under the refined "generalised monotone" criterion, a class of the form
     (determinant ground state) × ($n$-independent excitation) would count as monotone.
   - With fermionic letters many determinants vanish. Whether any of the 90 classes is of this form has not been
     analysed.

## 4. Singlet index at $(n,p)=(3,2)$: concentration evidence (2026-10-07)

**What it tests.** Concentration, not fortuity.
- A $U(n)^3$ singlet needs $N_A=N_B=N_C=m$, so all singlets sit at $k=3m$: a single $\mathbb Z_3$ class.
- The singlet index is $I_0=\sum_m(-1)^m n(3m)$, with $n(k)$ the singlet multiplicity. It lower-bounds the number
  of singlet BPS states, $h_{\rm sing}=\sum_kh^k\ge|I_0|$, with **equality iff all singlet BPS states sit in one
  degree**.
- So $|I_0|$ large means macroscopic singlet BPS degeneracy. Concentration itself needs the cohomology in addition.

**Method** (`scripts/quiver_singlet_index.py`).
- Dual Cauchy, $\Lambda(\mathbb C^p\otimes V\otimes\bar W)=\bigotimes_f\bigoplus_{\lambda\subset n\times n}S_\lambda V\otimes S_{\lambda'}\bar W$, gives
  $$n(t)=\mathrm{Tr}\big[(D_tGT)^3\big],\qquad G_{PQ}=\int_{U(n)}s_P\,\overline{s_Q}.$$
  Here $P,Q$ run over $p$-tuples of box partitions, $T$ transposes them, and $D_t$ grades by degree.
- $G$ is computed by exact torus quadrature (Weyl measure) and then rounded; the rounding error is
  $\le3\times10^{-14}$.
- This avoids the $3pn^2$-mode weight expansion that exceeded 6.6 GB before. At $(3,2)$ it takes 10 s and 0.6 GB.
- **Validation:** matches the direct Weyl-alternating weight count exactly at $(2,1),(2,2),(2,3)$, and reproduces
  the known singlet indices $0,90,1680$ (§4x).

**Results** (`results/data/quiver_singlet_index_n3_p2.json`).

| $(n,p)$ | modes | singlet states | $\vert I_0\vert$ | ratio |
|---|---|---|---|---|
| (1,2) | 6 | 10 | 6 | 0.60 |
| (1,3) | 9 | 56 | 0 | 0 |
| (1,4) | 12 | 346 | 90 | 0.26 |
| (2,1) | 12 | 4 | 0 | 0 |
| (2,2) | 24 | 1 274 | 90 | 0.071 |
| (2,3) | 36 | 889 040 | 1 680 | 0.0019 |
| (3,1) | 27 | 8 | 0 | 0 |
| **(3,2)** | **54** | **2 115 632** | **1 680** | **0.00079** |

- **At $(3,2)$ there are at least 1680 singlet BPS states**, against 90 at $(2,2)$. The singlet multiplicities are
  palindromic, $n(27)=558\,304$ at half filling.
- **Coincidences (observed, unexplained):** $|I_0|(2,3)=|I_0|(3,2)=1680$ and $|I_0|(1,4)=|I_0|(2,2)=90$. It is not
  an $n\leftrightarrow p$ symmetry, since $(1,2)$ gives 6 and $(2,1)$ gives 0. Products $np$ agree in both
  coincident pairs, but not for $(1,2)$/$(2,1)$. Not investigated.
- **Growth.** $\log|I_0|/n^2=1.12$ at $n=2$ and $0.82$ at $n=3$, both at $p=2$. Two points do not establish
  $e^{cn^2}$ growth.
- **The index fraction falls fast** (7% → 0.08%). The index is exponentially small compared with the singlet Hilbert
  space, as it should be for BPS states. It says nothing about saturation.

**What is not established.**
- Whether the $\ge1680$ singlet BPS states at $(3,2)$ sit in one degree (concentration) or the index undercounts.
  That needs the singlet cohomology at $n=3$. Its weight spaces are about $10^6$–$10^7$ states per degree near
  half filling, so exact ranks need block Wiedemann or a singlet-adapted basis.
- A cheaper partial check: if $h^k_{\rm sing}$ vanishes away from $k=27$, by Lanczos on a few singlet blocks, the
  index is saturated.
- `--p 4 --n 2` hit the 4 GB guard: the degree-pair products are stored, 289 × 13 MB. Fix by computing them on
  the fly.

## 5. Concentration test at $(n,p)=(2,3)$: all 1680 singlet BPS states at half filling (2026-10-07)

**Why $(2,3)$ and not $(3,2)$.**
- The zero-weight sectors at $(3,2)$ hold $4.8\times10^{10}$ states at $k=24$ and $6.3\times10^{10}$ at $k=27$. These
  are exact counts by a torus trace, computed in the scratchpad (`zw_dims.py`). Any Fock-vector method (ED, Lanczos,
  exact ranks) is out of reach.
- The singlet subspaces are only about $5\times10^5$, but there is no cheap explicit basis for them.
- $(2,3)$ has the same singlet index (1680), 36 modes, and zero-weight sectors up to $2.2\times10^7$.
- At $n=2$ the singlets span 13 possible degrees ($k=0,3,\dots,36$), so concentration is not kinematically forced.

**Method** (`scripts/quiver_singlet_zero_modes.py`).
- Restrict to the zero-weight sector invariant under the $(S_2)^3$ Weyl group. It is a genuine gauge action, it
  commutes with $Q$, and singlets are invariant. This is about 8× smaller.
- Diagonalise $K=H_{\rm sym}+\mu\sum_v\|E^{(v)}_{12}\psi\|^2$ ($\mu=1$). $K\ge0$, and $\ker K$ is the singlet BPS
  space in degree $k$.
- Lanczos (`eigsh`, `which='SA'`), or dense diagonalisation for small sectors.

**Validation at $(2,2)$.** Exactly 90 zero modes at $k=12$, with the next eigenvalue 6.0. None at $k=3,6,9,15$
(lowest 546, 98.3, 31.2, 31.2). The $k=9$ and $k=15$ spectra coincide, as particle–hole symmetry requires.

**Results at $(2,3)$, seed-3 couplings.**

| $k$ | zero-weight → W-symmetric dim | lowest eigenvalues of $K$ | residuals | peak |
|---|---|---|---|---|
| 0 | 1 → 1 | 2008 (dense) | — | 0.7 GB |
| 3 | 216 → 27 | 1506.0, 1506.0, 1507.25 (dense) | — | 0.7 GB |
| 6 | 17 712 → 2 187 | 238.86, 240.47, … (dense) | — | 1.0 GB |
| 9 | 445 184 → 55 648 | 96.53, 97.05, 97.08, 99.86 | $\le2\times10^{-9}$ | 1.6 GB |
| 12 | 4 057 749 → 511 569 | 0.33484, 0.90961, 2.170, 2.474 | $\le1.2\times10^{-11}$ | 4.2 GB |
| 15 | 14 703 336 → 1 837 917 | 0.05386, 0.13493, 0.20233, 0.24990 | $\le1.2\times10^{-11}$ | 8.5 GB |

- **Robustness.**
  - $k=12$, rerun with a different start vector and 8 eigenvalues: identical to 9 digits.
  - $k=15$, rerun with a different start vector and 6 eigenvalues: identical to 9 digits (0.053859138, 0.134927966,
    …), residuals $\le9\times10^{-12}$.
- **Conclusion (numerical).**
  - There are no singlet BPS states at $k=0,3,6,9,12,15$. By particle–hole symmetry ($h^k=h^{36-k}$, from the
    top-form pairing) there are none at $k=21,\dots,36$ either.
  - Hence all singlet BPS states sit at $k=18$, and by the index $h^{18}_{\rm sing}=|I_0|=1680$.
  - **The singlet sector of the $(2,3)$ quiver is R-charge concentrated**, at half filling, with 1680 states. This
    is the second concentrated data point after $(2,2)$, at a different flavour number.
- **Near-BPS singlets.** The lowest singlet-ish eigenvalues fall sharply towards the window (0.33 at $k=12$, 0.054
  at $k=15$). They are eigenvalues of $K$, not pure singlet energies, so they are only upper bounds on the singlet
  gaps when the penalty term is small. They suggest small near-BPS gaps next to the window: a target for the
  near-BPS programme.
- **Caveats.**
  - Lanczos is a numerical method; the robustness reruns guard against a missed eigenvalue.
  - One coupling seed.
  - $n=2$ only. The $(3,2)$ concentration remains untested directly.

**Indirect hint at $(3,2)$.** The sign of the index is consistent with concentration at half filling there:
$I_0=-1680$, and half filling $k=27$ has $m=9$ odd, so $(-1)^9h^{27}=-h^{27}$. Concentration at $k=24$ or $30$ (even
$m$) would give a positive index. This is consistency only, not evidence of concentration. A direct test needs a
singlet-adapted algorithm, e.g. the dual-Cauchy decomposition used for the index with $Q$ written as Pieri maps
between $(P_A,P_B,P_C)$ components, at a dimension of about $5\times10^5$.

## 6. Step 1: near-BPS singlet spectrum and chaos at $n=2$ (2026-10-07)

Code:
- `src/quiver_n2.py` (W-symmetric zero-weight sectors, general operator builder, flavour map $F$);
- `scripts/quiver_singlet_spectrum.py` (`pairs`, `edge`, `lmrs`);
- `scripts/quiver_p2_n4.py`;
- `scripts/quiver_singlet_pair_large.py`;
- `scripts/quiver_commutant.py`;
- `scripts/plot_quiver_step1.py`.

### 6.1 A fact used throughout: the Casimir penalty is exact

On W-invariant zero-weight vectors, $\sum_v\|E^{(v)}_{12}\psi\|^2=\langle\psi|\sum_vJ_v^2|\psi\rangle$, because
$J_-J_+=J^2$ at $J_z=0$. So the penalty $P$ of §5 is the sum of the three nodes' $SU(2)$ Casimirs. It commutes with
$H$, vanishes on singlets, and is $\ge2$ on every non-singlet.

**Consequence.** Eigenvalues of $K=H+\mu P$ below $2\mu$ are *exact* singlet energies. §5 called them upper bounds;
that was too cautious. Each level is further classified as lower ($Q^\dagger\psi=0$) or upper ($Q\psi=0$) member,
which assigns it to a $Q$-pair.

### 6.2 Exact singlet bases and $Q$-pair spectra

- **Singlet bases.** $U_k=\ker P$ in the W-symmetric sector, dense, wherever the dimension is at most 4000: all
  degrees at $(2,2)$; $k\le6$ and $k\ge30$ at $(2,3)$. In every case the singlet count equals the dual-Cauchy $n(k)$.
- **Checks.** $QU_k$ stays in the singlets to $10^{-15}$. The $(2,2)$ BPS count is 90 at $k=12$.
- **Multiplet energies.** For each pair $(k,k+3)$ they are the nonzero eigenvalues of $(QU_k)^T(QU_k)$.

### 6.3 Hidden $\mathcal N=4$ supersymmetry of the $p=2$ quiver on gauge singlets (verified numerically; not derived)

**Discovery path.** The $(2,2)$ pair spectra were full of exact coincidences:
- every level of pair $(6,9)$ also appears in pair $(9,12)$ (43 of 49);
- pair $(9,12)$ has 95 exact doublets.

Searching the 8-dimensional family $\tilde Q=Q(\tilde C)$ for $[H,\tilde Q]=0$ on singlets gives a clean
2-dimensional null space: $C$ and $\tilde C=(\varepsilon\otimes\varepsilon\otimes\varepsilon)\bar C$, with
$\varepsilon=\begin{psmallmatrix}0&1\\-1&0\end{psmallmatrix}$. For the seed-3 couplings,
$\tilde C_{abc}=(-1)^{a+b+c}C_{\bar a\bar b\bar c}$.

**Relations on every singlet sector** (to $10^{-15}$–$10^{-14}$; seed-3 integer, two Gaussian real and two Gaussian
complex coupling sets):
$$\tilde Q^2=0,\quad\{Q,\tilde Q\}=0,\quad\{\tilde Q,\tilde Q^\dagger\}=H,\quad\{Q^\dagger,\tilde Q\}=0 .$$
- This is an $\mathcal N=4$ supersymmetry algebra.
- The flavour map $F$ ($\varepsilon$ on every edge) is a unitary symmetry of $H$ on singlets for real couplings,
  with $FQF^{-1}=\tilde Q$ and $F^2=(-1)^k$.
- For complex couplings only the antiunitary $\Theta=KF$ survives, with $\Theta^2=(-1)^k$. This is seen as exact
  Kramers doubling at odd $k$.

**Only on singlets.** On the full (non-singlet) zero-weight sectors, $\{\tilde Q,\tilde Q^\dagger\}-H$ and
$\{Q^\dagger,\tilde Q\}$ are about 30% of $|H|$. The extended algebra closes only up to gauge transformations, as in
gauged supersymmetric QM. It is a property of the **gauged** model.

**Multiplet structure.** The singlet spectrum is entirely organised into long $\mathcal N=4$ multiplets, each
occupying degrees $k,k+3,k+3,k+6$. The bottom spaces $B_k=\ker Q^\dagger\cap\ker\tilde Q^\dagger$ (non-BPS) have
dimensions $1,6,43,196,43,6,1$ at $k=0,3,\dots,18$, and
$$n(k)=b(k)+2b(k-3)+b(k-6)+h(k)$$
holds in every degree, for all coupling sets.

**$p=3$ is different.** Searching the 27-dimensional cubic family at $(2,3)$ on singlets $0\to3\to6$ gives only $C$:
no hidden second supercharge. So $p=3$ is a genuine $\mathcal N=2$ model, and $p=2$ is special. This is plausibly
tied to the pseudoreality of the $SU(2)$ flavour doublet; that link is not established.

**Not done.**
- An analytic proof of the singlet-sector identities, and the $n$-dependence ($n=2$ only).
- Whether the $(1,4)$/$(2,2)$ and $(2,3)$/$(3,2)$ index coincidences of §4 are related.
- Whether an $SU(2)_R$ (needed for a small $\mathcal N=4$ super-Schwarzian) exists. No continuous flavour generator
  commutes with $H$: only the trivial $N_A,N_B,N_C$ do.

### 6.4 Level statistics

The statistic is the mean ratio of consecutive spacings $\langle r\rangle$, over the central 80% of each ensemble,
against size-matched references sampled here (Gaussian matrices of the same size; Poisson).

| ensemble | levels | $\langle r\rangle$ | random-matrix reference | Poisson |
|---|---|---|---|---|
| $(2,3)$ pair $(6,9)$, $\vert q\vert=10.5$ | 676 | $0.507\pm0.011$ | GOE-type $0.533\pm0.013$ | $0.386\pm0.013$ |
| **$(2,3)$ pair $(9,12)$, $\vert q\vert=7.5$** | **10 024** | **$0.5268\pm0.0028$** | GOE-type $0.5343\pm0.0026$ | $0.3869\pm0.0042$ |
| $(2,2)$ bottom $k=9$, $F=+i$, seed-3 integer | 98 | $0.402\pm0.035$ | GUE $0.598\pm0.032$ | $0.389\pm0.036$ |
| $(2,2)$ bottom $k=9$, $F=+i$, Gaussian real (two seeds) | 98 | $0.331$, $0.445$ ($\pm0.033$) | GUE $0.598$ | $0.389$ |
| $(2,2)$ bottom $k=9$, Kramers-reduced, Gaussian complex (two seeds) | 98 | $0.479$, $0.497$ ($\pm0.034$) | GSE $0.675\pm0.025$ | $0.389$ |
| $(2,2)$ bottom $k=6$, $F=\pm1$ blocks | 25 / 18 | $0.31$–$0.50$ | GOE $0.53\pm0.07$ | $0.39$ |

- **$(2,3)$: chaotic.**
  - Pair $(9,12)$ (10,024 multiplets, $|q|=7.5$): $\langle r\rangle$ within 2σ of the GOE class and about 30σ from
    Poisson, with no degeneracies. The locally unfolded spacing distribution follows the Wigner surmise (figure (c)).
  - The far pair $(6,9)$ agrees: 676 levels, about 10σ from Poisson.
  - Both pairs sit about 2σ below GOE. That could be a weak residual structure or a statistical fluctuation.
  - **Method for $(9,12)$** (`scripts/quiver_singlet_pair_large.py`):
    - project random vectors onto singlets with the exact Casimir polynomial
      $\prod_v\prod_{j=1}^6(1-J_v^2/j(j+1))$, which leaves a non-singlet residual of $1.4\times10^{-15}$;
    - compress $Q^TQ$ onto their span.
  - **Checks:** rank 10,700 equals $n(9)$; the kernel of 676 equals the rank of $Q_6$; 10,024 multiplets equals
    the predicted rank of $Q_9$; $E_0=151.34$ equals the exact Lanczos value. Validated against the dense pair
    $(6,9)$ (every level within $2.5\times10^{-6}$ relative).
  - **Precision.** Float32 storage of the projected vectors was rejected at this size: the Gram noise was 0.06, and
    $E_0$ was off by 0.04, about 10% of the local spacing. The reported numbers use float64 (peak 8.6 GB).
- **$(2,2)$: not random-matrix** in any symmetry class, after resolving $F$ (real couplings) or Kramers (complex).
  - Two exact degeneracies per $k=9$ block survive in every coupling class.
  - **Conserved-charge search:** within the 109-dimensional span of all gauge-invariant $k$-preserving operators
    with at most four fermion operators (bilinears, 288 single-trace quartic loops, double traces), the commutant
    of $H$ on singlets at $k=9$ is exactly $\{1,H\}$. $H$ lies in the span to $7\times10^{-15}$.
  - So no simple conserved charge explains the statistics. The options are a genuinely non-chaotic small sector
    (98 levels), higher-degree charges, or further discrete symmetries outside the span; the degeneracies point to
    the last. **Open.**

### 6.5 BPS chaos at $(2,2)$ (projected-operator statistic)

$P_{\rm BPS}OP_{\rm BPS}$ on the 90 singlet BPS states, which split 48/42 under $F$:
- generic simple operators: $\langle r\rangle=0.505\pm0.037$ ($N_{A0}$) and $0.462\pm0.030$
  ($\mathrm{Tr}A^0B^0\bar B^0\bar A^0$+h.c.), against GOE $0.530\pm0.038$;
- $F$-anticommuting operators (positive half, 42 levels): $0.32\pm0.05$ and $0.54\pm0.05$, against GOE $0.53\pm0.06$;
- an $F$-commuting operator, per $F$ block: $0.47\pm0.04$ and $0.45\pm0.05$.

Mostly within 1–2σ of GOE, with one outlier at 3.6σ low. **Inconclusive**: 42–90 levels, and the non-BPS spectrum of
the same model is not random-matrix-like. Not attempted at $(2,3)$, where the 1680 BPS vectors at $k=18$ are out of
reach.

### 6.6 Near-BPS singlet edges and the Turiaci–Witten law

Lowest singlet multiplet energy $E_0$ of each $Q$-pair. Sources:
- $(2,2)$: dense, exact;
- $(2,3)$, $|q|\ge10.5$: dense;
- $(2,3)$, $|q|=1.5,4.5,7.5$: Lanczos on $K$ with every level classified (§6.1). The $k=9,12,15$ runs found only
  singlets among the 6, 8 and 12 lowest eigenvalues of $K$, so these level lists are complete below the quoted
  energies.

| $\vert q\vert$ | pair at $(2,3)$ | $E_0$ at $(2,3)$ | next levels at $(2,3)$ | $E_0$ at $(2,2)$ |
|---|---|---|---|---|
| 1.5 | (15,18) | **0.05386** | 0.1349, 0.2023, 0.2499, 0.2988, 0.3163, 0.3605, 0.3894, 0.4548, 0.4813, 0.5403 | 95.56 |
| 4.5 | (12,15) | **0.3348** | 0.9096, 2.170, 2.474, 3.109, 3.557, 5.229, 5.381 | 199.95 |
| 7.5 | (9,12) | 151.34 | 151.63, 154.57, 156.80, 157.50, 159.50 | 546 |
| 10.5 | (6,9) | 575.23 | (676 levels, §6.4) | 728 |
| 13.5 | (3,6) | 1506 | | |
| 16.5 | (0,3) | 2008 | | |

**Cross-check.** The $|q|=4.5$ edge 0.3348 appears as a lower member at $k=12$ and as an upper member at $k=15$.

**Findings.**
- **Small near-BPS gaps exist at $(2,3)$, and only next to the BPS degree.**
  - The two innermost multiplet sectors have $E_0=0.054$ and $0.335$. Their partners sit at $E\approx150$ and above,
    on the microscopic scale set by the couplings.
  - At $(2,2)$ there is no small gap at all: the first multiplet is at 95.6, against total scales of order $10^3$.
  - So the near-BPS scale drops by more than three orders of magnitude from $(2,2)$ to $(2,3)$, as the singlet
    sector grows from 568 to 300,584 states at half filling. Two data points with different supersymmetry, so no
    scaling law is inferred.
- **Turiaci–Witten tested quantitatively at $(2,3)$:**
  - Fix $E_s=0.862$ from $E_0(1.5)$, using $E_0(q)=q^2E_s/(4\hat q^2)$ with $\hat q=3$. The predicted
    $E_0(4.5)=0.485$; the measured value is 0.335. The ratio is 6.2 against 9, a local exponent of 1.66 against 2.
  - Fix the normalisation of $\rho_q\propto\sinh(2\pi\sqrt{(E-E_0)/E_s})/E$ from the eleven $|q|=1.5$ levels. The
    predicted numbers of $|q|=4.5$ multiplets below $E=1.0$, 2.17, 3.11 and 5.38 are 5.4, 258, 2000 and
    $8.7\times10^4$; the measured numbers are 2, 3, 5 and 7.
  - **A single Schwarzian scale does not describe the two sectors.**
- **Why this is unsurprising.** Turiaci–Witten is a large-$e^{S_0}$ statement in which neighbouring charge sectors
  have the same entropy. At $(2,3)$ the singlet sectors at $k=12$, 15 and 18 have 71,685, 211,113 and 300,584 states,
  changing by factors of 3–4 per step $\Delta k=\hat q$. The test needs larger $n$; the bootstrap is the candidate
  route.

**Figure.** `results/figures/quiver_step1.{pdf,png}`:
- (a) the edges, with the $q^2$ lines drawn through $|q|=1.5$;
- (b) all $\langle r\rangle$ values with size-matched references;
- (c) the spacing distribution with local unfolding over ±8 levels. A global polynomial unfolding gave a spurious
  Poisson-like shape, so the global method is not used;
- (d) the counting-function test above.

## 7. Is the singlet BPS entropy macroscopic? Index growth with $n$ (2026-10-08)

**Motivation.** The first singlet indices all matched the $n=1$ quiver with $np$ flavours, i.e. Dixon's identity
$\sum_m(-1)^m\binom{np}m^3=(-1)^{np/2}(3np/2)!/((np/2)!)^3$:
$$|I_0(2,2)|=90,\qquad|I_0(2,3)|=|I_0(3,2)|=1680 .$$
If that held generally, the singlet BPS entropy would grow only like $1.65\,np$, which is not macroscopic in $n^2$.

**Method** (`src/quiver_index.py`, `scripts/quiver_index_dixon.py`).
- Only the index is needed, and $I_0=n(t=-1)$. Regrouping the dual-Cauchy sum by the irreps at each node gives
  $$I_0=\mathrm{Tr}\big[R^{(p)}(-1)^3\big],\qquad R^{(q)}=\sum_\lambda(-1)^{|\lambda|}L_\lambda^TR^{(q-1)}L_{\lambda^T},$$
  where $L_{\mu\lambda\nu}=c^\nu_{\mu\lambda}$ are Littlewood–Richardson coefficients (exact torus quadrature,
  rounding $\le2\times10^{-13}$). $R$ is the edge transfer matrix over $U(n)$ irreps, adding one flavour at a time.
- Exact integer arithmetic throughout.
- **Validation.** Reproduces every earlier value (dual-Cauchy route) and the direct weight counts at $n=2$. The total
  singlet counts from the same recursion at $t=+1$ reproduce 1,274, 889,040 and 2,115,632.

**Results** (`results/data/quiver_singlet_index_dixon.json`).

| $(n,p)$ | $I_0$ | Dixon$(np)$ | total singlets | $\ln\vert I_0\vert$ per mode | ln(singlets) per mode |
|---|---|---|---|---|---|
| (2,2) | 90 | 90 | 1 274 | 0.187 | 0.298 |
| (3,2) | −1 680 | −1 680 | 2 115 632 | 0.138 | 0.270 |
| **(4,2)** | **34 650** | **34 650** | $5.59\times10^{10}$ | 0.109 | 0.258 |
| (2,3) | 1 680 | −1 680 | 889 040 | 0.206 | 0.380 |
| (3,3) | 0 | 0 | $4.90\times10^{12}$ | — | 0.361 |
| **(4,3)** | **9 434 197 872** | 17 153 136 | $1.08\times10^{22}$ | 0.159 | 0.352 |
| (1,4) | 90 | 90 | 346 | 0.375 | 0.487 |
| (2,4) | 519 750 | 34 650 | $9.70\times10^{8}$ | 0.274 | 0.431 |
| (3,4) | 199 090 719 360 | 17 153 136 | $2.99\times10^{19}$ | 0.241 | 0.415 |
| (2,5) | 95 866 056 | −756 756 | $1.37\times10^{12}$ | 0.306 | 0.466 |
| (2,6) | 28 686 212 100 | 17 153 136 | $2.29\times10^{15}$ | 0.334 | 0.491 |
| (2,1), (3,1), (4,1) | 0 | −6, 0, 90 | | | |

**Findings.**
1. **$p=2$ (the hidden-$\mathcal N=4$ case) has a closed form, verified for $n=1,\dots,4$.**
   $$I_0(n,2)=(-1)^n\frac{(3n)!}{(n!)^3}=I_0(1,2n).$$
   - Growth $\sim27^n$: $\ln|I_0|$ per mode falls like $1/n$ (0.187 → 0.138 → 0.109). **Not macroscopic.**
   - Since $(2,2)$ is concentrated, the count there equals the index. If concentration persisted, the $p=2$ singlet
     BPS entropy would be $\approx3.3n$, against $6n^2$ fermion modes.
   - Unproved, and $n\ge5$ is not computed (the quadrature for $U(5)$ is too heavy). **[Proved for all $n$ on 2026-10-08: at $p=2$ the three Vandermondes and the edge factors combine into $|\Delta_{3n}|^2$; §8, D18.]**
2. **The general "depends only on $np$" pattern is false.** It holds only at $p=2$ and $n=1$; the $(2,3)$ match was a
   coincidence. Odd $np$ gives $I_0=0$ trivially: particle–hole symmetry maps $m\to pn^2-m$ with sign $(-1)^{pn^2}$.
3. **$p\ge3$ grows faster, consistent with macroscopic entropy; not established.**
   - $p=4$: $\ln|I_0|=4.50,\ 13.16,\ 26.02$ ($n=1,2,3$). The per-mode index entropy flattens, 0.375 → 0.274 → 0.241;
     the quadratic fit $\ln|I_0|\approx2.1n^2+2.4n$ gives an asymptote of 0.175 per mode. **[Superseded (§8): the true leading coefficient is $F^*(4)/12=0.213$ per mode; a fit through $n=1,2,3$ cannot see the $\frac12\ln n$ and $O(1)$ terms.]**
   - $p=3$: only even $n$ are informative. $n=2,4$ give 0.206 and 0.159 per mode. Two points cannot distinguish $n^2$
     growth with a large linear correction from slower growth; $(6,3)$ is out of reach.
   - At fixed $n=2$ the per-mode index entropy rises with $p$ (0.19 → 0.33), approaching the singlet-sector entropy
     (ratio 0.63 → 0.68): **large $p$ is macroscopic** (the SYK-like direction).
4. **Large-$n$ heuristic (derived here, not rigorous).**
   - **Setup.** At $t=-1$ the singlet index is the integral over $U(n)^3$ of $\prod_e\det(1-U_v\otimes U_{v+1}^\dagger)^p$.
     In the moments $u_{v,m}=\mathrm{Tr}\,U_v^m$, the Haar measure contributes $e^{-\sum_m|u_{v,m}|^2/m}$
     (Diaconis–Shahshahani), and each edge contributes $e^{-p\sum_m u_{v,m}\bar u_{v+1,m}/m}$. The quadratic kernel is
     $K_m=\tfrac1m(I+pP)$, with $P$ the 3-cycle permutation.
   - **Criterion.** Its Hermitian part has eigenvalues $1+p$, $1-p/2$, $1-p/2$. So the uniform (confined) eigenvalue
     saddle is **stable for $p<2$, marginal at $p=2$, unstable for $p>2$**. There is no fugacity suppression at
     $t=-1$, so the instability affects every $m$ equally.
   - **Analogy.** This is the Hagedorn/deconfinement criterion of Aharony–Marsano–Minwalla–Papadodimas–Van Raamsdonk
     (hep-th/0310285; Chen 2025 ref. [29]). An unstable uniform saddle signals $\ln I_0\propto n^2$.
   - **Consistency with the data.** It matches $p=1$ (index 0), $p=2$ (marginal, $e^{O(n)}$ with an exact closed
     form) and $p\ge3$ (super-linear).
   - **Not done.** Solving the deconfined saddle for the coefficient of $n^2$, to compare with $\approx2.1$ at $p=4$ and
     $\approx1.0$–$1.3$ at $p=3$. **[Done in §8: $F^*(4)=2.5550$, $F^*(3)=1.2373$.]**

**Consequences for the project.**
- The $p=2$ quiver, which has hidden $\mathcal N=4$, concentration at $n=2$, and an exact index, has sub-macroscopic
  BPS entropy at large $n$. It is a structured, solvable-looking corner, not a black-hole model.
- **$p=3$ remains the candidate.** It has $\mathcal N=2$, chaotic multiplet statistics, concentration at $n=2$, small
  near-BPS gaps, and index growth consistent with $e^{cn^2}$.
- Concentration at $n\ge3$ for $p=3$ is beyond exact methods; the bootstrap (or a saddle-point analysis of
  refined indices) is the remaining route.

## 8. The large-$n$ saddle of the singlet index: the coefficient of $n^2$ (2026-10-08)

**Question.** What is the coefficient of $n^2$ in $\ln|I_0(n,p)|$? Does the saddle that produces it also account for
the exact finite-$n$ values, in particular $(4,3)$? Full derivation: `docs/derivations.md` D18. Code:
`src/quiver_saddle.py`, `scripts/quiver_index_saddle.py`, `scripts/plot_quiver_index_saddle.py`. Data:
`results/data/quiver_index_saddle.json`. Figure: `results/figures/quiver_index_saddle.{pdf,png}`.

**Exact results from the eigenvalue form (proved).**
- *Representation.* Writing $1-e^{i\varphi}=2\sin(\varphi/2)e^{i(\varphi-\pi)/2}$, the phases telescope around the
  3-cycle. $I_0(n,p)$ becomes $(-i)^{3pn^2}/(n!)^3$ times the integral of a three-species log gas on the circle: like
  eigenvalues repel with weight 2, unlike ones with weight $p$. An extra sign $\sigma^p$ appears for odd $p$.
- *$p=2$.* All weights are 2, so the integrand is $|\Delta_{3n}|^2$ and $I_0(n,2)=(-1)^n(3n)!/(n!)^3$ for every $n$.
  This was conjectured in §7 from $n\le4$.
- *Signs and zeros.* For even $p$, $I_0\ne0$ with sign $(-1)^{pn/2}$. For odd $p$ and odd $n$, $I_0=0$.
- *Hölder.* For odd $p$, $|I_0(n,p)|^2\le|I_0(n,p-1)||I_0(n,p+1)|$.
- *CUE form.* $I_0(n,p)=(-1)^n\frac{(3n)!}{(n!)^3}\,\mathbb E_{\mathrm{CUE}(3n)}\prod(1-x_{vi}/x_{v+1,j})^{p-2}$.
- *Against the data.* All 15 exact indices satisfy these rules, as do the $n=1$ closed forms (Dixon).

**Large-$n$ saddle.**
- *Structure.* For $p>2$ the maximum of the log-gas energy has the three nodes' eigenvalues on three disjoint arcs,
  rotated by $2\pi/3$ (figure, panel d). The holonomies sit near $U_v\approx e^{2\pi iv/3}$, which turns the $(-1)^F$
  weight of each bifundamental into a phase $e^{\mp2\pi i/3}$ instead of cancelling it. The singlet count's saddle
  has all three nodes coincident on one arc.
- *Solution.* The continuum problem is solved to 12 digits: Chebyshev Euler–Lagrange solver with a soft-edge root.
- *Cross-check.* Independently, Newton maximisation of the discrete problem up to $n=2048$ extrapolates to the same
  numbers ($10^{-9}$–$10^{-13}$). Its $O(n)$ term matches the predicted $\beta=2$ self-energy to six digits. Random
  unconstrained starts never beat the symmetric maximum ($n\le12$).

| $p$ | $F^*$: $\ln\vert I_0\vert\sim F^*n^2$ | $G^*$: $\ln(\text{singlets})\sim G^*n^2$ | $F^*/G^*$ | $F^*$ per mode ($/3p$) |
|---|---|---|---|---|
| 2 | 0 | 1.4244 | 0 | 0 |
| **3** | **1.2373** | 3.0427 | 0.41 | 0.137 |
| 4 | 2.5550 | 4.7688 | 0.54 | 0.213 |
| 5 | 3.9297 | 6.5623 | 0.60 | 0.262 |
| 6 | 5.3455 | 8.4015 | 0.64 | 0.297 |
| $p\to\infty$ | $\frac{3p}2\ln3-\frac32\ln\frac{2p}3-\frac94$ | $3p\ln2-\frac32\ln\frac p2-\frac94$ | 0.79 | $\frac12\ln3$ |

**Limiting checks.** Near $p=2$, $F^*\ge(p/2-1)\,39\zeta(3)/2\pi^2$ (uniform arcs of length $2\pi/3$), tight as
$p\to2^+$. The large-$p$ forms come from semicircular arcs of radius $\sqrt{6/p}$ (index) and $\sqrt{8/p}$ (count).
The continuum values approach both.

**Finite $n$ and the $(4,3)$ test.**
- *The estimate.* Gaussian fluctuations about the discrete saddle, minus the universal $\beta=2$ crystal excess
  $3n\ln(e/\sqrt{2\pi})$. That correction is exact on the circle; with it the $p=2$ formula is reproduced to
  $-1/(36n)$. The estimate has no free parameters.
- *Against the data.* It reproduces all 27 nonzero exact values (indices and counts, $n\le4$, $p\le6$) to within
  $-0.007$ to $-0.60$ in the logarithm. Every residual is negative.
- *At $(4,3)$:* $F^*n^2=19.80$; corrected Laplace $22.69$; exact $22.97$. So $|I_0(4,3)|=9.43\times10^9$ is
  reproduced to 25%. The leading term alone misses by a factor of 24.

**Large-$n$ form and predictions** (index, $p=3$; corrected Laplace, so heuristic).
- $\ln|I_0(n,3)|\approx1.2373\,n^2+\frac12\ln n+2.4$. The fit gives $0.54\ln n$ at $p=3$ and $0.50$ for $p\ge4$.
- Predicted $\ln|I_0(n,3)|$ is $47.7$, $82.6$, $127.3$, $181.8$ at $n=6,8,10,12$. Going by the residuals at
  $n\le4$, the exact values should be about $0.1$–$0.6$ higher.
- These are targets for a Monte Carlo or improved exact method. The current recursion cannot reach $(6,3)$.

**What this does and does not establish.**
- It establishes that the singlet index, and therefore the number of singlet BPS states, grows like $e^{F^*n^2}$ with
  $F^*>0$ for $p\ge3$.
  - Even $p$: rigorous up to standard log-gas large deviations.
  - Odd $p$: also needs the sign $\sigma^p$ not to cancel the leading saddle. The arc gap at $p=3$ is only 0.039 rad,
    so sign-changing configurations are suppressed only by $e^{-O(n)}$ with a modest rate. The exact $(2,3)$ and
    $(4,3)$ show no cancellation.
- It does not address R-charge concentration. $|I_0|$ equals the BPS count only if there are no cancellations between
  degrees, for example under concentration. That is verified at $(2,2)$ and $(2,3)$ but open for $n\ge3$.
- Earlier estimates superseded: the §7 fit asymptote 0.175 per mode at $p=4$ (true 0.213); "$\alpha\approx1.0$–$1.3$"
  at $p=3$ (true 1.2373).

**Next.**
- *Refined-index saddle.* Add an R-charge fugacity $y$, giving a complex saddle. Concentration at half filling
  predicts $\ln|I(y)|=\ln|I_0|+\frac{pn^2}2\ln|y|$ exactly, i.e. zero curvature in $\ln|y|$. A nonzero curvature at
  order $n^2$ would show that the index-weighted R-charge distribution has width of order $n$. This is the large-$n$
  concentration test.
- *Monte Carlo.* Thermodynamic integration in $p$ from the exact $p=2$ point (positive integrand for even $p$; sign
  average for odd $p$). It would test the predictions at $n=6$–$12$.
