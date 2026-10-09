# The $U(n)^3$ fermionic quiver as a project: literature check and fortuity test (2026-10-07)

**Progress report (2026-10-08):** `research/tex/quiver_progress_report.tex` → `research/pdfs/quiver_progress_report.pdf` (12 pp.; model, literature, all exact and numerical results of §§1–8, prospects and next steps). New figure for it: `results/figures/quiver_concentration.{pdf,png}` (`scripts/plot_quiver_concentration.py`): singlet dimensions and BPS states per degree, and the lowest non-BPS singlet energy per degree, at $(2,2)$ and $(2,3)$.

## Status summary (updated 2026-10-08, evening; details in §§1–12)

| property | status | evidence |
|---|---|---|
| Escapes the adjoint-model obstruction | **derived** | bifundamentals have no zero weights ($V_0=0$), so the window-widening mechanism of Chen-type models is absent; this does not by itself prove concentration |
| Concentration (singlet sector) | **verified** at $(2,2)$ (exact ranks); **numerical** at $(2,3)$ | 90 at $k=12$; 1680 at $k=18$ (§5) |
| Concentration at $n=3$ | **$p=2$: numerical** (all 1680 at $k=27$); **$p=3$: partial** (no BPS singlets at $k\le21$ or $k\ge60$; open for $24\le k\le36$ and mirror) | $(3,2)$: Schur–Weyl Lanczos, none at $k\le24$ (§11). $(3,3)$: exact for $k\le12$; level-2 bootstrap with Farkas certificates at $k=18,21$ (§12, D21) |
| $p=1$ | trivial | every per-irrep index is $\pm1$; singlet index 0 |
| Macroscopic BPS entropy | **$p=2$: no, proved** ($I_0=(-1)^n(3n)!/(n!)^3$ for all $n$, D18); **$p\ge3$: yes, at the level of the index** ($\ln|I_0|=F^*(p)n^2+O(\ln n)$, $F^*(3)=1.2373$, $F^*(4)=2.5550$; 41% of the singlet log-count at $p=3$; rigorous up to standard log-gas large deviations for even $p$, plus a sign proviso for odd $p$) | §§7–8, D18 |
| Fortuity | **$(2,2)$: all 90 classes, all orders**; **$(2,3)$: all 1680 classes (numerically certified)**; ordinary sense | $(2,2)$: first-order obstruction (§3) and $h^{12}_{\rm sing}(3,2)=0$ (§11); $(2,3)$: $h^{18}_{\rm sing}(3,3)=0$ by bootstrap certificate (§12); Choi–Choi–Kim refinement not checked |
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
- **[Withdrawn 2026-10-08, §9: no R-charge refinement is protected in the singlet sector, so this does not test
  concentration.]** *Refined-index saddle.* Add an R-charge fugacity $y$, giving a complex saddle. Concentration at half filling
  predicts $\ln|I(y)|=\ln|I_0|+\frac{pn^2}2\ln|y|$ exactly, i.e. zero curvature in $\ln|y|$. A nonzero curvature at
  order $n^2$ would show that the index-weighted R-charge distribution has width of order $n$. This is the large-$n$
  concentration test.
- *Monte Carlo.* Thermodynamic integration in $p$ from the exact $p=2$ point (positive integrand for even $p$; sign
  average for odd $p$). It would test the predictions at $n=6$–$12$.

## 9. Next steps: what can test the remaining properties (2026-10-08)

**Correction to §8.** There is no protected R-charge refinement of the singlet index.
- Non-BPS pairs $(\psi,Q\psi)$ in degrees $3m$ and $3m+3$ contribute $(-1)^my^m(1-y)$ to $\mathrm{Tr}_{\rm sing}(-1)^ky^{k/3}$,
  so the trace is protected only at $y=1$. The refined indices $\mathrm{Tr}(-1)^Fe^{2\pi ijk/3}$ of $\mathcal N=2$ SYK all
  reduce to $I_0$ on singlets, since $k\equiv0$ mod 3 there.
- With generic couplings there is no flavour symmetry either, so $I_0$ is the only protected quantity of the gauged
  model.
- What the index can say about where concentration sits is its sign, $(-1)^{m_*}$, and that always agrees with half
  filling (Table in §7, D18). This is the singlet-sector version of Chang–Chen–Sia–Yang's argument that the index
  phase fixes the concentrated charge.
- The "refined-index saddle" of §8, the report draft, `docs/todo.md` and `docs/research_questions.md` (2026-10-08c) is
  therefore withdrawn as a concentration test. A protected refinement exists only for non-generic couplings with a
  flavour $U(1)$ (as in the two-flavour SYK model of Chang–Chen–Sia–Yang), which changes the model.

**Clarification of §3.** The $(2,2)$ fortuity result is complete in the ordinary sense, not just "leading order".
Non-vanishing of the first-order obstruction already shows $[z_0]\notin\mathrm{Im}\,\pi_{3\to2*}$. The remaining caveats
are floating-point ranks, one coupling, and the Choi–Choi–Kim refinement.

**Fortuity follows from concentration (derived here).**
- $\pi_{M\to N}$ preserves the degree, and $k_*=3pn^2/2$ moves by $3p(2n+1)/2\ (>3)$ per unit of rank.
- So if the rank-$(n+1)$ singlet cohomology vanishes in degree $k_*(n)$, no rank-$n$ class has a preimage and all are
  fortuitous. Concentration at rank $n+1$, or any window narrower than the shift, suffices.
- A sufficient computable condition is that $Q$ is injective on rank-$(n+1)$ singlets of degree $k_*(n)$.

**Reach of exact methods: singlet dimensions $n(k)$.** Computed with the new
`src/quiver_index.singlet_series_recursive`, via `scripts/quiver_singlet_series.py`, saved to
`results/data/quiver_singlet_series.json`. All checks pass: totals, indices, palindromy, and agreement with the
dual-Cauchy series at $n=2$.

| target | singlet dimension | feasible? |
|---|---|---|
| $(3,2)$, $k=12$: injectivity test for fortuity of the $(2,2)$ classes | 3 494 → 20 168 | trivially, with a singlet basis |
| $(3,2)$, half filling $k=27$ (neighbours $k=24,30$): concentration at the next rank | 558 304 (443 657) | yes, with a singlet basis and sparse $Q$ |
| $(2,3)$, $k=18$: the 1680 BPS vectors (BPS chaos, direct fortuity test) | 300 584 | yes, with a singlet basis (about 4 GB for the kernel) |
| $(3,3)$, $k=18$: injectivity test for fortuity of the $(2,3)$ classes | $1.2\times10^8$ → $1.2\times10^9$ | no (exactly) |
| $(3,3)$, $k=39,42$ (half filling 40.5) | $1.06\times10^{12}$ | no |
| $(2,4)$, half filling | $3.0\times10^8$ | no |
| $(4,2)$, $k=27$: fortuity of $(3,2)$ classes | $4.5\times10^7$ | hard |

So $(3,2)$ is the only next-rank point within exact reach, and only in the $p=2$ corner. At $p=3$, $n=2$ is the last
exact rank.

**Revised plan, in priority order.**
1. **Singlet-adapted toolkit.** An explicit singlet basis from dual-Cauchy components, with $Q$ as Pieri/LR maps, or
   multi-trace words with a Gram rank. It unlocks:
   - $(3,2)$: concentration at the next rank, the near-BPS gap at $n=3$ against $n=2$, multiplet statistics, hidden
     $\mathcal N=4$ at $n=3$, and an all-orders fortuity check;
   - $(2,3)$: the 1680 BPS vectors, giving BPS chaos (projected-operator statistics) and a direct fortuity test at
     $p=3$.
2. **Quiver singlet bootstrap.** Closed walks on the quiver as trace words, sector-resolved, with BPS constraints.
   Targets, easiest first:
   - a certificate of absence of BPS singlets at $(3,3)$, $k=18$, which would make all 1680 $(2,3)$ classes
     fortuitous by the argument above;
   - certificates for $k\le36$ at $(3,3)$, which would confine the BPS states to $k=39,42$ (concentration);
   - near-window gaps and their $n$-scaling, the super-Schwarzian test.
3. **Large-$p$ analytics: the quiver as gauged $\mathcal N=2$ SYK.**
   - Expectation: at large $p$ the holonomy saddle is the $\mathbb Z_3$ twist that makes the SYK index large, which
     explains $F^*/3p\to\frac12\ln3$.
   - With random $C_{abc}$ at fixed $n$, the flavour-melonic limit should give FGMS Schwinger–Dyson equations for
     three species, with a gauge correction of relative order $1/p$ ($\dim G/\dim V=1/p$).
   - Deliverables: the super-Schwarzian coupling, the near-BPS density and the chaos exponent at large $p$.
   - Check: reproduce the $1/p$ terms of $F^*(p)$ (D18).
   - Supporting fact: $F^*(p)$ is smooth for $p>2$, so there is no sign of a transition between $p=3$ and $\infty$.
   - Risk: melonic dominance with only $p^3$ independent couplings shared across colour indices needs checking.
4. **Cheap checks.** More coupling seeds at $(2,3)$ (concentration, chaos; about 8.5 GB peak at $k=15$); Monte Carlo
   tests of the finite-$n$ index predictions; the Choi–Choi–Kim generalised-monotone analysis.
5. **$p=2$ hidden $\mathcal N=4$** at $n=3$ (comes with item 1) and its derivation.

## 10. Step 3: the $n=1$ member as a large-$p$ test (2026-10-08)

**Model.**
- $n=1$ is three-species $\mathcal N=2$ SYK: $3p$ fermions, $Q=\sum C_{fgh}a^\dagger_fb^\dagger_gc^\dagger_h$, restricted to the
  charge-balanced blocks $N_a=N_b=N_c=m$ of dimension $\binom pm^3$.
- By D19 it is the large-$p$ limit of the whole family at fixed $n$, so it tests the large-$p$ expectation directly.
- Couplings: Gaussian, $\langle C^2\rangle=1/p^2$ (seed 1), so the effective SYK coupling is $O(1)$.
- Code: `src/quiver_n1.py` (single-species creation matrices; explicit sparse or matrix-free $Q$, agreeing to
  $10^{-15}$) and `scripts/quiver_n1_spectrum.py`.
- Data: `results/data/quiver_n1_spectrum_gauss_seed1*.json`.

**Results.**

| $p$ | $I_0$ | BPS states by block $m$ | edges $E_0$ by $\vert q\vert$ | largest pair: levels, $\langle r\rangle$ |
|---|---|---|---|---|
| 2 | $-6$ | 6 at $m=1$ only | 1.5: 1.062 | — |
| 3 | 0 | 2 at $m=1$, 2 at $m=2$ | 0 (central): 0.0097; 3: 2.267 | — |
| 4 | 90 | 90 at $m=2$ only | 1.5: 0.2464; 4.5: 2.946 | 63, $0.549\pm0.032$ |
| 5 | 0 | none | 0 (central): $5\times10^{-8}$; 3: 0.581; 6: 4.14 | 876 (central), $0.428\pm0.011$ |
| 6 | $-1680$ | 1680 at $m=3$ only | 1.5: 0.1118; 4.5: 1.683; 7.5: 5.27 | 3160, $0.5244\pm0.0051$ |
| 7 | 0 | none in $m\le2$ (central blocks not computed) | 3: 0.369; 6: 2.103; 9: 5.96 | 8919, $0.5291\pm0.0030$ |
| 8 | 34650 | none in $m\le3$ | 1.5: 0.0494; 4.5: 0.847; 7.5: 3.05; 10.5: 6.77 | 511, $0.515\pm0.012$ |
| 10 | $-756756$ | none in $m\le3$ ($m=4$ running) | 4.5: 0.731; 7.5: 2.388 | — |

**Findings.**
1. **Even $p$ concentrates exactly.** The BPS states sit only in the central block, and their number equals $|I_0|$
   at $p=2,4,6$. At $p=8$ there are none outside the centre.
2. **Odd $p$ has structure the index cannot see.**
   - $p=3$ has $2+2$ BPS states in the two central blocks (index $0$); $p=5$ has none.
   - *Structural reason (derived here).* For odd $pn^2$ the two central degrees are dual under the top-form pairing.
     $B(a,a')=\langle Qa,a'\rangle$ on the lower one satisfies $B(a,a')=(-1)^{k_1}B(a',a)$ with $k_1=3(pn^2-1)/2$, so it
     is antisymmetric for $pn^2\equiv3$ mod 4 and symmetric for $pn^2\equiv1$ mod 4. The top form is gauge invariant
     (the determinants cancel around the triangle), so this holds on singlets.
   - *Antisymmetric case:* the middle rank is even, giving the mod-2 relation
     $\sum_{k\le k_1}h^k\equiv\sum_{k\le k_1}n(k)$. At $(1,3)$, $Q_1$ is an antisymmetric $27\times27$ matrix of rank 24.
     The parity sums are even for $(1,3)$, $(1,7)$ and $(3,3)$, so no BPS state is forced in those cases.
   - *Symmetric case ($p=5$):* the central-pair energies are squares of the eigenvalues of a symmetric matrix, i.e.
     two superposed orthogonal sequences. That explains $\langle r\rangle=0.428$ (two superposed GOE sequences give
     about $0.42$).
   - The central $q=0$ pair is nearly gapless, consistent with the Turiaci–Witten $q=0$ sector ($\rho\propto1/\sqrt E$).
   - **Correction:** the expectation stated earlier in this session, that a zero index means no BPS states, is false
     ($p=3$).
3. **Chaos.** Every non-central multiplet sector with enough levels is in the orthogonal class:
   $0.5244\pm0.0051$ ($p=6$, 3160 levels) and $0.5291\pm0.0030$ ($p=7$, 8919), against 0.5307. Couplings are real.
4. **Edges do not yet follow Turiaci–Witten.**
   - $E_0(1.5)=1.062,\ 0.246,\ 0.112,\ 0.049$ at $p=2,4,6,8$ falls faster than the asymptotic Schwarzian $1/p$
     (local exponents 2.1, 1.9, 2.8).
   - The ratio $E_0(4.5)/E_0(1.5)=12.0,\ 15.1,\ 17.1$ drifts away from 9.
   - The $|q|=4.5$ edge, $2.95,\ 1.68,\ 0.85,\ 0.73$ at $p=4,6,8,10$, slows to a local exponent of about $0.7$ between
     $p=8$ and $10$, possibly approaching the asymptotic $1/p$. $E_0(1.5)$ at $p=10$ (block $m=4$, $9.3\times10^6$ states)
     is running.
   - With $3p\le24$ fermions this is pre-asymptotic. The $q^2$ law is not yet testable, as at $(2,3)$.

**Reading.** At $n=1$ the large-$p$ direction behaves qualitatively like SYK: exact concentration for even $p$,
BPS count equal to the index, and orthogonal-class chaos in every multiplet sector. The Schwarzian charge dependence
needs larger $N_{\rm eff}$ than exact diagonalisation reaches (the $p=10$ run is the last practical point).

## 11. Step 1: Schur–Weyl singlet basis and the $(3,2)$ next-rank laboratory (2026-10-08)

**Engine** (D20; `src/quiver_singlet_basis.py`; `scripts/quiver_singlet_spectrum_sw.py`).
- A singlet is a triple of $S_m$ group-algebra elements, one per node, in Young's orthogonal form. The finite-$n$
  relations are the condition $\ell(\lambda)\le n$, and fermion statistics are sign projections per edge.
- $Q$ is a Kronecker product of Gelfand–Tsetlin edge maps, and the Fock norm is diagonal.
- **Validation:**
  - singlet dimensions equal $n(k)$ in every degree at $(2,2)$, $(2,3)$ and $(3,2)$;
  - $Q^2=0$ to $10^{-13}$;
  - all 11 stored $(2,2)$ and $(2,3)$ pair spectra are reproduced level by level to $3\times10^{-15}$;
  - the 90 BPS states at $(2,2)$ are recovered.
- **Cost at $(3,2)$:** 191 sectors at half filling; one application of $Q$ takes 0.5 s at under 2 GB. Zero-weight
  Fock methods would need $5\times10^{10}$ states.

**$(3,2)$ results** (seed-3 integer couplings, as at $n=2$; `results/data/quiver_singlet_sw_n3_p2_seed3.json`).

| $k$ | singlets | method | BPS | lowest singlet level |
|---|---|---|---|---|
| 0 | 1 | dense | 0 | 2457 |
| 3 | 8 | dense | 0 | 2184 |
| 6 | 56 | dense | 0 | 1671.6 |
| 9 | 456 | dense | 0 | 1384.4 (doublet) |
| 12 | 3 494 | dense | 0 | 873.1 |
| 15 | 20 168 | Lanczos | 0 | 667.4 (doublet) |
| 18 | 81 280 | Lanczos | 0 | 328.0 |
| 21 | 229 544 | Lanczos | 0 | 207.4 (doublet) |
| 24 | 443 657 | Lanczos | 0 | 99.53 (lower member of the pair $(24,27)$, $\vert q\vert=1.5$) |

**Findings so far.**
1. **Fortuity of the $(2,2)$ classes to all orders.** The rank-3 singlet cohomology vanishes at $k=12$ (lowest level
   873). By the moving-window argument (§9), all 90 BPS classes at $(2,2)$ are fortuitous in the ordinary sense, to all
   orders. This upgrades §3, which was already complete at first order, with an independent route.
2. **Concentration at the next rank (numerical, established).** There are no BPS singlets at $k\le24$, and so none at
   $k\ge30$ by particle–hole symmetry.
   - Hence every singlet BPS state at $(3,2)$ sits at half filling, $k=27$, and there are exactly $|I_0|=1680$ of them.
   - This is the first concentration result at rank 3, and it holds where the adjoint models' window grew with the rank.
   - The near-BPS gap at $|q|=1.5$, relative to the vacuum energy $E_{\rm vac}=n^3\sum|C|^2$, falls from $0.131$ at
     $n=2$ ($95.56/728$) to $0.041$ at $n=3$ ($99.53/2457$). The $|q|=4.5$ edge falls from $0.275$ to $0.084$.
3. **The $p=2$ structure persists at $n=3$.**
   - Levels come in exact doublets at odd $k$ and not at even $k$, as expected from the flavour map $F$ with
     $F^2=(-1)^k$.
   - The 3087-level pair $(12,15)$ has $\langle r\rangle=0.386\pm0.006$, the Poisson value: further from random-matrix
     behaviour than at $n=2$.
   - **Hidden $\mathcal N=4$ at $n=3$: confirmed.** $\tilde Q=Q((\varepsilon\otimes\varepsilon\otimes\varepsilon)\bar C)$ satisfies
     $\{\tilde Q,\tilde Q^\dagger\}=H$, $\{Q^\dagger,\tilde Q\}=0$ and $\{Q,\tilde Q\}=0$ on singlets at every degree $k=3$–$24$ (all degrees through half filling), to
     $10^{-16}$–$4\times10^{-15}$. A structure seen at one rank now holds exactly at two, so a derivation for all $n$
     should exist (todo).

## 12. Concentration at $p=3$: the quiver singlet bootstrap (2026-10-08)

**Does concentration at $(3,3)$ need the bootstrap?** Yes, for degrees beyond the reach of exact methods.
- *Exact reach.* At $(3,3)$ the singlet dimensions are $n(k)=1,27,702,18\,981,460\,539$ for $k\le12$, then
  $8.7\times10^6$, $1.2\times10^8$, ..., up to $1.06\times10^{12}$ at $k=39,42$. The Schur–Weyl engine (§11) covers
  $k\le12$; $k=12$ takes 46 min of Lanczos at 0.67 GB. In that range there are no BPS singlets. The lowest levels are
  $6777,\ 6024,\ 4634.9,\ 3546.0,\ 2472.0$ (seed 3, `results/data/quiver_singlet_sw_n3_p3_seed3.json`).
- *The algebraic shortcut fails at $p=3$.* Specialising to diagonal couplings and using Künneth plus
  semicontinuity would bound the generic cohomology by the diagonal model's. But at $(2,3)$ the diagonal model is
  not concentrated: it has $\{12{:}15,\ 15{:}123,\ 18{:}1896,\ 21{:}123,\ 24{:}15\}$. So it gives no information off
  the middle.

**Method** (D21; `src/quiver_bootstrap.py`, `scripts/quiver_bootstrap_run.py`).
- The functional is real on the singlet sector of degree $k=3m$ (valid for real couplings), with $N_e=m$ on every
  edge.
- Cones:
  - adjoint(-projected) Gram matrices over open quiver words, blocked by endpoints and edge charges;
  - singlet Gram matrices over traced closed walks.
- Rows:
  - sector rows;
  - node-Casimir rows;
  - EOM rows;
  - BPS rows $\phi(XQ)=\phi(QX)=\phi(X\bar Q)=\phi(\bar QX)=\phi(XH)=\phi(HX)=0$.
- *Level 2* means open words of length $\le2$ and closed walks of length $\le3$, which gives cones of size
  $37,27,27,18,\dots$. Rows are generated by operators of length $\le4$. The resulting problem has
  $1.6\times10^5$ variables and $1.2\times10^5$ equality rows, needs about 3 GB, and solves in 1–9 min with Clarabel.
- *Verdicts:*
  - margin $t^*$ (max $t$ with all cones $\succeq t\mathbf 1$): clearly negative means excluded;
  - pure feasibility returning a verified Farkas certificate ($A^{\!\top}y\approx0$, $y\in K^*$, $b\cdot y<0$);
  - the energy lower bound without BPS rows: positive means excluded.
- *Validity:* every row family holds on exact functionals to $\le3\times10^{-15}$ (`tests/test_quiver_bootstrap.py`).

**Validation against exact singlet spectra** (seed-3 couplings, level 2;
`results/data/quiver_bootstrap_n2_p{2,3}_L2.jsonl`).

| $(n,p)$ | $k$ | exact $E_0$ | energy bound | margin $t^*$ | verdict |
|---|---|---|---|---|---|
| $(2,2)$ | 3 | 546 | 455.1 | $-1.334$ | excluded |
| | 6 | 199.95 | 182.6 | $-0.508$ | excluded |
| | 9 | 95.56 | 0 | $0$ | not excluded |
| | 12 | 0 (90 BPS) | 0 | $0$ | correct (BPS) |
| $(2,3)$ | 3 | 1506 | 1129.3 | $-1.342$ | excluded |
| | 6 | 575.23 | 253.7 | $-0.732$ | excluded |
| | 9 | 151.34 | 14.4 | $-0.163$ | excluded |
| | 12 | 0.3348 | 0 | $0$ | not excluded |
| | 15 | 0.0539 | 0 | $0$ | not excluded |
| | 18 | 0 (1680 BPS) | 0 | $0$ | correct (BPS) |

There are no contradictions: no BPS degree was excluded, and every energy bound lies below the exact value. Level 2
reaches half-way to the middle at both points.

**$(3,3)$ results** (seed-3 couplings
$C=[[[5,1,1],[2,1,5],[5,3,1]],[[1,2,3],[4,3,2],[1,4,4]],[[1,1,3],[2,5,3],[3,3,4]]]$, level 2;
`results/data/quiver_bootstrap_n3_p3_L2.jsonl`).

| $k$ | singlets | exact $E_0$ | energy bound | margin $t^*$ | Farkas certificate | verdict |
|---|---|---|---|---|---|---|
| 3 | 27 | 6024 | 5457.3 | $-5.894$ | | excluded |
| 6 | 702 | 4634.9 | 4139.9 | $-4.912$ | | excluded |
| 9 | 18 981 | 3546.0 | 2824.4 | $-3.956$ | | excluded |
| 12 | 460 539 | 2472.0 | 1511.2 | $-3.012$ | | excluded |
| 15 | $8.7\times10^6$ | — | 598.2 | $-2.102$ | | excluded |
| 18 | $1.2\times10^8$ | — | 218.6 | $-1.219$ | $\|A^{\!\top}y\|/\|y\|=2.3\times10^{-11}$, $b\cdot y/\|y\|=-4.0\times10^{-3}$ | **excluded** |
| 21 | $1.2\times10^9$ | — | 21.8 | $-0.374$ | $6.1\times10^{-12}$, $-1.2\times10^{-3}$ | **excluded** |
| 24 | $8.4\times10^9$ | — | 0 | $0$ | none | not excluded (also with the D21.4b rows) |
| 27–39 | up to $1.06\times10^{12}$ | — | 0 | $0$ | — | not excluded (controls near the middle) |

- The margin falls almost linearly in $k$, with steps 0.98, 0.96, 0.94, 0.91, 0.88, 0.85, and crosses zero between
  $k=21$ and 24.
- The two certified cells were also obtained independently by the margin optimisation; the energy bound is positive
  there too.
- The dual cone blocks of both certificates are PSD (smallest eigenvalue $\ge0$).

**Figure.** `results/figures/quiver_bootstrap.{pdf,png}` (`scripts/plot_quiver_bootstrap.py`).
- Panel (a): the level-2 BPS-exclusion margin $t^*$ against the filling fraction $k/(3pn^2)$ at $(2,2)$, $(2,3)$ and
  $(3,3)$. Filled markers are excluded degrees and hollow markers are not.
- Panel (b): the lowest singlet energy per degree in units of $E_{\rm vac}=n^3\sum C^2$, comparing level-2 lower
  bounds (dashed) with the exact values (solid) where known.

**Consequences.**
1. **Fortuity of all 1680 $(2,3)$ classes (numerically certified).**
   - The rank-3 singlet cohomology vanishes at $k=18$ for the seed-3 couplings.
   - By upper semicontinuity of $h^k$ in $C$ (D21.6) it also vanishes for generic couplings.
   - By the moving-window argument of §9 ($\pi_{3\to2}$ preserves the degree), no $(2,3)$ BPS class at $k=18$ has a
     rank-3 preimage. All are fortuitous in the ordinary sense, completing the $p=3$ analogue of §11's $(2,2)$ result.
   - *How strong the certificate is.* Suppose a BPS singlet existed. Its real functional $x^*$ would satisfy
     $Ax^*+s=b$ with $s\in K$, and $|\phi(m)|\le\|m\|\le n^{L(m)}\le3^8$ for every monomial. So
     $\|x^*\|_2\le3^8\sqrt{1.6\times10^5}\approx2.6\times10^6$. Then
     $0\le\langle y,s\rangle=b\cdot y-(A^{\!\top}y)\cdot x^*\le b\cdot y+\|A^{\!\top}y\|\,\|x^*\|$. Numerically the
     right-hand side is $(-4.0\times10^{-3}+6\times10^{-5})\|y\|<0$ at $k=18$ and
     $(-1.16\times10^{-3}+1.6\times10^{-5})\|y\|<0$ at $k=21$, a contradiction with a factor of about 70 to spare.
     This is a computer-assisted argument up to floating-point evaluation of $A^{\!\top}y$, $b\cdot y$ and the
     dual-block eigenvalues.
   - **Exact verification (2026-10-09, D21.8).** Both certificates were re-checked in exact integer arithmetic.
     Coefficients were placed on the grid $1/(3\cdot2^{20})$ (largest deviation $3\times10^{-8}$). The certificate
     was rounded to a $2^{-59}$ ($2^{-57}$) grid, and each dual block was shifted by at most $9\times10^{-13}$
     relative until positive definite, verified exactly by Bareiss elimination. The componentwise bound
     $\sum_j|r_j|n^{L_j}$ then comes to $1.03\times10^{-7}$ ($2.38\times10^{-7}$) of $|\mu\cdot b|$ at $k=18$ ($21$),
     so the floating-point estimate above was conservative by a factor of about $10^5$. An independent SCS solve at
     $k=18$ gives a different certificate, also verified exactly (ratio $4.4\times10^{-4}$). Certificates and records:
     `results/data/quiver_certificate_n3_p3_k{18,21}_clarabel.npz`, `..._k18_scs.npz`, `quiver_certificate_checks.jsonl`.
   - Caveats:
     - the proof rests on the correctness of the constraint generation (trace algebra and row families), which is
       tested against explicit Fock-space operators but not proved;
     - the $(2,3)$ count of 1680 is itself numerical (§5);
     - the Choi–Choi–Kim refinement is not addressed.
2. **Partial concentration at $(3,3)$.**
   - There are no BPS singlets at $k\le21$, and so, by particle–hole symmetry $h^k=h^{81-k}$, none at $k\ge60$.
   - Open: $k=24,\dots,36$ and their mirrors $45,\dots,57$, besides the middle pair $39,42$.
   - Since $I_0(3,3)=0$, the middle degrees carry equal numbers of states, possibly zero.
   - Concentration at $(3,3)$ is therefore **not** established. The bootstrap at this level covers $k\le21$ of the
     $k\le36$ needed.

**Why the frontier sits half-way: a diagnosis at $(2,2)$** (`results/data/quiver_bootstrap_frontier_diagnostics.txt`).
- *Non-singlet BPS states.* In the full Fock sector $N_e=(m,m,m)$, with all gauge irreps and no singlet projection
  (`scripts/quiver_sector_zero_modes.py`):
  - $m=2$ has no BPS states at all (lowest level 70.94);
  - $m=3$ has at least 12 zero modes, i.e. non-singlet BPS states, while the singlet gap there is 95.56.
  - So the level-2 frontier, excluded at $k=6$ and not at $k=9$, coincides with the onset of non-singlet BPS states
    at fixed $N_e$.
- *Gauge invariance barely enters the exclusions.* Removing the Casimir rows changes the margins only slightly
  ($-1.334\to-1.327$ at $k=3$, $-0.508\to-0.505$ at $k=6$), though it changes the energy bounds a lot (455 to 293,
  183 to 63). The exclusions are effectively whole-sector exclusions.
- *Nothing tried moves the $(2,2)$ frontier.* All of the following leave $k=9$ at $t^*\approx0$:
  - level 3 (cones up to 26);
  - Gauss rows;
  - finite-$n$ relations up to length 6;
  - sandwiched-Casimir rows (D21.4b, which state that $W\psi$ lies in a definite irrep for singlet $\psi$; they
    strengthen $k=6$ from $-0.51$ to $-0.77$);
  - these combined with level 3;
  - Gauss rows with longer words (open words $W$ up to length 4, prefactors $Y$ up to length 4; 2026-10-09).
  - **Level 4** (open words up to length 4, Gram blocks up to 104; 28 273 monomials, rows as at level 3; 2026-10-09).
    Clarabel reports Solved in 713 s with $b\cdot y/\|y\|=-8\times10^{-12}$, i.e. feasible, so $k=9$ is still not
    excluded. The control at the BPS degree $k=12$ is also feasible, as it must be (peak 14.0 GB under a 16 GB guard). SCS did not converge within 2.5 h, and its time limit was not honoured. This settles the yes-or-no
    question: raising the level from 2 to 3 to 4 does not move the frontier past the onset of non-singlet BPS states.
    In this formulation the trace bootstrap's exclusion power is effectively capped near filling $f_*\approx0.27$, and
    concentration near half filling at $n\ge3$ needs a different tool (§13 and the analytic plan).
  - A singlet cone over closed walks of length 4 exceeded the 6 GB guard and remains untested.
- *Reading.*
  - At low level the singlet constraint cannot overcome functionals that mimic non-singlet BPS states.
  - The reach of level 2 is about half-way to the middle at all three ranks: $6/12$, $9/18$, $21/40.5$. The
    single-node model behaved the same way.
  - Pushing concentration at $(3,3)$ past $k=21$ needs constraints that resolve the gauge representation of the
    state much more sharply than Casimir rows do, not merely more operators.

**The frontier as a function of $n$ at $p=3$ (2026-10-09; `scripts/quiver_bootstrap_frontier.py`,
`results/data/quiver_bootstrap_frontier_n{4,5,6}_p3.jsonl`).**
- *Method.* The level-2 SDP has the same size at every $n$; only its coefficients change. Each point is a scaled
  feasibility solve (D21, variable scaling) taking about 1 min. A coarse scan is followed by bisection and an exact
  certificate check at the last excluded degree (D21.8).

| $n$ | half filling | excluded through $k$ ($f=k/3pn^2$) | not excluded at $k$ ($f$) | exact check at the frontier |
|---|---|---|---|---|
| 2 | 18 | 9 (0.250) | 12 (0.333) | VERIFIED, ratio $4.5\times10^{-8}$ |
| 3 | 40.5 | 21 (0.259) | 24 (0.296) | VERIFIED, ratio $2.4\times10^{-7}$ |
| 4 | 72 | 39 (0.271) | 42 (0.292) | VERIFIED, ratio $2.0\times10^{-6}$ |
| 5 | 112.5 | 60 (0.267) | 63 (0.280) | VERIFIED, ratio $2.7\times10^{-6}$ |
| 6 | 162 | 87 (0.269) | 90 (0.278) | VERIFIED, ratio $1.1\times10^{-5}$ |

- *Reading.*
  - The frontier brackets shrink and are consistent with a constant filling fraction $f_*\approx0.27$
    ($0.269\le f_*<0.278$ at $n=6$).
  - At level 2, for every rank tested, singlet BPS states are excluded from filling fractions below about $0.27$, and
    above about $0.73$ by particle–hole symmetry. The BPS states, if any, occupy at most the central 46% of the
    filling range.
  - The large-$n$ convergence is an empirical extrapolation from five ranks, not a proof. Each finite-$n$ statement
    is certified.
- *Fortuity one rank up (rank tower).*
  - $(4,3)$ is excluded at $k=39$ (verified), so any $(3,3)$ BPS classes at $k=39$ are rank-fortuitous.
  - $k=42$ at $(4,3)$ is not excluded, so their partners at $k=42$ remain undecided.
  - In general the needed filling, $n^2/2(n+1)^2$, is $0.222$ for $2\to3$, $0.281$ for $3\to4$, $0.32$ for $4\to5$,
    and tends to $1/2$. With $f_*\approx0.27$ the level-2 bootstrap proves fortuity of half-filling classes only for
    $2\to3$, plus the lower half of the middle pair for $3\to4$.
- *$p=4$.* $(3,4)$ is excluded at $k=24$; the certificate is verified exactly (ratio $3.5\times10^{-7}$; a problem of
  $1.8\times10^6$ monomials at a peak of 9.7 GB under a 16 GB guard). So every $(2,4)$ singlet BPS class at half
  filling is rank-fortuitous. Concentration at $(2,4)$ is not established, so classes at other degrees, if any, are
  not covered. A full $p=4$ frontier scan needs one fresh process per point (`--subprocess`); an in-process scan was
  killed at the 10 GB guard.
- *Cost.* Each $n$ takes 12–19 min (8–10 solves) at a peak of 3.1–3.2 GB.

**Cost and reproducibility.**
- Level 2 at $(3,3)$: build 20 s, peak 3.0–3.2 GB, margin solves 6–9 min (with other jobs running), feasibility
  solves about 1 min.
- Runs used the uncommitted working tree on top of `bb78329`.
- A $\mathbb Z_3$ quiver-rotation reduction for cyclic couplings is implemented and verified to give identical
  optima. It is not needed at level 2.

## 13. Fortuity in which direction? Rank versus arrows (comparison with Chryssanthacopoulos–Vegh 2026) (2026-10-09)

**The comparison.** J. Chryssanthacopoulos and D. Vegh, *Fortuity and fragility in supersymmetric SYK* (2026;
`~/src/syk_gravity/papers/chryssanthacopoulos_vegh_2026.pdf`), test fortuity in a different direction from §§3, 9,
11, 12.
- Their main text is N = 2 SYK: the system grows by adding fermions, $N\to N+1$, at fixed charge.
- In the Discussion ("Fortuity and fragility beyond SYK", Table 2) they treat the Abelian three-node cyclic quiver with
  $m$ arrows per edge. The gauge group is U(1) at each node, and $m$ varies at fixed cohomological degree.
  - Ambient classes persist as $m$ grows (monotonous).
  - Variable (pure-Higgs) classes occur only at $m=n+2$ ("a fortuity window of width one").
  - They propose a "diagonal uplift" shifting the degree by one per unit of $m$.
- Notation: their degree $n$ ↔ our degree $k$ (or $m=k/3$); their arrows $m$ ↔ our flavours $p$. Our rank $n$ has no
  counterpart there (it is fixed at 1).

**Both are fortuity, relative to different families.** Fortuity is defined only relative to a family with
degree-preserving cochain maps.
- Our rank direction, $n\to n+1$ at fixed $p$, uses Tierz's projections. It is the planar large-$n$ limit, and here
  "supersymmetry relies on finite-$n$ trace relations" is literally the mechanism (Chang–Lin).
- The arrow direction, $p\to p+1$ at fixed $n$, uses the projection that deletes the new flavour, which is also a
  cochain map for nested couplings. It is the SYK direction, with the system growing by fermion modes, and its limit
  is the large-$p$ limit of D19.
- The two verdicts can differ for the same class.

**Why the directions behave differently: window speed against degree spacing.** Singlet degrees are spaced by 3.

| step | shift of half filling | consequence of concentration at the larger size |
|---|---|---|
| rank $n\to n+1$, fixed $p$ | $3p(2n+1)/2\ge9$ | one-step fortuity of every class (§9) |
| arrows $p\to p+1$, $n\ge2$ | $3n^2/2\ge6$ | also one-step |
| arrows $p\to p+1$, $n=1$ | $3/2$ | window two flavours wide; needs the actual maps |

The paper's SYK case and its quiver example are in the slow regime, which is why walks and diagonal uplifts arise
there. In our rank direction the window jumps past several degrees at every step.

**Our $n=1$ member is the fermionic version of their quiver example:** the same Abelian cyclic quiver with
$\mathrm{Tr}(ABC)$ and generic couplings, but with fermionic matter. Its singlet BPS counts by degree (§10) are the
counterpart of their Table 2:
- $k=3$: 6 at $p=2$, 2 at $p=3$;
- $k=6$: 2 at $p=3$, 90 at $p=4$;
- $k=9$: 1680 at $p=6$;
- $k=12$: 34650 at $p=8$;
- zero in every other cell computed.

The window is two flavours wide, as the table predicts. Unlike the bosonic quiver, there is no ambient
(monotonous) row: every BPS state sits at the moving middle, as in generic N = 2 SYK. The same holds in the rank
direction at all sizes computed. In the fermionic model fortuity covers the whole singlet BPS sector rather than
singling out a subsector.

**Resolution (adopted).**
1. *Rank stays the defining test for this project.* The claims are about the planar large-$n$ limit, and every
   fortuity statement so far is rank fortuity at fixed $p$:
   - (2,2) via (3,2);
   - (2,3) via the exact (3,3) certificate;
   - (3,3) classes at $k=39$ via (4,3);
   - (2,4) half-filling classes via (3,4).

   They are labelled as such.
2. *The arrow direction is the paper's notion, and is computed where it carries new content: $n=1$.* At $n\ge2$ it
   follows in one step from concentration. For example, the 90 (2,2) classes are arrow-fortuitous because (2,3) has no
   singlet cohomology at $k=12$.
3. *At $n=1$ we compute the paper's two uplift channels exactly* (`scripts/quiver_n1_arrow_tower.py`, nested
   Gaussian couplings $C^{(p)}=C[:p,:p,:p]$; results below):
   - horizontal: $\pi_*:H^k(p{+}1)\to H^k(p)$, deleting the new flavour;
   - diagonal: $\iota_*:H^k(p)\to H^{k+3}(p{+}1)$, multiplying by $a_p^\dagger b_p^\dagger c_p^\dagger$. The fully
     occupied states form a sub-complex, so this is a cochain map.

   The SYK statement "every class uplifts through at least one channel" becomes $\ker\iota_*\subseteq\mathrm{im}\,\pi_*$.

**$n=1$ results, $p\le5$** (`results/data/quiver_n1_arrow_tower_pmax5.json`). Couplings are nested Gaussian, seed 1.
The BPS counts reproduce §10. $\pi Q=Q\pi$ exactly (residual 0.0). $\iota$ anticommutes,
$Q^{(p+1)}\iota=-\iota Q^{(p)}$ exactly, because $Q$ passes the three fermions $a_pb_pc_p$; it therefore still maps
cohomology to cohomology.

| $p$ | $k$ | $\dim H^k(p)$ | horizontal: $\dim H^k(p{+}1)$, rank $\pi_*$ | diagonal: $\dim H^{k+3}(p{+}1)$, rank $\iota_*$ | lift through neither |
|---|---|---|---|---|---|
| 2 | 3 | 6 | 2, 2 | 2, 2 | **2** |
| 3 | 3 | 2 | 0, 0 | 90, 2 | 0 |
| 3 | 6 | 2 | 90, 2 | 0, 0 | 0 |
| 4 | 6 | 90 | 0, 0 | 0, 0 | **90** |
| 5 | — | 0 | — | — | — |

- **At odd $p$ the SYK pattern holds.** At $p=3$ every class uplifts through exactly one channel: the lower central
  degree diagonally, the upper horizontally. These are the paper's walks.
- **From even $p$ it fails.** At $p=2$, 2 of the 6 classes uplift through neither channel; the dimension count forces
  this, since $6>2+2$. At $p=4$ none of the 90 do, because $p=5$ has no BPS states at all.
- **Reason.** Adding a flavour adds three modes coupled by the singlet constraint. The partially occupied pieces of
  the filtration by new-flavour occupation (one or two of $a_p,b_p,c_p$ occupied) leave non-singlet remainders. The
  SYK long exact sequence becomes a spectral sequence, and its two extreme channels are not exhaustive.
- **In the arrow direction at $n=1$, every class is fortuitous within one step.** At fixed $k$ the BPS states live at
  $p=2k/3-1$ and $2k/3$ only.
- **So the fermionic quiver at $n=1$ is neither like the bosonic quiver nor like SYK.** It has no monotonous ambient
  classes, and SYK's guaranteed one-step uplift fails.
- The $p=6\to7$ row (the 1680 classes at $k=9$) needs the $p=7$ central cohomology. Those blocks have 42 875 states,
  and their lowest multiplet may sit below floating-point resolution (the central level was 0.0097 at $p=3$ and
  $5\times10^{-8}$ at $p=5$). That row is pending; exact modular ranks may be required.
