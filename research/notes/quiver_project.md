# The $U(n)^3$ fermionic quiver as a project: literature check and fortuity test (2026-10-07)

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
