# Bena, Berkooz, de Boer, El-Showk, Van den Bleeken (2012) — Scaling BPS Solutions and pure-Higgs States

**File:** `papers/bena_berkooz_de_boer_el-showk_van_den_bleeken_2012.pdf` (37 pp.; §§1–3 read in full, §4 and
appendices skimmed, 2026-10-07)

## Full citation

Iosif Bena, Micha Berkooz, Jan de Boer, Sheer El-Showk, Dieter Van den Bleeken, *Scaling BPS Solutions and
pure-Higgs States*, arXiv:1205.5023v1 [hep-th] (22 May 2012); JHEP 11 (2012) 171.

## Main question

In type II on a Calabi–Yau, BPS bound states are multicentre supergravity solutions at strong coupling and
Coulomb- or Higgs-branch states of a quiver quantum mechanics at weak coupling. For a three-centre quiver with a
closed loop, which Higgs-branch states map to Coulomb-branch (supergravity) states, and what are the rest?

## Physical system

- $\mathcal N=4$ quiver QM (Denef's), with three abelian nodes ($U(1)$; primitive charges, so $m=1$, p. 9).
- $a=\Gamma_{12}$, $b=\Gamma_{23}$, $c=\Gamma_{31}$ **bosonic** chiral multiplets (with fermion partners) on the
  three edges. These multiplicities are the "flavour" numbers.
- Generic cubic superpotential $W=w_{\alpha\beta\gamma}\phi^\alpha_{12}\phi^\beta_{23}\phi^\gamma_{31}$ (3.1).
  D-terms (3.2), F-terms (3.3).
- For FI parameters $\theta_1,\theta_2<0$ the Higgs branch is $\phi_{31}=0$ with $c$ quadratic constraints on
  $\mathbb{CP}^{a-1}\times\mathbb{CP}^{b-1}$: a complete intersection $\mathcal M^c_{ab}$ (§3.1).

## Important definitions

- BPS states on the Higgs branch are the cohomology of $\mathcal M^c_{ab}$. Lefschetz $SU(2)$ is identified with
  spatial angular momentum.
- **Pure-Higgs states**: the extra middle-cohomology classes $\beta(a,b,c)$ from the Lefschetz hyperplane theorem
  (3.4)–(3.5), which are not inherited from the ambient space. All are Lefschetz singlets, i.e. zero angular
  momentum.
- **Coulomb states** $N(a,b;c)$ (3.15): the inherited classes, in one-to-one correspondence with
  Coulomb-branch/supergravity states.
- Index $\Omega=\mathrm{Tr}(-1)^{2J_z}=(-1)^{a+b+c}\chi$ (3.16)–(3.17); $\Omega=(-1)^{a+b+c}N+\beta$ (3.18).

## Main assumptions

Generic $w_{\alpha\beta\gamma}$; abelian nodes; a decoupling limit in §4 (appendix C).

## Main analytical results

1. **Euler characteristic** generating function: (3.20), (3.26).
2. **Pure-Higgs generating function** (3.27):
   $$Z_\beta=\frac{x^2y^2z^2}{(1-xy)(1-xz)(1-yz)(1-xy-yz-zx-2xyz)} .$$
   It is symmetric in $a,b,c$ although the Higgs branch is not.
3. $\beta\ne0$ iff $a+b-2\ge c\ge2$ and cyclically (3.32), i.e. iff the Coulomb branch has a scaling point.
4. Combinatorics (3.28)–(3.31): $\beta$ is a convolution of 3-derangement numbers $D(a,b,c)$ with an
   even-triangle delta function.
5. **Exponential growth** (3.43): $\beta(a,b,c)\sim\frac2\pi\sqrt{\frac{abc(ABC)^3}{(aA+bB+cC)^7}}\frac{a^ab^bc^c}{A^AB^BC^C}2^{a+b+c}$
   with $A=-a+b+c$ etc.; $\beta(a,a,a)\sim\frac{2}{3^{7/2}\pi}\frac{2^{3a}}{a}$ (3.45).
6. §4: Higgs–Coulomb map; the pure-Higgs states lie in its kernel. The scaling point is "unreliable", and the Higgs
   branch has a mass gap that violates superconformal invariance at low energy.

## Important equations

(1.3)–(1.5), (3.1)–(3.5), (3.15)–(3.20), (3.26)–(3.32), (3.43)–(3.45).

## Numerical methods

None (generating functions; asymptotics of multivariate meromorphic generating functions, (3.33)–(3.34)).

## Relevant figures/results

Figs. 2–4 (Lefschetz structure of the cohomology).

## Limitations

- Three abelian nodes only, so no large-rank gauge group.
- Bosonic Higgs-branch description; the pure-Higgs interpretation (single-centre black hole versus fuzzball) is
  left open.

## Relationship to our project

1. **Same quiver, different theory.** Their quiver is our three-node quiver with cyclic cubic superpotential, but
   with bosonic chirals and abelian nodes. Their $a,b,c$ correspond to our flavour numbers per edge; we take
   $a=b=c=p$ and non-abelian rank $n$. Our model is purely fermionic and has no Higgs or Coulomb branch, and its BPS
   states are $Q$-cohomology of a finite Fock space. **No equivalence is claimed.**
2. **A concrete hook.** They write (p. 8) that the combinatorial form of $Z_\beta$ "strongly suggests some elegant
   combinatorial origin, perhaps related to fermionic degrees of freedom on strings stretched between the
   centers."
   - **Cheap test (proposed, not done):** our quiver at $n=1$ with $(a,b,c)$ fermionic flavours,
     $Q=\sum w_{\alpha\beta\gamma}a^\alpha b^\beta c^\gamma$ on $\Lambda^\bullet(\mathbb C^{a+b+c})$, has at most
     about 18 modes.
   - **Comparison:** set its BPS counts by sector against $\beta$, $N$ and $\Omega$.
   - **A first check already fails.** The naive singlet Euler characteristic $\sum_m(-1)^m\binom am\binom bm\binom cm$
     is $-6$ at $(2,2,2)$, against $\Omega=2$ and $\beta=1$. Any match must be subtler than "singlet sector equals
     pure-Higgs".
3. **What it supports.** Exponentially many zero-angular-momentum BPS states in a three-node cyclic quiver are the
   precedent for reading our concentrated singlets as black-hole-like.
