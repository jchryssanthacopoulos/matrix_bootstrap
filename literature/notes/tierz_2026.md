# Tierz (2026) — BPS spectra of $\mathrm{Tr}[\Psi^p]$ matrix models for odd $p$

**File:** `papers/tierz_2026.pdf`

> **NOTATION WARNING.** Tierz's $p$ is the **degree of the supercharge**, i.e. our $q$. He has a *single*
> fermion matrix, i.e. our $p=1$. Throughout this note we use **his** convention inside quoted results and
> flag the translation explicitly where it matters. In our variables this paper is the
> $(\text{flavours}=1,\ q\ \text{odd},\ N)$ axis of our family.

## Full citation

Miguel Tierz, *BPS spectra of $\mathrm{Tr}[\Psi^p]$ matrix models for odd $p$*, arXiv:2604.27164 (submitted 29
April 2026, revised 25 August 2026). Shanghai Institute for Mathematics and Interdisciplinary Sciences (SIMIS).
Same author as `papers/tierz_2017.pdf`. Found 2026-09-26; postdates the rest of the project's literature review.

## Main question

What is the BPS cohomology of a single $N\times N$ matrix of complex fermions with supercharge
$Q_p=\mathrm{Tr}(\Psi^p)$ for **odd $p\ge3$**, generalising Chen's cubic case? Specifically: the BPS
multiplicities at finite $N$, the factorisation structure of the BPS generating function, which classes are
fortuitous, the large-$N$ growth rate, and which parts of the Hamiltonian carry information beyond the quadratic
Casimir.

## Physical system

- One $N\times N$ matrix $\Psi_{ij}$ of fermion creation operators, $\{\Psi_{ij},\Psi^\dagger_{kl}\}=\delta_{ik}\delta_{jl}$.
- $Q_p=\mathrm{Tr}(\Psi^p)$, $p\ge3$ odd; $H=\{Q_p,Q_p^\dagger\}$ (eq. 1).
- **$SU(N)$ is a global symmetry and is explicitly *not* gauged**, so the full Fock space
  $\mathcal H=\Lambda^\bullet(\mathbb C^{N\times N})$ is retained. (Same convention as ours.)
- BPS space $=\ker Q_p/\mathrm{im}\,Q_p$ in each fermion-number sector $R$.

## Important definitions

- $h_R=\dim V_R-\mathrm{rank}\,Q_R-\mathrm{rank}\,Q_{R-p}$ (13); $r_R:=\mathrm{rank}\,Q_R$ (14).
- $Z^{(p,N)}_{\rm BPS}(x)=\sum_R h_Rx^R$; $q_{\min}$ = minimum BPS charge.
- Integral reduction $Q_p=p\widetilde Q_p$ (11), so the complex is defined over $\mathbb Z$ and $H=p^2\widetilde H$ (12).
- Orthogonal split $\Psi=A+\chi\mathbb 1/\sqrt N$ into traceless part $A$ and $U(1)$ trace fermion $\chi$.
- Projection $\pi_{M\to N}:\mathcal H_M\to\mathcal H_N$ deleting every monomial containing $\Psi_{ij}$ with $i>N$
  or $j>N$; induced maps $\Pi^{(R)}_{M\to N}$ on cohomology.
- Monotone subspace $B^{(p,N)}_{R,\rm mon}:=\bigcap_{M>N}\mathrm{Im}\,\Pi^{(R)}_{M\to N}$; a class is **fortuitous**
  if its coset in $B_R/B_{R,\rm mon}$ is non-zero.
- $K_k$ = the normal-ordered piece of $H/p^2$ with $k$ creation and $k$ annihilation operators.

## Main assumptions

- Odd $p$ throughout (even $p$ gives $\mathrm{Tr}\,\Psi^p\equiv0$).
- $p\le 2N-1$, stated as "the only range in which $Q_p$ can be nonzero".
- Positive BPS multiplicities are **numerical**, certified by agreement of modular ranks over three
  characteristics but not proved over $\mathbb Q$ (his §4.1 is explicit about this).

## Main analytical results

1. **Prop. 1** (traces of odd matrices): $\mathrm{Tr}(X^p)=0$ for even $p$ and $=p\widetilde T_p(X)$ for odd $p$,
   by a cyclic-orbit argument; every non-zero orbit has size exactly $p$. No invariant theory needed.
2. **Prop. 2**: $Q_p^2=0$ for every odd $p$, "using only the odd degree of $Q_p$, not its trace form."
   (This is our $\omega\wedge\omega=0$ observation, independently.)
3. **Thm. 1 (rank symmetry)**: the top-form pairing $\langle a,b\rangle=[a\wedge b]_{\Lambda^{N^2}}$ gives
   $\langle Q_Ra,b\rangle=(-1)^{pR}\langle a,Q_{N^2-p-R}b\rangle$, hence $r_R=r_{N^2-p-R}$ (15) and
   $R_{p,N}(x)=\sum_Rr_Rx^R$ is **palindromic** of degree $N^2-p$.
4. **Thm. 2 (exact $U(1)$ decoupling)**: $\mathrm{Tr}(\Psi^p)=\mathrm{Tr}(A^p)$ (16) for every odd $p\ge3$, so
   the complex splits as $\Lambda^\bullet(\mathbb C\chi)\otimes(\Lambda^\bullet W_0,Q_p)$ and
   $Z_{\rm BPS}$ acquires **one exact factor of $(1+x)$** from the free trace mode.
   *(Single-flavour only: with several flavours $\chi^a\chi^b\ne0$ and the argument fails.)*
5. **Prop. 6**: the **projection** $\pi_{M\to N}$ is a cochain map, $\pi Q^{(M)}_p=Q^{(N)}_p\pi$, and the
   $\Pi^{(R)}$ compose. He states explicitly that the **inclusion** $\mathcal H_N\hookrightarrow\mathcal H_{N+1}$
   is *generally not* a cochain map: "Even on the vacuum, $Q^{(N+1)}_p$ creates terms involving the new index that
   are absent from the image of $Q^{(N)}_p$." Fortuity is therefore defined by failure to lie in the image of the
   projections from above, not by an obstruction to lifting.
6. **Conjecture 1 (additional factors)**: for odd $p\ge5$ with $p\le2N-1$, $(1+x)^N\mid Z^{(p,N)}_{\rm BPS}(x)$
   (58), and after removing $x^{q_{\min}}$ and that factor the quotient lies in $p\mathbb Z_{\ge0}[x]$ (59), so
   $Z=p^{b_{p,N}}x^{q_{\min}}(1+x)^NT^{(p)}_N(x)$ (60). Only $(1+x)^2$ (even $N$) or $(1+x)^3$ (odd $N$) is
   proved, from Thm. 2 plus the Euler characteristic. **The rest is conjectural.**
7. **Thm. 4 (large $N$)**: $\lim_N N^{-2}\log\sum_a|I_a|=\log(2\cos(\pi/2p))$ (71), hence
   $\liminf_N N^{-2}\log Z^{(p,N)}_{\rm BPS}\ge\log(2\cos(\pi/2p))$ (72); any convergent subsequence limit lies in
   $[\log(2\cos\frac{\pi}{2p}),\log2]$.
8. **Hamiltonian structure (§3)**: $H/p^2$ normal-orders into an alternating sum of $K_k$; exact formulas for
   $K_0,K_1,K_2$ at $p=5$, $N\ge3$. Only $K_3,K_4$ can carry spectral information beyond the fermion number and
   the quadratic Casimir. At $N=3$ the Hodge star maps the invariant cubic form to the quintic form up to scale,
   so all five terms commute; at $N=4$ explicit integer computation gives non-zero commutators.

## Important equations

(1), (10)–(16), (39)–(40), (57)–(60), (71)–(72).

## Numerical methods

Exact integer sector maps in the ordered monomial basis; ranks via $A_R=[Q_R;Q^\dagger_{R-p}]$ and
$H_R=A_R^\dagger A_R$ (39)–(40), diagonalising $A_RA_R^\dagger$ when the target is smaller (what makes the middle
sectors of $(5,4)$ feasible). **Modular ranks over characteristics $2^{31}-1$, $2147483629$ and $p$**, all three
agreeing, plus floating-point singular-value gaps of $>4$ orders of magnitude ($N\le4$) and $>2$ ($N=5$).
Ancillary `rank_data.json`. Notably these are *the same two primes our own `cohomology.py` uses*.

Complete data at $(p,N)=(5,3),(5,4),(5,5),(7,4)$; partial at $(7,5)$ ($R\le6$ only).

## Relevant figures/results

Table 3 (summary): $(5,3)$ dim 512, $Z_{\rm BPS}=440$, index bound 400, $q_{\min}=2$, $b=1$, $\deg T=2$,
$T(1)=11$; $(5,4)$ dim 65 536, $Z=44\,000$, bound 38 000, $q_{\min}=4$, $b=3$, $T(1)=22$; $(7,4)$ $Z=59\,136$,
bound 55 468, $q_{\min}=3$, $b=1$, $T(1)=528$; $(5,5)$ $Z=15\,480\,000$, bound 11 781 250, $q_{\min}=8$, $b=4$,
$T(1)=774$. Table 5: BPS fraction $85.9\%,67.1\%,46.1\%$ and $N^{-2}\log Z=0.676,0.668,0.662$ for $N=3,4,5$.
$b_{5,N}=1,3,4$, no known formula.

**Figure 7 and §6.1 (the part closest to our project).** Charge support of the matrix model against $\mathcal N=2$
SYK at the *same* fermion count $N_f=N^2=16$:

| model | $Z_{\rm BPS}$ | index bound | difference | occupied sectors |
|---|---|---|---|---|
| SYK $\hat q=3$, $N_f=16$ | 8 748 | 8 748 | 0 | 3 |
| SYK $\hat q=5$, $N_f=16$ | 38 000 | 38 000 | 0 | 5 |
| SYK $\hat q=7$, $N_f=16$ | 55 468 | 55 468 | 0 | 7 |
| matrix $(5,4)$ | 44 000 | 38 000 | 6 000 | 9 of 17 |
| matrix $(7,4)$ | 59 136 | 55 468 | 3 668 | 11 of 17 |

"The matrix model has nonzero multiplicities in more sectors and lies strictly above the index bound." At
$\hat q=7$ the two central SYK multiplicities $h_7=11\,319$, $h_8=12\,838$ coincide with the matrix values;
differences occur away from the centre and in the extra outer sectors.

## Limitations

- $(1+x)^N$ divisibility is **conjectural** beyond $(1+x)^{2}$ / $(1+x)^{3}$; the power $b_{p,N}$ is unexplained.
- Positive multiplicities are numerical, not proved over $\mathbb Q$ (vanishings *are* certified modularly).
- Only four complete data points; large-$N$ rate undetermined between the two bounds.
- No gauging, no melonic limit, no chaos test. He is explicit that Fig. 7 "does not by itself determine fortuity
  under the projections at fixed $R$ used here or provide a test of chaos".
- Single flavour only; no multi-trace; nothing on several flavours or on general gauge algebras.

## Relationship to our project

This is the **closest published work to our report**, and it was found only after the report was drafted. It
occupies the flavour-$=1$, general-odd-$q$ axis of our $(p,q,N)$ family, which our own analysis treats only as a
degenerate corner.

**Cross-validation (2026-09-26).** Our machinery independently reproduces his $(p,N)=(5,3)$ result exactly:
$Z_{\rm BPS}=440$ on a Fock space of 512, with profile $h_{2..7}=20,75,125,125,75,20$. This matches his Table 3
row in every entry, since $440=5^1\cdot x^2(1+x)^3T$ with $T=4+3x+4x^2$, i.e. $b=1$, $\deg T=2$, $T(1)=11$.
Two independently written codes, agreeing on a published number.

**Points of genuine convergence.**
- His Prop. 2 ($Q^2=0$ from odd degree alone) is our $\omega\wedge\omega=0$.
- His §6.3 observes that traceless diagonal creation operators anticommute with $Q_p$ and make $B_0^\bullet$ a
  graded module over $\Lambda^\bullet(\mathbb C^{N-1})$, and that **freeness** of that module would explain the
  $(1+x)^N$ factor. This is exactly our D10 module structure, and his open question 1 is a close relative of our
  D10.5. Our slot factorisation *proves* that freeness **at the maximal weight** (Künneth over sites gives
  $\Lambda^\bullet\mathcal Z$-freeness there, hence the $(1+x)^N$ he sees); whether it holds on the full complex
  is open on both sides.
- His §6.1 measures the same phenomenon our no-go is about: at matched fermion number, SYK occupies exactly
  $\hat q$ sectors and saturates the index, while the matrix model occupies more and exceeds it. He reports this
  without drawing a concentration no-go from it.

**What he has that we do not.** The rank symmetry theorem (our symmetric profiles are observed, not proved); the
exact $U(1)$ decoupling; the Witten-index large-$N$ bound; the normal-ordered $K_k$ decomposition; and the
correct functorial treatment of rank comparison (projections, not inclusions).

**What we have that he does not.** Several flavours; general reductive $\g$; multi-trace; the weight/Kostant
decomposition into irreps; the rank law $W=r(w_s-1)+1$; the Weyl-group mechanism; and the identification of SYK
as the $N=1$ member.

**Two corrections to our report forced by this paper.**
1. **Our §5 fortuity framing is wrong.** We assert that the block embedding $\gl(N)\hookrightarrow\gl(N+1)$ "is a
   chain map and induces $H_N\to H_{N+1}$", then prove that map is zero. His §4.5 shows the inclusion is *not* a
   cochain map, and he is right: $\omega_{N+1}\ne\iota(\omega_N)$, since $\mathrm{Tr}_{N+1}$ contains terms in the
   new index. So no induced map exists and "the stabilisation map vanishes" is not a well-formed statement. Our
   *computational* content survives intact and in fact independently corroborates his remark: the embedded
   representative is never closed at rank $N+1$ (verified 9/9 and 27/27). The right statement is that the
   inclusion fails to descend to cohomology, and the correct functorial comparison is his projection
   $\pi_{M\to N}$.
2. **Our §2.4 "loophole" paragraph overclaims.** We say that at one flavour, $w_s=2$ gives $W=r+1$ for any $q$,
   so $q\ge N+1$ satisfies $W\le q$ at every $N$. That is true **at the maximal weight only**: we measure $W=N+1$
   with profile $\binom Nj$ at $(q,N)=(3,3),(3,4),(5,3),(5,4),(7,4)$. But the *total* support at $(5,3)$ is six
   sectors $R=2,\dots,7$, whose residues mod 5 are $2,3,4,0,1,2$, so $R=2$ and $R=7$ collide and concentration
   **fails**. Irrep-independence of the window, which we verified only at $q=3$, is therefore false at $q=5$.
   This does not threaten the main no-go (which is at $q=3$, where the Casimir descents verify
   irrep-independence) and it strengthens rather than weakens "raising $q$ hurts", but the loophole paragraph as
   written is incorrect and the rank law must be stated as a maximal-weight result.
