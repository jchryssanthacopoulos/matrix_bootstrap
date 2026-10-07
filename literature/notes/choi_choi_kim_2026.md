# Choi, Choi, Kim (2026) — Fortuity on the conifold

**File:** `papers/choi_choi_kim_2026.pdf` (94 pp.; §1, §2.1 opening, §3.1–3.2.1 and §6 read; §§4–5 and
appendices skimmed, 2026-10-07)

## Full citation

Jaehyeok Choi, Sunjin Choi, Seok Kim, *Fortuity on the conifold*, arXiv:2609.19297v1 [hep-th] (16 Sep 2026).

## Main question

How should the fortuity programme (BPS black holes = fortuitous cohomology classes) be set up in a typical
$\mathcal N=1$ holographic SCFT? Such theories are strongly coupled fixed points and contain baryonic operators
that depend explicitly on $N$.

## Physical system

- Klebanov–Witten theory: $SU(N)\times SU(N)$, chirals $A_{1,2}\in(N,\bar N)$ and $B_{1,2}\in(\bar N,N)$,
  $W=\lambda\,\mathrm{Tr}(A_1B_1A_2B_2-A_1B_2A_2B_1)$ (2.1).
- Defined as the IR limit of the mass-deformed $\mathcal N=2$ orbifold theory. The classical cohomology is
  formulated in the UV, and a field redefinition reduces it to the light conifold letters (§2).

## Important definitions

- **Ordinary covering** (§3.1, following Chang–Lin): the multi-trace space $\tilde{\mathcal H}$ with rank-$N$
  relations $I_N$. The projection $\pi_N$ is a cochain map (their footnote 3).
  - *Monotone* = in $\mathrm{Im}\,\pi^*_N$.
  - *Fortuitous* = a non-trivial finite-$N$ class independent of the monotones. For a fortuitous class,
    $\tilde Q\tilde O\in I_N$ represents a non-trivial class of $H^*(I_N)$.
- **Baryons break the ordinary covering.** $\epsilon\epsilon$ contractions of two different gauge groups, such as
  $\det A$ (3.1), cannot be written as multi-traces.
- **Baryonic covering space** (3.5): multi-traces with formal inverse letters $a^{-1}=\mathrm{Cof}(a)/\det a$ (3.4).
  Baryonic excitations factor as an $N$-dependent ground state $(\det a)^{|B|}$ times an $N$-independent multi-trace
  (3.2)–(3.3). It is realised by fermionic integrals (3.8).
- **Generalised monotone**: $Q$-closed in the baryonic covering. **Strongly fortuitous**: the rest, proposed as
  the black-hole states.

## Main assumptions

Classical (one-loop-level) cohomology in the UV description; quantum corrections are not addressed (§2.1).

## Main analytical results

1. The cohomology problem of the KW theory is reduced to the light letters (§2).
2. The refined fortuity programme with baryonic covering spaces (§3).
3. Mesonic and baryonic monotones (§4), plus a digression on ABJM/BLG: the recently found ABJM fortuitous states
   exhibit generalised monotonicity (§4.3).
4. **Infinitely many strongly fortuitous cohomologies at $N=2$** (§5): detected by an index with the monotones
   subtracted, constructed explicitly, and some interpreted as hairy black holes with baryonic-condensate hair.
5. A noted partial failure of the "moduli space method" for counting multi-gravitons (§5).

## Important equations

(2.1), (2.6), (3.1)–(3.8), §5 indices.

## Numerical methods

Computer algebra for cohomology at $N=2$ and index expansions.

## Relevant figures/results

§5.3 (the explicit strongly fortuitous representatives).

## Limitations

- Classical cohomology only.
- $N=2$ for the black-hole states.
- The authors note that the refined criterion may reclassify states and leaves an ambiguity in the fortuity
  programme, "to be fixed by (black hole) physics".

## Relationship to our project

1. **Fortuity in a bifundamental quiver.** This is the closest field-theory analogue of our $U(n)^3$ quiver, with
   cyclic bifundamental superpotential and fortuity at $N=2$.
2. **Caveat for our quiver fortuity test.**
   - We used the ordinary covering, via Tierz's projection $\pi_{3\to2}$. For $U(n)^3$ gauge singlets every
     $\epsilon$ pairs with an $\bar\epsilon$ of the same node, so singlets are multi-trace and the ordinary covering
     is complete. For example, $\det A\,\det B\,\det C=\det(ABC)$.
   - But such $\det$-type singlets have $n$-dependent multi-trace expansions. In the spirit of this paper,
     classes of the form ($\det$-type ground state) × ($n$-independent excitation) could count as **generalised
     monotones**.
   - For fermionic letters many such determinants vanish or need flavour antisymmetrisation. **Not analysed.**
   - So our "90/90 fortuitous" verdict is in the ordinary (Chang–Lin/Tierz) sense, and the refined criterion
     could reclassify some classes.
3. **Methodological support.** Their footnote 3 states the cochain-map property of $\pi_N$ that our test relies on.
