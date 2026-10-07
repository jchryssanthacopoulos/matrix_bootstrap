# Chang, Colin-Ellerin, Rangamani (2018) — On melonic supertensor models

**File:** `papers/chang_collin-ellerin_rangamani_2018.pdf` (57 pp.; §1, §3, §4 overview, §6 and App. A read;
§5 and Apps. B–C, the four-point and $SL(2,\mathbb R)$ technology, skimmed, 2026-10-07). The file name misspells
the second author, who is Sean **Colin**-Ellerin.

## Full citation

Chi-Ming Chang, Sean Colin-Ellerin, Mukund Rangamani, *On melonic supertensor models*, JHEP 10 (2018) 157,
arXiv:1806.09903v3 [hep-th] (10 Dec 2018).

## Main question

Does supersymmetry simplify, or survive in, melonic tensor quantum mechanics?

## Physical system

- $\mathcal N=2$ (two supercharges) QM of real tensor **superfields** $\Psi_{a_1\cdots a_{q-1}}$ in the
  $(q-1)$-fundamental of $O(N)^{q-1}$, $q\ge4$ even.
- Melonic $q$-body superpotential. In components, fermions appear at most bilinearly, coupled to $q-2$
  **dynamical bosons** (a melonic Yukawa term).
- Global $O(N)^{q-1}$; no gauging.

## Important definitions

The melonic contraction $[\Phi^q]$, in which each pair of fields shares one index; super-Schwinger–Dyson equations;
twisted versus standard matching conditions for fermionic $SL(2,\mathbb R)$ eigenfunctions (App. B).

## Main assumptions

Large $N$ with melonic scaling; a regularisation that fine-tunes a bare mass (§2), needed because the bosonic sector
has UV divergences even in QM.

## Main analytical results

1. **At finite $N$ supersymmetry is unbroken.** The Witten index is $(-1)^N$ (3.7), with a unique supersymmetric
   ground state (3.10) whose norm vanishes as $N\to\infty$.
2. **At large $N$ the melonic IR fixed point breaks supersymmetry.** Singlet excitations do not form
   supermultiplets, and there is no goldstino. The breaking is traced to regularising bosons and fermions
   independently, i.e. it is explicit along the RG flow. There are $O(N^2)$ light modes from the $O(N)^{q-1}$
   rotations.
3. **Four-point functions and spectrum** of singlet composites (§5); the chaos exponent is maximal, from the
   reparametrisation mode.
4. **App. A, §6: no other supersymmetric melonic model works.**
   - $\mathcal N=1$ with fermionic superfields gives an odd-fermion potential.
   - $\mathcal N=4$ produces derivative couplings.
   - The FGMS $\mathcal N=2$ SYK uses **Fermi superfields with odd $q$** (their eq. 6.4), "and it is therefore
     unclear how to promote this to a melonic tensor model."
   - The authors conclude that melonic dominance is "intrinsically at tension with supersymmetry".

## Important equations

(3.7), (3.10), (6.1), (6.4); App. A (A.1)–(A.3).

## Numerical methods

Numerical solution of conformal eigenvalue equations for the operator spectrum.

## Relevant figures/results

Fig. 3 (supergraph Schwinger–Dyson equation); the spectrum tables of §5.3.

## Limitations

- Only models with dynamical bosons are analysed.
- Purely fermionic $\mathcal N=2$ models with a cubic supercharge are explicitly left open (their eq. 6.4).

## Relationship to our project

1. **Their open case is our setting.** Their obstruction is to purely fermionic $\mathcal N=2$ models with an
   odd-degree supercharge and *melonic* contractions. Our matrix models (Chen's adjoint models, the $U(n)^3$
   quiver) are purely fermionic, $\mathcal N=2$, with a cubic $Q$ and *no dynamical bosons*. That sidesteps their
   supersymmetry-breaking mechanism, at the price of being planar (quiver) or non-melonic (adjoint) rather than
   melonic.
2. **Who has filled the gap since.** Biggs–Lin–Maldacena 2026 found a melonic, purely fermionic, cubic $\mathcal N=2$
   model by replacing tensor contractions with $SU(2)$ 3j symbols. The fermionic coloured quiver is the planar
   ($D=2$) analogue.
3. **A caution for our bootstrap and large-$N$ work.** Unbroken supersymmetry at finite $N$ does not guarantee it at
   large $N$: here the finite-$N$ ground state's norm vanishes. In our purely fermionic, finite-dimensional models
   the Witten index and the BPS cohomology are exact at each $N$, but their large-$N$ continuation should not be
   assumed.
