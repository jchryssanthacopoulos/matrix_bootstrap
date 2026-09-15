# Klebanov, Milekhin, Popov, Tarnopolsky (2018) — Spectra of Eigenstates in Fermionic Tensor Quantum Mechanics

**File:** `papers/klebanov_milekhin_popov_tarnopolsky_2018.pdf`

## Full citation

Igor R. Klebanov, Alexey Milekhin, Fedor Popov, Grigory Tarnopolsky, *Spectra of Eigenstates in Fermionic Tensor Quantum Mechanics*, Phys. Rev. D 97 (2018) 106023, arXiv:1802.10263v4 [hep-th] (31 Jul 2018). Cited by Chen 2025 as ref. [53].

## Main question

Understand the energy spectrum of the $O(N_1)\times O(N_2)\times O(N_3)$ Majorana tensor QM (a disorder-free SYK-like model): rigorous energy bounds per representation, singlet counting, gaps between singlet and non-singlet sectors, and the exactly solvable fermionic *matrix* models obtained for $N_3=1,2$.

## Physical system

Majorana $\psi^{abc}$, $\{\psi^{abc},\psi^{a'b'c'}\}=\delta\delta\delta$ (1.1), tetrahedral Hamiltonian
$$H=\tfrac g4\psi^{abc}\psi^{ab'c'}\psi^{a'bc'}\psi^{a'b'c}-\tfrac g{16}N_1N_2N_3(N_1-N_2+N_3)\quad(1.2)$$
Melonic large-$N$ limit with $gN^{3/2}=J$ fixed. Special cases: $N_3=1$ ⇒ real $O(N_1)\times O(N_2)$ matrix model (6.1); $N_3=2$ ⇒ complex $O(N_1)\times O(N_2)\times U(1)$ matrix model (6.16)–(6.18); also the complex $SU(N_1)\times SU(N_2)\times U(1)$ matrix model (6.8) with 't Hooft limit $gN=\lambda$.

## Important definitions

- $SO(N_i)$ charges (2.8); discrete symmetries (interchange $P_{23}$, parities (2.9), time reversal (2.11)); $H$ is odd under interchange when two ranks coincide ⇒ spectrum symmetric under $E\to-E$.
- Composite $O(N_2N_3)$ generators $A^{bc,b'c'}$ (3.1); Casimir identity $C_2^{O(N_2N_3)}+C_2^{O(N_1)}=\frac{N_1N_2N_3}{8}(N_1+N_2N_3-2)$ (3.4).
- Singlet density matrix / projector on an irrep $\rho_R=\frac1{\dim R}\sum|e_i\rangle\langle e_i|$ (3.15).

## Main assumptions

- Melonic dominance for equal large ranks; 't Hooft scaling for the matrix cases.
- Gauging = truncation to singlets.

## Main analytical results

1. **Basic bounds** from $C_2\pm\frac12A\!\cdot\!A=\frac14(A\pm A)^2\ge0$ (3.7): $\frac Hg\le\frac18N_1N_2N_3(\dots)-C_2^{O(N_1)}$ etc. (3.8)–(3.9); for $N_3=2$ singlets $|H|\le\frac g8N_1N_2(N_1+N_2)$ (3.10), saturated by the exact solution.
2. **Refined bounds** (§3.2) via Cauchy–Schwarz on an invariant density matrix: $\mathrm{tr}[\rho_Rh]^2\le\mathrm{tr}[\rho_Rh^\dagger h]$ (3.18) rotated back with the symmetry, giving
   $$|E_R|\le\tfrac g{16}N_1N_2N_3\Big[N_1N_2N_3+N_1^2+N_2^2+N_3^2-4-\tfrac8{N_1N_2N_3}\sum_i(N_i+2)C_i^R\Big]^{1/2}\quad(3.21)$$
   Singlets (3.22); $N_3=2$ (3.23) = (3.10); equal ranks $|E|\le\frac g{16}N^3(N+2)\sqrt{N-1}\to JN^3/16$ (3.24); Casimir correction $\propto\sum_iC_i^R/N^3$ (3.25). Estimates of how close the true ground state is to the bound via an angle $\cos^2\theta=E_0^2/E_{\rm bound}^2$ (3.26)–(3.29).
3. **Sigma model** (§4): in the strong-coupling limit the global symmetry becomes an emergent gauge symmetry of the $G$–$\Sigma$ action; the leading breaking gives $H_{\rm gauge}=\frac{J}{N^2A}\sum_iC_2(O_i(N))$ (4.9) ⇒ singlet/non-singlet gap $\sim J/N$; singlet-sector level spacings expected $\sim c^{-N}$.
4. **Singlet counting** (§5): integral formulae; singlets exist only when all $N_i$ are even (anomaly interpretation §5.2); $O(N)^3$: 2, 36, 595354780 singlets for $N=2,4,6$ (Table 1); large-$N$ growth $\exp(N^3\log2/2-3N^2\log N/2)$.
5. **Fermionic matrix models** (§6):
   - $O(N_1)\times O(N_2)$: $H=-\frac g2C_2^{SO(N_2)}+\frac g{16}N_1N_2(N_2-1)=\frac g2C_2^{SO(N_1)}-\frac g{16}N_1N_2(N_1-1)$ (6.3); $C_2^{SO(N_1)}+C_2^{SO(N_2)}=\frac18N_1N_2(N_1+N_2-2)$ (6.4) ⇒ no doubly-singlet states; for even $N$ the ground state is an $O(N)_1$ singlet in the $\frac N2\times\frac N2$ square Young diagram of $SO(N)_2$, $E_0=-gN^2(N-1)/16$, gap $E_1-E_0=\frac g2(N-1)$ (6.5), finite in the 't Hooft limit ⇒ no quasi-conformal behaviour. Spectra Table 4.
   - $SU(N_1)\times SU(N_2)\times U(1)$ complex matrix model (6.8): Casimirs (6.12), constraint $C_2^{SU(N_1)}+C_2^{SU(N_2)}=\frac{N_1+N_2}{2N_1N_2}\big(\frac{(N_1N_2)^2}4-Q^2\big)$ (6.13) ⇒ only two doubly-singlet states (empty and filled); $H=g\big(2C_2^{SU(N_1)}+\frac1{N_1}Q^2-\frac14N_1N_2\big)$ (6.15); all states in the same irrep degenerate. Spectra Table 5.
   - $O(N_1)\times O(N_2)\times U(1)$ ($N_3=2$): $H=-2g\big(4C_2^{SU(N_1)}-C_2^{SO(N_1)}+C_2^{SO(N_2)}+\frac2{N_1}Q^2+(N_2-N_1)Q-\frac14N_1N_2(N_1+N_2)\big)$ (6.23) (quoted as (4.2) in Gaitan et al.); integer spectrum in units of $g$; energy *not* fixed by the $O\times O\times U(1)$ Casimirs alone; for even $N_1,N_2$ the ground state is a singlet saturating (3.23).

## Important equations

(1.1)–(1.2), (3.4)–(3.10), (3.15)–(3.25), (4.9), (6.1)–(6.5), (6.8), (6.12)–(6.15), (6.16)–(6.18), (6.23).

## Numerical methods

Exact diagonalisation of small tensor/matrix models ($N_1,N_2\le4$ for the $N_3=2$ model); character-integral evaluation of singlet counts; representation-theoretic enumeration (Appendix A, B).

## Relevant figures/results

Table 1 (singlet counts), Table 4–5 (matrix spectra), Appendix B explicit singlet states.

## Limitations

- Non-supersymmetric; the solvable matrix models are Casimir-type (integrable) and gapped in the 't Hooft limit.
- Bounds are not tight for large equal ranks.

## Relationship to our project

Two direct connections. (1) **Fermionic matrix models whose Hamiltonians are Casimirs** are the norm, not the exception: (6.3), (6.15), (6.23) are the non-SUSY cousins of Chen's $H=3N(N^2-1)-9\hat C_2$, and the structural lessons (no singlet states except trivial ones for $SU\times SU$; $O(1)$ gaps in the 't Hooft limit; integer spectra) are all directly relevant to expectations for Chen's single-matrix model and to what must *change* in the 3-matrix model to be SYK-like. (2) **Representation-resolved energy bounds** (§3.1–3.2) are a hand-made version of exactly the bootstrap we plan: positivity of $(A\pm A)^2$ and Cauchy–Schwarz with a symmetry-averaged density matrix $\rho_R$ give rigorous bounds on $E$ *within each irrep*. Their $\rho_R$ is precisely the "multiplet-averaged" functional of our covariant formulation, and (3.21)'s Casimir-dependent correction shows how a bootstrap organised by $SU(N)$ irrep naturally produces $E$-vs-$C_2$ bounds. The sigma-model gap $\sim J/N$ between singlet and adjoint sectors is the tensor analogue of the adjoint-gap bounds in Lin–Zheng 2025 App. D.
