# Witten (2016) — An SYK-Like Model Without Disorder

**File:** `papers/witten_2016.pdf` (13 pp., read in full 2026-10-07)

## Full citation

Edward Witten, *An SYK-Like Model Without Disorder*, arXiv:1610.09758v2 [hep-th] (3 Nov 2016).

## Main question

Can the large-$N$ physics of SYK (melonic dominance, hence the same correlators and thermodynamics) be obtained
from an ordinary quantum system without quenched disorder?

## Physical system

- $q=D+1$ real fermion fields $\psi_0,\dots,\psi_D$, each with $n^D$ components, so $N=(D+1)n^D$ real fermions.
- Symmetry $G_0=\prod_{a<b}G_{ab}\cong O(n)^{D(D+1)/2}$ (2.1), modulo a discrete subgroup acting trivially (2.2).
  $\psi_a$ is the tensor product of the vector representations of the $G_{ab}$, $b\ne a$.
- Action (2.3): $H=i^{q/2}j\,\psi_0\psi_1\cdots\psi_D$, with each pair $(a,b)$ sharing exactly one contracted
  $G_{ab}$ index (the Gurau "coloured" contraction; for $q=4$ the vertex is a tetrahedron, fig. 2).
- Global symmetry, which "could be gauged". Since $\dim G\ll N$ for $D\ge3$, gauging hardly affects the
  thermodynamics (p. 2).

## Important definitions

- Faces $F_{ab}$: closed index loops of type $ab$; a diagram scales as $n^F$ (2.4)–(2.5).
- Coupling scaling $j=J/n^{D(D-1)/4}$ (3.6); for $D=3$, $j=J/n^{3/2}$ (2.6).
- Degree $\omega(G)=\sum_{\mathcal J}g_{\mathcal J}$ over the cyclic orderings $\mathcal J$, $g_{\mathcal J}$ the
  genus of the ribbon graph obtained by keeping adjacent strands (3.1)–(3.2).

## Main assumptions

Large $n$ at fixed $J$; perturbation theory in $j$.

## Main analytical results

1. A diagram scales as $n^{D-2\omega(G)/(D-1)!}=N^{1-2\omega/D!}$ (3.7), so the leading graphs are exactly those
   with $\omega=0$.
2. $\omega=0$ graphs are the melonic (SYK) graphs. The proof follows Bonzom–Gurau–Riello–Rivasseau. Its key step
   is $F_1=2D+\sum_{s\ge3}(s-2)F_s+\tfrac{D(D-3)}4v_0$ (3.10), which forces a two-vertex face, **valid for
   $D\ge3$**.
3. The $1/n$ corrections differ from SYK's $1/N$ expansion. For example, fig. 3(b) is $O(N^{-2/3})$ at $D=3$.

## Important equations

(2.1)–(2.6), (3.1)–(3.10).

## Numerical methods

None.

## Relevant figures/results

Fig. 2 (the tetrahedral vertex and strands); fig. 3 (leading versus subleading two-point diagrams).

## Limitations

- Leading order only; no systematic $1/n$ expansion.
- **The melonic argument needs $D\ge3$** (eq. 3.10). Witten remarks that the cubic random-tensor variants
  (FGMS supersymmetric SYK, his ref. [16]) are not obviously expressible without disorder (p. 9).

## Relationship to our project

1. **Our $U(n)^3$ quiver is the $D=2$ member of this coloured family.** With $q=3$ fields of rank $D=2$, each pair
   of fields shares one index, so the fields are bifundamentals of $U(n)^3$ and the interaction is
   $\mathrm{Tr}(ABC)$. We use complex fermions, $p$ flavours and the $\mathcal N=2$ supercharge
   $Q=\sum C_{abc}\mathrm{Tr}(A^aB^bC^c)$ rather than a Hamiltonian vertex.
2. **At $D=2$ the melonic proof fails.** The last term of (3.10) is negative, and $\omega$ reduces to the genus of
   ribbon graphs, so the large-$n$ limit is **planar, not melonic**. The quiver is a matrix model, and its large-$n$
   physics is not that of SYK by this argument. Whether it is SYK-like is open.
3. **Gauging is not innocuous at $D=2$.** $\dim G=3n^2$ against $3pn^2$ complex fermion modes is a ratio of $p$,
   not a large one. This is why the singlet sector is a severe restriction at $p=1$ (cohomology notes §4x: every
   per-irrep index is $\pm1$) and becomes macroscopic only at $p\ge2$.
4. Witten's remark that the cubic supersymmetric case is not obviously disorder-free is the gap addressed
   melonically by Biggs–Lin–Maldacena 2026 ($SU(2)$ 3j couplings) and, at the planar level, by the quiver.
