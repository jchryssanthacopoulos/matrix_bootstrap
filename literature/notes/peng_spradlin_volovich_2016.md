# Peng, Spradlin, Volovich (2016) — A Supersymmetric SYK-like Tensor Model

**File:** `papers/peng_spradlin_volovich_2016.pdf`

## Full citation

Cheng Peng, Marcus Spradlin, Anastasia Volovich, *A Supersymmetric SYK-like Tensor Model*, JHEP 05 (2017) 062,
arXiv:1612.03851v2 [hep-th] (v2: 28 Apr 2017). Brown University. 32 pp.

## Main question

Is there a disorder-free tensor model with the same large-$N$ physics as the supersymmetric SYK model of Fu,
Gaiotto, Maldacena and Sachdev (FGMS), in the way the Gurau–Witten tensor model reproduces ordinary SYK?

## Physical system

- **$\mathcal N=1$** superspace in $0+1$ dimensions: one real Grassmann $\theta$, off-shell $Q=\partial_\theta-i\theta\partial_t$,
  $Q^2=-i\partial_t$ (eqs. 1.8–1.14). On shell the supercharge is **Hermitian and $H=Q^2$**. It mirrors the
  $\mathcal N=1$ FGMS model, not the $\mathcal N=2$ one.
- **Quarks**: four real tensor superfields $\Psi_a=\psi_a+\theta\beta_a$, $a=0,\dots,3$, the Gurau–Witten
  colouring: each $\psi_a$ carries three vector indices $i_{ab}$ of three of six groups $G_{ab}\cong O(n)$ (eq. 1.4).
- **Mesons**: six superfields $\Pi_{ab}=\chi_{ab}+\theta\pi_{ab}$, rank-4 tensors carrying every index of $\Psi_a$
  and $\Psi_b$ except the shared one (eqs. 2.2–2.4).
- **Global** symmetry $O(n)^6/\mathbb Z_2^2$; no gauging, no singlet projection. Degrees of freedom
  $\mathcal N=4n^3+6n^4$ (eq. 3.1).
- Lagrangian (2.5): superpotential $\frac{ig}{2}\Pi_{ab}(\Psi_a\Psi_b+\frac12\epsilon_{abcd}\Psi_c\Psi_d)+\frac{ih}{6}\Pi_{ab}\Pi_{bc}\Pi_{ca}$.

## Important definitions

- The motivating obstruction, stated explicitly on p. 8: **"it is not possible to construct a scalar under all
  of the $O(n)$ groups that is cubic in the tensorial fields $\Psi_a$"** of the $q=4$ model. The mesons exist
  solely to make a cubic superpotential possible.
- On-shell supercharge after eliminating auxiliaries at $h=4g$ (eq. 3.19):
  $Q\propto g\big(\chi_{ab}\psi_a\psi_b+\tfrac12\epsilon_{abcd}\chi_{ab}\psi_c\psi_d\big)+\tfrac{\sqrt2g}{3}\chi_{ab}\chi_{bc}\chi_{ca}$,
  cubic in fermions, with $Q^2=H$ (eq. 3.20).
- "Pillow" operator $(\psi_a\psi_b)(\psi_a\psi_b)$ (eq. 2.13), the only other $O(n)^6$-invariant quartic besides
  the Gurau–Witten tetrahedron.

## Main assumptions

Large $n$ with $g\sim n^{-1}$ (eq. 3.2); IR conformal ansatz; the precise ratio $h/g$ taken to be immaterial to
the large-$n$ limit (asserted, p. 9).

## Main analytical results

1. **Well-defined large-$n$ limit** at $g\sim1/n$, proved in the appendix ("mesonic melon patch"): quark loops
   are $1/n$-suppressed and quarks propagate through a background of meson melons.
2. **Non-commuting limits**: integrating out the auxiliary bosons first (with $h=4g$) decouples mesons entirely
   (eqs. 2.11–2.12) and leaves the Gurau–Witten quartic *plus* the pillow operator, whose presence changes the
   large-$n$ scaling. Taking large $n$ first gives the result below.
3. **IR dimensions** (eq. 3.24): $\Delta_\chi=\Delta_\psi=\frac16$, $\Delta_\pi=\Delta_\beta=\frac23$, the unique
   solution of the transcendental equation (3.44) compatible with supersymmetry. These equal the $\mathcal N=1$
   FGMS values, from which the authors infer that chaos and the OPE spectrum agree with FGMS.
4. **4-point functions** dominated by meson-exchange ladders (§3.2).
5. **Purely mesonic model** (set $\Psi_a=0$) conjectured to have the same IR physics (§4).

## Important equations

(1.5), (1.16)–(1.19), (2.5)–(2.16), (3.1)–(3.2), (3.19)–(3.24), (3.44), (3.46).

## Numerical methods

None.

## Relevant figures/results

Fig. 4–5 (dominant $\psi\psi$ and $\chi\chi$ diagrams, each block $\sim g^2n^2$); Fig. 11 and appendix (proof of
the large-$n$ expansion by resolving vertices into index strands).

## Limitations

$\mathcal N=1$ only; no finite-$n$ spectrum or exact diagonalisation; large-$n$ proof in the appendix is the only
rigorous content; the mesons dominate, so the model is "better regarded as a tensor version of the FGMS model"
(their words, p. 3) than a supersymmetrisation of Gurau–Witten.

## Relationship to our project

**Is our model covered? No.** Different supersymmetry, matter and symmetry:

| | Peng–Spradlin–Volovich | our $(p,q,N)$ family |
|---|---|---|
| supersymmetry | $\mathcal N=1$: Majorana, $Q=Q^\dagger$, $H=Q^2$ | $\mathcal N=2$: complex, $Q^2=0$, $H=\{Q,\bar Q\}$ |
| R-symmetry / grading | none | $U(1)_R$, $[N_\Psi,Q]=qQ$ |
| BPS states | none exact: in $\mathcal N=1$ FGMS SUSY is broken non-perturbatively, $E_0\sim e^{-\alpha N}$ (our FGMS note) | exponentially many, $\ker Q/\operatorname{im}Q$ |
| matter | rank-3 quarks + rank-4 mesons, $O(n)^6$ | adjoint matrices, $U(N)$ |
| large-$N$ limit | melonic | planar |

So R-charge concentration is not even definable in PSV, and none of their results overlap with ours.

**The structural analogy is strong, though, and it clarifies where our model sits.** FGMS's supercharge is
$Q=i\sum C_{ijk}\psi^i\psi^j\psi^k$ with random $C$. Both PSV and our family are *de-randomisations* of it,
replacing $C$ by a symmetry-invariant contraction:

| | random $C$ | invariant contraction, melonic | invariant contraction, planar |
|---|---|---|---|
| $\mathcal N=1$ | FGMS $\mathcal N=1$ | **PSV** (tensor quarks + mesons) | |
| $\mathcal N=2$ | FGMS $\mathcal N=2$ = our $N=1$ member | **empty** | **our family** (adjoint matrices, trace) |

PSV's on-shell $Q$ (eq. 3.19) is exactly our kind of object: in our language its triangles are
$\{\chi_{ab},\psi_a,\psi_b\}$, $\{\chi_{ab},\psi_c,\psi_d\}$ and $\{\chi_{ab},\chi_{bc},\chi_{ca}\}$, with the coupling
fixed by $O(n)^6$ invariance the way ours is fixed by the trace.

**Confirmation of our parity observation.** In the Kim 2018 review we derived that a supercharge built from
$O(N)^3$ trifundamentals alone must have even degree, since $O(N)$ invariants pair indices, while nilpotency
needs odd degree. PSV state the same obstruction (p. 8) and resolve it with mesons. The observation should now
be attributed to them.

**The empty cell is the natural next model.** Complexify PSV: complex quarks and mesons, holomorphic cubic
supercharge $Q=\omega\wedge-$ with $\omega=g\sum(\chi_{ab}\psi_a\psi_b+\frac12\epsilon_{abcd}\chi_{ab}\psi_c\psi_d)+h\sum\chi_{ab}\chi_{bc}\chi_{ca}$.
It has odd degree, so $Q^2=0$ automatically; a $U(1)_R$ under which every field has charge one; our cohomology
framework applies verbatim. And every field carries three or four $O(n)$ vector indices, so **at even $n$ there are
no zero weights**, which removes both mechanisms of our no-go (Cartan self-loops at the maximal weight; Cartan
excess in the critical-set bound). Caveats: (a) we have **not** checked that PSV's melonic proof survives
complexification; (b) $4n^3+6n^4=128$ modes already at $n=2$ puts exact cohomology out of reach except at the
top weights, though the refined index is a product of characters and cheap at $n=2,3$.

**A smaller testbed for the same hypothesis (our suggestion).** The obstruction to a cubic invariant is specific
to rank-3 tensors. With rank-2 fields a cubic *does* exist without mesons: the three-node quiver
$Q=\Tr(\Psi_1\Psi_2\Psi_3)$ with $\Psi_1\in(n,\bar n,1)$, $\Psi_2\in(1,n,\bar n)$, $\Psi_3\in(\bar n,1,n)$ of $U(n)^3$.
Bifundamentals have weights $e^{(A)}_i-e^{(B)}_j$ with $A\ne B$, so **no zero weights at any $n$**; $3n^2$ modes, so
complete spectra are within reach at $n=2,3$ and top weights at $n=4$ with existing code. The planar rather than
melonic limit is the price. This isolates the one variable our analysis says controls the obstruction.

**Directly relevant find while checking this paper** (abstract level only, not in `papers/`): Biggs, Lin,
Maldacena, *A melonic quantum mechanical model without disorder*, arXiv:2601.08908 (Jan 2026). $SU(2)$-invariant
fermions with a supercharge built from the $SU(2)$ $3j$ symbol, same low-energy physics as supersymmetric SYK, a
melonic expansion, exact diagonalisation at small $N$, and BPS states with a simple description near maximal
angular momentum. In our language this is a **rank-one** symmetry with matter in a large irrep, exactly the
regime our rank law allows to concentrate. If the supercharge is nilpotent (the search snippet says
$\mathcal N=2$, with $Q\sim(3j)\psi^3$ and its conjugate), it is the cleanest external data point for the thesis
that concentration is a rank effect. Needs to be pulled in and read.
