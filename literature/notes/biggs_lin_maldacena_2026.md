# Biggs, Lin, Maldacena (2026) — A melonic quantum mechanical model without disorder

**File:** `papers/biggs_lin_maldacena_2026.pdf`

## Full citation

Anna Biggs, Loki L. Lin, Juan Maldacena, *A melonic quantum mechanical model without disorder*,
arXiv:2601.08908v1 [hep-th], 13 Jan 2026 (LITP-26-01). 45 pp. Princeton / Michigan / IAS.

## Main question

Is there a disorder-free quantum mechanics with the low-energy physics of $\mathcal N=2$ SYK, and what can be
said about its BPS states? The mechanism is new: melonic dominance from $SU(2)$ angular-momentum constraints
rather than from tensor indices (following Amit–Roginsky and Benedetti–Delporte).

## Physical system

- $N=2j+1$ **complex** fermions $\psi_m$, $m=-j,\dots,j$, forming one spin-$j$ irrep of $SU(2)$.
- **$\mathcal N=2$**: $Q=\frac{1}{3!}\sqrt{2\mathsf J N}\sum C^j_{m_1m_2m_3}\psi_{m_1}\psi_{m_2}\psi_{m_3}$ with $C^j$ the
  Wigner $3j$ symbol $(j\,j\,j;m_1m_2m_3)$ (eqs. 2–5); $Q^2=0$, $H=\{Q,Q^\dagger\}$. Cubic, built from creation
  operators only.
- $C^j$ is symmetric for even $j$, antisymmetric for odd $j$, zero for half-integer $j$, so **$j$ must be odd** and
  the symmetry is really $SO(3)$ (p. 4). Global, not gauged.
- $U(1)_R$ with $R=N_\psi/3$, $N_\psi=k-\frac{2j+1}{2}$ (eqs. 8–9); a $\mathbb Z_3$ subgroup commutes with $Q$.
- Normal-ordered $H=\frac{\mathsf J}{3}\big[(2j+1)-3(N_\psi+j+\frac12)+3\sum_mO^\dagger_{j,m}O_{j,m}\big]$ (eq. 14),
  a single Haldane pseudopotential.
- Picture: fermions in the lowest Landau level on a sphere with $N$ flux quanta ("fuzzy sphere"); $H$ is non-local
  on the sphere.

## Important definitions

- $d_n$ states with $J_3=n$ (eqs. 16–17); spin-$\ell$ multiplicity $D_\ell=d_\ell-d_{\ell+1}$ (eq. 18).
- Maximal spin $\ell_{\max}=j(j+1)/2$ (eq. 19), reached by filling the northern hemisphere; exactly **two** such
  irreps, and both are BPS (p. 6, §7).
- Refined Witten index $W_r=\operatorname{tr}[(-1)^Fe^{2\pi irR}]=\omega^{-(2j+1)/2}(1-\omega^r)^{2j+1}$ (eq. 23),
  and with a $J_3$ fugacity (eq. 24).

## Main assumptions

Large $j$ with $\mathsf J$ fixed for the melonic analysis; for the BPS counts (25)–(26), the **assumption** that
every BPS state has $R=\pm\frac16$, justified numerically up to $O(1)$ exceptions.

## Main analytical results

1. **BPS counts from the index under concentration** (eqs. 25–26): $D_{\rm BPS}(j,\pm\frac16)=3^j$ and
   $Z_{\rm BPS}(\theta)=2\,e^{-i\theta j(j+1)/2}\prod_{m=1}^j(1+e^{i\theta m}+e^{2i\theta m})$, with Gaussian and
   Cardy-like asymptotics (27).
2. **Melonic dominance** from a kinematic enhancement of melon diagrams by the degeneracy of angular-momentum
   configurations (§3.1); leading non-melonic corrections suppressed only by $(\log j)/j$.
3. **Nearly conformal IR** matching $\mathcal N=2$ SYK for $SU(2)$-singlet observables, plus operators with
   nontrivial spin absent in SYK (§4); Schwarzian gap prediction (eq. 40).
4. **Near maximal spin the model is a $c=1$-type 2d CFT** with two decoupled sectors, one with zero Hamiltonian
   that generates the BPS states explicitly as oscillator modes $\alpha_{-n}=\sum_m\psi^\dagger_{-n+m}\psi_m$,
   $n\equiv\pm1$ mod 3 (§7, App. C).
5. Largest-$R$ and largest-spin-at-fixed-$R$ sectors solved (§5–6).

## Important equations

(2)–(9), (14), (16)–(19), (23)–(27), (40), (56)–(59), (77)–(79), (94)–(97).

## Numerical methods

Exact diagonalisation up to $j=11$ ($2^{23}$ states), block-diagonalised by $N_\psi$ and $J_3=0$ with $H$ and $J^2$
diagonalised together (largest block 32,540 at $j=11$); floating-point; Gaussian-filtered spectral form factor.

## Relevant figures/results

**Sporadic BPS multiplets at $|R|\ne\frac16$** (table 79): none at $j=5,11$; $j=7$: $R=\pm\frac12$, $\ell=7$, one
each; $j=9$: $R=\pm\frac12$, $\ell=9$, two each, and $R=\pm\frac56$, $\ell=0$, two each. Fig. 11–13: gaps approach
the Schwarzian prediction slowly and non-uniformly; the gap is **not** a function of the $SU(2)$ Casimir.
Fig. 14: ramp and plateau in the SFF at $R=\frac16$, $\ell=16$.

## Limitations

$j\le11$ only; the sporadic states are unexplained; the BPS counts (25)–(26) assume concentration; no refined index
separating the exceptional charges; no description of BPS states away from maximal spin.

## Relationship to our project

**This is the most relevant paper so far, and it fits our framework exactly.** It is an $\mathcal N=2$ model with
$Q=\omega\wedge-$ for a cubic $\omega\in\Lambda^3\mathbb C^N$ invariant under a compact group, i.e. precisely the
object our cohomology machinery computes, with $SU(2)$ spin-$j$ matter in place of $U(N)$ adjoint matter.

**Independent reproduction (2026-09-28).** `scripts/run_su2_3j_model.py` computes $\ker Q/\operatorname{im}Q$ per
block of fixed $(k,J_3)$ with **exact** ranks over two primes (the 3j symbols are made rational by the diagonal
rescaling $\psi_m\to\psi_m/\sqrt{(j+m)!(j-m)!}$, which preserves $J_3$ and all cohomology dimensions), and resolves
spins by $h_\ell=h(J_3{=}\ell)-h(J_3{=}\ell+1)$. It reproduces **every entry of their table (79)** and their counts:

| $j$ | $R=\pm\frac16$ | $R=\pm\frac12$ | $R=\pm\frac56$ |
|---|---|---|---|
| 3 | $27=3^3$ each | none | none |
| 5 | $243=3^5$ each | none | none |
| 7 | $2187=3^7$ each | one $\ell=7$ multiplet each | none |
| 9 | $19685=3^9+2$ each | two $\ell=9$ multiplets each | two $\ell=0$ multiplets each |

Their numerics are floating-point; ours are exact, so this confirms theirs and validates our code on a third
published dataset (after Chen and Tierz). Data in `results/data/su2_3j_bps_j*.txt`.

**Their maximal-spin states are our maximal-weight complex.** Filling the northern hemisphere is our greedy
maximal-weight state; the two maximal-spin irreps are the zero-weight mode $\psi_0$ empty or filled. In our
language: $Z_{\max}=t^{j}(1+t)$, so $W=2$ at the top, and $W\ge2\alpha_0-n_0+1=2$ from the critical-set bound with
one zero-weight mode, both tight.

**Their eq. (26) is our eq. (20) mechanism.** Pair the modes $\pm m$: each pair's Fock factor
$(1+yx^m)(1+yx^{-m})$ becomes $y(x^m+1+x^{-m})$ at a primitive sixth root of unity, and the single zero-weight mode
contributes $(1+y)$. That is exactly our index-congruence derivation, with $SO(3)$ weights $m$ in place of roots
$e_i-e_j$. Their argument (assume the BPS states sit in two adjacent degrees, then solve the index) is the
special case of ours where the band has width $\le3$, in which the congruence mod $1+t^3$ determines $Z_{\rm BPS}$
by itself. The sporadic states are precisely the failure of our hypothesis (B), the band condition, which our
report already flags as the unproved input.

**The sporadic states are our trichotomy, made non-generic by symmetry.** $N=2j+1\equiv3$ mod 4 for odd $j$. For a
*generic* 3-form on $\mathbb C^N$ with $N\equiv3$ mod 4 our §4r finds $w_s=2=q-1$: every nonzero-index class
saturated and the zero-index class empty. BLM's $SO(3)$-invariant form departs from that in two ways, both at $O(1)$:
- $j=7,9$: the zero-index class ($R\equiv\frac12$ mod 1) is occupied, exactly the "$q+1$" branch of our trichotomy;
- $j=9$: the $R=-\frac16$ class carries BPS states in two degrees ($R=-\frac16$ and $\frac56$), so a
  *nonzero-index* class fails to saturate, which generic forms never do in our 27 checks.
This is the rank-one version of our Weyl-invariance finding: symmetric forms are non-generic, but at rank one the
damage is $O(1)$ multiplets, while for the adjoint at rank $\ge3$ it widens the window in every irrep. Per irrep,
strict concentration fails only in the sporadic complexes (e.g. $\ell=9$ at $j=9$ has degrees $k=8$ and $11$, both
$\equiv2$ mod 3), and those sit in zero-index classes, so they are not the macroscopic-index complexes the CCSY
conjecture is about.

**What it means for our thesis.** A rank-one symmetry with a huge matter irrep concentrates (up to $O(1)$
exceptions) and is SYK-like at large $j$; our rank-$\ge3$ adjoint models do not concentrate. That is strong external
support for "concentration is controlled by the zero-weight/rank sector". It also sharpens the variable: see the
generalised maximal-weight lemma in `research/notes/cohomology_results.md` §4w. The relevant quantity is the
**zero-weight subspace $V_0$ of the matter**, which is $p\cdot r$-dimensional for adjoint matter but
**one-dimensional** for $SO(3)$ spin $j$, however large $j$ is. That is how BLM evade the obstruction.

**Footnote 14 addresses Chen's model directly**: $Q\sim\Tr\psi^3$ with adjoint $\psi$ "can also be written in terms of
3j symbols for the adjoint. However, for the adjoint case, we do not expect melon diagrams, since the quantum
numbers are not large." So they place our class outside theirs for the same reason our reports call it planar
rather than melonic.

**Their proposed generalisation is a sharp test of ours.** §9.2 suggests $SU(N)$ with large Dynkin labels. Our
lemma predicts the answer at the top: $W_{\max}$ is the width of $\ker/\operatorname{im}$ of
$\omega|_{V_0}$ on $\Lambda^\bullet V_0$. For $SU(3)$ irreps the zero-weight multiplicity grows with the Dynkin
labels, and $S_3$ does not contain $-\mathrm{id}$, so $\omega|_{V_0}$ can be non-zero and Weyl-invariant, the regime
where our $n=12$ experiment found widening. So large-label $SU(N)$ generalisations may *not* inherit BLM's
concentration, and that is computable with our code.
