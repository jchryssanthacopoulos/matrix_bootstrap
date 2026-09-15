# Lin, Zheng (2025) — High-Precision Bootstrap of Multimatrix Quantum Mechanics

**File:** `papers/lin_zheng_2025.pdf`

## Full citation

Henry W. Lin, Zechuan Zheng, *High-Precision Bootstrap of Multimatrix Quantum Mechanics*, arXiv:2507.21007v3 [hep-th] (v3: 11 Jun 2026). Cited by Chen 2025 as ref. [55].

## Main question

Can the Hamiltonian matrix bootstrap, pushed to high level with nonlinear relaxation and full use of $O(D)$ symmetry, produce rigorous, high-precision (eight significant digits) values for ground-state observables of large-$N$ bosonic multi-matrix quantum mechanics?

## Physical system

$D$ traceless Hermitian $N\times N$ bosonic matrices $X_I$, $[(X_I)_{ij},(P_J)_{kl}]=i\delta_{il}\delta_{jk}\delta_{IJ}$,
$$H=\tfrac12\sum_I\mathrm{Tr}P_IP_I+M^2\mathrm{Tr}X_IX_I-\tfrac{g_{YM}^2}4\sum_{I,J}\mathrm{Tr}[X_I,X_J]^2\quad(1)$$
$SU(N)$ gauge generators $C=\sum_I(-i[X_I,P_I]-N\mathbf 1)$; ground state is a gauge singlet; global $O(D)$. Cases: $D=2$ and $D=9$ ("bosonic BFSS"), massless $M=0$ and $M^2=1$. Complex letters $Z,\bar Z,P,\bar P$ for $D=2$ (2). Normalisation (3): $\mathrm{tr}=\frac1N\mathrm{Tr}$, $X\to X/\sqrt N$, $P\to P/\sqrt N$ so single-trace correlators are $O(1)$.

## Important definitions

- **Dynamical constraints**: $\langle[H,O]\rangle=0$ for all single-trace $O$ (4).
- **Kinematic constraints**: cyclicity of the trace (generating double traces), gauge invariance, Hermiticity/time reversal ($\langle\mathrm{tr}O_1\cdots O_n\rangle=\pm\langle\mathrm{tr}O_n\cdots O_1\rangle$, sign by number of $P$'s).
- **Inner-product positivity** $M_{ij}=\langle\mathrm{tr}\bar O_iO_j\rangle\succeq0$ (5).
- **Ground-state positivity** $N_{ij}=\langle\mathrm{tr}\bar O_i[H,O_j]\rangle\succeq0$ (6), the $T\to0$ limit of the thermal bootstrap; uses adjoint operators with two open matrix indices (justified since the ground state is a gauge singlet).
- **Level / hierarchy**: $\ell(X_I)=1$, $\ell(P_I)=2$, so $\ell([H,O])=\ell(O)+1$. Free variables = search space after quotienting equality constraints (Table I: e.g. $D=2$ level 14 → 1569 free of 192374; $D=9$ level 11).
- **Large-$N$ factorisation** inserted via cyclicity: double traces → products of single traces, making constraints quadratic.
- **Nonlinear relaxation** (App. B): replace $Q_{ij}=x_ix_j$ by a matrix $\mathcal Q$ with $\begin{pmatrix}1&\vec x^T\\\vec x&\mathcal Q\end{pmatrix}\succeq0$ (B6) **and** $M\succeq\mathcal Q$ (B7) (covariance positivity); stronger than earlier relaxations.
- **Irrep decomposition of positivity** (App. F): $M_{ij}=\sum_{R,\bar R,r}a_{\bar R,R}(C_{\bar R})^r_i(C_R)^r_j$ (F4); only $a_{\bar R,R}\succeq0$ per irrep need be imposed (F8); construction of $O(D)$ projectors from Young symmetrisers with trace subtraction (F18)–(F30); multiplicity handled by intertwiners $\iota_{m,n}$.

## Main assumptions

- Infinite-$N$ ('t Hooft) limit via factorisation; single-trace observables suffice.
- Ground state is a gauge singlet and $O(D)$ invariant.
- Comparison to finite-temperature Monte Carlo assumes a confined phase where leading-order observables are $T$-independent.

## Main analytical results

- Level-5 analytic bounds (7)/(E39):
  $$\frac ED\ge\max\Big\{\frac{3}{16\langle\mathrm{tr}X^2\rangle}+\frac{M^2\langle\mathrm{tr}X^2\rangle}4,\ M^2\langle\mathrm{tr}X^2\rangle\Big\},\qquad \frac ED\le\frac18\big[2M^2\langle\mathrm{tr}X^2\rangle+3(2(D-1)\langle\mathrm{tr}X^2\rangle+M^2)^{1/2}\big]$$
  For $M=0$: $\langle\mathrm{tr}X^2\rangle\ge(8D-8)^{-1/3}$; (7) is optimal through level 8 for the massless case; massive case gives an island already at level 5.
- Worked level-5 example (App. E): 19 variables, constraints (E4)–(E32), three free variables; first quadratic (double-trace) equation appears at level 7 (E33). Positivity blocks by irrep (E35)–(E38). Virial derivation (E40)–(E42).
- Adjoint gap bound (App. D): $\Delta E_{\rm adj}\le\langle\mathrm{tr}X_I[H,X_I]\rangle/\langle\mathrm{tr}X_IX_I\rangle=1/(2\langle\mathrm{tr}X^2\rangle)$ (D1); improved (D3). $D=9$: $\Delta E_{\rm adj}\le1.96339$ vs MC $2.043(76)$; for BFSS $\le1.41$.
- Large-$D$ expansion of $\langle\mathrm{tr}X_IX_IX_JX_J\rangle$ (C1), (App. G).

## Important equations

(1), (3)–(7), Table I–II, (B6)–(B7), (D1), (E39), (F4), (F8).

## Numerical methods

SDPA-GMP (arbitrary precision) because the primal SDP lacks strict feasibility and standard solvers are unstable (App. A); Mathematica notebook to generate the SDP (GitHub). Levels up to 14 ($D=2$) and 11 ($D=9$). Results (Table II): $M=0,D=2$: $E\in[0.707832,0.707868]$, $\langle\mathrm{tr}X_IX_I\rangle\in[1.15420,1.15460]$; $M^2=1,D=2$: $E\in[1.172098376,1.172098408]$; $M=0,D=9$: $E\in[6.69946,6.69968]$, $\langle\mathrm{tr}X_IX_I\rangle\in[2.29195,2.29218]$ vs MC $6.695(5)$, $2.291(1)$.

## Relevant figures/results

Fig. 1: peninsula → island at level 10 for massless $D=2$; level-10 island for $M^2=1$ already excludes the Monte Carlo point. Fig. 2: $D=9$ level 10/11 islands. Table I (variable counts), Table II (bounds).

## Limitations

- Bosonic only; flat directions make upper bounds on $\langle\mathrm{tr}X^2\rangle$ hard (island only at level 10 for $M=0$).
- Requires large-$N$ factorisation and a singlet ground state.
- Numerical instability requires arbitrary-precision solvers.
- Discussion notes the $\mathcal N=2$ SUSY $SO(D)$ models [54,55] as an interesting future target.

## Relationship to our project

This is the state-of-the-art template for the *technology*: hierarchy/level truncation, symmetry-blocked positivity (App. F is directly reusable for an $O(p)$/$S_p$ flavor symmetry and for $SU(N)$ covariant blocks), ground-state positivity (6) for two-sided bounds, nonlinear relaxation (B6)–(B7) if we ever impose factorisation, and the practical warnings about solver precision. Two differences for us: our model is *purely fermionic* (finite Hilbert space, no flat directions, the operator algebra truncates by Pauli exclusion), and for the single-matrix warm-up the interesting states are *non-singlet*, so the singlet-ground-state justification for using adjoint operators in (6) must be re-examined (see synthesis). The adjoint-gap bound (D1) is a model for bounding gaps between $SU(N)$ sectors in our models.
