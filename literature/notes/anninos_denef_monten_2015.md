# Anninos, Denef, Monten (2015) — Grassmann Matrix Quantum Mechanics

**File:** `papers/anninos_denef_monten_2015.pdf`

## Full citation

Dionysios Anninos, Frederik Denef, Ruben Monten, *Grassmann Matrix Quantum Mechanics*, JHEP 04 (2016) 138, arXiv:1512.03803v1 [hep-th] (11 Dec 2015). Cited by Chen 2025 as ref. [50].

## Main question

Can quantum mechanical models whose only degrees of freedom are Grassmann-valued (rectangular) matrices, hence with a finite-dimensional Hilbert space, develop at low energies an effective description in terms of *bosonic* Hermitian matrix quantum mechanics with an emergent semiclassical phase space? Motivation: finite-dimensional holography (de Sitter static patch).

## Physical system

- Vector model (§2): $\psi^\alpha_A$, $A=1..M$, $\alpha=1,2$ an $SU(2)$ spinor index; $\{\bar\psi^\alpha_A,\psi^\beta_B\}=\delta^{\alpha\beta}\delta_{AB}$; action (2.1) with quartic $g(\bar\psi\sigma^a\psi)(\bar\psi\sigma^a\psi)$; Hilbert space $2^{2M}$.
- Matrix model (§3): $\psi^\alpha_{Ai}$, $A=1..M$, $i=1..N$, symmetry $U(M)\times SU(N)\times SU(2)$; $S=\int dt\,[i\bar\psi_{iA}\partial_t\psi_{Ai}+g(\bar\psi_{iA}\sigma^a\psi_{Aj})(\bar\psi_{jB}\sigma^a\psi_{Bi})]$ (3.1); $H=-g\sum:\bar\psi_{iA}\sigma\psi_{Aj}\bar\psi_{jB}\sigma\psi_{Bi}:$ (3.2); Hilbert space $2^{2NM}$.

## Important definitions

- Spin operators $\hat J_a=\bar\psi\sigma^a\psi/2$; $U(M)$ generators (2.2) with normal-ordering constant $c$.
- Matrix spin operators $\hat S^a_{ij}=\sum_A\bar\psi_{iA}\sigma^a\psi_{Aj}/2$, $\hat S^0_{ij}$; they close into $u(2N)$ (3.3)–(3.5); $H$ commutes only with the diagonal $U(N)$.
- Auxiliary bosonic fields: $x$ (vector case), three Hermitian $N\times N$ matrices $\Sigma^a_{ij}$ (on-shell $\Sigma^a=2\hat S^a$), $R\equiv\Sigma^a\otimes\sigma^a$ (3.6).
- Bloch coherent states (2.24)–(2.25) on $CP^1$ with Fubini–Study metric $ds^2=2M\,dz d\bar z/(1+z\bar z)^2$ (2.29); Berezin coherent states $|Z^\dagger_{ij}\rangle$ (3.33) with Kähler potential $K=M\log\det(I+ZZ^\dagger)$ (3.31), metric (3.30), $U(2N)$ isometry (3.29).

## Main assumptions

- Large $M$ (many "flavours") for the small-velocity expansion; the long rectangular limit is needed for semiclassical bosonic dynamics.
- $SO(3)$-invariant potential; quartic interactions.

## Main analytical results

1. Vector model exactly solvable: $\hat H=-4\hat J\cdot\hat J+3\hat n$ (2.3), $E=-4J(J+1)+3n$; ground states = maximally spinning spin-$M/2$ multiplet, $E_g=-M(M-1)$; $J=0$ singlet has *maximal* energy. Degeneracies $d_J=\binom{2M}{M+2J}-\binom{2M}{M+2J+2}$ (2.5), $d_{J,n}$ (2.7). Gauging $U(M)$ with $c=-M$ leaves only the $M+1$ ground states.
2. Effective action $S_{\rm eff}=-M\,\mathrm{Tr}\log(-\partial_\tau+\sigma\cdot x)+\int r^2/4$ (2.10); reparametrisation symmetry (2.13) ⇒ no kinetic term for $r$; $V_{\rm eff}=-Mr+r^2/4$ (2.15), minimum $r=2M$, $V=-M^2$; Berry-phase term = monopole of strength $M/2$ (2.16); quadratic term (2.17). Small-velocity expansion valid for $\omega_c\ll r\sim M$ (2.18). Thermal $V_{\rm eff}(\beta)=-\frac{2M}\beta\log\cosh\frac{r\beta}2+\frac{r^2}4$ (2.21); bosonic description breaks down for $\beta\lesssim1/M$.
3. Matrix model: generating function (3.6)–(3.7) with $e^{M\mathrm{Tr}\log(-\partial_\tau+R)}$; linear-velocity term depends on $2N^2$ of the $3N^2$ variables; emergent phase space = Berezin's compact Kähler manifold; classical Hamiltonian $H[Z,Z^\dagger]=-NM^2+M^2\mathrm{tr}(S^0)^2$ (3.38) minimised when $[Z,Z^\dagger]=0$, $E_{\min}=-NM^2$; in the rescaled large-$M$, $M\gg N$ limit the dynamics reduces to $H\propto\mathrm{tr}[\tilde Z,\tilde Z^\dagger]^2$ (3.42).
4. Counting: $\dim\mathcal H_K=\prod_{j=1}^N\frac{\Gamma[N+M+j]\Gamma[j]}{\Gamma[N+j]\Gamma[M+j]}$ (3.35) = number of $U(M)$-invariant states (Appendix); limits $N\gg M$: $\sim2^{2MN}$; $M\gg N$: $\sim M^{N^2}$; $M=\alpha N$: $\log\dim=f(\alpha)N^2$ (3.36)–(3.37).

## Important equations

(2.3), (2.5), (2.10), (2.15)–(2.17), (2.29), (3.1)–(3.6), (3.30)–(3.31), (3.35)–(3.38), (3.42).

## Numerical methods

None beyond plotting exact formulas (Figs. 1–2).

## Relevant figures/results

Fig. 1 ($d_J$ vs $J$ for $M=70$, peak at $J\approx\sqrt{M/8}$); Fig. 2 ($r_{\min}$ vs $\beta$).

## Limitations

- Matrix-model spectrum not solved exactly ("combinatorial problem ... rather difficult"); only the low-energy effective description.
- Requires the long-rectangular limit $M\gg N$; not a square-matrix 't Hooft limit.
- Not supersymmetric; no chaos diagnostics.

## Relationship to our project

Background for *purely fermionic matrix QM*: it establishes the generic structural features we also meet — finite Hilbert space, energy fixed by a Casimir in the solvable vector case (cf. Chen's $H=3N(N^2-1)-9\hat C_2$, where likewise the singlet is the *highest* state and a maximal-spin/maximal-Casimir multiplet is the ground state), and the bosonisation via Hubbard–Stratonovich into a bosonic matrix with an emergent Berezin/Kähler phase space. The counting of $U(M)$-invariant states (3.35) is the kind of representation-theoretic bookkeeping we need for resolving our Hilbert space into $(N_\Psi,SU(N)\text{-irrep})$ blocks. The paper's lesson that semiclassical bosonic matrix dynamics emerges only for $M\gg N$ is a warning: for square adjoint fermion matrices there is no obvious small parameter making the fortuitous sector semiclassical (consistent with Chen's remark that fortuitous states are not captured by a classical saddle).
