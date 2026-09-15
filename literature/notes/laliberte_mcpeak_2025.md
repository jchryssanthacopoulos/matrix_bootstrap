# Laliberte, McPeak (2025) — Bootstrapping supersymmetric (matrix) quantum mechanics

**File:** `papers/laliberte_mcpeak_2025.pdf`

## Full citation

Samuel Laliberte, Brian McPeak, *Bootstrapping supersymmetric (matrix) quantum mechanics*, arXiv:2510.01356v2 [hep-th] (v2: 3 Dec 2025). Cited by Chen 2025 as ref. [56].

## Main question

What is the full set of constraints needed to bootstrap supersymmetric quantum mechanics, and in particular the Marinari–Parisi supersymmetric matrix model with a cubic superpotential? Can one obtain rigorous bounds on ground-state data, including when SUSY is spontaneously broken?

## Physical system

1. $\mathcal N=1$ SUSY QM on a line: $Q=(p+iW')\psi$, $Q^\dagger=(p-iW')\psi^\dagger$ (2.3), $\{Q,Q^\dagger\}=2H$, $Q^2=0$ (2.1), $H=\frac12p^2+\frac12W'^2+\frac12[\psi^\dagger,\psi]W''$ (2.7); sectors $\epsilon=2F-1=\pm1$ (2.9). Cases: $W=\frac12\omega x^2$, $W=\frac{x}{\sqrt{2}g}+\frac{g}{3\sqrt2}x^3$, $W=\frac12x^2+\frac g3x^3$.
2. Marinari–Parisi matrix model: Hermitian $X_{ij}$, complex fermion matrix $\Psi_{ij}$ with $[X_{ij},P_{kl}]=i\delta_{il}\delta_{jk}$, $\{\Psi^\dagger_{ij},\Psi_{kl}\}=\delta_{il}\delta_{jk}$ (3.4);
   $$H=\tfrac12P_{ij}P_{ji}+\tfrac12\frac{\partial W}{\partial X_{ij}}\frac{\partial W}{\partial X_{ji}}+\tfrac12[\Psi^\dagger_{ij},\Psi_{kl}]\frac{\partial^2W}{\partial X_{ij}\partial X_{kl}}\quad(3.5)$$
   supercharges $Q=(P_{ij}+i\,dW/dX_{ji})\Psi_{ji}$ (3.9); gauge generator $G=i[X,P]-(\Psi\Psi^\dagger+\Psi^\dagger\Psi)+2N\mathbf 1$ (3.13); fermion number $F=\mathrm{tr}\Psi^\dagger\Psi\in[0,N^2]$ (3.14). Superpotentials $W=\frac a2\mathrm{Tr}X^2$ (3.21) and $W=\frac12X^2+\frac g3X^3$ (3.30) with $H$ (3.31).

## Important definitions

- Bootstrap matrix $M_{ij}=\langle O_i^\dagger O_j\rangle\succ0$ (1.4) in an eigenstate; Heisenberg equation $\langle[H,O]\rangle=0$ (1.5); stronger $\langle OH\rangle=E\langle O\rangle$ (non-convex, used only for $N=1$).
- Constraints for the matrix model: reality $\langle O_1\cdots O_n\rangle^*=\langle O_n^\dagger\cdots O_1^\dagger\rangle$ (3.34); gauge $\langle\mathrm{tr}GO\rangle=0$ (3.12)/(3.17); fermion number $\langle[F,O]\rangle=0$ (3.18) ⇒ block-diagonal $M$ by fermion number; supercharge constraint $\langle[Q^\dagger Q,O]\rangle=0$ (3.19) — found **redundant** at level 8; ground-state constraint $\langle O_i^\dagger[H,O_j]\rangle\succeq0$ (3.20)/(3.35) from the zero-temperature thermal bootstrap.
- Level counting as in Lin–Zheng: $\ell(X)=1$, $\ell(P)=2$, $\ell(\Psi)=\ell(\Psi^\dagger)=3/2$; "level-8 bootstrap" = all strings of level $\le4$ in the basis (44 strings ⇒ $44\times44$).

## Main assumptions

- Normalisability $\langle1\rangle=1$ (or $N$) selects normalisable states; for odd superpotentials this excludes the non-normalisable $E=0$ state, so bounds are for the lowest normalisable state.
- Large $N$ enters as a parameter ($N=100,1000$); no explicit factorisation — SDP is linear; the large-$N$ scaling is read off from the results ($E\propto N^2$ to $10^{-11}$).
- Rescaling $g\to g/\sqrt N$ to compare with Marinari–Parisi.

## Main analytical results

- General recursion for polynomial $W$ in $N=1$: (2.13)–(2.15).
- Quadratic matrix superpotential: from a $3\times3$ block, $\mathrm{tr}X^2\ge N^2/(2a)$ (3.24); from $\{\Psi,\Psi^\dagger\}$ block $0\le F\le N^2$ (3.27); $\langle H\rangle=a^2\mathrm{tr}X^2+\frac a2(2F-N^2)\ge0$ with $E=0$ requiring $F=0$ (3.28).
- Cubic MP at large $g$: $H=g^{2/3}[\frac12P^2+X^4+([\Psi^\dagger,\Psi]X+X[\Psi^\dagger,\Psi])]+O(g^{-2/3})$ (3.36) ⇒ $E_0/N^2=\kappa_0g^{2/3}+\dots$ (3.37); bootstrap $\kappa_0\simeq0.196$ vs WKB (bosonic free fermions) $0.242$ (3.39).
- Instanton estimate for $N=1$: $E_0=\frac1{2\pi}e^{-2S}$, $S=1/(6g^2)$ (2.30), matches bounds for $g\lesssim0.4$.

## Important equations

(1.4)–(1.5), (2.1), (2.3), (2.7), (2.15), (2.28), (3.4)–(3.5), (3.9), (3.12)–(3.13), (3.19)–(3.20), (3.23)–(3.28), (3.31)–(3.32), (3.34)–(3.37).

## Numerical methods

Mathematica `SemidefiniteOptimization`, SDPB, and a nonlinear solver (found impractical). Level-8 matrix model: blocks by fermion number $20\times20$ ($F=0$), $8\times8$ ($F=\pm1$), $4\times4$ ($F=\pm2$); 143 raw variables → 9 after reality → 44 after gauge → 27 after EOM; ground-state constraint adds only a $5\times5$ block and did not yield upper bounds. Constraint generation took ~8 h on a laptop; each SDPB point 5–10 min. Reality constraints needed words up to length 11, gauge up to length 8 to exhaust "useful" constraints.

## Relevant figures/results

Fig. 1 (SUSY HO spectrum recovered, non-convex islands); Fig. 3 (SUSY-breaking $E_{\min}>0$ converging with $K$); Fig. 4–7 ($N=1$ MP: upper and lower bounds, large-$g$ $g^{2/3}$ scaling, small-$g$ instanton match); Fig. 10 (matrix: $E/N^2\approx0.196g^{2/3}$); Fig. 11–13 (spurious kink at $g_0=\sqrt2g_c$ with $(g-g_0)^2$ scaling; levels 7 and 8 agree to $10^{-9}$; level 6 was extremely weak, $\approx-2800$).

## Limitations

- The ungauged matrix model could not be solved (too many free variables) — gauge constraints were essential.
- No upper bounds for the matrix model; convergence in level unclear (level 7 ≈ level 8 but $N=1$ experience shows plateaus followed by jumps).
- Kink at $\sqrt2g_c$ is unexplained (attributed to truncation).
- Only one bosonic matrix; fermions enter through $\Psi^\dagger\Psi$ bilinears in the basis.

## Relationship to our project

The closest published bootstrap of a **supersymmetric matrix QM with fermionic matrices** and therefore the most direct precedent for constraint engineering: block-diagonalisation by fermion number (exactly our $N_\Psi$ blocks), reality constraints, gauge constraints, and the observation that $\langle[Q^\dagger Q,O]\rangle=0$ is redundant given EOM. Important cautionary lessons: (i) the ungauged model was intractable for them — our single-matrix warm-up is intrinsically ungauged/non-singlet, so we need a smarter parametrisation (e.g. $SU(N)$-covariant blocks / Wigner–Eckart) rather than brute force; (ii) upper bounds required more than their small ground-state-positivity block; (iii) level-to-level convergence can be non-monotone-looking (plateaus). Their $E\ge0$ argument from positivity of $\{\Psi,\Psi^\dagger\}$ and $\{Q,Q^\dagger\}$ blocks is the same sum-of-squares mechanism that pins $E_0=0$ in our BPS sectors.
