# Klebanov, Lin, Meshcheriakov (2026) — Regge trajectories from the adjoint sector of Matrix Quantum Mechanics

**File:** `papers/klebanov_lin_meshcheriakov_2026.pdf`

## Full citation

Igor R. Klebanov, Henry W. Lin, Pavel Meshcheriakov, *Regge trajectories from the adjoint sector of Matrix Quantum Mechanics*, arXiv:2603.04522v2 [hep-th] (11 Apr 2026).

## Main question

What is the large-$N$ spectrum of the *adjoint* (first non-singlet) sector of $SU(N)$-symmetric Hermitian one-matrix quantum mechanics near the critical coupling ($c=1$ / 2D string limit), and how should it be interpreted in the dual 2D string theory?

## Physical system

Bosonic Hermitian $N\times N$ matrix QM with $V=\mathrm{tr}(\frac12X^2+gX^4)$ (1), $\hbar=1/N$; also cubic (13) and double-well potentials. Singlet sector = $N$ free fermions in $V$; adjoint sector governed by the Marchesini–Onofri (MO) integral equation
$$\Delta_n\Phi_n(x)=\fint_{x_1}^{x_2}dy\,\rho(y)\frac{\Phi_n(x)-\Phi_n(y)}{(x-y)^2}\quad(6)$$
with $\rho(x)=\frac1\pi\sqrt{2(\mu_F-V)}$ (4), constraint $\int\rho\Phi_n=0$; adjoint energies $E^{\rm adj}_{(s,n)}=E^{\rm singlet}_s+\Delta_n$ (5).

## Important definitions

Critical coupling $g_c=-\frac{\sqrt2}{6\pi}$ (corrects BIPZ); $\mu=V_{\max}-\mu_F\to0$; double scaling $N\mu$ fixed; effective Hamiltonian $H\phi=-\rho\fint\frac{\phi(y)}{(x-y)^2}+\eta(g,x)\phi$ (7); "short" folded strings (Regge regime) vs "long" strings (WKB regime).

## Main assumptions

Large $N$; MO equation valid for the adjoint; semiclassical approximations for $n\gg1$.

## Main analytical results

- Regge trajectories at criticality: $\Delta_n^{\rm Regge}\approx\sqrt{2n+\frac23}-\frac{2\sqrt2}\pi$ for $n\lesssim n_{\max}\approx\frac1{4\pi^2}\log^2\mu$ (3); crossover to WKB $\Delta_n^{\rm high}=(n+1)\omega+\eta$ (2). Adjoint gap $\Delta_1\approx0.7416$ (agreeing with two-point-function bootstrap [28]) — finite, so low-energy dynamics is singlet-dominated.
- Two trajectories (even/odd $n$) from the $Z_2$ of the quartic potential; cubic potential has one. Higher representations ($B_2,C_2,A_2$, dims $\sim N^4$): energies are sums $\Delta_i+\Delta_j$ (Discussion). Allowed representations must contain a zero weight; Young diagrams with box number divisible by $N$.
- Remarks on non-singlet contributions to thermodynamics and the BKT/deconfinement transition; factor-of-2 correction to the earlier WKB estimate.

## Important equations

(1)–(7), (13), (A1), (A8)–(A9).

## Numerical methods

Discretised MO equation with $M$ up to 36000 grid points and extrapolation; $\mu$ down to $10^{-14}$.

## Relevant figures/results

Fig. 2 (Regge vs WKB spectra, quartic and cubic), Fig. 3 (eigenvalues vs $\mu$, short/long string crossover).

## Limitations

- Bosonic, single matrix, large $N$ only; adjoint sector only (higher reps sketched); no SUSY, no fermions.

## Relationship to our project

Indirect but conceptually important: it is a modern, quantitative study of **non-singlet sectors of matrix quantum mechanics at large $N$** — the regime our single-matrix warm-up lives in (its fortuitous states are in the maximal-Casimir representation, as far from singlet as possible). It demonstrates (i) that non-singlet sectors have their own tractable large-$N$ description (MO equation ≈ 't Hooft equation) rather than being simply "heavy", (ii) that bootstrap methods (the two-point bootstrap of ref. [28]) can access the adjoint gap, and (iii) which $SU(N)$ representations are admissible (zero-weight / box-count divisible by $N$). For the gauged version of the 3-matrix model, the singlet-vs-adjoint gap is exactly the type of quantity Lin–Zheng 2025 App. D bounds, and this paper gives the exact answer in a solvable bosonic example to calibrate against.
