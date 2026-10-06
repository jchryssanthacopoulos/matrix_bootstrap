# Adams (2025) — Thermal Bootstrap of Large-$N$ Matrix Models via Conic Optimization

**File:** `papers/adams_2025.pdf` (15 pp., read in full 2026-10-06)

## Full citation

Sophia M. Adams, *Thermal Bootstrap of Large-$N$ Matrix Models via Conic Optimization*, arXiv:2511.01209v2
[hep-th] (15 Feb 2026). Code: github.com/smadams821/Thermal-Bootstrap-of-Large-N-MQM, Zenodo
10.5281/zenodo.17497230. Cited as ref. [30] by Klebanov–Lin–Meshcheriakov 2026.

## Main question

Can the thermal (KMS) bootstrap of large-$N$ matrix QM avoid the semidefinite relaxation of the matrix logarithm
used by Cho–Gabai–Sandor–Yin 2024, by imposing the operator relative-entropy cone directly with a nonlinear conic
solver? Does that reach higher levels and tighter bounds, and can low-temperature long-string data (gap,
couplings) be extracted from the bounds?

## Physical system

- Ungauged one-matrix anharmonic oscillator $H=\mathrm{Tr}(\frac12P^2+\frac12X^2+\frac gNX^4)$ (2.1), at $g=2$.
- Two-matrix model $H=\mathrm{Tr}(\frac12(P_1^2+P_2^2)+\frac12(X_1^2+X_2^2)-\frac gN[X_1,X_2]^2)$ (2.5), at
  $g=0.1$, rewritten with complex $X,\bar X,P,\bar P$ (2.6)–(2.7) to make the $SO(2)=U(1)$ charge manifest.
- 't Hooft limit throughout.

## Important definitions

- **Constraints** (§2.1), on single-trace words up to length $L/2$:
  - moment matrix $\mathcal M_{ij}=\langle\mathrm{Tr}(O_i^\dagger O_j)\rangle_\beta\succeq0$;
  - parity, time reversal, reality;
  - stationarity $\langle[H,\mathrm{Tr}O_i]\rangle=0$;
  - cyclicity with large-$N$ factorisation;
  - normalisation, and the large-$N$ commutator (2.4).
- **KMS as a conic constraint** (3.1)–(3.5): $\beta C\succeq A^{1/2}\log(A^{1/2}B^{-1}A^{1/2})A^{1/2}$, where
  - $A_{ij}=\frac1{N^2}(\langle\mathrm{Tr}O_i^\dagger O_j\rangle-\frac1N\langle\mathrm{Tr}O_i^\dagger\rangle\langle\mathrm{Tr}O_j\rangle)$;
  - $B$ is the same with the order reversed;
  - $C_{ij}=\frac1{N^2}\langle\mathrm{Tr}(O_i^\dagger[H,O_j])\rangle$;
  - operators have length $(L-4)/2$.
  This is the operator relative-entropy cone, equivalent to KMS. At zero temperature it is replaced by
  $\langle HO_i\rangle=E\langle O_i\rangle$.
- **Two ways to impose it.**
  - Linear relaxation: Gauss–Radau rational upper bounds $r_{m,k}$ on $\log$ (3.6)–(3.8), with $(m,k)=(3,3)$, solved
    by MOSEK.
  - Directly: QICS (Quantum Information Conic Solver; Skajaa–Ye homogeneous interior point with Hypatia's stepping),
    which handles the non-symmetric cone using only the primal barrier.

## Main assumptions

Thermal state, ungauged model, large-$N$ factorisation. The spectral extraction assumes the low-temperature
long-string effective theory (4.5)–(4.12).

## Main analytical results

- **Low-temperature long-string expansion.** In the planar limit (4.5),
  $E_{\rm L.T.}(\beta)/N^2=e_0+\Delta_1e^{-\beta\Delta_1}+\dots$ (4.10). To second order (4.12) it includes the
  coupling term $[\Delta_1+h_{1111}(1-2\beta\Delta_1)]e^{-2\beta\Delta_1}$, with $h_{abcd}$ defined by the integral
  (4.6) over the Marchesini–Onofri eigenfunctions (4.4).
- **Benchmark values** at $g=2$ (from Marchesini–Onofri): $e_0=0.8654577$, $\Delta_1=2.1281936$,
  $\Delta_2=4.6201131$, $\Delta_3=7.0716122$, $\Delta_4=9.5258038$. A Nyström computation gives
  $h_{1111}\approx0.3278$ (Table 1).

## Important equations

(2.1), (2.4), (2.5)–(2.7), (3.1)–(3.8), (4.4)–(4.6), (4.10)–(4.12), Tables 1–3.

## Numerical methods

- **One-matrix model.**
  - At $L=8$, MOSEK with the log relaxation and QICS agree, both in seconds.
  - At $L=10$, MOSEK is numerically unstable unless the tolerance is relaxed from $10^{-8}$ to $10^{-7}$; QICS
    solves it in about a minute.
  - At $L=12$ only QICS works, taking 10–100 minutes per optimisation.
  - Hardware: a 10-core laptop.
- **Two-matrix model.** The $U(1)$-charge basis cuts the variables from 25 to 12 at $L=4$, and from 220 to 81 at
  $L=6$. $L=6$ is unstable for both solvers, and QICS fails above $T=0.56$.

## Relevant figures/results

- Figs. 1–2: thermal energy bounds at $L=8,10,12$, tighter than Cho–Gabai–Sandor–Yin.
- **Fit to the $L=12$ bounds, low-$T$ part** (Table 2):
  - $e_0=0.865457750210\pm3\times10^{-7}$;
  - $\Delta_1=2.1283360\pm2\times10^{-4}$, against $2.1281936$;
  - $h_{1111}=0.32731\pm7\times10^{-2}$, against $0.3278$.
- **Fit to the full curve** (Table 3): $\Delta_1=2.1281758\pm5\times10^{-3}$, which is the "within 0.001%" of the
  abstract; and $h_{1111}=0.49983\pm4\times10^{-2}$.
- Fig. 6: two-matrix bounds at $L=4,6$.

## Limitations

- The extraction of $\Delta_1$ and $h_{1111}$ is by fitting rigorous bounds to an assumed effective form, so the
  extracted values are not themselves rigorous.
- Double precision limits accuracy, and the two-matrix case is numerically unstable at $L=6$.
- Bosonic, ungauged, planar only; no fermions, no supersymmetry.

## Relationship to our project

1. **If we take the thermal route to the near-BPS density, impose KMS exactly.** QICS (a Python package) handles
   the relative-entropy cone without the log relaxation and stays stable at higher level, where MOSEK fails.
2. **A template for extracting spectral data from rigorous thermal bounds:** fit to a low-temperature effective
   form. Our analogue would fit low-$T$ bounds to the super-Schwarzian form
   $Z=N_{\rm BPS}+\sum_q\int dE\,\rho_q(E)e^{-\beta E}$, with $\rho_q\propto\sinh(2\pi\sqrt{E-E_0(q)})/E$
   (Turiaci–Witten (3.11)), and treat the Schwarzian scale as the fit parameter. As here, the extracted numbers
   would be estimates.
3. **Caveat specific to our model.** A Gibbs state restricted to one fermion-number sector satisfies KMS only for
   charge-neutral operators. The grand-canonical $e^{-\beta(H-\mu N_\Psi)}$ satisfies it for every operator, with
   $H\to H-\mu N_\Psi$ as the modular Hamiltonian.
4. **Numerical lessons that transfer.** Exploit charge conservation (for us $N_\Psi$ and $\mathbb Z_3$) to cut
   variables, and expect to need higher precision at higher level.
