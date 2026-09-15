# Derivations

*Substantial derivations and theoretical notes. Each entry states assumptions, defines variables, derives step by step, checks limits, and says what is proven vs. conjectured.*

---

## D1. Single-trace form of $H=\{Q,\bar Q\}$ for the single-matrix model (2026-09-15)

**Status: proven (analytic derivation below) and verified numerically at $N=2,3,4$ to machine precision** (`scratchpad/check_traceH.py`, run 2026-09-15). This supersedes the identity quoted in `research/notes/matrix_syk_largeN_engine.md` §2 and in the checkpoint PDF §6.2, whose constant term $-9(N-1)^2$ was obtained by a numerical fit at $N=2,3$ only and is **wrong for $N\ge4$** (relative residual $0.16$ at $N=4$). The two constants agree at $N=2,3$, which is why the fit did not detect the error.

### Assumptions and conventions

Chen 2025 conventions: $\Psi_{ij}$ creation, $\bar\Psi_{ij}\equiv\Psi^\dagger_{ji}$ annihilation,
$$\{\Psi_{ij},\bar\Psi_{kl}\}=\delta_{il}\delta_{jk},\qquad Q=\mathrm{Tr}\Psi^3=\Psi_{ij}\Psi_{jk}\Psi_{ki},\qquad \bar Q=Q^\dagger=\mathrm{Tr}\bar\Psi^3 .$$
Repeated indices summed. $N_\Psi=\mathrm{Tr}[\Psi\bar\Psi]=\sum_{ij}\Psi_{ij}\Psi^\dagger_{ij}$ is the number operator; $\mathrm{Tr}[\bar\Psi\Psi]=N^2-N_\Psi$. We take as input Chen's identity
$$H=3N(N^2-1)-9\hat C_2,\qquad \hat C_2=\sum_{a=1}^{N^2-1}J^aJ^a,\quad J^a=\mathrm{Tr}[\Psi[T^a,\bar\Psi]],\quad \mathrm{Tr}[T^aT^b]=\tfrac12\delta^{ab},\tag{D1.1}$$
(Chen (2.6)–(2.8); verified as an operator identity at $N=2,3$ by `mm_model.py`). The task is to express $\hat C_2$ through traces.

### Step 1: $\hat C_2$ as a trace over an operator-valued matrix

For any $N\times N$ c-number matrix $A$, $J(A)\equiv\mathrm{Tr}[\Psi[A,\bar\Psi]]=\mathrm{Tr}[A\,\mathcal K]$ with the operator-valued matrix
$$\mathcal K_{kj}\equiv\mathcal L_{kj}-\mathcal R_{kj},\qquad \mathcal L_{kj}=\sum_i\Psi_{ij}\bar\Psi_{ki},\qquad \mathcal R_{kj}=\sum_i\Psi_{ki}\bar\Psi_{ij}.$$
($\mathrm{Tr}\,\mathcal K=N_\Psi-N_\Psi=0$.) With the completeness relation $\sum_aT^a_{ij}T^a_{kl}=\tfrac12(\delta_{il}\delta_{jk}-\tfrac1N\delta_{ij}\delta_{kl})$,
$$\hat C_2=\sum_a\mathrm{Tr}[T^a\mathcal K]\,\mathrm{Tr}[T^a\mathcal K]=\tfrac12\sum_{ij}\mathcal K_{ji}\mathcal K_{ij}-\tfrac1{2N}(\mathrm{Tr}\mathcal K)^2=\tfrac12\sum_{ij}\mathcal K_{ji}\mathcal K_{ij}.\tag{D1.2}$$

### Step 2: normal-order the four quartic pieces

Write $c\equiv\Psi$, $a\equiv\bar\Psi$ and use $c\,a\,c'a'=c\{a,c'\}a'-c\,c'\,a\,a'$. Define the two normal-ordered quartics
$$\mathcal N_1\equiv\sum c_{mi}c_{nj}a_{jm}a_{in}\ (\text{index chain of }\mathrm{Tr}[\Psi\bar\Psi\Psi\bar\Psi]),\qquad \mathcal N_2\equiv\sum c_{mi}c_{in}a_{nj}a_{jm}=\mathrm{Tr}[\Psi\Psi\bar\Psi\bar\Psi].$$
Then, using $\{c_{nj},a_{jm}\}=\delta_{nm}$ (shared index $j$), $\{c_{in},a_{jm}\}=\delta_{im}\delta_{nj}$, etc.:
$$\begin{aligned}
\sum\mathcal L_{ji}\mathcal L_{ij}&=\sum c_{mi}a_{jm}c_{nj}a_{in}=N\,N_\Psi-\mathcal N_1,\\
\sum\mathcal L_{ji}\mathcal R_{ij}&=\sum c_{mi}a_{jm}c_{in}a_{nj}=\mathrm{Tr}\Psi\,\mathrm{Tr}\bar\Psi+\mathcal N_2,\\
\sum\mathcal R_{ji}\mathcal L_{ij}&=\sum c_{jm}a_{mi}c_{nj}a_{in}=\mathrm{Tr}\Psi\,\mathrm{Tr}\bar\Psi+\mathcal N_2,\\
\sum\mathcal R_{ji}\mathcal R_{ij}&=\sum c_{jm}a_{mi}c_{in}a_{nj}=N\,N_\Psi+\mathcal N_1 .
\end{aligned}$$
(The last quartic is $-\sum c_{jm}c_{in}a_{mi}a_{nj}$; relabelling shows it equals $+\mathcal N_1$.) Inserting into (D1.2), $\mathcal N_1$ cancels:
$$\boxed{\ \hat C_2=N\,N_\Psi-\mathrm{Tr}\Psi\,\mathrm{Tr}\bar\Psi-\mathrm{Tr}[\Psi\Psi\bar\Psi\bar\Psi]\ }\tag{D1.3}$$
Here $\mathrm{Tr}\Psi\,\mathrm{Tr}\bar\Psi=N\,\psi^0\bar\psi^0$ is $N$ times the occupation of the decoupled trace mode, so $N N_\Psi-\mathrm{Tr}\Psi\mathrm{Tr}\bar\Psi=N\sum_{a\ge1}\psi^a\bar\psi^a$ involves only $su(N)$ modes, as it must.

**Checks of (D1.3).** Empty state: all terms vanish, $C_2=0$ ✓ (singlet). Filled state $|F\rangle$: $N N_\Psi\to N^3$, $\mathrm{Tr}\Psi\mathrm{Tr}\bar\Psi\to N$, $\mathrm{Tr}[\Psi\Psi\bar\Psi\bar\Psi]\to N^3-N$ (count of triples $(m,i,n)$ with $(m,i)\ne(i,n)$), total $0$ ✓. One-particle state $\Psi_{12}|0\rangle$: $C_2=N$ ✓ (adjoint). $\Psi_{11}|0\rangle$: $C_2=N-1=N(1-\tfrac1N)$ ✓ (adjoint fraction $1-1/N$).

### Step 3: the Hamiltonian

From (D1.1) and (D1.3),
$$\boxed{\ H=9\,\mathrm{Tr}[\Psi\Psi\bar\Psi\bar\Psi]-9N\,N_\Psi+9\,\mathrm{Tr}\Psi\,\mathrm{Tr}\bar\Psi+3N(N^2-1)\ }\tag{D1.4}$$
(normal-ordered form). Normal-ordering the reversed word gives
$$\mathrm{Tr}[\bar\Psi\bar\Psi\Psi\Psi]=\mathrm{Tr}[\Psi\Psi\bar\Psi\bar\Psi]-2N\,N_\Psi+2\,\mathrm{Tr}\Psi\,\mathrm{Tr}\bar\Psi+N(N^2-1),\tag{D1.5}$$
so the symmetric form, in which the trace-mode term cancels, is
$$\boxed{\ H=\tfrac92\big(\mathrm{Tr}[\Psi\Psi\bar\Psi\bar\Psi]+\mathrm{Tr}[\bar\Psi\bar\Psi\Psi\Psi]\big)-\tfrac32N(N^2-1)\ }\tag{D1.6}$$
The earlier constant $-9(N-1)^2$ equals $-\tfrac32N(N^2-1)$ only at $N=2$ ($9$) and $N=3$ ($36$); at $N=4$ the correct value is $90$, not $81$.

**Limiting checks.** $\langle0|H|0\rangle=\tfrac92N(N^2-1)-\tfrac32N(N^2-1)=3N(N^2-1)$ ✓ (top of spectrum, singlet). Coefficient $9/2$ is $N$-independent ✓ (as found numerically). Numerical check on random vectors at $N=2,3,4$: residuals $\le2\times10^{-16}$ for (D1.4) and (D1.6); $0.16$ for the old constant at $N=4$.

### Remarks

- The dynamics is single-trace plus quadratic/constant terms; the only "double-trace" piece, $\mathrm{Tr}\Psi\,\mathrm{Tr}\bar\Psi$, is the decoupled trace mode and cancels in (D1.6).
- Lesson recorded for the project: **do not infer $N$-dependence from fits at two values of $N$.** Polynomials in $N$ of degree $\le3$ appear naturally in constants (cf. $N^3-N$ above) and cannot be pinned by two points.
- The analogous computation for the 3-matrix model is not yet done; it should be carried out analytically (same technique) and verified at $N=2$ and, if feasible, $N=3$ on random vectors with sparse operators.
