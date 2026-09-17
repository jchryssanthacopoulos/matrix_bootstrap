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

---

## D2. Single-trace form of $H=\{Q,\bar Q\}$ for cubic multi-matrix supercharges, and the three-matrix model (2026-09-16)

**Status: proven for an arbitrary cubic supercharge (D2.6) and verified numerically to machine precision** on random vectors for $(N,p)=(2,1),(3,1),(2,2),(3,2),(2,3)$ with random complex couplings and with Chen's couplings (`scripts/verify_trace_hamiltonian.py`). The Casimir statements of D2.9 are proven analytically for the structure and verified numerically at $N=2$ by an exhaustive symmetry-sector test.

### D2.1 Setup and conventions

$p$ complex fermion matrices $\Psi^a_{ij}$, $a=1,\dots,p$, $i,j=1,\dots,N$, with
$$\{\Psi^a_{ij},\bar\Psi^b_{kl}\}=\delta^{ab}\delta_{il}\delta_{jk},\qquad\bar\Psi^a_{ij}\equiv(\Psi^a_{ji})^\dagger,\qquad\{\Psi,\Psi\}=\{\bar\Psi,\bar\Psi\}=0 .\tag{D2.1}$$
Note that $\bar\Psi^a_{kl}$ is the conjugate of $\Psi^a_{lk}$ (transposed pairing). Fermion number $N_\Psi=\sum_a\mathrm{Tr}[\Psi^a\bar\Psi^a]$; flavor bilinears $F^{cd}\equiv\mathrm{Tr}[\Psi^c\bar\Psi^d]$, so $N_\Psi=\sum_cF^{cc}$.

**General cubic supercharge.** For a complex tensor $C_{abc}$,
$$Q=C_{abc}\,\mathrm{Tr}[\Psi^a\Psi^b\Psi^c]=C_{abc}\,\Psi^a_{ij}\Psi^b_{jk}\Psi^c_{ki}\qquad(\text{sums over all indices}).\tag{D2.2}$$
Since a cyclic rotation of the three anticommuting letters costs $(-1)^2=+1$, $\mathrm{Tr}[\Psi^a\Psi^b\Psi^c]$ is cyclically symmetric and only the cyclically symmetric part of $C$ contributes; **we assume $C_{abc}=C_{bca}=C_{cab}$ throughout.** No other symmetry is assumed. $Q^2=0$ automatically (only creation operators). The conjugate is
$$\bar Q=Q^\dagger=\bar C_{abc}\,(\Psi^c_{ki})^\dagger(\Psi^b_{jk})^\dagger(\Psi^a_{ij})^\dagger=\bar C_{abc}\,\bar\Psi^c_{ik}\bar\Psi^b_{kj}\bar\Psi^a_{ji}=\bar C_{abc}\,\mathrm{Tr}[\bar\Psi^c\bar\Psi^b\bar\Psi^a].\tag{D2.3}$$

**Chen's three-matrix model** (Chen 2025, eq. 4.3): $Q=\sum_{1\le a\le b\le c\le3}\mathrm{Tr}[\Psi^a\Psi^b\Psi^c]$, i.e. $C=\mathcal S_{\rm cyc}[\tilde C]$ with $\tilde C_{abc}=1$ if $a\le b\le c$, else $0$, and $\mathcal S_{\rm cyc}[T]_{abc}=\tfrac13(T_{abc}+T_{bca}+T_{cab})$. Explicitly:
$C_{aaa}=1$; $C_{abc}=\tfrac13$ for every cyclic rotation of $(1,1,2),(1,1,3),(1,2,2),(1,3,3),(2,2,3),(2,3,3)$ and of $(1,2,3)$; $C_{132}=C_{213}=C_{321}=0$.
Two facts follow immediately and are used below: (i) $C$ is invariant under **cyclic** flavor relabelings $1\to2\to3\to1$ but **not** under transpositions (the class of $\mathrm{Tr}[\Psi^1\Psi^2\Psi^3]$ is present, that of $\mathrm{Tr}[\Psi^1\Psi^3\Psi^2]$ is absent) — the flavor symmetry of Chen's supercharge is $\mathbb Z_3$, not $S_3$ (confirmed numerically in D2.9; the earlier project note claiming "$S_3$" is incorrect); (ii) $C$ is real.

**Flavor tensors.** Define
$$K_{bc;ed}\equiv\sum_aC_{abc}\bar C_{dea},\qquad M_{cd}\equiv\sum_{ab}C_{abc}\bar C_{abd},\qquad L_{cd}\equiv\sum_{ab}C_{abc}\bar C_{bad},\qquad \|C\|^2\equiv\sum C_{abc}\bar C_{abc},\qquad \langle C,C^{\rm rev}\rangle\equiv\sum C_{abc}\bar C_{cba}.\tag{D2.4}$$
$M$ and $L$ are Hermitian; $\langle C,C^{\rm rev}\rangle$ is real; $K_{bc;ed}=\overline{K_{de;cb}}$.

### D2.2 Wick structure of the anticommutator

Write $Q=C_{abc}\,c_1c_2c_3$ with $c_1=\Psi^a_{ij},\ c_2=\Psi^b_{jk},\ c_3=\Psi^c_{ki}$, and $\bar Q=\bar C_{def}\,a_1a_2a_3$ with $a_1=\bar\Psi^f_{lm},\ a_2=\bar\Psi^e_{mn},\ a_3=\bar\Psi^d_{nl}$ (this is (D2.3) with $(a,b,c)\to(d,e,f)$). Let $D(a_p,c_q)\equiv\{a_p,c_q\}$ denote the c-number contraction; from (D2.1), $\bar\Psi^x_{uv}$ contracts with $\Psi^y_{rs}$ iff $x=y$, $u=s$, $v=r$.

**Lemma.** For $A=c_1c_2c_3$ and $B=a_1a_2a_3$, $\{A,B\}$ equals the sum of all Wick terms of $a_1a_2a_3c_1c_2c_3$ with at least one contraction. *Proof:* Wick's theorem gives $BA=\,:\!BA\!:\,+(\text{contracted terms})$, and $:\!a_1a_2a_3c_1c_2c_3\!:\,=(-1)^{3\cdot3}c_1c_2c_3a_1a_2a_3=-AB$, which cancels $AB$. $\square$

Hence
$$H=\{Q,\bar Q\}=H_4+H_2+H_0,\tag{D2.5}$$
with $H_4$ (one contraction: quartic, normal-ordered), $H_2$ (two contractions: quadratic) and $H_0$ (three contractions: c-number). Wick signs: contracting $a_p$ with $c_q$ in $a_1a_2a_3c_1c_2c_3$ and normal-ordering the remainder gives the sign $(-1)^{p+q}$ times $D(a_p,c_q)\,c_{q'}c_{q''}a_{p'}a_{p''}$ with the survivors in their original relative order (move $a_p$ rightward past $3-p$ annihilators and $q-1$ creators, then the remaining $aacc\to ccaa$ costs $(-1)^4$). Multiple contractions are handled below by explicit permutation counting.

**Cyclic reduction.** Because $C$ and $\bar C$ are cyclic and cyclic rotation of either word is a symmetry of the *summed* expression (sign $+1$), the Wick term obtained by contracting "letter $q$ of $Q$ with letter $p$ of $\bar Q$" is, after summation, independent of $(p,q)$. So the 9 single contractions are all equal, the 18 double contractions fall into two $\mathbb Z_3\times\mathbb Z_3$ orbits of 9, and the 6 triple contractions into two orbits of 3. It suffices to compute one representative of each class; the numerical verification in D2.7 tests the whole statement.

### D2.3 Quartic part (single contractions)

Representative $(p,q)=(1,1)$, sign $+$: $D(\bar\Psi^f_{lm},\Psi^a_{ij})=\delta^{fa}\delta_{jl}\delta_{im}$; remainder $c_2c_3a_2a_3=\Psi^b_{jk}\Psi^c_{ki}\bar\Psi^e_{mn}\bar\Psi^d_{nl}\to\Psi^b_{jk}\Psi^c_{ki}\bar\Psi^e_{in}\bar\Psi^d_{nj}=\mathrm{Tr}[\Psi^b\Psi^c\bar\Psi^e\bar\Psi^d]$ (indices $j\to k\to i\to n\to j$; the operator order is already creation-first). Therefore
$$H_4=9\,K_{bc;ed}\ \mathrm{Tr}[\Psi^b\Psi^c\bar\Psi^e\bar\Psi^d]\equiv9\,K_{bc;ed}\,V^{bced},\qquad V^{bced}\equiv\Psi^b_{ij}\Psi^c_{jk}\bar\Psi^e_{kl}\bar\Psi^d_{li}.\tag{D2.6}$$
(As a check of the cyclic-reduction claim, the representative $(2,1)$ gives $-\,C_{abc}\bar C_{daf}\,\Psi^b_{jk}\Psi^c_{ki}\bar\Psi^f_{lj}\bar\Psi^d_{il}=+C_{abc}\bar C_{daf}\mathrm{Tr}[\Psi^b\Psi^c\bar\Psi^d\bar\Psi^f]$, which equals the $(1,1)$ term after relabeling and using cyclicity of $\bar C$.)

**Compact form.** Define the matrix-valued operators
$$X_a\equiv C_{abc}\,\Psi^b\Psi^c\ \ \big((X_a)_{ik}=C_{abc}\Psi^b_{ij}\Psi^c_{jk}\big),\qquad \bar X_a\equiv\bar C_{abc}\,\bar\Psi^c\bar\Psi^b,\tag{D2.7}$$
so that $Q=\mathrm{Tr}[\Psi^aX_a]$, $(\bar X_a)_{ki}=\big((X_a)_{ik}\big)^\dagger$, and $\partial Q/\partial\Psi^a_{ji}=3(X_a)_{ij}$. Then
$$H_4=9\sum_{a=1}^p\mathrm{Tr}[X_a\bar X_a]\tag{D2.8}$$
— the normal-ordered "$|\partial Q|^2$", as expected for $\{Q,\bar Q\}$ with $Q$ built from creation operators only. In particular $K$ has rank $\le p$ as a $p^2\times p^2$ matrix (rank 3 for Chen's $C$).

### D2.4 Quadratic part (double contractions)

*Parallel class*, representative $(a_1c_1)(a_2c_2)$, remainder $c_3a_3$. Permutation $a_1a_2a_3c_1c_2c_3\to a_1c_1a_2c_2a_3c_3$ has sign $-1$; the normal-ordered remainder is $:\!a_3c_3\!:=-c_3a_3$; net $+D_{11}D_{22}\,c_3a_3$. The deltas force $l=j$, $m=i=k$, $n=j$, leaving $\sum_{i,j}\Psi^c_{ii}\bar\Psi^d_{jj}$ with flavor factor $\delta^{fa}\delta^{eb}$:
$$T_\parallel=C_{abc}\bar C_{dba}\,\mathrm{Tr}\Psi^c\,\mathrm{Tr}\bar\Psi^d=L_{cd}\,\mathrm{Tr}\Psi^c\,\mathrm{Tr}\bar\Psi^d.$$
*Crossed class*, representative $(a_1c_2)(a_2c_1)$, remainder $c_3a_3$. Permutation sign $+1$, normal ordering $-1$; net $-D(a_1,c_2)D(a_2,c_1)\,c_3a_3$. The deltas force $l=k$, $m=j$, $n=i$ with $j$ free (factor $N$), leaving $\sum_{i,k}\Psi^c_{ki}\bar\Psi^d_{ik}$ with flavor $\delta^{fb}\delta^{ea}$:
$$T_\times=-N\,C_{abc}\bar C_{dab}\,\mathrm{Tr}[\Psi^c\bar\Psi^d]=-N\,M_{cd}\,F^{cd}.$$
Each class has 9 members, so
$$H_2=9\,L_{cd}\,\mathrm{Tr}\Psi^c\,\mathrm{Tr}\bar\Psi^d-9N\,M_{cd}\,\mathrm{Tr}[\Psi^c\bar\Psi^d].\tag{D2.9}$$

### D2.5 Constant part (triple contractions)

*Even class*, representative $(a_1c_1)(a_2c_2)(a_3c_3)$: sign $-1$; all six indices are forced equal (one free sum, factor $N$); flavor $\delta^{fa}\delta^{eb}\delta^{dc}$: $T_{\rm id}=-N\,C_{abc}\bar C_{cba}$.
*Odd class*, representative $(a_1c_1)(a_2c_3)(a_3c_2)$: sign $+1$; the deltas force $l=j$, $m=i$, $n=k$ with $i,j,k$ free (factor $N^3$); flavor $\delta^{fa}\delta^{ec}\delta^{db}$: $T_{\rm odd}=N^3\,C_{abc}\bar C_{bca}=N^3\|C\|^2$.
Each class has 3 members:
$$H_0=3N^3\,\|C\|^2-3N\,\langle C,C^{\rm rev}\rangle.\tag{D2.10}$$

### D2.6 Result for a general cyclic $C$

$$\boxed{\ H=\{Q,\bar Q\}=9\sum_a\mathrm{Tr}[X_a\bar X_a]\;-\;9N\,M_{cd}\,\mathrm{Tr}[\Psi^c\bar\Psi^d]\;+\;9\,L_{cd}\,\mathrm{Tr}\Psi^c\,\mathrm{Tr}\bar\Psi^d\;+\;3N^3\|C\|^2-3N\langle C,C^{\rm rev}\rangle\ }\tag{D2.11}$$
with $X_a$, $\bar X_a$ from (D2.7) and the tensors from (D2.4). All terms are normal-ordered (creation operators to the left), so (D2.11) is a canonical form. $H$ is manifestly Hermitian term by term (using $(V^{bced})^\dagger=V^{decb}$ and $K_{bc;ed}=\overline{K_{de;cb}}$).

Structure: a single-trace quartic, a single-trace quadratic, a product of two length-one traces, and a constant. The only "double-trace" piece is $\mathrm{Tr}\Psi^c\,\mathrm{Tr}\bar\Psi^d=N\,\psi^{0c}\bar\psi^{0d}$, a bilinear in the $p$ trace-mode fermions $\psi^{0a}\equiv\mathrm{Tr}\Psi^a/\sqrt N$. Large-$N$ scaling: $H_4\sim N^3$, $-9NM\cdot F\sim N\cdot N^2$, $H_0\sim N^3$, but the trace-mode term is $O(N)$, i.e. suppressed by $1/N^2$; it is kept exactly at finite $N$.

### D2.7 Checks

1. **Single matrix ($p=1$, $C=1$):** $X=\Psi^2$, $M=L=\|C\|^2=\langle C,C^{\rm rev}\rangle=1$, so $H=9\mathrm{Tr}[\Psi\Psi\bar\Psi\bar\Psi]-9NN_\Psi+9\mathrm{Tr}\Psi\mathrm{Tr}\bar\Psi+3N(N^2-1)$, which is (D1.4). ✓
2. **Vacuum energy.** All normal-ordered terms annihilate $|0\rangle$, so $E_0(N_\Psi{=}0)=\langle0|H|0\rangle=\|Q|0\rangle\|^2=H_0$. The particle–hole map $\Psi^a_{ij}\to(\Psi^a_{ij})^\dagger$ sends $Q\to-\bar C_{abc}\to C_{abc}$-conjugated $\bar Q$ (using $\mathrm{Tr}[A^TB^TC^T]=-\mathrm{Tr}[CBA]$ for anticommuting entries), so the filled state has the same energy: $E_0(pN^2)=E_0(0)=H_0$.
3. **One-particle sector (exactly solvable).** $H_4$ annihilates one-particle states. On $\Psi^e_{ij}|0\rangle$: $F^{cd}\to\delta^{de}\Psi^c_{ij}|0\rangle$ and $\mathrm{Tr}\Psi^c\mathrm{Tr}\bar\Psi^d\to\delta^{de}\delta_{ij}\mathrm{Tr}\Psi^c|0\rangle$. Hence on gauge-traceless one-particle states $H=H_0-9N\,M$ (acting on the flavor index), and on the $p$ trace-mode states $H=H_0-9N(M-L)$. Levels are $H_0-9N\mu_M$ and $H_0-9N(\mu_M-\mu_L)$ with $\mu$'s the eigenvalues in a common eigenbasis when $[M,L]=0$.
4. **Numerical verification** of (D2.11) as an operator identity on random vectors (relative residual $\le8\times10^{-16}$): $(N,p)=(2,1),(3,1)$ [$C=1$]; $(2,2),(3,2)$ [random complex cyclic $C$; $(3,2)$ has $\dim=2^{18}$]; $(2,3)$ [Chen's $C$ and two random complex $C$]. $\langle0|H|0\rangle=H_0$ verified in all cases.

### D2.8 Specialisation to Chen's three-matrix model

With the $C$ of D2.1 (computed and cross-checked by code):
$$\|C\|^2=\tfrac{16}3,\quad\langle C,C^{\rm rev}\rangle=5,\quad M=\tfrac59\,\mathbb 1+\tfrac{11}9\,\mathbb J,\quad L=\tfrac49\,\mathbb 1+\tfrac{11}9\,\mathbb J,\qquad(\mathbb J_{cd}=1\ \forall c,d),\tag{D2.12}$$
so $M-L=\tfrac19\mathbb 1$, and in the $\mathbb Z_3$ Fourier flavor basis $\{s,\omega,\bar\omega\}$ (with $\Psi^s=(\Psi^1+\Psi^2+\Psi^3)/\sqrt3$) both are diagonal: $\mu_M=(\tfrac{38}9,\tfrac59,\tfrac59)$, $\mu_L=(\tfrac{37}9,\tfrac49,\tfrac49)$. The Hamiltonian is
$$\boxed{\ H=9\sum_{a=1}^3\mathrm{Tr}[X_a\bar X_a]-5N\,N_\Psi-11N\sum_{c,d=1}^3\mathrm{Tr}[\Psi^c\bar\Psi^d]+4\sum_c\mathrm{Tr}\Psi^c\,\mathrm{Tr}\bar\Psi^c+11\sum_{c,d}\mathrm{Tr}\Psi^c\,\mathrm{Tr}\bar\Psi^d+16N^3-15N\ }\tag{D2.13}$$
The two "double-trace" terms are products of length-one traces, $\mathrm{Tr}\Psi^c=\sqrt N\,\psi^{0c}$, i.e. together they are the *one-body* operator $4N\,n_0+33N\,n_{0s}$ ($n_0=\sum_c\psi^{0c}\bar\psi^{0c}$, $n_{0s}$ the occupation of the flavor-symmetric trace mode) — a chemical potential for the three trace-mode fermions, of size $O(N)$ against the $O(N^3)$ of the other terms. Here, e.g., $X_1=\Psi^1\Psi^1+\tfrac13\big(\Psi^1\Psi^2+\Psi^2\Psi^1+\Psi^1\Psi^3+\Psi^3\Psi^1+\Psi^2\Psi^2+\Psi^3\Psi^3+\Psi^2\Psi^3\big)$ and $X_2,X_3$ by cyclic relabeling (note the asymmetric last term: $\Psi^2\Psi^3$ appears, $\Psi^3\Psi^2$ does not). Equivalently $-9NM\cdot F=-38N\,N_s-5N(N_\omega+N_{\bar\omega})$ in the Fourier basis, where $N_s=\mathrm{Tr}[\Psi^s\bar\Psi^s]$ etc. (these individual numbers are *not* conserved; only $N_\Psi$ and the $\mathbb Z_3$ charge are).

**Exact consequences (all confirmed by ED at $N=2$, and $E_0(0),E_0(1)$ at $N=3$ against the values recorded in `research/notes/analytics_and_plots.md`):**
- $E_0(N_\Psi{=}0)=E_0(3N^2)=16N^3-15N$: $98$ ($N{=}2$), $387$ ($N{=}3$).
- One-particle sector: gauge-traceless states at $16N^3-15N-38N=16N^3-53N$ (flavor $s$) and $16N^3-20N$ (flavors $\omega,\bar\omega$); trace-mode states at $16N^3-16N$ (all three, since $M-L\propto\mathbb 1$). At $N=2$: $22,\ 88,\ 96$; at $N=3$: $E_0(1)=273$. ✓
- The trace modes **do not decouple** for $p=3$ (unlike $p=1$): expanding $\Psi^a=\psi^{0a}/\sqrt N+\hat\Psi^a$, $Q$ contains $\tfrac{3}{\sqrt N}C_{a[bc]}\,\psi^{0a}\,\hat\psi^{b x}\hat\psi^{c x}$, and $C_{a[bc]}\ne0$ for Chen's $C$ (e.g. $C_{123}-C_{132}=\tfrac13$). The traceless ($su(N)$-valued) version of the model is therefore a *different* model (differing at relative order $1/N^2$ in $H$ but with a different Hilbert space); the project's ED uses the full $U(N)$ model.

### D2.9 Which parts are Casimirs?

**Symmetry algebra.** By construction $[H,N_\Psi]=0$ and $[H,J^x]=0$ for the gauge $su(N)$ generators $J^x=\sum_a\mathrm{Tr}[\Psi^a[T^x,\bar\Psi^a]]$; Chen's $C$ also gives the discrete cyclic flavor symmetry $\mathbb Z_3$. Numerically at $N=2$: the space of one-body operators $\sum x_{mn}c^\dagger_mc_n$ commuting with $H$ is exactly 4-dimensional and equals $\mathrm{span}(J^1,J^2,J^3,N_\Psi)$ — **no continuous flavor symmetry survives** (in particular no $U(1)$ beyond $N_\Psi$, and no conserved trace-mode number); the cyclic flavor permutation commutes with $H$ exactly, the three transpositions do not ($\|[H,U_{\rm swap}]\|_\infty=6$).

**Gauge Casimir in trace form (general $p$; derived as in D1 and verified numerically for $(N,p)=(2,2),(2,3),(3,2)$):**
$$\hat C_2^{\rm gauge}=\frac{N(p+1)}2\,N_\Psi-\sum_a\mathrm{Tr}\Psi^a\mathrm{Tr}\bar\Psi^a-\sum_{a,b}\mathrm{Tr}[\Psi^a\Psi^b\bar\Psi^b\bar\Psi^a]+\frac12\sum_{a,b}\Big(\mathrm{Tr}[\Psi^a\bar\Psi^a\Psi^b\bar\Psi^b]-\mathrm{Tr}[\Psi^a\bar\Psi^b\Psi^b\bar\Psi^a]\Big).\tag{D2.14}$$
For $p=1$ the last bracket vanishes and (D1.3) is recovered. For $p\ge2$ the Casimir contains *alternating* words $\mathrm{Tr}[\Psi\bar\Psi\Psi\bar\Psi]$ (index topology different from $\mathrm{Tr}[\Psi\Psi\bar\Psi\bar\Psi]$; the two are independent normal-ordered operators at generic $N$), and its $\mathrm{Tr}[\Psi\Psi\bar\Psi\bar\Psi]$ part has the specific flavor structure $\delta^{bd}\delta^{ce}$ in $V^{bced}$.

**Comparison with $H$.** $H_4=9K_{bc;ed}V^{bced}$ contains no alternating words, and $K$ is not proportional to $\delta^{bd}\delta^{ce}$ (for Chen's $C$, $K$ has rank 3 out of 9). Projecting $K$ onto the Casimir structure gives $\kappa\,\delta^{bd}\delta^{ce}$ with $\kappa=\|C\|^2/p^2=\tfrac{16}{27}$, i.e. a component $-\tfrac{16}{3}\hat C_2^{\rm gauge}$ (versus $-9\hat C_2$ for $p=1$), but the remainder — the traceless part of $K$ plus the alternating words that must be added back — is non-zero. The flavor Casimir of $U(p)$, $\sum_{cd}F^{cd}F^{dc}$, is a genuine double trace and cannot help. So, unlike $p=1$, **$H$ is not a polynomial in the Casimirs of its symmetry algebra**; only the quadratic and constant parts are "Casimir-like" ($N_\Psi$ and flavor-rotated number operators).

**Direct finite-$N$ proof ($N=2$).** Simultaneously resolving $N_\Psi$, $\hat C_2^{\rm gauge}$, $J_z$ and the $\mathbb Z_3$ charge splits the 4096-dimensional space into 353 joint sectors (largest dimension 42). In 245 of them the restriction of $H$ has more than one distinct eigenvalue (e.g. $N_\Psi=4$, singlet, $\mathbb Z_3$-fixed: 15 states, 15 distinct energies). Since a function of Casimirs is constant on such sectors, $H$ is not one. Combined with the commutant computation (no hidden continuous symmetry), this settles the question at $N=2$; the structural argument above shows the same for all $N$ at the level of independent trace words.

### D2.10 Consequences for the single-trace bootstrap

1. **Linear EOM.** $[H,\mathrm{Tr}\,w]$ is a linear combination of single-trace words for the $H_4$ and $-9NM\cdot F$ parts; the trace-mode term contributes $\mathrm{Tr}\Psi^c\,[\mathrm{Tr}\bar\Psi^d,\mathrm{Tr}\,w]+\dots$, a product of a length-one trace with a single trace. In an $su(N)$-truncated (traceless) engine this term is absent; in the $U(N)$ model it must be kept as $p$ extra fermionic "letters" $\psi^{0a}$, or dropped at leading large $N$ where it is $O(N^{-2})$ relative.
2. **Operator content.** The words that appear in $H$ are $V^{bced}=\mathrm{Tr}[\Psi^b\Psi^c\bar\Psi^e\bar\Psi^d]$ with the rank-3 flavor tensor $K$, i.e. effectively the three composite matrices $X_a$ and their conjugates; $F^{cd}$; and $\mathrm{Tr}\Psi^c$. A natural minimal letter set for the bootstrap is $\{\Psi^a,\bar\Psi^a\}$ with the flavor structure imposed through $K,M,L$; the $\mathbb Z_3$ Fourier basis diagonalises $M$ and $L$ and block-diagonalises everything by $\mathbb Z_3$ charge.
3. **No Casimir shortcut.** The sector-resolved ground energy $E_0(N_\Psi,\mathbf r)$ is a genuine dynamical quantity; the SUSY floor $\phi(H)\ge0$ must be lifted by constraints that see the sector (open-index operators), exactly as diagnosed for $p=1$.
4. **Analytic anchors** for validating a bootstrap at any $N$: $E_0(0)=16N^3-15N$ and the one-particle levels of D2.8 (exact), plus $E_0(N_\Psi)=0$ inside the BPS window.

### Open items

- Derive $E_0$ in the two-particle sector analytically (a $\binom{3N^2}{2}$-dimensional problem with $H_4$ acting non-trivially); compare with $N=2$ ED ($E_0(2)=5.165\ldots$, irrational) — likely the first genuinely interacting level.
- Decide whether the traceless ($su(N)$) version of the model is the better bootstrap target; check whether it also concentrates at $N=2$ ($2^9=512$ states).
- Generalise the $\mathbb Z_3$-Fourier block structure to the $U(N)$-covariant operator basis.

---

## D3. Sector-resolved bootstrap with adjoint-valued operators (formulation) (2026-09-16)

**Status: formulation used by `src/sector_bootstrap.py`; validated exactly on the single-matrix model at $N=3$ (all sectors) and $N=4$ (15 of 17 sectors); results in `research/notes/sector_bootstrap_results.md`.**

### D3.1 The functional

Fix a fermion-number sector $k$ with projector $P_k$ and dimension $d_k$. The bootstrap variable is a linear functional
$$\phi(O)=\mathrm{Tr}\big[\rho_k\,P_kOP_k\big],\qquad \rho_k=\rho_k^\dagger\ \text{on the $k$-block, not assumed positive},\qquad \phi(\mathbf 1)=1.\tag{D3.1}$$
Restricting to the block is legitimate for bounding $E_0(k)\equiv\min_{|\psi\rangle\in\mathcal H_k}\langle\psi|H|\psi\rangle$: the minimiser is a state in the block, hence a member of this class, so $\min\phi(H)$ over any set of *valid* constraints is a rigorous lower bound. It is strictly stronger than imposing $\phi(N_\Psi)=k$, $\phi(N_\Psi^2)=k^2$ on a functional over the whole space (as `sm_covariant_bootstrap.py` did): it sets $\phi(O)=0$ for every charged $O$ and identifies $\phi(O)=\phi(P_kOP_k)$. Because $[H,N_\Psi]=0$, $P_k[H,O]P_k=[H_k,O_k]$ for charge-neutral $O$, with $H_k=P_kHP_k$.

### D3.2 Operators

*Open words.* For letters $\ell\in\{\Psi^a,\bar\Psi^a\}$ an open word $w=\ell_1\cdots\ell_L$ is the operator-valued $N\times N$ matrix $w_{ij}=\sum(\ell_1)_{im}(\ell_2)_{mn}\cdots(\ell_L)_{rj}$, with charge $q(w)=\#\Psi-\#\bar\Psi$; $w_{ij}$ maps $\mathcal H_k\to\mathcal H_{k+q}$. Its conjugate word is $(w^\dagger)_{ji}\equiv(w_{ij})^\dagger$ (so $\Psi^\dagger=\bar\Psi$ as matrices). **All orderings** of letters are enumerated. (An earlier version of this note claimed that normal-ordered words suffice; that is false: reordering letters with (D2.1) keeps the index contraction pattern of the original word, e.g. $(\Psi\bar\Psi\Psi)_{ij}=N\Psi_{ij}-\Psi_{im}\Psi_{nj}\bar\Psi_{mn}$, and the last term is not a matrix product of letters. The restriction demonstrably weakened the bounds — see `research/notes/sector_bootstrap_results.md` §5.)

*Channels.* For a gauge-invariant functional the Gram matrix of the vectors $\{w_{ij}|\psi\rangle\}$ decomposes into $U(N)$ tensor structures (Cho et al. 2024, eq. 3.5). We impose the two consequences that are valid for **any** positive $\rho_k$, gauge-invariant or not:
$$\text{adjoint channel:}\quad A_{ww'}=\phi\big(\mathrm{Tr}[w^\dagger w']\big)=\sum_{ij}\phi\big((w_{ij})^\dagger w'_{ij}\big)\succeq0\ \ (q(w)=q(w')),\tag{D3.2}$$
$$\text{singlet channel:}\quad S_{ww'}=\phi\big((\mathrm{Tr}\,w)^\dagger\,\mathrm{Tr}\,w'\big)\succeq0\ \ (q(w)=q(w')),\tag{D3.3}$$
(D3.2) is the sum of the diagonal $ij$-blocks of the full Gram matrix, hence PSD; (D3.3) is the Gram matrix of traced words and includes $Q=\mathrm{Tr}[\Psi^aX_a]$ and $\bar Q$ when $L_{\rm sing}\ge3$ (they are added explicitly otherwise), which yields the SUSY floor $\phi(H)=\phi(\bar QQ)+\phi(Q\bar Q)\ge0$. Blocks with $q(w)\ne q(w')$ vanish identically by (D3.1). (D3.2) is what was missing in the earlier gauge-invariant-only formulation: it contains entries such as $\phi(\mathrm{Tr}[\Psi\bar\Psi])$ and $\phi(\mathrm{Tr}[\bar\Psi^\dagger X])$ whose positivity/Cauchy–Schwarz relations encode that annihilators act first on the sector.

*Equations of motion.* $\phi([H,X])=0$ for every charge-neutral operator $X$ that appears as an entry of (D3.2)–(D3.3), plus the charge-neutral traced words up to length $L_{\rm eom}$. (For the single-matrix model all these $X$ are gauge invariant and commute with $H=3N(N^2-1)-9\hat C_2$, so the EOM are empty there; they are non-trivial for $p\ge2$.)

*Ground-state positivity (optional).* $\phi(\mathrm{Tr}[w^\dagger[H,w']])\succeq0$ and its singlet analogue, **only for $q(w)=q(w')=0$**: for charge-changing $w$ the state $w|\Omega_k\rangle$ leaves the sector, where the energy may be lower, so the inequality is not valid.

*Gauge Ward identities (optional).* $\phi([J^x,X])=0$; legitimate because the sector minimiser may be taken gauge invariant (multiplet-averaged). Not needed in the runs so far.

### D3.3 Reduction and SDP

Let $V$ be the span of the Hermitian and anti-Hermitian parts of all operators appearing in the constraints (closed under $\dagger$ by construction). $\phi$ restricted to these operators depends only on the projection of $\rho_k$ onto $V$, and every projection is realised by some Hermitian $\rho_k$, so we may parametrise $\rho_k=\sum_\alpha y_\alpha E_\alpha$ with $\{E_\alpha\}$ an orthonormal Hermitian basis of $V$ and $y\in\mathbb R^r$, $r=\dim V$. Then $\phi(X)=\sum_\alpha y_\alpha\mathrm{Tr}[E_\alpha X]$ is linear in $y$ with computable complex coefficients, and the problem is a standard SDP: minimise $\phi(H_k)$ subject to $\phi(\mathbf 1)=1$, the PSD blocks (D3.2)–(D3.3) (and GS blocks), and the linear EOM. $r$ is the number of SDP variables and is independent of $d_k$ in principle; in the single-matrix model $r\le5$ at level 2 because the touched operators are functions of $\hat C_2$ and $N_\Psi$ on the block.

The basis is obtained from the Gram matrix $\mathcal G_{ab}=\langle O_a,O_b\rangle$ of the touched operators; this costs $O(n_h^2)$ memory and is the step that must be budgeted (see the results note).

### D3.4 Validity summary

Every constraint above is satisfied by the exact sector ground state: (D3.2)–(D3.3) by positivity of the Hilbert-space inner product, EOM by stationarity within the block, GS positivity by the variational principle restricted to sector-preserving perturbations, Ward identities by gauge invariance of the multiplet average. Hence $\min\phi(H)$ is a rigorous lower bound on $E_0(k)$ at every truncation level, and it can only increase with $L_{\rm adj},L_{\rm sing},L_{\rm eom}$.


## D4. Symmetry reduction of the sector bootstrap: invariant functionals and isotypic blocks (2026-09-17)

*Implemented in `src/symmetry_reduction.py` (`IsotypicReducer`) and the `symmetry=True` path of `src/sector_bootstrap.py`. Numerical checks are listed at the end.*

### D4.1 Assumptions and the group

Let $G=SU(2)_{\rm gauge}\times\mathbb Z_p$ act unitarily on the $k$-block $\mathcal H_k$ (dimension $d$), with $[U(g),H_k]=0$. The gauge $SU(2)$ is generated by $J^a=\mathrm{Tr}[\bar\Psi_b\,T^a\Psi_b]$-type charges already in `model['Ja']` (for $N=2$ only; for $N\ge3$ the code currently falls back to the flavor group alone). $\mathbb Z_p$ is the cyclic flavor rotation $\Psi^c\to\Psi^{c+1}$, a symmetry of Chen's $C_{abc}$ (which is cyclically symmetric but **not** $S_3$-symmetric; see the ED note). The Fock-space unitary $U_f$ implementing $\mathbb Z_p$ is built with the fermionic sign of the permutation and verified against $U_fc^\dagger_mU_f^\dagger=c^\dagger_{\pi(m)}$.

**Claim.** $E_0(k)=\min\{\phi(H_k)\}$ may be computed over $G$-invariant functionals only. *Proof.* If $\rho$ is any density matrix on $\mathcal H_k$ (in particular the projector on a ground state), its group average $\bar\rho=\frac1{|G|}\sum_gU(g)\rho U(g)^\dagger$ is again a density matrix with $\mathrm{Tr}[\bar\rho H_k]=\mathrm{Tr}[\rho H_k]$ (and $\bar\rho$ is the multiplet average of the ground state when $\rho$ is a ground-state projector). All constraints of D3 are of the form "$\phi$ is the expectation in a state" and are therefore satisfied by $\bar\rho$ as well; the ones that are not manifestly $G$-covariant (the singlet/adjoint blocks are built from $G$-covariant word sets, so they are) are not used. Hence restricting the SDP variable to invariant $\rho_k$ loses nothing and the minimum is unchanged. $\square$ (The same holds for any subgroup; the compact $SU(2)$ average is over Haar measure.)

The restriction does two things: (i) it **removes variables** — non-invariant components of $\rho_k$ never enter $\phi$ of an invariant operator, and after the reduction every operator that appears is invariant (see D4.3); (ii) it makes the **adjoint-projected channel** (Cho et al. 2024, eq. 3.5) available: for invariant $\rho$ the Gram matrix $\phi\big((w_{ij})^\dagger w'_{kl}\big)$ is a $U(N)$-invariant tensor, so its singlet and adjoint components are separately PSD, giving $\phi(\mathrm{Tr}[w^\dagger w'])-\frac1N\phi\big((\mathrm{Tr}\,w)^\dagger\mathrm{Tr}\,w'\big)\succeq0$ in addition to (D3.2)–(D3.3) (flag `adjoint_projected`).

### D4.2 Isotypic decomposition and the reduced vector

Decompose $\mathcal H_k=\bigoplus_R\,\mathbb C^{m_R}\otimes V_R$ into isotypic components ($R$ = irrep label $(j,\omega)$ with $j$ the $SU(2)$ spin and $\omega$ the $\mathbb Z_p$ eigenvalue; $m_R$ = multiplicity, $\dim V_R=2j+1$). By Schur's lemma every invariant operator has the form
$$\rho=\bigoplus_R\sigma_R\otimes\mathbf 1_{V_R},\qquad \sigma_R\in\mathbb C^{m_R\times m_R},\tag{D4.1}$$
so the real dimension of the space of invariant Hermitian operators is $D=\sum_Rm_R^2$ (versus $d^2$ for all Hermitian operators). With an orthonormal basis $|R,a,\mu\rangle$ ($a=1..m_R$ multiplicity index, $\mu=1..\dim V_R$ magnetic index, chosen so that the group acts only on $\mu$ and identically for all $a$), for **any** operator $X$
$$\phi(X)=\mathrm{Tr}[\rho X]=\sum_R\mathrm{Tr}\big[\sigma_R\,\tilde X^R\big],\qquad \tilde X^R_{ba}\equiv\sum_\mu\langle R,b,\mu|X|R,a,\mu\rangle .\tag{D4.2}$$
The concatenation of the $\tilde X^R$ is the **reduced vector** $\tilde X\in\mathbb C^D$ (`reduce`). Two facts follow: $\phi$ depends on $X$ only through $\tilde X$; and $\widetilde{X^\dagger}=(\tilde X)^\dagger$ blockwise (`adjoint`), so the Hermitian/anti-Hermitian parts can be taken in the reduced space. The inverse on invariant operators is $X=\sum_R\frac1{\dim V_R}\sum_\mu V_\mu\tilde X^RV_\mu^\dagger$ (`lift`; checked to round-trip to $10^{-15}$).

*Construction of the basis (N=2).* $J_z$, $\hat C_2=\vec J^2$ and the Hermitian and anti-Hermitian parts of $U_f$ commute; a generic real combination is diagonalised once and the eigenvectors are labelled by $(J_z,\hat C_2,\omega)$. Highest-weight vectors ($J_z=j$, $\hat C_2=j(j+1)$) with the same $(j,\omega)$ span the multiplicity space; they are orthonormalised (QR) and the lower components are generated with $J_-$, normalised by the $a$-independent factor $\sqrt{j(j+1)-m(m-1)}$, so that $J^a$ acts identically on every copy — this is what makes (D4.2) true with the *same* $\sigma_R$ for all $\mu$. Completeness $\sum_Rm_R\dim V_R=d$ is asserted.

### D4.3 The reduced SDP

Parametrise the invariant $\rho$ by real coordinates in an orthonormal basis $\{e_\alpha\}$ of the real span (inside $\mathbb R^{2D}\ni[\mathrm{Re}\,\tilde X,\mathrm{Im}\,\tilde X]$) of the Hermitian and anti-Hermitian parts of all reduced vectors that occur (identity, $H_k$, Gram entries, commutators, observables). Then $y_\alpha=\phi(e_\alpha)$ and for any touched operator $\phi(X)=\sum_\alpha\big(e_\alpha\!\cdot\!\widetilde{X_h}+i\,e_\alpha\!\cdot\!\widetilde{X_a}\big)y_\alpha$, $X=X_h+iX_a$. The number of SDP variables is $r=\dim{\rm span}\le D$, and once $r=D$ the functional is *completely general* within the invariant class: further improvement of the bound can then come only from additional constraints, not from resolving $\rho$ better. (This is what happens for the three-matrix model at $k=2$, $D=90$, already at level 3.)

Because only the Gram matrix of the touched span is needed, it is **accumulated in streaming fashion** as a $2D\times2D$ matrix $M=\sum e\,e^{\!\top}$ over unit-normalised rows, and likewise for the span of EOM rows; nothing of size $n_h\times D$ is stored. An orthonormal basis is read off from the eigenvectors of $M$ with $\lambda>10^{-10}\lambda_{\max}$, and the EOM constraints become "$\phi$ vanishes on an orthonormal basis of the EOM span" — this replaces the earlier SVD pruning of $10^4$ redundant rows. Only the upper triangle of each Gram cone is kept, since $\phi(X^\dagger)=\phi(X)^*$ implies $c(X^\dagger)=\overline{c(X)}$ for the coefficient rows.

*A numerical trap, recorded because it produced an infeasible SDP.* Commutators $[H_k,X]$ with $X$ Hermitian are anti-Hermitian up to round-off; their "Hermitian part" has norm $\sim10^{-12}$, and an absolute cutoff of $10^{-12}$ let some of these through, after which unit-normalisation turned round-off into a spurious constraint direction that the exact ground state violates ($\phi\ne0$ at the $4\times10^{-2}$ level). The cut is now **relative**, $\|X_{h,a}\|\le10^{-9}\|X\|$, in all three code paths. The diagnostic that caught it, `SectorSDP.exact_coordinates` + `check_feasibility` (evaluate every constraint on the exact multiplet-averaged ground state), is kept and should be run whenever a new constraint type is added.

### D4.4 What the symmetry does and does not buy

* The touched operators of D3 are already gauge singlets (traces), so $SU(2)$ invariance alone does not change $r$; it changes the *cost*: the linear algebra is done in $D=\sum m_R^2$ components instead of $d^2$ ($4356\to90$ at $k=2$; $48\,400\to666$ at $k=3$; $245\,025\to2727$ at $k=4$ for $N=2$, $p=3$).
* The $\mathbb Z_3$ grading of the Fourier letters $\Psi^{(m)}=3^{-1/2}\sum_c\omega^{mc}\Psi^c$ splits every Gram cone into three ($z$-charge) blocks and removes the $z\neq0$ components of $\rho$ from the problem ($r$: $267\to87$ at $k=2$, level 2). This is the same bound with a three times smaller SDP.
* Gauge Ward identities $\phi([J^a,X])=0$ are automatically satisfied and add no rows (verified: identical EOM count with `ward=True`).
* Q-descendant words $[Q,w]$ at fixed length are linear combinations of the words already present, so they cannot change the span of the Gram cones (the cone is invariant under a change of basis of the word set); they would matter only as a *sparser* basis at a higher length, which is not pursued.

### D4.5 Checks performed

1. Unreduced and reduced SDPs agree on every regression value: single-matrix $N=2$ (all sectors) and $N=3$ (all ten sectors, $72,45,18,0,\dots$); three-matrix $N=2$: $k=1\to22.000$, $k=2$ level 2 $\to0$ with and without GS blocks, with and without the adjoint-projected channel.
2. `lift∘reduce` and `adjoint` consistency to $10^{-13}$ on all $8058$ commutators of the $k=2$ level-3 set; the exact ground-state functional satisfies all EOM rows to $10^{-14}$ and all cones to $-10^{-13}$ after the fix of D4.3.
3. Eigenstate constraints in reduced form: feasible at $E=E_0$, infeasible at $E_0+0.05$ ($k=2$, level 2).
