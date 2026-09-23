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


## D6. Refined Witten index of the fermionic matrix models (2026-09-19)

*Implemented in `src/bps_index.py`; results in `research/notes/bps_index_results.md`.*

**Setup.** Fock space $\mathcal F=\wedge^\bullet V$, $V=\mathbb C^p\otimes\mathbf{adj}_{U(N)}$, spanned by the modes $\Psi^{(m)}_{ij}$ ($m$ the $\mathbb Z_p$ Fourier flavor charge, $(i,j)$ the $U(N)$ weight $u_i/u_j$). The gradings commuting with $Q=C_{abc}\mathrm{Tr}\Psi^a\Psi^b\Psi^c$ are $U(N)$, the flavor $\mathbb Z_p$ (Chen's $C$ is cyclic), and $N_\Psi$ **mod 3** (since $[N_\Psi,Q]=3Q$). The Fock character is
$$Z(t,z;u)=\prod_{m=0}^{p-1}\prod_{i,j=1}^N\big(1+t\,z^m\,u_i/u_j\big),\tag{D6.1}$$
and $n(k,\lambda,\omega)$, the multiplicity of the $U(N)$ irrep $\lambda$ ($\sum\lambda_i=0$; all states are $U(1)$-neutral) in the $(N_\Psi{=}k,\ \omega)$ subspace, is obtained from the weight multiplicities $W(k,\omega,\mu)=[t^kz^\omega u^\mu]Z$ by the Weyl alternating sum
$$n(k,\lambda,\omega)=\sum_{\sigma\in S_N}\mathrm{sgn}(\sigma)\,W\big(k,\omega,\lambda+\delta-\sigma(\delta)\big),\qquad\delta=(N-1,\dots,1,0).\tag{D6.2}$$
(D6.1)–(D6.2) are evaluated exactly: $W$ by multiplying the $pN^2$ linear factors as integer shifts of an array indexed by $(k,\omega,e)$, $e$ the $SU(N)$-torus exponents ($u_N=1/u_1\cdots u_{N-1}$; each $e_i\in[-2p(N-1),2p(N-1)]$, stored modulo $M=4p(N-1)+1$ without aliasing).

**Index.** The complex $(c,\lambda,\omega)$ is $\cdots\to\mathcal H_{k}\xrightarrow{Q}\mathcal H_{k+3}\to\cdots$ restricted to $\lambda,\omega$ with $k\equiv c$; its Euler characteristic is
$$I_{c,\lambda,\omega}=\sum_{k\equiv c\,(3)}(-1)^{k}\,n(k,\lambda,\omega)=\pm\sum_{j}(-1)^j\dim H_Q^{c+3j}\big|_{\lambda,\omega},\tag{D6.3}$$
(the overall sign per complex is conventional). Hence $|I|\le\#\{\text{BPS multiplets in the complex}\}$ with equality iff the cohomology is concentrated in a single degree.

**Closed form of the class totals.** Summing over all irreps and flavors, $\sum_k(-1)^k\binom{n}{k}z^k=(1-z)^n$ with $n=pN^2$ gives, with $\varpi=e^{2\pi i/3}$ and $1-\varpi=\sqrt3\,e^{-i\pi/6}$,
$$I_c(N)=\frac13\sum_{r=0}^{2}\varpi^{-rc}(1-\varpi^r)^{n}=\frac23\,3^{n/2}\cos\!\Big(\frac{2\pi c}{3}+\frac{\pi n}{6}\Big)\;\xrightarrow{\,p=3\,}\;\frac23\,3^{3N^2/2}\cos\!\Big(\frac{2\pi c}{3}+\frac{\pi N^2}{2}\Big).\tag{D6.4}$$
Writing $k=c+3j$ and $k_*=n/2$, (D6.4) is $\propto\cos\big(\pi(k-k_*)/3\big)$: the Turiaci–Witten profile about half filling. Checks: $N=2$: $486,-243,-243$ (= ED: $972$ BPS states at $k=5,6,7$ with $243{:}486{:}243$); $N=3$: $0,-3^{13},3^{13}$; $N=4$: $2\cdot3^{23},-3^{23},-3^{23}$. The refined tables reproduce these totals exactly.

**Verified.** At $N=2$ the refined index equals the BPS multiplet count in every complex (36 complexes, spins $0..3$, three flavor charges) computed by exact diagonalisation with the isotypic reducer of D4 — index saturation, as concentration requires and as the $\{5,6,7\}$ window forces.

**Not derived here.** The degree in which the cohomology of each complex sits; the per-irrep large-$N$ asymptotics (a saddle point of (D6.1) against Schur functions — open); the analogous formula for the traceless ($su(N)$) variant, which would differ by the trace-mode factors $(1+tz^m)^{N}$.


## D7. The quadratic gauge Casimir as a trace polynomial, and irrep-resolved constraints (2026-09-21)

*Implemented as `trace_algebra.casimir(p)` / `casimir_value(N, lam)`; used by the `casimir=` option of `TraceSDP`.*

**Generators.** Under $\Psi^a\to U\Psi^aU^\dagger$ the creation operators $c^\dagger_{(a,i,j)}=\Psi^a_{ij}$ transform as $c^\dagger_{(a,i,j)}\to U_{ii'}\bar U_{jj'}c^\dagger_{(a,i',j')}$. For $U=e^{i\theta T}$ the one-body generator is $E(T)=\sum c^\dagger_\alpha M_{\alpha\beta}c_\beta$ with $M=T\otimes\mathbf 1-\mathbf 1\otimes T^{T}$, i.e.
$$E(T)=\sum_{a}\Big[\sum_{i,i',j}T_{ii'}\,\Psi^a_{ij}\bar\Psi^a_{ji'}-\sum_{i,j,j'}T_{j'j}\,\Psi^a_{ij}\bar\Psi^a_{j'i}\Big]
=\sum_{ii'}T_{ii'}\Big[(\Psi^a\bar\Psi^a)_{ii'}+(\bar\Psi^a\Psi^a)_{ii'}\Big]-N\,\mathrm{Tr}\,T\cdot p ,$$
where the second form uses $\sum_i\Psi_{ij}\bar\Psi_{j'i}=-(\bar\Psi\Psi)_{j'j}+N\delta_{jj'}$ (the anticommutator (D2.1)) and matrix words are in operator order. For traceless $T$ the constant drops: $E(T)=\sum_{ii'}T_{ii'}M_{ii'}$ with the matrix word
$$M=\sum_{a=1}^p\big(\Psi^a\bar\Psi^a+\bar\Psi^a\Psi^a\big),\qquad \mathrm{Tr}\,M=\sum_a\big(N_\Psi^a+(N^2-N_\Psi^a)\big)=pN^2 .$$

**Casimir.** With $\mathrm{Tr}\,T^aT^b=\tfrac12\delta^{ab}$ and the completeness relation $\sum_aT^a_{ii'}T^a_{kk'}=\tfrac12(\delta_{ik'}\delta_{i'k}-\tfrac1N\delta_{ii'}\delta_{kk'})$,
$$\boxed{\ \hat C_2=\sum_aJ^aJ^a=\tfrac12\,\mathrm{Tr}[M^2]-\frac{(\mathrm{Tr}M)^2}{2N}=\tfrac12\sum_{a,b}\mathrm{Tr}\big[(\Psi^a\bar\Psi^a+\bar\Psi^a\Psi^a)(\Psi^b\bar\Psi^b+\bar\Psi^b\Psi^b)\big]-\tfrac12p^2N^3\ }\tag{D7.1}$$
a trace polynomial of degree four, flavor-$U(p)$ invariant (hence identical in the $\mathbb Z_p$ Fourier letter basis). The algebra of D5/M1 canonicalises it to $\sum_a\big(N\,\mathrm{Tr}[\bar\Psi^a\Psi^a]-\mathrm{Tr}\bar\Psi^a\mathrm{Tr}\Psi^a-\mathrm{Tr}[\bar\Psi^a\bar\Psi^a\Psi^a\Psi^a]\big)+\sum_{a\ne b}(\dots)$, which for $p=1$ is (D1.3). **Verified** against the explicit $\sum_aJ^aJ^a$ of `fermion_matrix_model` to $10^{-16}$ at $(N,p)=(2,3),(2,1),(3,1)$, in both letter bases. Eigenvalue on the $U(N)$ irrep $\lambda$ ($\sum\lambda_i=0$): $c_\lambda=\tfrac12\sum_i\lambda_i(\lambda_i+N+1-2i)$ ($=j(j+1)$ at $N=2$; $3$ for the adjoint at $N=3$).

**Irrep-resolved functional.** A functional supported on the $\hat C_2=c_\lambda$ eigenspace satisfies $\phi((\hat C_2-c_\lambda)X)=\phi(X(\hat C_2-c_\lambda))=0$ for all $X$ — linear rows of exactly the same form as the sector rows $\phi((N_\Psi-k)X)=0$. Combined with the BPS rows this gives a **BPS-exclusion test per $(k,\hat C_2)$**. Legitimacy: if a BPS state exists in sector $k$ and irrep $\lambda$, its multiplet average is a functional obeying all rows, so infeasibility certifies absence in that $(k,c_\lambda)$ cell; different irreps sharing a Casimir value are not separated (they could be, with the cubic Casimir, also a trace polynomial).

**First results (2026-09-21, `research/notes/bps_exclusion_results.md` §4).** At $N=2$ the Casimir rows do not extend the reach ($k=4$ not excluded in any spin at level 3). At $N=3$, where the plain sector $k=7$ is not excludable, the level-2 test excludes its large-Casimir cells $C_2=36,38,42,48$ — an exclusion frontier in the $(k,C_2)$ plane that penetrates the window at high Casimir, which is where the index says the largest complexes live.

### D7.4 Which irreps occur, and the closed form for $C_2^{\max}$ (2026-09-23)

`scripts/list_casimirs.py` enumerates the $U(N)$ content of the $p$-flavour Fock space and its Casimirs.

**Which highest weights occur.** Three conditions, the third being the subtle one:

1. *Traceless.* Each mode $\Psi^a_{ij}$ has weight $e_i-e_j$, of total zero, so every Fock state has total weight
   zero and only $\sum_i\lambda_i=0$ appears.
2. *Dominant*, $\lambda_1\ge\cdots\ge\lambda_N$, integers.
3. *Dominated by the maximal weight*: $\lambda\le\lambda^*$, i.e. $\sum_{i\le m}\lambda_i\le\sum_{i\le m}\lambda^*_i$
   for every $m$, where $\lambda^*_i=p(N+1-2i)$.

The box bound $|\lambda_i|\le p(N-1)$ — from $\lambda_i=(\text{row }i)-(\text{col }i)$ with $n_{ii}$ cancelling,
leaving at most the $p(N-1)$ off-diagonal modes of row $i$ — is **necessary but not sufficient**. For example
$(12,12,0,-12,-12)$ at $N=5,p=3$ respects the box yet cannot occur: row 1 full forces $n_{12}=p$, contradicting
column 2 empty. Dominance is the exact criterion, verified against the Weyl alternating sum of D6: it reproduces
the occurring set exactly, $25/25$ at $N=3$ and $213/241$ at $N=4$. *(Verified at $N=3,4$; not proved here.)*

**Casimir.** $c_2(\lambda)=\tfrac12\sum_i\lambda_i(\lambda_i+N+1-2i)$ (D7.1).

**Maximum.** At $\lambda^*$ one has $\lambda^*_i+N+1-2i=(p+1)(N+1-2i)$, so

$$C_2^{\max}=\tfrac12\,p(p+1)\sum_{i=1}^N(N+1-2i)^2=\frac{p(p+1)\,N(N^2-1)}{6}, \tag{D7.6}$$

using $\sum_i(N+1-2i)^2=N(N^2-1)/3$. For $p=1$ this is $N(N^2-1)/3$, which is **Chen 2025 (2.18)**, so (D7.6)
is its multi-flavour generalisation; for $p=3$ it gives $C_2^{\max}=2N(N^2-1)$, i.e. $48,120,240$ at $N=3,4,5$
(all three confirmed numerically).

**Counts.** $N=3$: 25 irreps, 16 distinct Casimirs. $N=4$: 213 irreps, 69 distinct. $N=5$: 2131 irreps, 208
distinct.

## D8. Exact $Q$-cohomology by weight-space decomposition, and Kostant inversion to per-irrep BPS counts (2026-09-22)

**Goal.** Compute the exact number of BPS states of the three-matrix model at finite $N$, per degree $k=N_\Psi$,
per $U(N)$ irrep and per $\mathbb Z_p$ flavour charge — without the bootstrap and without diagonalising anything.

**Assumptions.** Only that $H=\{Q,\bar Q\}$ with $Q^2=0$ (verified numerically in `tests/test_cohomology.py` and
established in D2), so that BPS states ($H|\psi\rangle=0$) are in bijection with $Q$-cohomology classes:
$\mathcal H_{\rm BPS}\cong\ker Q/\operatorname{im}Q$. No large-$N$, no truncation, no conjecture.

**Variables.** $N$ = gauge rank, $p$ = number of flavours (here $p=3$), modes $\Psi^a_{ij}$ with
$a\in\{0,\dots,p-1\}$, $i,j\in\{1,\dots,N\}$; $k=N_\Psi$ = fermion number; $\lambda$ = a $U(N)$ weight (integer
$N$-vector, here always of total charge $0$); $\mu$ = a dominant weight, i.e. an irrep highest weight;
$K(\mu,\lambda)$ = multiplicity of the weight $\lambda$ in the irrep $\mu$; $h^k_\mu$ = number of BPS multiplets
of irrep $\mu$ at degree $k$; $W_\lambda(k)$ = the weight-$\lambda$ subspace of $\mathcal H_k$.

### D8.1 The complex splits over weights

$Q=\sum_{abc}C_{abc}\sum_{ijk}\Psi^a_{ij}\Psi^b_{jk}\Psi^c_{ki}$ is a $U(N)$ singlet: every term is a closed index
loop, so each raised index $i$ is matched by a lowered $i$. Hence $[\,\mathfrak h,Q\,]=0$ for the Cartan
$\mathfrak h$, and $Q$ preserves the weight. Since $[N_\Psi,Q]=3Q$, for each weight $\lambda$ we get a complex

$$\cdots\xrightarrow{\;Q\;}W_\lambda(k)\xrightarrow{\;Q\;}W_\lambda(k+3)\xrightarrow{\;Q\;}\cdots,\qquad
\dim H^k(W_\lambda)=\dim W_\lambda(k)-\operatorname{rank}Q_k-\operatorname{rank}Q_{k-3}. \tag{D8.1}$$

The mode $\Psi^a_{ij}$ carries weight $e_i-e_j$, so a state's weight depends only on the *occupation matrix*
$n_{ij}\in\{0,\dots,p\}$ counting filled flavours of $(\cdot,i,j)$: the weight is $\lambda_i=\sum_j n_{ij}-\sum_j n_{ji}$
(row sums minus column sums). Enumerating $W_\lambda(k)$ is therefore a small search over $n_{ij}$ followed by
$\prod_{ij}\binom{p}{n_{ij}}$ flavour choices, never a scan of the $2^{pN^2}$ Fock space.

**Why this is the decisive gain.** At $N=3$, $\dim\mathcal H_{13}=2\times10^7$, hopeless for a rank. But
$\dim W_{(6,0,-6)}(13)=126$. The cost is set by the weight, not by the sector.

Ranks are computed by Gaussian elimination over $\mathbb F_P$, $P=2^{31}-1$. A mod-$P$ rank can only *under*-estimate
the rational rank (never over-estimate), which would inflate $\dim H^k$; all reported ranks were re-checked against
a second prime $2147483629$ and against floating point, with exact agreement.

### D8.2 Kostant inversion: from weights to irreps

A weight space collects contributions from every irrep containing that weight:

$$\dim H^k(W_\lambda)=\sum_{\mu\ \ge\ \lambda}K(\mu,\lambda)\,h^k_\mu, \tag{D8.2}$$

where $\mu\ge\lambda$ means $\mu-\lambda$ is a non-negative integer combination of positive roots. For the
**maximal** weight of the Fock space there is no $\mu>\lambda$, so $H^k(W_\lambda)$ *is* the multiplicity space and
$h^k_\lambda=\dim H^k(W_\lambda)$ directly. For lower weights (D8.2) is triangular and is solved from the top down.

The weight multiplicities come from Kostant's formula

$$K(\mu,\lambda)=\sum_{w\in W}\operatorname{sgn}(w)\,\mathcal P\!\big(w(\mu+\rho)-(\lambda+\rho)\big),\qquad \rho=(N-1,N-2,\dots,0), \tag{D8.3}$$

with $W=S_N$ and $\mathcal P$ the Kostant partition function (the number of ways to write a vector as a
non-negative integer combination of the positive roots $e_i-e_j$, $i<j$), evaluated by a recursion pruned with the
fact that such combinations have non-negative partial sums. *(Note $\rho$ may be shifted by any multiple of
$(1,\dots,1)$: $w$ permutes, the constant is $W$-invariant, and it cancels between $w(\mu+\rho)$ and $\lambda+\rho$.)*
Implementation `cohomology.kostant`; verified by $\sum_\lambda K(\mu,\lambda)=\dim\mu$ (Weyl) for every irrep used.

**Closure.** The inversion needs every $\mu>\lambda$ with $h_\mu\neq0$. Higher weights are more extreme and so have
*smaller* weight spaces: the set of affordable weights is automatically closed upward, and the peel over the eight
$N=3$ weights with $C_2\ge35$ is exact and complete, not a truncation. This was checked exhaustively — those eight
are all the irreps of the $N=3$ Fock space with $C_2\ge35$.

### D8.3 Checks performed

1. **$N=2$ against ED.** Cohomology supported only at $k=5,6,7$; per-weight peeling reproduces the ED multiplet
   counts exactly ($j=3\to9$, $j=2\to18$, $j=1\to27$, $j=0\to9$ at $k=5$; total $243$).
2. **Weight dimensions.** $\dim W_\lambda(k)=\sum_\mu m(k,\mu)K(\mu,\lambda)$ for every $(\lambda,k)$, against the
   independently computed multiplicities of D6 — exact for all eight weights.
3. **Euler characteristic.** For every peeled irrep and every class $c$,
   $\sum_{k\equiv c}(-1)^k h^k_\mu = I_{c,\mu}$, the refined index of D6 — exact in all cases, and in all $27$
   flavour-refined complexes $(\lambda,w,c)$ of the top three irreps.
4. **Non-negativity.** Every peeled $h^k_\mu\ge0$ (a strong constraint: the peel subtracts large numbers).
5. **Particle–hole.** $h^k=h^{pN^2-k}$ emerges, unimposed.
6. **Rank robustness.** Two primes plus floating point agree on every rank.

### D8.4 Result and its interpretation

At $N=3$, for all eight irreps with $C_2\ge35$, the cohomology is supported on exactly
$k\in\{12,13,14,15\}$. Since the complexes are graded by $c=k\bmod3$:

$$c=1:\ \{13\}\quad\text{(single degree)},\qquad c=2:\ \{14\}\quad\text{(single degree)},\qquad c=0:\ \{12,15\}.$$

**Proven** (by exact computation, for these irreps at $N=3$): the $c=1$ and $c=2$ complexes — the ones carrying
macroscopic index — are *index-saturated*, so R-charge concentration holds exactly there and the refined index of D6
counts the BPS states without error. The $c=0$ complexes have vanishing index yet non-zero cohomology, in
cancelling pairs at $k=12$ and $k=15$: $149{,}526$ index-invisible BPS states, $19.7\%$ of the $758{,}808$ total.

**Not proven, and not claimed:** that this persists at larger $N$; that the pattern $1:3:3:1$ observed for the
maximal-Casimir family is general (it already fails at $(5,0,-5)$, where $h^{12}:h^{13}=38:486$); or anything about
the low-Casimir bulk of the spectrum, which this computation has not reached.

**Limiting cases.** At $N=2$ the method reproduces ED exactly (check 1). For the maximal weight the Kostant step is
the identity and (D8.2) degenerates correctly to $h=\dim H(W)$.

### D8.5 $N=4$ (2026-09-23)

The enumeration of $W_\lambda(k)$ must be pruned to reach $N=4$: the direct product over the $N^2$ cells of the
occupation matrix is $(p+1)^{N^2}=4^{16}\approx4.3\times10^9$. `cohomology.occupation_matrices` instead recurses
over rows, pruning on (i) the remaining $k$ budget against the remaining capacity $pN(N-m)$, (ii) for each fixed
row $i'<m$, the requirement that its column finish at $r_{i'}-\lambda_{i'}$, reachable only if the shortfall lies
in $[0,\,p(N-m)]$, and (iii) for each unfixed row $t\ge m$, that $r_t=c_t+\lambda_t$ remain a legal row sum.
Validated exhaustively: $\sum_\lambda\dim W_\lambda(k)=\binom{pN^2}{k}$ for every $k$ at $N=3$, and agreement with
the D6 generating function for every $k$ at the four highest $N=4$ weights.

For the maximal weight $\lambda=(9,3,-3,-9)$ there is nothing above it, so (D8.2) is trivial and the weight-space
cohomology is already the multiplet count. In 6 s at 0.33 GB (against $\dim\mathcal H_{24}=3.2\times10^{13}$):

$$h^{22},\dots,h^{26} \;=\; 81,\;324,\;486,\;324,\;81 \;=\; 81\binom{4}{j},$$

with all other degrees empty. Together with $N=2$ ($1{:}2{:}1$) and $N=3$ ($1{:}3{:}3{:}1$ at the top of the
spectrum) this supports, but does not prove, a window of width $N+1$ centred at $k=pN^2/2$ with profile
$\binom{N}{j}$.

**Proven at $N=4$ for this irrep** (exact integers, Euler characteristics matching the refined index in all nine
$(c,w)$ complexes): the $c=0$ complex is saturated (cohomology at $k=24$ alone), while the $c=1$ and $c=2$
complexes have cohomology in *two* degrees, $\{22,25\}$ and $\{23,26\}$, with **non-vanishing** index
($-81$ per flavour charge) and $135$ multiplets against $|I|=81$. This is stronger than the $N=3$ statement, where
saturation failed only in the zero-index complexes: at $N=4$ it fails where the index is non-zero. Over the irrep,
$1296$ multiplets versus $972$ counted by the index, a $25\%$ excess.

**Not proven:** that the excess persists or grows with $N$ (the two data points, $19.7\%$ at $N=3$ and $25\%$ at
$N=4$, are over different portions of the spectrum and are not directly comparable), or anything about the
low-Casimir bulk at $N=4$.

### D8.6 Blocked and flavour-blocked exact ranks; thirteen irreps at $N=4$ (2026-09-23)

Two devices move the rank wall far enough to reach weight spaces of $3\times10^4$ (from $924$), at $\sim30$ min
and $6.1$ GB for the whole $N=4$ run.

**(a) Blocked elimination with an exact BLAS update** (`rank_mod_p_blocked`). Right-looking LU with partial
pivoting, skipping columns in which no pivot is found (the matrices are rank-deficient by construction — that is
what is being measured). Per panel of $b$ columns: factor the panel, solve $U_{12}=L_{11}^{-1}A_{12}$, then update
the trailing block with the single matmul $A_{22}\mathrel{-}=L_{21}U_{12}$, reduced mod $P$. Exactness comes from
choosing $P$ small enough that float64 holds every partial sum of the matmul as an exact integer,

$$b\,(P-1)^2 < 2^{53},\qquad P\approx2^{20},$$

so `dgemm` performs integer arithmetic. No rescaling of pivot rows is done anywhere: normalising the pivot row
across the panel only, or normalising before rather than after the elimination (the two operations do not
commute), both give wrong ranks, and both were caught by the plain routine on random matrices of known rank.

**(b) $\mathbb Z_p$ flavour blocking** (`complex_cohomology_blocked`). $Q$ commutes with the flavour rotation
$\sigma$, so each weight-space complex splits into $p$ blocks of $\approx n/p$: a $p^2$-fold saving in both time
and memory. To keep it exact, work over $\mathbb F_P$ with $P\equiv1\pmod p$, where a primitive $p$-th root of
unity $\omega$ exists. Orbits of $\sigma$ on a weight basis have size $1$ or $p$, with signs $d_t$ defined by
$\sigma^t|s_0\rangle=d_t|s_t\rangle$; the charge-$w$ eigenvector on an orbit is $v_w=\sum_t\omega^{-wt}d_t|s_t\rangle$,
a fixed point necessarily has $d=+1$ and contributes only to $w=0$, and from $|u_r\rangle=(d_r/p)\sum_w\omega^{wr}v'_w$
together with $Qv_w=\sum_t\omega^{-wt}\sigma^t(Q|s_0\rangle)$ one gets

$$(Q_w)_{O',O}=\sum_r a_{u_r}\,d_r\,\omega^{wr},\qquad Q|s_0\rangle=\sum_u a_u|u\rangle. \tag{D8.4}$$

Validated at $N=2,3,4$ against the dense computation: identical totals, identical per-charge counts, and the block
ranks sum to the dense rank for every $k$.

**Result.** All thirteen $N=4$ irreps with $C_2\ge107$ (`results/data/cohomology_N4_top.json`), totalling
$3\,008\,825\,568$ BPS states exactly, of which $19.0\%$ are invisible to the index. Uniformly across all
thirteen: the window is exactly $k\in\{22,\dots,26\}$; the $c=0$ complex is saturated at $k=24$; and the $c=1,2$
complexes are never saturated, occupying $\{22,25\}$ and $\{23,26\}$ with non-vanishing index.

**Proven / not proven.** The $\binom Nj$ profile holds for eleven of the thirteen and **fails** for the two at
$C_2=107$ ($114,1572,2916,1572,114$), exactly as it fails at $(5,0,-5)$ for $N=3$; in both cases the centre is
unchanged and the edges are suppressed. So the profile is a feature of the high-Casimir top, not a theorem. The
*width* $N+1$ centred at $pN^2/2$ has held for all 25 irreps computed across $N=2,3,4$ without exception, but
remains a conjecture: nothing here reaches the low-Casimir bulk, where the weight spaces run to $\sim10^9$.


## D9. The maximal-weight BPS sector in closed form, and a rank-vs-degree no-go for concentration (2026-09-23)

**Setting.** $u(N)$ adjoint complex fermions $\Psi^a_{ij}$, $a=1..p$, with
$Q=\sum C_{a_1\ldots a_q}\mathrm{Tr}[\Psi^{a_1}\cdots\Psi^{a_q}]$, $q$ odd (even $q$ gives
$\mathrm{Tr}[\Psi^q]\equiv0$ by cyclicity, and $Q^2=0$ is automatic precisely for odd $q$). Write
$r=\mathrm{rank}\,G$ ($=N$ for $u(N)$).

### D9.1 The maximal-weight complex factorises over Cartan slots

At the maximal weight $\lambda^*_i=p(N+1-2i)$ every off-diagonal mode is saturated, so the only unoccupied modes
are the diagonal ones. $Q$ creates a closed index loop; any loop entering an index cannot leave it (its outgoing
edges are occupied), so the only surviving loops are **self-loops** at a single index $i$, creating $q$ diagonal
modes of distinct flavours. Hence $Q=\sum_{i=1}^{r}Q_i$ with $Q_i$ acting only on slot $i$'s flavour space, and
$Q_i$ is wedging with the totally antisymmetric part of $C$,

$$\omega=\mathrm{Alt}(C)\in\Lambda^q\mathbb C^p,\qquad Q_i=\omega\wedge-\ \text{ on }\ \Lambda^\bullet\mathbb C^p .$$

The complex is a tensor product of $r$ identical factors, so by Künneth

$$Z_{\max}(t)=t^{k_0}\big(z_{\rm slot}(t)\big)^{r},\qquad
z_{\rm slot}(t)=\sum_j \dim H^j\!\left(\Lambda^\bullet\mathbb C^p,\ \omega\wedge\right) t^j . \tag{D9.1}$$

**Verified** against every measured $u(N)$ case: $(p,q)=(1,3),(2,3),(3,3),(4,3),(3,5)$ at $N=2,3,4,5$ — ten cells,
exact agreement including $27,81,81,27$; $81,324,486,324,81$; $243,1215,2430,2430,1215,243$;
$27,162,405,540,405,162,27$; $1,6,15,20,15,6,1$.

### D9.2 Width

If $z_{\rm slot}$ is supported on $w_s$ degrees then (D9.1) gives window width

$$W=r\,(w_s-1)+1. \tag{D9.2}$$

This reproduces all the empirical patterns: $p=1,3$ with $q=3$ have $w_s=2$, giving $W=r+1$ (Chen's $(1+q)^N$);
$p=2,4$ have $w_s=3$, giving $W=2r+1$; $q>p$ forces $\omega=0$, $z_{\rm slot}=(1+t)^p$, $w_s=p+1$, $W=rp+1$.

### D9.3 $w_s\ge2$, hence a no-go

Let $E_c=\sum_{j\equiv c\ (\mathrm{mod}\ q)}(-1)^j\binom pj$ be the Euler characteristic of the class-$c$
subcomplex; $E_c\ne0$ forces $H^j\ne0$ for some $j\equiv c$. Two facts:

* $\sum_{c}E_c=\sum_j(-1)^j\binom pj=(1-1)^p=0$ for $p\ge1$ — so **exactly one** non-zero $E_c$ is impossible.
* $E_c=\frac1q\sum_{m=1}^{q-1}\zeta^{-mc}(1-\zeta^m)^p$ with $\zeta=e^{2\pi i/q}$ and $1-\zeta^m\ne0$ — so they
  cannot **all** vanish.

Therefore at least two classes carry non-zero Euler characteristic, so $z_{\rm slot}$ has support on at least two
degrees: $w_s\ge2$. (Checked exhaustively for $1\le p\le25$, odd $3\le q\le15$: no exceptions, and no
single-degree slot cohomology found in a direct scan of $p\le10$.) With (D9.2),

$$W\ \ge\ r+1 .$$

**No-go.** Concentration requires the window to fit inside one grading period, $W\le q$, hence

$$\boxed{\ \mathrm{rank}(G)\ \le\ q-1\ } \tag{D9.3}$$

So for any **fixed** supercharge degree $q$ the rank is bounded and concentration fails at large $N$. Increasing
$q$ does not help at fixed $p$: once $q>p$ one has $\omega=0$ and $W=rp+1$, the widest possible. Since $\omega\ne0$
requires $q\le p$, keeping concentration as $N\to\infty$ demands

$$p\ \ge\ q\ \ge\ \mathrm{rank}+1 ,$$

i.e. both the flavour count and the interaction order must grow with $N$ — an interaction touching $O(N)$ of the
$N^2$ matrix entries, outside the class of $q$-local "matrix SYK" models.

**Scope.** Proved for the maximal weight of $u(N)$. The factorisation uses that the $N$ diagonal modes $E_{ii}$ are
independent, which fails for $su(N)$ (the Cartan directions $E_{mm}-E_{m+1,m+1}$ straddle slots) — consistent with
the measured $su(4)$ anomaly $13,81,81,13$ against the factorised $27,81,81,27$. Lower-Casimir irreps are not
covered by (D9.1); empirically (note section 4g) their window matches the maximal one at $N=3$, but that is not
derived here.


### D9.4 Solving $z_{\rm slot}$: the complete list of concentrating models (2026-09-23)

Concentration needs $W=r(w_s-1)+1\le q$, i.e. $r\le(q-1)/(w_s-1)$, so the question is how small $w_s$ can be.

**With $\omega\ne0$ (i.e. $q\le p$).** Direct computation of $z_{\rm slot}$ for generic $\omega$ over all
$q\le p\le14$:

* $w_s=2$ occurs **only at $q=3$**, and there exactly for $p\equiv3\pmod4$ ($p=3,7,11,\dots$), with support
  $\{(p-1)/2,\ (p+1)/2\}$.
* For $q\ge5$ the minimum is $w_s=q-1$ (observed at $q=3,5,7,9\Rightarrow w_s^{\min}=2,4,6,8$), so

$$r_{\max}=\Big\lfloor\tfrac{q-1}{w_s-1}\Big\rfloor=\Big\lfloor\tfrac{q-1}{q-2}\Big\rfloor=1\quad(q\ge5),
\qquad r_{\max}=2\quad(q=3).$$

The exact maximum of $r_{\max}$ over every cell with $q\le p\le14$ is **2**, attained at $(p,q)=(3,3)$ and its
relatives $p\equiv3\ (4)$. The necessary Euler condition agrees at much larger scale: scanning $p\le400$ and odd
$q\le31$, exactly two non-zero classes $E_c$ occurs for every odd $p$ when $q=3$, but for $q\ge5$ only at $p=1$.

**Conclusion.** The complete list of concentrating models in this family is $\mathrm{rank}\,G\le2$ — i.e. $u(2)$
and $su(3)$ (plus the trivial $u(1)$, $su(2)$) — with $q=3$ and $p\equiv3\pmod4$. This matches the numerics
exactly: $u(2)$ and $su(3)$ concentrate with identical content $9,18,9$, and nothing else does.

**Scaling $p$ does not help.** It produces infinitely many concentrating models ($p=3,7,11,15,\dots$ at $q=3$) but
never lifts the rank ceiling above 2, because only $q$ can raise the ceiling and $w_s^{\min}$ grows with $q$ at
least as fast as the ceiling would.

**The $\omega=0$ loophole.** When $q>p$, $\omega\in\Lambda^q\mathbb C^p=0$, giving $z_{\rm slot}=(1+t)^p$; at
$p=1$ this is $w_s=2$ for every $q$, hence $r\le q-1$, unbounded. But $\omega=0$ means $Q$ annihilates the
maximal-weight sector outright — every state there is BPS because nothing acts on it, not because of fortuity.
That is the single-matrix/solvable situation (Chen 2025 item 1: $H$ a Casimir function, non-chaotic). So the
counting criterion alone is satisfiable in the free branch; a useful model needs $Q$ to act non-trivially **and**
the window to fit inside one period, and D9.3–D9.4 say those two demands are incompatible beyond rank 2.


### D9.5 Random couplings do not help, and at $p=3$ there are none to draw (2026-09-23)

A natural objection to D9 is that SYK owes its behaviour to disorder, so perhaps random couplings rescue
concentration. They do not, for two reasons.

**(i) At $p=3$, $q=3$ the coupling space is a point.** By (D9.1) the maximal-weight cohomology depends on $C$ only
through $\omega=\mathrm{Alt}(C)\in\Lambda^q\mathbb C^p$, and only through its $GL(p)$ orbit. Since
$\dim\Lambda^3\mathbb C^3=1$, $\omega$ is unique up to scale, and rescaling does not change cohomology. **Chen's
tensor is the only non-degenerate three-flavour cubic supercharge up to $GL(3)$.** Drawing $C$ at random can
therefore only reproduce it ($\mathrm{Alt}(C)\ne0$) or destroy it ($\mathrm{Alt}(C)=0$, whence $\omega=0$,
$z_{\rm slot}=(1+t)^3$ and the window *widens* to $3N+1$). This accounts exactly for the earlier sample of 30
random cyclic tensors: 26 gave Chen's profile and 4 gave width 10 — not 26 models agreeing, but 26 draws of one
model. Likewise $\dim\Lambda^3\mathbb C^4=4$ with a single non-zero $GL(4)$ orbit (every 3-form in 4 variables is
decomposable), so $p=4$ is also rigid.

**(ii) Generic couplings are optimal.** Matrix ranks are lower semicontinuous in $\omega$, so the cohomology
dimensions are upper semicontinuous: specialising $\omega$ can only enlarge $H^j$ and hence widen the support.
Random draws sample the generic orbit, so they *minimise* $w_s$. Verified directly at $q=3$:

| $p$ | generic | $e_{123}$ | $e_{123}+e_{456}$ | $e_{123}+e_{145}$ |
|---|---|---|---|---|
| 6 | **3** | 5 | 3 | 5 |
| 7 | **2** | 6 | 4 | 6 |
| 8 | **3** | 7 | 5 | 7 |

No special orbit beats the generic one. Every width quoted in D9.2--D9.4 is therefore the best case available, and
the no-go (D9.3) applies at the optimum rather than at some arbitrary choice of couplings.

**Contrast with SYK.** There, disorder over $\binom Nq$ independent couplings is what produces chaos and the
$\cos(\pi k/q)$ BPS profile. Here the couplings collapse to a single point of $\Lambda^q\mathbb C^p/GL(p)$ at
small $p$, and to the generic orbit at larger $p$ — and the obstruction sits at the generic orbit.
