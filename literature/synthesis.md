# Literature synthesis: bootstrapping matrix SYK models

*Written 2026-09-15 after reading the 12 papers then in `papers/`; extended the same day with five further papers (Han–Hartnoll–Kruthoff 2020, Chang–Lin 2024, Fu–Gaiotto–Maldacena–Sachdev 2016, Cho–Gabai–Sandor–Yin 2024, Turiaci–Witten 2023) — see §8. Per-paper notes are in `literature/notes/`. Equation numbers refer to the cited paper. Statements marked **[project]** come from the earlier in-repo work in `research/notes/` and have not been independently re-verified in this session (no Python environment with numpy/scipy is currently installed); statements marked **[my estimate]** are derived here and should be checked.*

---

## 1. The three strands of the literature

The twelve papers fall into three groups that Chen 2025 ties together (its refs. [50]–[56] are exactly the other papers in the library, apart from Chang et al. 2024, Lin–Zheng 2024, Lin 2025 and Klebanov–Lin–Meshcheriakov 2026):

| Strand | Papers | Role for us |
|---|---|---|
| **A. Fortuity & R-charge concentration** | Chang–Lin 2024 (origin of the taxonomy); Fu–Gaiotto–Maldacena–Sachdev 2016 ($\mathcal N=2$ SYK, refined index, super-Schwarzian); Turiaci–Witten 2023 ($\mathcal N=2$ random-matrix ensemble); Chang–Chen–Sia–Yang 2024; Chen 2025 | Defines the physics question, the target models, and the quantitative predictions of the super-Schwarzian |
| **B. Matrix / SUSY quantum-mechanical bootstrap** | Han–Hartnoll–Kruthoff 2020 (founding paper); Lin–Zheng 2024; Lin–Zheng 2025; Laliberte–McPeak 2025; Cho–Gabai–Sandor–Yin 2024 (thermal/KMS, adjoint-valued operators in ungauged MQM); Lin 2025 (TASI §6) | Defines the method |
| **C. Purely fermionic matrix (and tensor) quantum mechanics** | Anninos–Denef–Monten 2015; Anninos–Silva 2016; Tierz 2017; Klebanov–Milekhin–Popov–Tarnopolsky 2018; Gaitan et al. 2020; Klebanov–Lin–Meshcheriakov 2026 (adjoint sector of bosonic MQM) | Prior art on the *kind of system* we bootstrap: finite Hilbert spaces, Casimir-type solvability, non-singlet sectors, thermodynamics |

```mermaid
graph TD
  CCSY[Chang-Chen-Sia-Yang 2024<br/>fortuity, R-charge concentration,<br/>supercharge chaos conjecture] --> CHEN[Chen 2025<br/>single-matrix Q=Tr Psi^3: solvable, fortuitous,<br/>NOT concentrated; proposes 3-matrix 'matrix SYK']
  CHEN --> US[(This project:<br/>bootstrap E_0 per R-charge sector<br/>of matrix SYK models)]
  LZ24[Lin-Zheng 2024<br/>BFSS SUSY bootstrap, SO(9) blocks] --> LZ25[Lin-Zheng 2025<br/>high-precision multimatrix bootstrap,<br/>ground-state positivity, relaxation]
  LZ25 --> US
  LM[Laliberte-McPeak 2025<br/>SUSY matrix QM bootstrap<br/>with fermion matrices] --> US
  TASI[Lin 2025 TASI §6<br/>method review] --> US
  ADM[Anninos-Denef-Monten 2015] --> AS[Anninos-Silva 2016] --> TZ[Tierz 2017]
  KMPT[Klebanov et al 2018<br/>fermionic tensor/matrix QM,<br/>Casimir Hamiltonians, rep-resolved bounds] --> GK[Gaitan et al 2020<br/>density of states, Hagedorn]
  AS --> CHEN
  KMPT --> CHEN
  GK --> CHEN
  KLM[Klebanov-Lin-Meshcheriakov 2026<br/>adjoint sector of MQM] -.non-singlet sectors.-> US
  ADM -.-> US
  KMPT -.rep-resolved positivity bounds.-> US
```

---

## 2. Strand A — what "fortuity" and "R-charge concentration" mean precisely

**Fortuity** (Chang et al. §2, Defs. 1–2): a BPS state is *monotonous* if its $Q$- (or $Q^\dagger$-) cohomology class can be uplifted, preserving fermion number, to a nontrivial class in the theory with every larger $N$; otherwise it is *fortuitous*. Operationally: BPS states whose allowed R-charge range depends explicitly on $N$ are fortuitous (Chen §2.2). Fortuitous states become null when $N$ is *decreased* (trace relations / Young diagrams exceeding $N-1$ rows, Chen Fig. 4) and lifted when $N$ is *increased*.

**R-charge concentration** (Chang et al. (1.1), p. 4, Conjecture 1): decompose the Hilbert space into *irreducible cochain complexes* of $Q$ — complexes that cannot be refined further by symmetries commuting with $Q$. For a cubic supercharge these labels include $N_\Psi$ mod 3 (Chang et al. (2.6)–(2.7)); **for the matrix models they also include the $SU(N)$ irrep and any flavor quantum numbers**, since $[J^a,Q]=0$. Concentration means: in every irreducible complex with a macroscopic index, all BPS states sit in a *single* degree. Consequences:

- A window of **three consecutive** $N_\Psi$ values is compatible with strict concentration (one degree per $\mathbb Z_3$ class); this is the situation in generic $\mathcal N=2$ SYK, $\frac{N-3}2\le p\le\frac{N+3}2$ (Chang et al. Fig. 2; Chen p. 8).
- A window of **four or more consecutive** values within one $SU(N)$ irrep is *not* concentrated. Chen's single-matrix model has window $[\tfrac{N(N-1)}2,\tfrac{N(N+1)}2]$, i.e. $N+1$ sectors, all inside $r_*$ (Chen (2.20)–(2.22)); it is concentrated for $N=2$ (trivially) and non-concentrated for all $N\ge3$. The earlier project comparison "window width 3 vs 3 at $N=2$" **[project]** is correct but should be read in this refined sense: the decisive test is whether any $(\mathbb Z_3\text{ class},\ SU(N)\text{ irrep},\ \text{flavor})$ complex acquires BPS states in two degrees as $N$ grows.
- Concentration ⇒ the (refined) index equals the BPS count in each complex with no cancellations (Chang et al. p. 10). Deconcentration shows up as *index no longer saturated* (§5.1). This gives a purely combinatorial diagnostic: compute $\mathrm{Tr}[(-1)^Fq^{N_\Psi}]$ per irrep (a unitary-matrix integral with a character insertion, Chen (3.2)) and compare with the actual BPS count.

**Why concentration matters physically**: it is the "smoking-gun" of the $\mathcal N=2$ super-Schwarzian (Chang et al. §1.1–1.2), the boundary avatar of the non-linear charge constraint of BPS black holes (§5.2), and, per Conjecture 1, a property of *generic* $q$-local supercharges (Turiaci–Witten ensemble near the BPS states). Chen's model fails it because the $SU(N)$ symmetry makes the supercharge extremely sparse and structured: "not all the states that can be lifted are lifted."

**Chen's diagnosis and proposal.** The single-matrix model is solvable because $H=3N(N^2-1)-9\hat C_2$ (2.8): energy is a Casimir, the singlet sector sits at the *top* of the spectrum, and all BPS states form the maximal-Casimir irrep $r_*$ with $C_2=N(N^2-1)/3$, $\dim r_*=3^{N(N-1)/2}$, $Z_{\rm BPS}(q)=3^{N(N-1)/2}(1+q)^Nq^{N(N-1)/2}$. Fortuity is present (staircase diagram with $N-1$ rows), concentration and chaos are absent. Chen conjectures (§4) that within purely fermionic cubic-supercharge matrix models the minimal chaotic, concentrating example needs **three** matrices, e.g. $Q=\sum_{i\le j\le k\le3}\mathrm{Tr}[\Psi_i\Psi_j\Psi_k]$ (4.3), and suggests the matrix bootstrap for it. The only evidence given is that an $N=2$ simulation "suggests R-charge concentration"; no numbers, no chaos test, no larger $N$.

**Sparseness heuristic [my estimate].** Chang et al. (§5.1, (5.8)) parametrise genericity by the fraction $p_s$ of nonzero 3-form components and find deconcentration below $p_s^{\rm crit}\sim N^{-\gamma}$, $\gamma\gtrsim1$ (their Fig. 12 overlays slopes $-1$ and $-3/2$); for $\mathcal N=4$ SYM they estimate $p_s\sim\tilde N^{-3/2}$. For the 3-matrix model with $\tilde N=3N^2$ fermions, $Q$ has $\lesssim10N^3$ monomials out of $\binom{\tilde N}{3}\approx4.5N^6$, so $p_s\sim2/N^3\sim10\,\tilde N^{-3/2}$ — the same $\tilde N^{-3/2}$ scaling as SYM. The single-matrix model has $p_s\sim N^3/(N^6/6)\sim6\tilde N^{-3/2}$ as well, yet does *not* concentrate, so sparseness alone cannot be the criterion; the algebraic structure (Casimir Hamiltonian) is what kills concentration. Whether the 3-matrix model escapes this is exactly the open question.

---

## 3. Strand B — the bootstrap method, and what changes for a purely fermionic model

All four bootstrap papers share one skeleton (Lin 2025 §6.2; Laliberte–McPeak §1; Lin–Zheng 2025 "Bootstrap ingredients"):

1. **Positivity** $M_{ij}=\langle O_i^\dagger O_j\rangle\succeq0$ for a finite operator set (Hilbert-space positivity only; fermions are no obstacle).
2. **Stationarity** $\langle[H,O]\rangle=0$ (any eigenstate or stationary $\rho$); at finite $N$ optionally the stronger, non-convex $\langle OH\rangle=E\langle O\rangle$ which resolves individual eigenvalues ("archipelago", Lin–Zheng 2024 App. B).
3. **Kinematic constraints**: algebra/normal ordering, reality/time reversal, symmetry Ward identities, gauge constraint $\langle\mathrm{tr}\,GO\rangle=0$ when gauged.
4. **Ground-state selection**: ground-state positivity $\langle O_i^\dagger[H,O_j]\rangle\succeq0$ (Lin–Zheng 2025 (6); $T\to0$ limit of the energy–entropy inequality, Lin 2025 (102)–(104)); for a SUSY vacuum instead the supercharge equations $\langle\{Q,O\}\rangle=0$ / $\langle QO\rangle=\langle OQ\rangle=0$ (Lin–Zheng 2024 (14), Lin 2025 (128)), which *imply* ground-state positivity (Lin–Zheng 2024 App. A).
5. **Large $N$** enters only through cyclicity + (anti)commutators producing double traces, replaced by products of single traces (factorisation); the resulting quadratic constraints are handled by scanning or **nonlinear relaxation** (Lin–Zheng 2025 App. B: $\begin{pmatrix}1&x^T\\x&\mathcal Q\end{pmatrix}\succeq0$ and $M\succeq\mathcal Q$).
6. **Symmetry blocking**: decompose operators into irreps of the global symmetry; positivity reduces to $a_{\bar R,R}\succeq0$ per irrep (Lin–Zheng 2024 (16)–(24), Lin–Zheng 2025 App. F), with crossing kernels/6j-symbols relating channels.
7. **Hierarchy**: levels $\ell(X)=1,\ell(P)=2,\ell(\psi)=3/2$; truncate by level; count free variables after solving equalities (Lin–Zheng 2024 Table 1; 2025 Table I).

Practicalities worth remembering: SDPA-GMP/SDPB rather than generic solvers because the primal problem is not strictly feasible (Lin–Zheng 2025 App. A; Laliberte–McPeak used SDPB at small $g$); constraint *generation* (algebra) dominates run time; the ungauged Marinari–Parisi SDP was intractable (Laliberte–McPeak §4); level-to-level convergence can plateau then jump (Laliberte–McPeak §3.2.2); ground-state positivity blocks were too small to give upper bounds at level 8 there but were decisive at level $\ge10$ in Lin–Zheng 2025.

**What is different for Chen's models (structural, not just technical):**

- *Finite-dimensional, purely fermionic Hilbert space.* Pauli exclusion means the operator algebra of words in $\Psi,\bar\Psi$ closes at finite length; there is no continuum of moments. At finite $N$ the bootstrap can in principle be exact (complete at length $2pN^2$), and the eigenstate-resolving constraint $\langle OH\rangle=E\langle O\rangle$ is available. No flat directions, so the "peninsula" pathology of BFSS/bosonic YM is not expected.
- *No $\hbar$, no semiclassical large-$N$ saddle for the fortuitous sector.* Chen's fortuitous states are non-singlet, maximal-Casimir states not captured by a classical saddle (Chen §1, §3.3; cf. Anninos–Denef–Monten's finding that semiclassical bosonic dynamics needs the long-rectangular limit). Large-$N$ factorisation of single-trace correlators — the standard input of strand B — is *not* established in such states. **This is the main conceptual gap between the literature and our problem.**
- *Which sector is the ground state.* In the single-matrix model the singlet sector is the top of the spectrum; a gauged (singlet) bootstrap would bound the *highest* energy. The bootstrap must be run in the ungauged theory, organised by $(N_\Psi,SU(N)\text{ irrep})$ blocks, with the Casimir as a variable — exactly the "multiplet-averaged density matrix" $\rho_R$ that Klebanov et al. 2018 §3.2 use to derive representation-resolved bounds by hand. For the 3-matrix model the earlier project ED **[project]** indicates that gauge-singlet BPS states *do* exist at $N=2$ (36 singlet multiplets among 972 BPS states), so a gauged bootstrap becomes meaningful there — but the standard justification for factorisation still needs an argument.
- *SUSY structure is directly usable.* $E=\langle H\rangle=\|Q|E\rangle\|^2+\|\bar Q|E\rangle\|^2\ge0$ is a sum-of-squares certificate: including $Q,\bar Q$ in the operator basis makes $E\ge0$ manifest; the BPS window is a linear feasibility problem ($\langle\bar QO\rangle=\langle OQ\rangle=0$). The hard direction is proving $E_0(N_\Psi)>0$ in *lifted* sectors, which needs ground-state positivity and/or higher level. Caution: ground-state positivity within a fixed-$N_\Psi$ (or fixed-irrep) sector is only valid for operators $O$ that do **not** leave the sector, since $O|\Omega_{\rm sector}\rangle$ could otherwise land in a sector with lower ground energy.

---

## 4. Strand C — what fermionic matrix QM has taught us so far

- **Casimir Hamiltonians are generic in fermionic matrix QM.** Anninos–Denef–Monten's vector model ($\hat H=-4\hat J^2+3\hat n$), Klebanov et al. 2018's $O(N_1)\times O(N_2)$ (6.3), $SU(N_1)\times SU(N_2)\times U(1)$ (6.15) and $O(N_1)\times O(N_2)\times U(1)$ (6.23) models, Gaitan et al.'s vector models (2.7), and Chen's $H=3N(N^2-1)-9\hat C_2$ all have spectra fixed (wholly or largely) by quadratic Casimirs; energies are integers; singlet states are rare or absent (Klebanov et al. (6.4), (6.13)); in the 't Hooft limit gaps stay $O(1)$ and the spectrum does not become dense ⇒ no quasi-conformal/Schwarzian regime (Klebanov et al. §6.1). Escaping this class is *necessary* for an SYK-like matrix model, and a Casimir shortcut is what the 3-matrix model must lack. The earlier project ED **[project]** reports 168 distinct, irrational energies for $p=3,N=2$, consistent with non-Casimir structure, but chaos could not be diagnosed at $N=2$.
- **Exact thermodynamics via bosonisation.** Anninos–Silva and Tierz show that quartic fermionic matrix models reduce to unitary/Stieltjes–Wigert eigenvalue integrals with $\sin^2$/$\sinh^2$ Vandermondes — the same class of integrals as Chen's $r_*$-multiplicity integral (3.2), (3.9). Character-insertion integrals are the natural tool for exact (refined) indices per irrep in our models.
- **Density of states of fermionic matrix models** (Gaitan et al. §4): $\log\rho(E)\approx N^2\log2-\frac12(E/\lambda N)^2$ from planar moments $\mathrm{tr}H^n$; vector-type models show Hagedorn growth $e^{-|E|/\lambda}$ instead. The moment method transfers directly to $H=\{Q,\bar Q\}$ for the 3-matrix model and provides analytic checks of bulk spectral data.
- **Representation-resolved positivity bounds** (Klebanov et al. 2018 §3): $C_2\pm\frac12A\cdot A=\frac14(A\pm A)^2\ge0$ and Cauchy–Schwarz with a symmetry-averaged $\rho_R$ yield rigorous $E$-vs-Casimir bounds (3.8), (3.21), saturated by exact solutions in the matrix cases — a proof of concept that low-level "bootstrap" reasoning is sharp in Casimir-type fermionic models.
- **Non-singlet sectors at large $N$** (Klebanov–Lin–Meshcheriakov 2026): the adjoint sector of bosonic MQM has a controlled large-$N$ description (Marchesini–Onofri) with finite gap $\Delta_1\approx0.7416$ at criticality and Regge towers, and the two-point bootstrap reproduces it; only representations with a zero weight (box count divisible by $N$) appear. This is the closest thing in the library to a large-$N$ theory of non-singlet sectors, which our single-matrix warm-up needs.

---

## 5. How the papers relate to one another (dependency summary)

- Chang et al. 2024 → Chen 2025: Chen adopts the definitions of fortuity/concentration, the "phantom saddle" language, and the SYK comparison window; Chen's model is a maximally structured (non-generic, sparse) special case of $\mathcal N=2$ SYK. Chang et al. (5.9)–(5.11) already wrote a Lie-algebra fermionic matrix supercharge $Q\propto\mathrm{Tr}\psi\psi\bar\psi$ (R-charge 1); Chen footnote 2 reports that this $\mathrm{Tr}[\Psi^2\bar\Psi]$ model has no fortuitous states, and instead takes $Q=\mathrm{Tr}\Psi^3$ (R-charge 3).
- Lin–Zheng 2024 → Lin–Zheng 2025 → Lin 2025 TASI: a single programme; 2024 introduces the SUSY (supercharge-EOM) bootstrap and $SO(9)$ blocking for BFSS; 2025 adds ground-state positivity, nonlinear relaxation and high level for bosonic models; TASI reviews both and states the general method. Klebanov–Lin–Meshcheriakov 2026 is by an overlapping author and connects the adjoint sector to a two-point bootstrap.
- Laliberte–McPeak 2025 is the independent SUSY-matrix-QM bootstrap (Marinari–Parisi), the only one with fermionic matrices in the operator basis and fermion-number block structure; it and Lin–Zheng 2025 are the two papers Chen recommends.
- Anninos–Denef–Monten → Anninos–Silva → Tierz: successive exact treatments of quartic Grassmann matrix QM (emergent bosonic matrix / Kähler phase space → exact eigenvalue integral → $q$-orthogonal-polynomial solution).
- Klebanov et al. 2018 → Gaitan et al. 2020: the Princeton fermionic tensor/matrix programme (Casimir solvability, singlet counting, rep-resolved bounds → densities of states, Hagedorn). Chen cites [50–54] collectively as "other fermionic matrix models that are non-supersymmetric".

---

## 6. Open questions the literature leaves (and that we can address)

1. **Does the 3-matrix model (Chen (4.3)) R-charge-concentrate for $N\ge3$?** Nothing beyond "suggests" at $N=2$ exists. Sharp version: in each irreducible complex (fixed $N_\Psi$ mod 3, $SU(N)$ irrep, $S_3$/flavor labels) are BPS states in one degree only? Equivalently: is $E_0(N_\Psi,\mathbf r)>0$ outside a width-$\le3$ window?
2. **Is it chaotic?** No test exists. Diagnostics from Chang et al. §3.4: LMRS projected-operator statistics, information/entanglement entropy of BPS states, level statistics of the non-BPS spectrum within symmetry-resolved sectors; from Klebanov et al.: does the spectrum become dense (gaps $\to0$) as $N\to\infty$, unlike Casimir models?
3. **Is "three" the magic number, and which cubic couplings?** Chen's (4.3) is "not particularly special". What about $p=2$ with mixed traces, or $p=3$ with only mixed-flavor traces, or the $O(p)$-symmetric $d_{ijk}$-type couplings? Which minimal structure breaks the Casimir shortcut?
4. **Large-$N$ description of the fortuitous sector.** Fortuitous states are not a classical saddle (Chen); which large-$N$ organising principle (if any) replaces factorisation for non-singlet, maximal-Casimir states? For the gauged 3-matrix model, is the singlet BPS sector amenable to standard planar factorisation?
5. **Genericity vs. structure.** Conjecture 1 (Chang et al.) says generic supercharges concentrate; matrix models are non-generic by construction. Where on the sparse-to-generic axis do cubic matrix supercharges sit, and is there a matrix analogue of $p_s^{\rm crit}$?
6. **Super-Schwarzian**: Chen expects it as the low-energy theory. What is the near-BPS gap scaling in a concentrated sector — $E_{\rm gap}\sim1/N$ or $1/N^2$ (Chang et al. §3.2 quote $\sim1/N$ for SYK with $N$ fermions; with $N^2$ fermions the natural guess is $1/N^2$)? Untested.
7. **Index saturation**: are refined indices per irrep saturated by BPS counts in the 3-matrix model? Computable exactly by unitary-matrix integrals and checkable against ED at $N=2$.
8. **Following $N$ in the matrix integral for a concentrating model** (Chen §3.3, §4): does a single complex saddle vanish upon integrating the last eigenvalue, or is cancellation among many saddles generic?

---

## 7. Implications for the research plan (see `docs/research_questions.md`)

- Use the single-matrix model as the exactly solvable **calibration** of a purely fermionic, $SU(N)$-covariant bootstrap: targets $E_0(N_\Psi)=3N(N^2-1)-9C_2^{\max}(N_\Psi)$, the BPS window (2.22), degeneracies (2.20), and the $C_2^{\max}=N(N^2-1)/3$ plateau. The earlier project work **[project]** reports a tight $N=2$ covariant bootstrap and tight BPS-window results at $N=3$, with lifted sectors at $N=3$ not yet tight — the natural next technical target.
- Formulate concentration in the refined sense (per irreducible complex) and compute **refined indices per $SU(N)$ irrep** as rigorous companions to bootstrap bounds.
- For the 3-matrix model, carry over the machinery, add $S_3$/$O(3)$-flavor blocking (Lin–Zheng App. F style), use $Q,\bar Q$ in the basis for the BPS sectors and ground-state positivity restricted to sector-preserving operators for lifted sectors, and treat any large-$N$ factorisation as a hypothesis to be tested against ED at $N=2$ (and against exact index data).
- Independently derive large-$N$ spectral moments (Gaitan et al. method) for $\{Q,\bar Q\}$ as analytic checks.

---

## 8. Addendum (second batch of papers, 2026-09-15)

### 8.1 What the new papers add to each strand

- **Chang–Lin 2024** is the origin of the monotone/fortuitous taxonomy via a *covering* $\mathcal H_N\simeq\tilde{\mathcal H}/I_N$ and the long exact sequence (2.5): a fortuitous state is one whose lift $\tilde O$ satisfies $\tilde Q\tilde O\neq0$ but $\tilde Q\tilde O\in I_N$ (a relation that holds only at rank $N$). For the matrix models the covering is concrete (formal polynomials in matrix entries; $I_N$ = rank-$N$ relations such as Young diagrams with $\ge N$ rows), which is exactly Chen's Fig. 4 mechanism. Their Conjecture 3 (fortuitous states exponentially outnumber monotone ones at black-hole charges) is what a "black-hole-like" matrix SYK model must satisfy; Chen's single-matrix model has only $2^N$ copies of $r_*$ — polynomially many multiplets.
- **Fu–Gaiotto–Maldacena–Sachdev 2016** is the parent $\mathcal N=2$ SYK model of which both of Chen's models are sparse, non-random special cases. It supplies the refined index $W_r=(1-e^{2\pi ir/\hat q})^N$ (5.5), the statement that the index is saturated (SUSY unbroken), the $\mathcal N=2$ super-Schwarzian thermodynamics (5.39)–(5.42), and the exact ED counts for $\hat q=3$: BPS states at three adjacent charges with multiplicities $3^{N/2-1}(1,2,1)$ for even $N$ (5.7).
- **Turiaci–Witten 2023** gives the random-matrix ensemble behind Conjecture 1 of Chang et al.: independent AZ $(1+2\nu,2)$ ensembles per $(k,k+\hat q)$ multiplet; BPS states only for $|k|<\hat q/2$ with $N_{\rm BPS}(k)\propto\cos(\pi k/\hat q)$ (3.10); non-BPS multiplet gap $E_0(q)=q^2/4\hat q^2$ and edge density $\sinh(2\pi\sqrt{E-E_0})$ (3.11); genericity forbids BPS states at charges differing by $\hat q$. Their SYK check uses the $r$-ratio of the *singular values of $Q_k$* per charge sector (§2.5).
- **Han–Hartnoll–Kruthoff 2020** is the founding bootstrap paper (positivity + $\langle[H,O]\rangle=0$ + gauge + cyclicity/factorisation + reality) and contains the convergence argument (footnote [20]) that becomes exact for our finite fermionic Hilbert spaces.
- **Cho–Gabai–Sandor–Yin 2024** formulates the KMS condition as a matrix inequality (1.6) (whose $\beta\to\infty$ limit is ground-state positivity), makes it SDP-representable, and — most relevant — bootstraps an **ungauged** matrix model with **adjoint-valued (open-index) operators** whose Gram matrix is decomposed into $U(N)$ tensor structures (3.5). This is the missing ingredient identified in the project's covariant bootstrap (addendum to `research/notes/matrix_syk_covariant_bootstrap.md`). They also show the KMS inequality trivialises for gauged models at large $N$ with traced operators only.

### 8.2 A quantitative prediction we can already compare with

Combining FGMS (5.7) and Turiaci–Witten (3.10): a concentrating $\hat q=3$ supercharge governed by the $\mathcal N=2$ super-Schwarzian has BPS states at three consecutive charges $k=-1,0,1$ (relative to the centre) with counts in ratio $\cos(\pi k/3)=(\tfrac12,1,\tfrac12)$, i.e. **1:2:1**. The project's three-matrix ED at $N=2$ gives $243{:}486{:}243$ at $N_\Psi=5,6,7$ — this pattern. The single-matrix model gives binomial weights $\binom{N}{N_\Psi-N(N-1)/2}$ (Chen (2.20)), which equal $1{:}2{:}1$ only at $N=2$ and become $1{:}3{:}3{:}1$ at $N=3$. So the $N=2$ data cannot distinguish the models, but at $N\ge3$ the prediction is sharp: **a concentrating, Schwarzian-governed three-matrix model should keep three charges in ratio $\approx1{:}2{:}1$**; per irreducible complex ($N_\Psi$ mod 3 and $SU(N)$ irrep) these three charges are one per $\mathbb Z_3$ class. (Caveat: (3.10) is the leading large-$e^{S_0}$ answer; FGMS (5.7) shows it is exact in SYK for even $N$; odd fermion number gives the modified pattern in (5.7).)

### 8.3 Consequences for the plan

1. **Chaos diagnostic**: replace $H$-level statistics by the $r$-ratio of singular values of $Q_k$ restricted to a $(k,k+3)$ pair of sectors and a fixed $SU(N)$ irrep (Turiaci–Witten §2.5); expect $\beta=2$ (or $\beta=1$ for the $CT$-self-conjugate multiplet when $pN^2$ is odd). Feasible by sparse SVD in sectors of dimension $10^4$–$10^5$ (e.g. $p=2$, $N=3$; $p=3$, $N=3$ edge sectors).
2. **Refined index per irrep** (FGMS (5.5) generalised with an $SU(N)$ character insertion) as the rigorous companion to bootstrap bounds; check saturation at $N=2$.
3. **Redesign of the covariant bootstrap** with adjoint-valued operators and $U(N)$/$SU(N)$ tensor-structure decomposition (Cho et al. (3.5); Lin–Zheng App. F), plus, later, the KMS inequality if thermal data (entropy, $E(\beta)$ near the BPS window) are wanted.
4. **Targets for lifted sectors**: gaps growing as (distance from window)$^2$ in Schwarzian units (Turiaci–Witten (3.11)); the overall scale with $N$ is unknown for the matrix model and is itself a question.

---

## 9. Addendum (2026-10-06): two-point and thermal bootstraps, and the near-BPS gap

Two papers read in full: Cho–Gabai–Lin–Yeh–Zheng 2025 (`notes/cho_gabai_lin_yeh_zheng_2025.md`) and Adams 2025
(`notes/adams_2025.md`). Both extend Strand B from one-point functions to dynamics. (The four papers added since
§8, Tierz 2026, Kim 2018, Peng–Spradlin–Volovich 2016 and Biggs–Lin–Maldacena 2026, have their own notes but are
not yet folded into this synthesis.)

### 9.1 What they add to Strand B

- **Two-point correlators with rigorous bounds** (CGLYZ). The variables are
  $\mathcal M(\tau)=\langle\bar{\mathcal O}_i(\tau)\mathcal O_j(0)\rangle$. The constraints are reflection
  positivity, Heisenberg equations ($\partial_\tau\mathcal M=-\mathcal MD$), and either ground-state positivity or
  KMS. The dual turns the equations of motion into "inequalities of motion", and finite spline or polynomial
  ansätze for the multipliers give rigorous bounds at continuous $\tau$.
- **Non-singlet gaps from the global ground state** (CGLYZ §4.4). Adjoint operators acting on the singlet ground
  state of the *ungauged* model give
  $\Delta_{\rm adj}\le\max\{\Delta:\tilde{\mathcal N}\succeq\Delta\tilde{\mathcal M}\}$, with $\mathcal M,\mathcal N$
  the zero-time two-point data. This is a rigorous upper bound, converging to $10^{-14}$ by level 16 in the
  one-matrix model. Excited levels and matrix elements come from extremal functionals and fits, which are
  estimates.
- **Log-convexity** of connected correlators gives $\Delta\le-G'(0)/G(0)$, and a new derivation of the
  energy–entropy balance inequality.
- **Exact KMS** (Adams). The relative-entropy cone is solved directly by QICS, which reaches $L=12$ where MOSEK
  with the log relaxation becomes unstable. The adjoint gap $\Delta_1$ is extracted to about $10^{-5}$ relative by
  fitting thermal bounds to the long-string low-$T$ expansion.

### 9.2 Consequences for this project

1. **Our near-BPS gap is the supersymmetric analogue of the adjoint gap.**
   - The global ground states are the BPS states $B$: $E=0$, $QB=\bar QB=0$.
   - For charged operators carrying $B$ into a sector $k$ next to the window,
     $E_0(k,\lambda)\le\max\{\Delta:\mathcal N-\Delta\mathcal M\succeq0\}$, with
     $\mathcal M=\langle B|O^\dagger O|B\rangle$ and
     $\mathcal N=\langle B|O^\dagger HO|B\rangle=\langle B|[Q,O\}^\dagger[Q,O\}|B\rangle+\langle B|[\bar Q,O\}^\dagger[\bar Q,O\}|B\rangle$.
   - The near-BPS gap is therefore the smallest relative size of the supersymmetry variation of an operator acting
     on a BPS state.
   - This supplies the **upper bounds** the project has lacked. The sector bootstrap supplies lower bounds, which
     currently vanish within five sectors of the window.
2. **Ground-state positivity with charged and non-singlet operators** is valid when the reference is a BPS state.
   It is not valid for sector ground states, which is why it added nothing to the sector bootstrap.
3. **The density route.** $\langle B|O^\dagger e^{-\tau H}O|B\rangle$ decomposes over near-BPS states, so the
   two-point bootstrap is the natural tool for the near-BPS density seen by simple probes, which is the
   super-Schwarzian prediction proper. The thermal route (Adams) is the alternative. With $N_\Psi$-sectors it needs
   the grand-canonical ensemble for KMS to hold for charged operators.
4. **Tooling.** SDPB / polynomial matrix programs, MOSEK, QICS, and arbitrary precision at high level.

The plan built on these points is `docs/nearbps_bootstrap_plan.md`.

---

## 10. Addendum (2026-10-07): coloured models, quivers, and a refined fortuity

Four papers read: Witten 2016, Bena–Berkooz–de Boer–El-Showk–Van den Bleeken 2012, Chang–Colin-Ellerin–Rangamani
2018, and Choi–Choi–Kim 2026 (notes in `literature/notes/`). They were prompted by the move to the $U(n)^3$
fermionic quiver (`research/notes/quiver_project.md`).

### 10.1 Where the quiver sits

- **It is the $D=2$ coloured model.** Witten's disorder-free SYK uses $q=D+1$ fields of rank $D$, each pair sharing
  one index.
  - At $q=3$, $D=2$ the fields are bifundamentals of a three-node quiver and the vertex is $\mathrm{Tr}(ABC)$.
  - The melonic proof needs $D\ge3$ (Witten eq. 3.10), so at $D=2$ the limit is planar.
  - The quiver is the supersymmetric, complex, flavoured version of that matrix-level member.
- **The purely fermionic cubic $\mathcal N=2$ case was an acknowledged gap.**
  - Witten (p. 9): the cubic supersymmetric SYK is not obviously disorder-free.
  - Chang–Colin-Ellerin–Rangamani (eq. 6.4): Fermi-superfield $\mathcal N=2$ SYK has odd $q$, and "it is
    therefore unclear how to promote this to a melonic tensor model". Their bosonic supertensor models break
    supersymmetry at large $N$.
  - Biggs–Lin–Maldacena 2026 fill the gap melonically. The quiver fills it at the planar level, with no
    supersymmetry breaking from dynamical bosons, since there are none.
- **Gauging matters at $D=2$.**
  - Witten's "gauging is harmless because $\dim G\ll N$" holds for $D\ge3$. For the quiver, $\dim G=3n^2$ against
    $3pn^2$ fermion modes, so the singlet sector is a genuine restriction.
  - This is consistent with §4x of the cohomology notes: concentration is vacuous at $p=1$ and appears at $p\ge2$.

### 10.2 Black-hole interpretation and fortuity

- **The bosonic precedent.** Bena et al. compute the Higgs-branch BPS spectrum of the bosonic, abelian three-node
  quiver with cyclic cubic superpotential.
  - There are exponentially many zero-angular-momentum "pure-Higgs" states, with closed-form generating function
    (3.27) and growth (3.43).
  - They exist iff the Coulomb branch has a scaling (AdS$_2$) point, and are read as candidate single-centre
    black-hole microstates.
  - Their own suggestion of a fermionic combinatorial origin (p. 8) motivates a cheap $n=1$ comparison with our
    model. That test is proposed; a naive "singlet Euler characteristic = $\Omega$" match already fails.
- **A refined fortuity.** Choi–Choi–Kim extend the Chang–Lin covering to bifundamental theories with baryons, using
  formal inverse letters. $N$-dependent baryon-dressed classes become *generalised monotones*, and the rest are
  *strongly fortuitous*.
  - For our $U(n)^3$ singlets the ordinary multi-trace covering is complete, since every $\epsilon$ pairs within a
    node.
  - But determinant-type singlets have $n$-dependent expansions. Our 90/90 fortuity verdict at $(n,p)=(2,2)$ is
    therefore in the ordinary sense, and the refined criterion could reclassify some classes (not analysed).

### 10.3 Consequences

1. Framing for a quiver paper: "the supersymmetric $D=2$ coloured model", a purely fermionic, cubic, $\mathcal N=2$
   planar model whose singlet sector concentrates. Cite Witten, Chang–Colin-Ellerin–Rangamani and
   Biggs–Lin–Maldacena for the gap it fills, and Bena et al. and Choi–Choi–Kim for the black-hole and fortuity
   context.
2. Two open checks follow:
   - the $n=1$, $(a,b,c)$ comparison with the pure-Higgs counts;
   - the baryonic (generalised-monotone) refinement of the fortuity test.
