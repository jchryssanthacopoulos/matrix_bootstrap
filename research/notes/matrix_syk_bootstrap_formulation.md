# Bootstrapping matrix SYK models — I. Formulation and the solvable warm-up $Q=\mathrm{Tr}[\Psi^3]$

*Working note. Goal: set up the quantum‑mechanical matrix bootstrap for the purely‑fermionic "matrix SYK" models proposed in Chen, "Fortuity with a Single Matrix" (arXiv:2511.00790), §4. This first note fixes the formulation and validates it on Chen's exactly‑solvable single‑matrix model $Q=\mathrm{Tr}[\Psi^3]$, whose entire spectrum we already know. Numerical implementation follows in a later note.*

---

## 0. Why this model, and why the bootstrap

Chen's discussion (§4, "Simplest chaotic matrix model with $R$‑charge concentration") conjectures that the minimal *chaotic* supersymmetric matrix model with the black‑hole‑like feature of $R$‑charge concentration is a **purely fermionic matrix model with a cubic supercharge** — a "matrix SYK model." The explicit example is

$$
Q \;=\; \sum_{1\le i\le j\le k\le 3}\mathrm{Tr}\!\left[\Psi_i\Psi_j\Psi_k\right],
\qquad H=\{Q,\bar Q\},
$$

with three $U(N)$‑adjoint complex fermion matrices $\Psi_i$. At $N=2$ (12 complex fermions) it already shows $R$‑charge concentration. Unlike the single‑matrix model it is **not** exactly solvable, so Chen suggests the quantum‑mechanical matrix bootstrap (Lin–Zheng arXiv:2507.21007; Laliberte–McPeak arXiv:2510.01356) as the tool of choice.

The plan of this project:

1. **Warm‑up / validation (this note).** Bootstrap Chen's *single*‑matrix model $Q=\mathrm{Tr}[\Psi^3]$, which is exactly solvable ($H=3N(N^2-1)-9\hat C_2$). Because we know the answer, this is the controlled setting to build and check a *purely fermionic* matrix bootstrap — a setting neither cited paper treats directly (both carry a bosonic matrix).
2. **Target.** Carry the same machinery over to the multi‑matrix cubic model above, where $H$ is no longer a Casimir and the bounds become genuine.

A one‑line summary of what makes the fermionic case special: **the Hilbert space is finite‑dimensional and every fermionic mode appears at most once (Pauli).** Consequently the operator algebra closes at finite size and the bootstrap is, in principle, *exact at finite level* — there is no continuum of moments to converge through as in bosonic matrix QM.

---

## 1. The single‑matrix model

### 1.1 Fields and algebra

A single $U(N)$‑adjoint complex fermion matrix $\Psi=(\Psi_{ij})$, $i,j=1,\dots,N$, i.e. $N^2$ complex fermions, with

$$
\{\Psi_{ij},\bar\Psi_{kl}\}=\delta_{il}\delta_{jk},\qquad
\{\Psi_{ij},\Psi_{kl}\}=\{\bar\Psi_{ij},\bar\Psi_{kl}\}=0,\qquad
\bar\Psi_{ij}\equiv\Psi^\dagger_{ji}.
$$

$\Psi_{ij}$ are creation operators, $\bar\Psi_{ij}$ annihilation operators, and $\bar\Psi_{ij}|0\rangle=0$.

**Supercharge and Hamiltonian** ($\mathcal N=2$):

$$
Q=\mathrm{Tr}[\Psi^3]=\Psi_{ij}\Psi_{jk}\Psi_{ki},\qquad
\bar Q=Q^\dagger=\mathrm{Tr}[\bar\Psi^3],\qquad
H=\{Q,\bar Q\}.
$$

$Q$ is nilpotent, $Q^2=0$ (verified below), so $E=\langle H\rangle=\|Q|E\rangle\|^2+\|\bar Q|E\rangle\|^2\ge 0$, with $E=0$ iff $Q|E\rangle=\bar Q|E\rangle=0$ (BPS).

**Adjoint basis.** Writing $\Psi=\tfrac{\psi^0}{\sqrt N}\mathbf 1+\sum_{a=1}^{N^2-1}\sqrt 2\,\psi^a T^a$ with $\mathrm{Tr}[T^aT^b]=\tfrac12\delta^{ab}$, $[T^a,T^b]=if^{abc}T^c$, one has $\{\psi^a,\bar\psi^b\}=\delta^{ab}$ and

$$
Q=\frac{i}{\sqrt2}\,f^{abc}\,\psi^a\psi^b\psi^c .
$$

The trace mode $\psi^0$ **decouples**; its only effect is to double every degeneracy. We work in the $\mathfrak{su}(N)$ sector ($N^2-1$ fermions) and restore the factor of two at the end.

### 1.2 Symmetries (the good quantum numbers for the bootstrap)

- **$R$‑charge / fermion number.** $N_\Psi=\mathrm{Tr}[\Psi\bar\Psi]$ with $[N_\Psi,Q]=3Q$. Energy eigenstates have definite $N_\Psi$. It ranges over $0\le N_\Psi\le N^2$.
- **$SU(N)$ global symmetry.** Generators $J^a=\mathrm{Tr}[\Psi[T^a,\bar\Psi]]=-if^{abc}\psi^b\bar\psi^c$, with $[J^a,J^b]=if^{abc}J^c$. Quadratic Casimir $\hat C_2=\sum_a J^aJ^a$.
- **SUSY.** $Q^2=0$, $\{Q,\bar Q\}=H$, $[Q,N_\Psi]=-3Q$, $[Q,J^a]=0$.

The single fact that makes this model solvable:

$$
\boxed{\,H=3N(N^2-1)-9\,\hat C_2\,}
$$

so the spectrum is fixed entirely by the $SU(N)$ representation content — the model is **non‑chaotic**. In the fermionic Hilbert space the representations are exactly those appearing in $\wedge^p\mathbf{adj}$, $p=N_\Psi$.

### 1.3 Exact solution = the bootstrap's validation targets

Numerically constructing $\psi^a$ (Jordan–Wigner), $Q$, $H$, $\hat C_2$ from the generalized Gell‑Mann $T^a$ and the structure constants $f^{abc}$, we confirm (this is our ground truth):

| Check | $N=2$ (3 fermions, $\dim=8$) | $N=3$ (8 fermions, $\dim=256$) |
|---|---|---|
| $Q^2=0$ | ✔ | ✔ |
| $H=\{Q,\bar Q\}$ Hermitian | ✔ | ✔ |
| $H=3N(N^2-1)-9\hat C_2$ (operator identity) | ✔ | ✔ |
| Spectrum $E$ (degeneracy, $\mathfrak{su}(N)$ sector) | $0\,(6),\ 18\,(2)$ | $0\,(108),\ 18\,(80),\ 45\,(64),\ 72\,(4)$ |
| BPS Casimir $C_2$ at $E=0$ | $2=N(N^2-1)/3$ | $8=N(N^2-1)/3$ |
| $\dim\mathbf r_*=3^{N(N-1)/2}$ | $3$ | $27$ |

The energies coincide with $E=3N(N^2-1)-9\,C_2(\mathbf r)$ evaluated on the irreps of Table 1 of the paper (e.g. $N=3$: singlet $C_2=0\Rightarrow72$, $\mathbf{8}\Rightarrow45$, $\mathbf{10}\oplus\overline{\mathbf{10}}\Rightarrow18$, $\mathbf{27}\Rightarrow0$).

Two remarks the bootstrap must respect:

- **Ground state = maximal representation $\mathbf r_*$.** BPS states saturate the *largest possible* Casimir, $C_2^{\max}=N(N^2-1)/3$, realized by the staircase Young diagram with $N-1$ rows and $2N-2i$ boxes in row $i$. Minimizing energy is *maximizing the Casimir*.
- **Degeneracy bookkeeping.** In the $\mathfrak{su}(N)$ sector the BPS count is $3^{N(N-1)/2}\cdot 2^{N-1}$ (e.g. $6=3\cdot2$, $108=27\cdot4$): $\dim\mathbf r_*=3^{N(N-1)/2}$ states per copy, and $2^{N-1}$ copies distributed over $\tfrac{N(N-1)}2\le N_\Psi\le\tfrac{N(N+1)}2$ from $Z_{\mathbf r_*}(q)=(1+q)^N q^{N(N-1)/2}$. The trace mode restores the full factor $2^N$.

---

## 2. The bootstrap: general formulation

The strategy is the standard eigenstate bootstrap (Han–Hartnoll–Kruthoff; Lin; and the matrix versions in the two cited papers), adapted to a purely fermionic, finite‑dimensional matrix model.

### 2.1 Objects and the master positivity condition

All expectation values $\langle\,\cdot\,\rangle\equiv\langle E|\,\cdot\,|E\rangle$ are taken in a single energy eigenstate. For any finite set of operators $\{O_I\}$, the Gram (moment) matrix

$$
M_{IJ}=\langle O_I^\dagger O_J\rangle \succeq 0
$$

must be positive semidefinite, because $M$ is the Gram matrix of the states $\{O_I|E\rangle\}$. Being in an eigenstate gives the **Heisenberg / equation‑of‑motion** constraints

$$
\langle[H,O]\rangle=0 \qquad\text{for all }O,
$$

and the stronger linear relation $\langle OH\rangle=E\langle O\rangle$. These relate the many entries of $M$ so that positivity has bite. Adding operators only tightens the bounds (principal‑minor monotonicity). Extremizing $E=\langle H\rangle$ subject to these constraints yields rigorous two‑sided bounds on the ground‑state energy and on any expectation value.

### 2.2 Operator basis: fermion words and why they truncate

Take the $O_I$ to be **words** in the elementary fermions,

$$
O_I \;=\; \Psi_{i_1 j_1}\cdots \bar\Psi_{k_1 l_1}\cdots ,
$$

organized (for the ungauged/covariant bootstrap) into definite $SU(N)$ tensors, or (for the gauged bootstrap) into gauge‑invariant traces $\mathrm{Tr}[\,\cdots]$. Three features specific to the fermionic matrix model:

1. **Pauli truncation / normal ordering.** Using $\{\Psi_{ij},\bar\Psi_{kl}\}=\delta_{il}\delta_{jk}$ and $\Psi_{ij}^2=\bar\Psi_{ij}^2=0$, every word reduces to a canonical normal‑ordered form, and each of the $N^2$ modes appears at most once. Hence the set of *independent* correlators is **finite**: the algebra closes exactly. A "level‑$L$" truncation (words with $\le L$ letters) reaches the full algebra at $L=2N^2$.
2. **$R$‑charge block‑diagonality.** $N_\Psi$ is conserved, so $\langle O_I^\dagger O_J\rangle=0$ unless $O_I$ and $O_J$ carry equal $N_\Psi$. $M$ splits into blocks labeled by charge.
3. **Fermion‑number/Hermiticity phases.** Reality of the eigenstate gives $\langle O_1\cdots O_n\rangle^*=\langle O_n^\dagger\cdots O_1^\dagger\rangle$, removing further variables (as in both cited papers).

### 2.3 Gauged vs. ungauged — a decision forced by fortuity

This is the central structural choice and differs from the cited (gauged, singlet‑sector) bootstraps.

- **Gauged / singlet sector.** Keep only $SU(N)$‑invariant traces; impose the gauge (Gauss‑law) constraint $\langle G\,O\rangle=0$ with $G$ the traceless generator. *But* in this model the singlets have $C_2=0$, hence all sit at the **maximum** energy $E=3N(N^2-1)$ — there are no BPS states in the singlet sector. Bootstrapping the gauged model therefore probes the *top* of the spectrum, not fortuity.
- **Ungauged, targeting $\mathbf r_*$.** Fortuity lives in the nontrivial representation $\mathbf r_*$. We must keep non‑singlet operators. Two equivalent ways to handle the ground‑state degeneracy:
  - **(a) Fixed reference state.** Evaluate in the highest‑weight state $|\lambda\rangle=\prod_{i>j}\Psi_{ij}|0\rangle$ (Eq. 2.15 of the paper). Concrete and computable; ideal for a first consistency check of the constraints against an explicit state.
  - **(b) Multiplet‑averaged (covariant) bootstrap.** Use the $SU(N)$‑invariant object $\langle O\rangle=\tfrac{1}{\dim}\mathrm{Tr}_{E}[O]$ (normalized trace over the degenerate eigenspace). Only singlet operators acquire nonzero value; Wigner–Eckart then expresses all correlators through a handful of reduced matrix elements. This is the natural state‑agnostic variable and is what generalizes cleanly to the SYK target.

For validation we will use **(b)** as the primary formulation and **(a)** as an independent cross‑check.

### 2.4 The full constraint set

Collecting everything, the SDP imposes:

1. **Positivity.** $M_{IJ}=\langle O_I^\dagger O_J\rangle\succeq0$ on each ($R$‑charge, irrep) block.
2. **Reality/Hermiticity.** $\langle O_1\cdots O_n\rangle^*=\langle O_n^\dagger\cdots O_1^\dagger\rangle$.
3. **Algebra.** Anticommutators + Pauli nilpotency → canonical form (removes most variables).
4. **Heisenberg / EOM.** $\langle[H,O]\rangle=0$ and $\langle OH\rangle=E\langle O\rangle$. Equivalently, since $H=3N(N^2-1)-9\hat C_2$ here, $\langle[\hat C_2,O]\rangle=0$ and $\langle O\hat C_2\rangle=\tfrac{1}{9}\big(3N(N^2-1)-E\big)\langle O\rangle$.
5. **Global Ward identities.** $\langle[N_\Psi,O]\rangle=0$ (charge), and $\langle[J^a,O]\rangle$ relations fixing $SU(N)$ covariance / selecting the target irrep. In the gauged formulation these become the Gauss‑law constraints $\langle G\,O\rangle=0$.
6. **SUSY constraints.** Nilpotency $\langle Q^2\,\cdots\rangle=0$; the algebra $\langle\{Q,\bar Q\}O\rangle=\langle HO\rangle$; and, for the BPS target, the annihilation conditions $Q|E{=}0\rangle=\bar Q|E{=}0\rangle=0$, which turn into linear constraints $\langle \bar Q\,O\rangle=\langle O\,Q\rangle=0$ and make the $E=0$ sector a *linear feasibility* problem.
7. **Ground‑state positivity (thermal bootstrap at $T=0$).** $\langle O_I^\dagger[H,O_J]\rangle\succeq0$, which supplies **upper** bounds on $E$ (lower bounds on $\hat C_2$). This is the ingredient that gave two‑sided control in both cited works.

### 2.5 What the bootstrap must reproduce (and the Casimir shortcut)

Because $H$ is an affine function of $\hat C_2$, the whole exercise reduces to bounding a Lie‑algebra Casimir realized on $N^2$ fermions:

- **Minimizing $E$ ⇔ maximizing $\langle\hat C_2\rangle$.** The lower bound on the ground‑state energy is $E\ge 3N(N^2-1)-9\,\langle\hat C_2\rangle_{\max}$. The bootstrap should return $\langle\hat C_2\rangle_{\max}=N(N^2-1)/3$, i.e. $E_{\min}=0$.
- **How positivity sees $C_2^{\max}$.** Build a moment matrix out of the $J^a$ (and their products) and impose the $\mathfrak{su}(N)$ algebra $[J^a,J^b]=if^{abc}J^c$ together with the relation that the $J^a$ are the specific fermion bilinears $-if^{abc}\psi^b\bar\psi^c$ (which encodes the $\wedge^p\mathbf{adj}$ realization). Positivity then caps the achievable Casimir at exactly the largest irrep hosted by the fermions. Reproducing $C_2^{\max}$ is the sharp, quantitative validation.
- **Spectrum on a line.** Every eigenvalue must satisfy $E=3N(N^2-1)-9C_2(\mathbf r)$ for an allowed $\mathbf r$; the bootstrap, scanned by fixing $E$ and testing feasibility, should find solutions only at the discrete $E$ values in the table above.

This Casimir shortcut is a feature of the *solvable* model; it is exactly what will be **absent** in the matrix SYK target, and is why the validation is meaningful — the machinery is exercised on a case with a known nontrivial answer.

### 2.6 SDP structure

- **Blocks.** Positivity is imposed block‑by‑block in $(N_\Psi,\ SU(N)\text{ irrep})$. Different irreps are orthogonal (Wigner–Eckart), exactly as Lin–Zheng decompose by $O(D)$ irrep and by $U(1)$ charge.
- **Finiteness.** Fermionic nilpotency bounds each block; there is a level $L=2N^2$ at which the bootstrap is complete and the "bounds" become the exact answer. For validation we track convergence in $L$ well below this, mirroring the level‑7/8 studies of Laliberte–McPeak.
- **Solver.** Small blocks → Mathematica/`SDPB`‑style solves; the algebra (normal ordering + constraint generation) is the main cost, as both cited papers emphasize.

---

## 3. Validation protocol for $Q=\mathrm{Tr}[\Psi^3]$

Concrete deliverables for the implementation note, each with a known target:

1. **Lower bound on $E$ → 0.** Extremize $\langle H\rangle$ with constraints 1–5; confirm the bound rises to $0$ as $L$ increases and never exceeds it.
2. **Upper bound on $E$ → 0.** Add ground‑state positivity (7); confirm two‑sided squeeze to $0$.
3. **$C_2^{\max}=N(N^2-1)/3$.** Independently bound $\langle\hat C_2\rangle$ from the $J^a$‑block; match $2$ (N=2), $8$ (N=3).
4. **Excited levels.** Scan feasibility at fixed $E$; recover the discrete ladder $\{0,18\}$ (N=2), $\{0,18,45,72\}$ (N=3).
5. **BPS linear system.** Impose $\langle\bar Q\,O\rangle=\langle O\,Q\rangle=0$; verify the $E=0$ moment data is consistent with the explicit $|\lambda\rangle$ correlators (cross‑check of formulations (a) vs (b)).
6. **Small‑$N$ sanity by ED.** All of the above are checkable against exact diagonalization ($\dim=2^{N^2}$: 16 for $N=2$, 512 for $N=3$), which we have already used to fix the targets in §1.3.

Passing 1–5 with these exact numbers certifies a working *purely fermionic* matrix bootstrap.

---

## 4. Generalization to the matrix SYK target

With the machinery validated, the multi‑matrix cubic model is a near drop‑in replacement.

### 4.1 The models

$p$ adjoint fermion matrices $\Psi_i$ ($i=1,\dots,p$), each obeying (1.1), with a cubic supercharge. Two natural families:

- **Symmetric cubic (Chen's example, $p=3$):** $\displaystyle Q=\sum_{1\le i\le j\le k\le 3}\mathrm{Tr}[\Psi_i\Psi_j\Psi_k]$, $H=\{Q,\bar Q\}$. At $N=2$: $pN^2=12$ complex fermions, $\dim=2^{12}=4096$ — **exactly diagonalizable**, so $R$‑charge concentration is directly checkable by ED (as Chen reports) and serves as a second validation layer for the bootstrap.
- **Two decoupled matrices (singlet fortuity):** $Q=\mathrm{Tr}[\Psi_1^3]+\mathrm{Tr}[\Psi_2^3]$ — gives fortuitous states in the *singlet* sector via $\mathbf 1\subset\mathbf r_{*,\Psi_1}\otimes\mathbf r_{*,\Psi_2}$, but no $R$‑charge concentration. A useful intermediate: it is solvable (product of two copies) yet has singlet BPS states, so the *gauged* bootstrap becomes nontrivial here.

### 4.2 What changes, what carries over

- **Carries over verbatim:** moment‑matrix positivity, reality, Pauli/normal‑ordering, $R$‑charge blocking, gauge/Ward constraints, SUSY nilpotency and $\{Q,\bar Q\}=H$, and ground‑state positivity. The $O(p)$ (flavor) global symmetry gives an *additional* block decomposition, directly analogous to Lin–Zheng's $O(D)$ organization that let them reach high level.
- **What is lost:** the identity $H=3N(N^2-1)-9\hat C_2$. $H$ is no longer a Casimir; the spectrum is chaotic and the ground‑state energy is a genuine unknown. The Casimir shortcut of §2.5 becomes an inequality the bootstrap must *discover* rather than an identity — this is where the method earns its keep.
- **Large $N$.** Purely fermionic ⇒ finite‑dimensional, so there is no semiclassical $\hbar\to0$; the large‑$N$ limit is the $N^2\to\infty$ fermion limit. Whether large‑$N$ factorization (single‑trace dominance) holds as cleanly as in the bosonic models of Lin–Zheng needs to be checked; at minimum, $N=2,3$ are ED‑accessible anchors, and the bootstrap targets the trend in $N$.

### 4.3 Physics questions the bootstrap can address

- **$R$‑charge concentration.** Bound the distribution of BPS states over $N_\Psi$; test whether it concentrates near $N_\Psi\approx pN^2/2$ (SYK‑like) as opposed to the spread‑out single‑matrix result.
- **Ground‑state energy and BPS degeneracy vs. $N$.** Two‑sided bounds on $E_0(N)$ and on the $E=0$ multiplicity; probe the sharp $N\to N-1$ sensitivity ("phantom" mechanism) that defines fortuity.
- **Chaos diagnostics** compatible with the bootstrap (e.g. bounds on simple correlators / level‑statistics proxies), and the conjectured super‑Schwarzian low‑energy behavior.

---

## 5. Immediate next steps

1. Fix conventions and enumerate the level‑$L$ fermion‑word basis for $N=2$, $Q=\mathrm{Tr}[\Psi^3]$, in both formulations (a) and (b); count independent variables after constraints 1–5.
2. Implement normal‑ordering + constraint generation (the expensive algebra step) and the $(N_\Psi,\text{irrep})$ block SDP.
3. Run the §3 validation protocol; confirm the exact targets in §1.3.
4. Swap in $p=3$, add the $O(p)$ blocks, and produce first bounds for the $N=2$ SYK model, cross‑checked against $2^{12}$ ED.

---

### Sources (project library)

- Y. Chen, *Fortuity with a Single Matrix*, arXiv:2511.00790 — model, exact solution, fortuity, §4 matrix‑SYK proposal.
- C.‑M. Chang, Y. Chen, B. S. Sia, Z. Yang, *Fortuity in SYK Models*, arXiv:2412.06902 — fortuity, $R$‑charge concentration background.
- S. Laliberte, B. McPeak, *Bootstrapping supersymmetric (matrix) quantum mechanics*, arXiv:2510.01356 — SUSY matrix bootstrap constraints, ground‑state positivity, 44×44 SDP.
- H. W. Lin, Z. Zheng, *High‑Precision Bootstrap of Multimatrix Quantum Mechanics*, arXiv:2507.21007 — irrep decomposition, ground‑state positivity, high‑level techniques.
