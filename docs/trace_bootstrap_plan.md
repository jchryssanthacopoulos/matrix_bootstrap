# Working paper: an exact-finite-$N$ trace bootstrap for the matrix SYK models

*Status: plan, 2026-09-18. Nothing here is implemented yet. Supersedes the "large-$N$ engine" outlook of `research/notes/matrix_syk_largeN_engine.md` §3–4 in the light of the sector-bootstrap results (`research/notes/sector_bootstrap_results.md` §8, derivations D3–D4).*

## 0. One-paragraph summary

The Hilbert-space sector bootstrap is now exact for the three-matrix model at $N=2$ in every sector tested ($k=2$ at level 3, $k=3$ at level 4), and we know *which* constraints do the work: the charge-sector restriction, open-index (adjoint-valued) words, the equations of motion, and the singlet/adjoint split of the open-word Gram matrix. Its cost, however, is set by the block dimension $d_k$ and it cannot go past $N=2$. The next engine re-expresses exactly the same constraints as identities among expectation values of **trace words**, with $N$ entering only as a numerical coefficient (no factorisation, no $1/N$ expansion, no computer algebra). Its cost is set by the operator level, not by $N$, so the same code runs at $N=2$ (validation against everything we have), $N=3$, and $N=10^2$–$10^3$. Realistic reach: sectors with $k\lesssim4$–$5$ particles at any $N$ — the fixed-$k$, large-$N$ regime, where $k=0,1$ are known analytically (D2) and $k\ge2$ is new — not the half-filling BPS window. The window question at $N=3$ is a sparse-Lanczos problem and is planned separately (§8).

## 1. Goal, deliverables, non-goals

**Goal.** Rigorous lower bounds $E_0(k;N)$ on the ground energy of the $N_\Psi=k$ sector of Chen's three-matrix model (and of the single-matrix model as a control), from a constraint set that is exact at every finite $N$, evaluated for $N=2,3,4,\dots$ up to $N\sim10^3$.

**Deliverables.**
1. `src/trace_algebra.py`: the fermionic trace-word algebra (canonical forms, reordering identities with double-trace terms, commutators with $H$, finite-$N$ trace relations), with every identity numerically verified against explicit operators at $N=2$ (and $N=3$, $p=1$).
2. `src/trace_bootstrap.py`: constraint assembly with $N$-polynomial coefficients, evaluation at numerical $N$, reduction to free variables, cone assembly, and a **direct** solver interface (SCS/Clarabel arrays or SDPA/SDPB file export; no cvxpy in the hot path).
3. `scripts/run_trace_bound.py --N --p --k --L_adj --L_sing --L_eom --finiteN_len --solver`, results as JSON lines in `results/data/`.
4. Physics output: $E_0(k;N)$ for $k\le4$ and a range of $N$; comparison with D2 anchors, with ED where available, and with the large-$N$ scaling.

**Non-goals (explicitly).** The half-filling window at $N\ge3$ (level $\sim k$ would be needed; §7); the planar/'t Hooft limit with factorisation (it is blind to concentration, which is $1/N$-subleading); symbolic linear algebra in $N$.

## 2. Formulation

### 2.1 Objects

Letters $\ell\in\{\Psi^{(m)},\bar\Psi^{(m)}\}$, $m=0,1,2$, in the $\mathbb Z_3$ Fourier flavor basis (charge $z$, D4), all Grassmann-odd; $\{\Psi^a_{ij},\bar\Psi^b_{kl}\}=\delta^{ab}\delta_{il}\delta_{jk}$, all other anticommutators zero.

- **Open word** $w=\ell_1\cdots\ell_L$: the operator-valued $N\times N$ matrix $w_{ij}=(\ell_1)_{ii_1}(\ell_2)_{i_1i_2}\cdots(\ell_L)_{i_{L-1}j}$, operator order = matrix order. Charges $q(w)=\#\Psi-\#\bar\Psi$, $z(w)$ mod 3, parity $|w|=L$ mod 2.
- **Trace word** $T=\mathrm{Tr}\,w$; it is a fermionic operator when $L$ is odd.
- **Monomial** $T_1T_2\cdots T_m$ (ordered product of trace words); **variable** $x_{\mathbf T}=\phi(T_1\cdots T_m)$.
- **Functional** $\phi(O)=\mathrm{Tr}[\rho O]$ for a density matrix $\rho$ on the full Fock space, later restricted to sector $k$ (2.3). Only $\phi(\mathbf 1)=1$, reality $\phi(O^\dagger)=\overline{\phi(O)}$, and Hilbert-space positivity are assumed.

Because $\rho$ may be taken $U(N)\times\mathbb Z_3\times U(1)_{N_\Psi}$-invariant (group averaging, D4.1), $\phi$ vanishes on every monomial with $q\ne0$ or $z\ne0$, and $\phi$ of an invariant operator is unchanged — so the variables are the **neutral** monomials only, and the Cho et al. adjoint-projected cone is legitimate without any Ward-identity rows.

### 2.2 The identities (all exact at finite $N$)

Two facts fix what is and is not an identity among trace words; both were learned the hard way in the Hilbert-space version (results note §5):

(i) **Cyclic moves preserve matrix routing; adjacent swaps do not.** Moving the first letter of a trace to the end past the other $L-1$ (all fermionic) letters gives
$$\mathrm{Tr}[\ell_1\ell_2\cdots\ell_L]=(-1)^{L-1}\,\mathrm{Tr}[\ell_2\cdots\ell_L\ell_1]+\sum_{m=2}^{L}(-1)^{m}\,\{\ell_1,\ell_m\}\text{-term},\qquad
\{\ell_1,\ell_m\}\text{-term}=\delta_{\ell_1\bar\ell_m}\ \mathrm{Tr}[\ell_2\cdots\ell_{m-1}]\ \mathrm{Tr}[\ell_{m+1}\cdots\ell_L],\tag{P.1}$$
where $\delta_{\ell_1\bar\ell_m}=1$ iff $\ell_m$ is the conjugate letter of $\ell_1$ (same flavor, opposite type), an empty trace is $N$, and the two factors keep their operator order. The index contraction $\delta_{il}\delta_{jk}$ splits the trace into two traces *with matrix routing preserved* — this is why cyclicity closes on trace words (it is the fermionic version of Han–Hartnoll–Kruthoff eq. 14). Swapping two adjacent letters inside a trace does **not** produce a trace word (the routing becomes $A_{mi}\bar\Psi_{jk}\Psi_{ij}B_{km}$), so **no normal ordering is ever performed**; words with all letter orders are kept, exactly as in the Hilbert-space version. **(P.1) with the sign $(-1)^m$ as written is verified numerically** (2026-09-18: 160 random words of length 2–5 at $N=2$, $p=3$, Fourier letters, residual $<10^{-8}$ in absolute norm).

(ii) **Trace products reorder with merging terms.** $T_AT_B=(-1)^{|A||B|}T_BT_A+\sum(\text{single traces})$: each anticommutator between a letter of $A$ and its conjugate in $B$ joins the two traces into one trace of length $|A|+|B|-2$ with routing $A_{i+1}\cdots A_{|A|}A_1\cdots A_{i-1}\,B_{j+1}\cdots B_{|B|}B_1\cdots B_{j-1}$. *Caution (learned from a first numerical attempt, 19/60 sign failures):* the merged letters are then in matrix order only up to cyclic rotations of the $A$- and $B$-blocks, and restoring matrix order costs further (P.1)-type moves with their own double-trace terms; so (ii) is not a one-line formula but the composition of elementary letter moves. The implementation will therefore have a single primitive — "move one letter one step, tracking operator order and routing" — from which cyclic rotation, trace reordering and commutators are all built, and the harness (§3) certifies each composite. This makes a canonical ordering of the factors in a monomial possible.

From (i)–(ii): a **canonical form** exists — cyclic representative for each trace (chosen by a fixed total order on letter strings), factors of a monomial sorted — and every non-canonical monomial rewrites, deterministically, into canonical monomials of no greater total length (terminating by induction on length). Different rewriting paths can give formally different results that agree as operators; the difference is a finite-$N$ relation and may be added as a constraint or ignored (still rigorous).

(iii) **Equations of motion.** $\phi([H,X])=0$ for every neutral canonical monomial $X$ up to level $L_{\rm eom}$, with $H$ from D2.11/D2.13 ($H$ is single-trace plus the length-one double traces $\mathrm{Tr}\Psi^c\mathrm{Tr}\bar\Psi^d$, both handled by (i)–(ii)); the commutator is computed letter by letter with the anticommutator, producing trace words with routing preserved.

(iv) **Sector restriction.** $\rho$ supported on $N_\Psi=k$ $\iff$ $\rho(N_\Psi-k)=0$ $\iff$ $\phi\big((N_\Psi-k)\,X\big)=0$ for all $X$ (and the adjoint statement). With $N_\Psi=\sum_a\mathrm{Tr}[\Psi^a\bar\Psi^a]$ this is a family of linear rows in double-trace variables. It is exact, needs no projector, and reproduces D3.1's $P_kOP_k$ logic at the trace level. ($\phi(N_\Psi)=k$, $\phi(N_\Psi^2)=k^2$ are the first two instances.)

(v) **Positivity cones**, both as in D3/D4 with trace-word entries:
$$\text{singlet:}\ \ \phi\big(T_a^\dagger T_b\big)\succeq0\ \text{over canonical trace words (and products) of equal }(q,z);\qquad
\text{adjoint-projected:}\ \ \phi\big(\mathrm{Tr}[w_a^\dagger w_b]\big)-\tfrac1N\,\phi\big((\mathrm{Tr}\,w_a)^\dagger\mathrm{Tr}\,w_b\big)\succeq0\ \text{over open words of equal }(q,z).$$
$Q$ and $\bar Q$ are length-3 traces and enter the singlet cone automatically at $L_{\rm sing}\ge3$ (giving the SUSY floor). Ground-state positivity is available in the same form but was never needed at $N=2$ (D4.4), so it is optional.

(vi) **Finite-$N$ trace relations (optional, needed at small $N$).** The antisymmetriser over $N+1$ matrix indices vanishes identically, $\delta^{[i_1}_{j_1}\cdots\delta^{i_{N+1}]}_{j_{N+1}}=0$; contracting it with any product of $N+1$ letter-entries *in a fixed operator order* gives, for each ordering, a linear relation among monomials of total length $N+1$ (and, with longer matrix blocks in place of letters, of any length $\ge N+1$). No reordering of operators is required, so the statistics play no role in the derivation. At $N=2$ these start at length 3 (e.g. relations among $\mathrm{Tr}(abc)$, $\mathrm{Tr}(ab)\mathrm{Tr}(c)$, $\mathrm{Tr}a\mathrm{Tr}b\mathrm{Tr}c$) and are what distinguishes $2\times2$ matrices from generic ones; at $N=100$ they are irrelevant at any reachable level. Parameter `finiteN_len` = maximum total length of relations included.

### 2.3 What is *not* imposed

No factorisation. No planar counting. No Cayley–Hamilton beyond the explicit family (vi). No assumption on the state beyond invariance (which is free). Every bound is therefore a rigorous lower bound on $E_0(k;N)$ at the chosen integer $N$.

### 2.4 Variables, coefficients, normalisation

- Coefficients are integers times powers of $N$ (from empty traces and from index contractions in $\mathrm{Tr}[w^\dagger w']$ expansions). They are stored as short integer arrays (polynomials in $N$ of degree $\le$ level); products are convolutions. The constraint set is built **once** and evaluated at any $N$ in $O(1)$ per coefficient.
- Variables are generated by closure ("touched" monomials), never pre-enumerated; the set is capped by total level only. A monomial with more traces is never truncated away — it becomes a variable (rigorous; possibly weaker).
- For $N\gg1$ variables are rescaled by $N^{-(\text{number of traces})}$ (and words by $N^{-L/2}$ if needed) to keep coefficients $O(1)$; this is the 't Hooft normalisation and affects conditioning only.

### 2.5 Reduction and SDP

At each numerical $N$: assemble the sparse linear system of rows (P.1)-rewriting residues if any, (iii), (iv), reality; find an independent set and a parametrisation of the free variables by sparse QR/SVD (as `_independent_rows` does now; exact fraction-free elimination as a fallback if conditioning at small $N$ demands it); express every cone entry in the free variables; hand the SDP to the solver directly. The solver stage is the one place where the present pipeline (cvxpy) is known to waste memory by $3\times$ (results note §8.3); this engine will write the cone data straight into SCS/Clarabel structures or an SDPA-format file for SDPA-GMP/SDPB (Lin–Zheng 2025 App. A recommend arbitrary precision because these SDPs are not strictly feasible).

## 3. Verification harness (the non-negotiable part)

At $N=2$ ($p=3$, $2^{12}$ states) and at $N=2,3$ ($p=1$) every trace word and monomial can be evaluated as an explicit sparse operator with the existing `fermion_matrix_model.word_matrix`. Therefore:

1. **Every identity the algebra emits is checked numerically** as an operator equation (rewriting (P.1)–(ii), commutators (iii), sector rows (iv), finite-$N$ relations (vi)) on random vectors, to $10^{-12}$, in the unit tests. Signs and $N$-powers are fixed by this, not by hand.
2. **Every SDP row is checked against the exact multiplet-averaged ground state**: $\phi_{\rm GS}$ evaluated on all monomials must satisfy all equalities and all cones (this is the `exact_coordinates`/`check_feasibility` diagnostic that caught the round-off bug in D4.3, transplanted).
3. **The bound must be $\le$ the exact value and $\le$ the Hilbert-space bound** at the same level, and must approach the Hilbert-space bound as `finiteN_len` grows.

## 4. Validation ladder and physics targets

| step | model | check | expected |
|---|---|---|---|
| V1 | $p=1$, $N=2,3,4$ | $E_0(k)$ all sectors (exact: Casimir ladder; $N=3$: $72,45,18,0,\dots$) | tight with finite-$N$ relations at low level; without them, rigorous but looser at $N=2$ |
| V2 | $p=3$, $N=2$ | $k=0,1$: $98,22$; $k=2$: $5.16536$ at level 3; $k=3$: $1.22706$ at level 4 | reproduce with `finiteN_len` large enough; quantify the loss without |
| V3 | $p=3$, $N=3$, $k=0,1$ | D2 anchors $387$, $273$ | exact at level $\le3$ |
| V4 | $p=3$, $N=3$, $k=2,3,4$ | sector ED ($d_k=351,\,2925,\,17550$ — trivial) | tightness vs level at $N=3$; first test of "level $\approx k+1$" beyond $N=2$ |
| V5 | $p=3$, $N=4,5,\dots,10^3$, $k\le4$ | $k=0,1$ analytic for all $N$; $k\ge2$ new | large-$N$ form of $E_0(k;N)$ at fixed $k$; is $E_0(k;N)-E_0(0;N)$ a polynomial in $N$ as for $k=0,1$? |

**Where the engine says something ED cannot.** For fixed small $k$, sector ED is cheap at any moderate $N$ ($d_k=\binom{3N^2}{k}$), so the *numbers* in V4 and much of V5 are checkable. The engine's specific value is (a) rigour and $N$-uniformity — one constraint set, a bound at every $N$, and a certificate; (b) the regime $k\ge3$, $N\gtrsim6$ where $d_k\gtrsim10^7$; (c) being the only route to a bootstrap statement at $N$ beyond ED at all. This should be said plainly in any write-up.

## 5. Cost estimates (to be replaced by measurements)

Six letters. Open words at length $L$: $6^L$; graded by $(q,z)$ into $\sim15$–$30$ cones: largest cone $\approx28$ ($L=3$), $170$ ($L=4$), $840$ ($L=5$), $\sim4500$ ($L=6$). Neutral single-trace canonical words of length $2\ell$: $\binom{2\ell}{\ell}3^{2\ell}/(2\ell\cdot3)$ — about $40$ ($\ell=2$), $800$ ($\ell=3$), $2\times10^4$ ($\ell=4$), $5\times10^5$ ($\ell=5$); the Gram entries $\mathrm{Tr}[w^\dagger w']$ have length $2L_{\rm adj}$, so the variable count is dominated by them: $\sim10^4$ at level 4, $\sim10^5$–$10^6$ at level 5 before reduction. (Lin–Zheng 2025 handle $2\times10^5$ constraints → $1.6\times10^3$ free variables at their level 14 with four letters; our level 4–5 is comparable in size.) **Level 4 is the realistic first production level, level 5 the stretch, level 6 the ceiling of any word-based approach with six letters.** Memory: with a direct solver interface, level 4 should fit in a few GB; level 5 will need SDPB-class tooling. All runs continue to be planned with an estimate first and a watchdog (project rule).

## 6. Milestones

- **M1 — algebra + harness** (`trace_algebra.py`, tests): canonical forms, (P.1), (ii), commutators, sector rows, finite-$N$ relations; all verified against explicit operators at $N=2,3$. *Acceptance:* zero identity failures on $10^3$ random words up to length 6. **Done 2026-09-18** (see todo for the summary). Design note: (P.1) and (ii) were not implemented as formulas at all — the only primitive is the adjacent swap of a *routed term*; both identities, and the finite-$N$ relations, fall out of bubble-sorting operator order into routing order. The self-symmetry reductions ($\mathrm{Tr}[\ell\ell]=0$ etc.) were added after the harness exposed $\mathrm{Tr}[\ell\ell]$ surviving as a dead variable.
- **M2 — assembly + reduction + direct solver** (`trace_bootstrap.py`): V1 at $N=2,3$ reproduces the Casimir ladder; $\phi_{\rm GS}$ feasibility check passes on every instance. *Acceptance:* V1 tight at $N=3$, level $\le3$.
- **M3 — three-matrix validation** (V2, V3): $k=2$ at level 3 and $k=3$ at level 4 reproduced at $N=2$; $387/273$ at $N=3$. *Acceptance:* agreement to $10^{-4}$ with `finiteN_len` reported; the gap without finite-$N$ relations recorded.
- **M4 — first new numbers** (V4, V5 up to $N\approx10$): $E_0(k;N)$, $k=2,3,4$, level 4; ED cross-checks at $N=3,4$. *Acceptance:* tightness table by $(k,N,\text{level})$; test of the level-$\approx k+1$ hypothesis.
- **M5 — large $N$** ($N=10^2,10^3$, 't Hooft normalisation): fixed-$k$ scaling; compare $k=0,1$ with D2 exactly; fit $k\ge2$.
- **M6 — write-up** in `research/notes/` + updates to `docs/derivations.md` (D5: the fermionic trace algebra and its identities, with the harness-fixed signs) and `docs/research_questions.md`.

## 7. Risks and what we do about them

1. **Level needed grows with $k$ (and possibly with $N$).** Evidence: $k=2\to$ level 3, $k=3\to$ level 4 at $N=2$; single-matrix $N=4$ edge sectors not tight at level 5. If confirmed in V4, the reach is $k\le4$–$5$, full stop; this is a result about fortuitous/near-BPS states resisting low-level truncation and is worth recording as such.
2. **Variable blow-up at level 5** ($10^5$–$10^6$ monomials before reduction). Mitigations: exploit $(q,z)$ neutrality, canonical forms, and the fact that only monomials touched by cones/EOM are variables; sparse elimination; if still too large, restrict length-5 open words to the $q=\pm1$ cones (the ones that carried the information at level 4).
3. **Finite-$N$ relations at $N=2$ may be numerous** (all lengths up to 8). If V2 cannot be reproduced exactly, that is acceptable *provided* the bound is rigorous and the deficit is understood; the $N=2$ Hilbert-space numbers remain the reference. At $N\ge3$ the relations start at length $\ge4$ and matter less.
4. **Conditioning.** Coefficients up to $N^{L}$; 't Hooft rescaling plus the diagonal cone congruence already in use; arbitrary precision (SDPA-GMP) as fallback, as Lin–Zheng advise.
5. **Sign errors in the fermionic rewriting.** Mitigated entirely by the harness (§3): no identity is used before it is verified numerically.

## 8. Parallel track (separate plan): concentration at $N=3$ by cohomology ranks, and the division of labour

**What the bootstrap does and does not buy for concentration.** Concentration (CCSY Conjecture 1) is a statement about the degree in which $Q$-cohomology sits in each complex $(k\bmod3,\,SU(N)\text{ irrep},\,\text{flavor})$. The direct tools are (a) exact cohomology ranks per sector, $\#\mathrm{BPS}(k)=\dim\ker Q_k-\operatorname{rank}Q_{k-3}$, i.e. sparse rank computations over a prime field on the $k$-particle Fock sectors (cost $\sim d_k$, exact integers) — how the fortuity literature decides monotone vs fortuitous at finite $N$; and (b) refined indices (Chen §3, Anninos–Silva §4, FGMS 5.7), polynomial cost at any $N$, with **index saturation $\Leftrightarrow$ concentration** (CCSY §5.1) — but saturation needs (a) or an exclusion of BPS states in the other degrees. A bootstrap bound $E_0(k,R)>0$ is such an exclusion, but the degrees to exclude are the window-adjacent, many-particle, small-gap sectors where the bootstrap is weakest (§7). Hence: **concentration $\Rightarrow$ cohomology ranks + indices; near-BPS energies and their $N$-scaling (Q5) $\Rightarrow$ bootstrap.** Cohomology and indices say nothing about energies; that is the bootstrap's genuine, non-redundant product.

**$N=3$ track.** Exact BPS counts per $(k,\text{irrep},z)$ for $n=27$ ($d_k\le2\times10^7$): matrix-free $Q_k$ in the Jordan–Wigner representation restricted by $N_\Psi$, $\mathbb Z_3$ and the $U(3)$ Cartan weights; ranks by sparse elimination mod a large prime (or by Lanczos on $Q_k^\dagger Q_k$ counting zero modes, which also yields the gaps and the Turiaci–Witten singular-value statistics for Q2). Output: window width at $N=3$ (3 = SYK-like vs 4 = single-matrix-like), per-irrep counts vs the refined index (Q1b), the $1{:}2{:}1$ test (Q1d). Memory estimate first ($2\times10^7$-dimensional vectors, $0.3$ GB each); this does not depend on the trace engine and should run alongside M1–M3.

## 9. Housekeeping decisions

- The Hilbert-space sector bootstrap (`sector_bootstrap.py`, D3–D4) is frozen as the $N=2$ reference; remaining item there is $k=4$ at level 4–5 after the direct-solver interface exists (it will be shared with the trace engine).
- The "$U(3)$ highest weights" and "never form $H_k$ densely" items are **retired**: they belonged to the Hilbert-space route past $N=2$, which is closed by design (§0).
- Dependencies: numpy/scipy for the algebra and reduction; the SDP backend to be decided at M2 (SCS/Clarabel via their Python APIs are already installed through cvxpy; SDPA/SDPB would be added to `requirements.txt`/README if used).
