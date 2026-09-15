# Lin, Zheng (2024) — Bootstrapping Ground State Correlators in Matrix Theory, Part I

**File:** `papers/lin_zheng_2024.pdf`

## Full citation

Henry W. Lin, Zechuan Zheng, *Bootstrapping Ground State Correlators in Matrix Theory, Part I*, arXiv:2410.14647v3 [hep-th] (v3: 27 Aug 2025).

## Main question

Can the bootstrap, working directly at infinite $N$ and zero energy, bound simple correlators of the BFSS (D0-brane) matrix model — a strongly coupled system with 9 bosonic and 16 fermionic matrices — with modest computational resources, by imposing the *supercharge* equations of motion and exploiting $SO(9)$?

## Physical system

BFSS: 9 traceless Hermitian bosonic matrices $X_I$ and 16 fermionic $\psi_\alpha$ with $[X^I_{ij},P^J_{kl}]=i\delta_{il}\delta_{jk}\delta^{IJ}$, $\{\psi_\alpha,\psi_\beta\}=\delta_{\alpha\beta}\delta_{il}\delta_{jk}$ (2);
$$H=\mathrm{Tr}\big(\tfrac{g^2}2P_I^2-\tfrac1{2g^2}[X_I,X_J]^2-\psi_\alpha\gamma^I_{\alpha\beta}[X_I,\psi_\beta]\big)\quad(3)$$
16 supercharges $Q_\alpha$ (4) with $\{Q_\alpha,Q_\beta\}=2\delta_{\alpha\beta}H+2\gamma^I_{\alpha\beta}\mathrm{Tr}X^IC$ (5), $C=-i[X^I,P^I]-\psi_\alpha\psi_\alpha-N\mathbf 1$ (6). Unique normalisable SUSY ground state assumed (refs [26–30]).

## Important definitions

- Variables: $SO(9)$-singlet single-trace words with invariant tensors from $\delta^{IJ},\epsilon^{I_1\cdots I_9},\gamma^I_{\alpha\beta}$ (8)–(10); normalisation (7).
- **Kinematic constraints** (§2.2): cyclicity with double-trace remainder, e.g. (12); Hermiticity/time reversal ($\langle\mathrm{tr}O_1\cdots O_n\rangle=\pm\langle\mathrm{tr}O_n\cdots O_1\rangle$); gauge invariance $\langle\mathrm{tr}\,OC\rangle=0$ (13).
- **Dynamical constraint**: $\langle\{Q_\alpha,O_\alpha\}\rangle=0$ for any single-trace $SO(9)$ spinor $O_\alpha$ (14) — the infinite-$N$ form of $Q_\alpha|\Omega\rangle=0$. Footnote 3: for bosonic $O$ one would use $\langle[Q,O_B]\rangle=0$, but $SO(9)$ already enforces it.
- **Positivity by $SO(9)$ blocks**: $M_{ij}=\langle\mathrm{tr}\bar O_iO_j\rangle\succeq0$ (15); decomposition $O_i=\sum_R\sum_r(C_R)^r_i(O_R)_r$ (16); $\langle\mathrm{tr}\bar O_{\bar R\bar r}(O_R)_r\rangle=\delta\,a_{\bar R,R}$ (19); reduces to $a_{\bar R,R}\succeq0$ per irrep (24).
- **Hierarchy**: $\ell(X)=1,\ell(P)=2,\ell(\psi)=3/2$; $\ell(\{Q_\alpha,O_\alpha\})=\ell(O_\alpha)+\tfrac12$ (25); cyclicity's double-trace term has level $\ell-3$; gauge constraint gives level $\ell+3$; positivity up to level $\ell_{\rm cut}$ needs operators up to $\ell_{\rm cut}/2$. Table 1: free variables 3, 4, 11, 18, 59, 149 at levels 4–9 (from 11, 38, 140, 569, 2528, 12077).
- **Crossing kernel** $F_{R_s,R_t}$ (32)–(35): a $6j$-symbol of $SO(9)$; for four spinors a $5\times5$ rational matrix with $F^2=1$ (Fierz identities). Analogy to conformal bootstrap (footnote 12).

## Main assumptions

- Large-$N$ factorisation (46) for double traces — argued to be subtle for BFSS because of flat directions / power-law tails (moments $\langle\mathrm{tr}X^\ell\rangle$ diverge for $2\ell\ge9$ at finite $N$; §4).
- Ground state is $SO(9)$-invariant, SUSY, gauge singlet.

## Main analytical results

- 4-fermion worked example (§2.6): $16\otimes16=1+9+36+84+126$ (26); solving cyclicity (29)–(31) with the crossing kernel eliminates 3 of 5 variables; positivity (36)–(39) gives $1\le\langle\mathrm{tr}OO\rangle\le2$ (40), improved to $\le1.53125$ with 6 fermions.
- Level-6 analytic bounds (43)–(44) and lower bound $\langle\mathrm{tr}X^2\rangle\ge\frac34(3/50)^{1/3}\approx0.2936$ (45).
- **Appendix A**: ground-state positivity $\langle O^\dagger[H,O]\rangle\ge0$ (47) is **redundant** once the supercharge constraints and inner-product positivity of $Q$-exact operators are imposed: $\langle O^\dagger_{ji}[H,O_{ij}]\rangle=\frac18\langle([Q_\alpha,O])^\dagger([Q_\alpha,O])\rangle\ge0$ (51), similarly (54). Holds for any SUSY QM with at least one supercharge annihilating the state.
- **Appendix B** (toy supermembrane $H=\frac12(p_x^2+p_y^2)+\frac12gx^2y^2$): with $\langle OH\rangle=E\langle O\rangle$ the archipelago resolves individual eigenstates (Fig. 4); with only $\langle[H,O]\rangle=0$ plus ground-state positivity (60) the region shrinks from a peninsula to an island by level 8 (Fig. 5). The equation $\langle OH\rangle=E\langle O\rangle$ becomes trivial under factorisation (58).

## Important equations

(2)–(6), (11)–(16), (19)–(25), (28)–(40), (43)–(47), (51), (54), (57)–(60).

## Numerical methods

Symbolic solution of kinematic/dynamical constraints to a set of free variables (level 9: ~1 s), then SDP with a scan over $\langle\mathrm{tr}X^2\rangle$ (constraints linear in the other variables at fixed $\langle\mathrm{tr}X^2\rangle$ up to level 9; at level 11 quadratic terms appear and nonlinear relaxation or more scanning is needed). Gamma algebra with the GAMMA package, LieART. Everything ran on a laptop within an hour.

## Relevant figures/results

Table 2: lower bounds on $\langle\mathrm{tr}X^2\rangle$: 0.260 (L5), 0.294 (L6), 0.329 (L7), 0.340 (L8+), 0.355 (L9) vs Monte Carlo $\approx0.378\pm0.04$ / $[0.346,0.430]$. Fig. 1: 3D "shoe" allowed region at level 9, tip $(\langle\mathrm{tr}X^2\rangle,\langle K\rangle,\langle\mathrm{tr}OO\rangle)\approx(0.355,0.495,1.02)$. Fig. 2: peninsulas levels 4–9. Fig. 3: kink in $\langle\mathrm{tr}X^2X^2\rangle$ bound near the tip.

## Limitations

- No upper bound on $\langle\mathrm{tr}X^2\rangle$ through level 9 (peninsula, not island) — attributed to flat directions.
- Large-$N$ factorisation may fail for sufficiently complicated operators in BFSS (§4).
- Group theory (crossing kernels for higher-point $SO(9)$ blocks) is the computational bottleneck.

## Relationship to our project

Three transferable ideas. (1) **SUSY ground-state bootstrap**: for a BPS state one can impose $\langle\{Q,O\}\rangle=0$ / $\langle QO\rangle=\langle OQ\rangle=0$ directly, which is stronger than $\langle[H,O]\rangle=0$ and (App. A) automatically implies ground-state positivity. For our BPS sectors ($E_0(N_\Psi)=0$) this is the natural formulation; for *lifted* sectors it does not apply and ground-state positivity must be imposed separately. (2) **Symmetry blocking with crossing kernels**: their $SO(9)$ machinery is what an $SU(N)$-covariant (Wigner–Eckart) treatment of our non-singlet fortuitous sector would look like; the "conformal-bootstrap analogy" (irrep ↔ dimension, $a_{R}$ ↔ OPE$^2$) is the right mental model. (3) The toy-model appendix shows how $\langle OH\rangle=E\langle O\rangle$ resolves *individual* eigenvalues (archipelago) — relevant because our finite fermionic Hilbert spaces make an eigenstate-resolved bootstrap conceivable at finite $N$.
