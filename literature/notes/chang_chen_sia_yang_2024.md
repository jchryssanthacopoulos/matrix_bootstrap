# Chang, Chen, Sia, Yang (2024) — Fortuity in SYK Models

**File:** `papers/chang_chen_sia_yang_2024.pdf`

## Full citation

Chi-Ming Chang, Yiming Chen, Bik Soon Sia, Zhenbin Yang, *Fortuity in SYK Models*, JHEP 08 (2025) 003, arXiv:2412.06902v3 [hep-th] (v3: 23 Nov 2025).

## Main question

How does the fortuity/monotony classification of BPS states (originally from $\mathcal N=4$ SYM and D1–D5) play out in supersymmetric SYK models, and what universal features does a *generic* $q$-local supercharge have near its fortuitous states? The paper proposes **R-charge concentration** as the smoking-gun signature and formulates the **supercharge chaos conjecture**.

## Physical system

- $\mathcal N=2$ SUSY SYK: $N$ complex fermions, $Q=\sum C_{i_1\cdots i_q}\psi_{i_1}\cdots\psi_{i_q}$ (1.2), $\{Q,Q^\dagger\}=H$, $Q^2=0$ (1.3), R-charge $R=N_\psi=\sum\psi_i\bar\psi_i$ (1.4) with $[R,Q]=\hat qQ$, $(-1)^F=e^{\pi iR}$ (1.5). $\mathcal H^p\cong\Lambda^p(\mathbb C^N)$, $D_p=\binom Np$ (1.7).
- Two-flavor model of Heydeman–Turiaci–Zhao: $Q=\sum C_{ijk}\psi_i\psi_j\bar\chi_k$ (3.1); conserved $J=N_\psi+2N_\chi$, $R=-N_\chi$ (3.5).
- Modified two-flavor models (3.18), (3.26) with symmetrised couplings that admit monotonous states.
- Comments on $\mathcal N=4$ SYM and on a Lie-algebra model $Q=\sum f_{ijk}\psi^i\psi^j\bar\psi_k$ (5.9)–(5.11).

## Important definitions

- **Cochain complex / grading**: $N_\psi=nq+q_f$ (2.6); $q_f\in\mathbb Z_q$ commutes with $Q$, so each fixed-$q_f$ sector is a complex $\mathcal H^{\bullet}_{N,q_f}$ (2.7); $Q$-cohomology (2.8) ↔ BPS states (Hodge).
- **Irreducible cochain complex**: one that cannot be further refined by flavor symmetries commuting with $Q$ (p. 9).
- **R-charge concentration** (1.1) and p. 4: along an irreducible cochain complex with non-zero (macroscopic) index, *all* BPS states sit in a **single** space $\mathcal H^{n_c}$.
- **Projection between $N$ and $N+1$** (2.9)–(2.12), short exact sequence (2.11), long exact sequence (2.15) with connecting map $\delta_N$ (2.13).
- **Monotonous cohomology** (Def. 1): a class in $H^\bullet(\mathcal H_N)$ that can be iteratively uplifted to a nontrivial class for all $N'>N$ via (2.16), (2.18). The uplift preserves fermion number.
- **Monotonous / fortuitous state** (Def. 2): BPS state is monotonous if it corresponds to a monotonous $Q$- or $Q^\dagger$-class; otherwise fortuitous.
- **Supercharge chaos (Conjecture 1, p. 9)**: for a generic $q$-local complex supercharge on an $L$-dimensional Hilbert space with an R-symmetry, in an irreducible cochain complex with macroscopic index $\sim L^\nu$, BPS states are sharply concentrated at a single R-charge, and near the BPS states $Q$ is approximated by the Turiaci–Witten $\mathcal N=2$ random matrix ensemble (Altland–Zirnbauer $(\alpha,\beta)=(1+2n_{\rm BPS},2)$ in adjacent spaces). Corollary: index = number of BPS states (no cancellations).
- **Following $N$** (§4.1): continuous interpolation $C_{ijk}(N)$, e.g. (4.1) (decouple one fermion) or the analytic (4.2)–(4.4) with weights $w_a(N)$.
- **LMRS criterion** (3.48): random-matrix statistics of a simple operator projected into the BPS subspace; **information entropy** (3.52); **entanglement entropy**.

## Main assumptions

- Genericity of couplings (measure-zero exclusion, (1.8)–(1.9)).
- No ensemble average unless stated.
- For the large-$N$ two-flavor analysis: no spontaneous breaking of $U(K)$ / $U(1)_\chi$ (§3.3, App. A).
- Conjecture 1 excludes complexes with vanishing index (footnote 5: e.g. odd $N$, $q=3$, where $p=(N\pm3)/2$ can each hold $O(1)$ BPS states).

## Main analytical results

1. **All BPS states in generic $\mathcal N=2$ SYK are fortuitous**: no cohomology for $p<\frac{N-q}2$ (rank argument, (1.8)–(1.9)); by charge conjugation none for $p>\frac{N+q}2$; hence BPS charges lie in $\frac{N-q}2\le p\le\frac{N+q}2$ (Fig. 2). Since $q_f=N_\psi$ mod $q$ splits this into distinct complexes, each complex has BPS states in one degree ⇒ concentration.
2. **Index determines the concentrated charge**: (1.11)–(1.12) $r_c=\mathrm{Im}[ir_*+\frac1\pi S(r_*)]$; for SYK with refined index $I(\mu)$, $\mu=2\pi k/q$, saddle $r_*=1+1/(e^{i\mu}-1)$ gives $r_c=1/2$ (1.13)–(1.14). Footnote 11: exact refined indices for $q=3$, even $N$: $I_0=2\cdot3^{N/2-1}\cos\frac{N\pi}6$, $I_1=-2\cdot3^{N/2-1}\sin\frac{(N+1)\pi}6$, $I_2=2\cdot3^{N/2-1}\sin\frac{(N-1)\pi}6$.
3. **Genericity argument** (1.10): among 3-local supercharges the $q=3$ SYK form is the most generic; $Q'=\sum f_{ij}^k\psi^i\psi^j\bar\psi_k+\sum c_i\psi^i$ has a much smaller solution space of $Q'^2=0$.
4. **Two-flavor model**: non-linear charge constraint $f(n_\psi,n_\chi)=0$ (3.11) derived from the index saddle (3.13)–(3.14) or from fugacity saddle (3.15)–(3.16); ED at $N=9$ confirms concentration along fixed-$J$ complexes (Fig. 3b: $\{0,0,3231,0,0\}$ at $J=10$). BPS states are *not* at the maximal-dimension space in general (Fig. 4).
5. **Monotonous states** from symmetrised couplings: $[Q,V]=0$, $V=\sum\psi_i\chi_i$ (3.19); $V,V^\dagger,N_V$ form $su(2)$ (3.22); $|V^k\rangle$ are BPS and monotonous (3.24)–(3.25); $K$-flavor version gives $sp(2K)$ (App. D) and $D_{\rm mono}\lesssim\binom{n+2}2$ (3.46) vs. $e^{O(N)}$ fortuitous (3.29). Large-$N$ $G$–$\Sigma$ action unchanged (3.42).
6. **Fortuitous ≫ monotonous in chaos**: LMRS spectrum of $\hat O_f$ shows GUE spacing and long linear ramp (Fig. 7); $S_{\rm info}$: 8.06 (fortuitous) vs 5.19 (monotonous) vs 8.074 (random) (3.54); $S_{EE}$: 6.97 vs 2.70 (3.56).
7. **Following $N$ / chaos invasion** (§4.2, Fig. 8, Table 1): in SYK $q=3$, sector $R=4$, $3^{N/2-1}=81$ states become BPS at $N=10$; information entropy stays near random throughout (Fig. 9). Gapless sectors at $R=\frac N2\pm\frac32$ for odd $N$ (4.5).
8. **Sparse SYK** (§5.1): concentration robust down to $p_s\sim10^{-2}$ at $N=14$; working definition $p_s^{\rm crit}(N,\alpha)$ (5.1) expected to scale as $c(\alpha)/N^\gamma$ with $\gamma\gtrsim1$ (5.2, Fig. 12). Deconcentration = index no longer saturated.
9. **$\mathcal N=4$ SYM**: the non-linear charge constraint of BPS black holes is reinterpreted as R-charge concentration of *core* fortuitous states (5.6)–(5.7). Schematic $Q=\mathrm{Tr}\,\Psi\Psi\frac{\delta}{\delta\Psi}$ (5.8) is SYK-like with sparseness $p_s\sim\tilde N^{-3/2}$ ($\tilde N=N^2$).
10. **Lie-algebra / matrix reformulation** (5.9)–(5.11): $Q=\sum f_{ijk}\psi^i\psi^j\bar\psi_k$ with $f$ structure constants; $Q=-\frac i2\mathrm{Tr}\,\psi\psi\bar\psi$, $H=\mathrm{Tr}\,\psi\psi\bar\psi\bar\psi+\mathrm{Tr}\{\psi,\bar\psi\}^2$ (their normalisation $\mathrm{Tr}T_iT_j=\delta_{ij}$). Left open whether it concentrates. (Chen 2025 footnote 2 states this model has no fortuitous states.)
11. **Discussion** (§6): other solutions of the Turiaci–Witten $Q^2=0$ constraint are *not* concentrated; naive counting of independent components of $Q$ would favour *spread-out* BPS distributions — an open puzzle. Fortuity ↔ null states ("phantom black holes"); monotony replaces the gauge principle in disordered models.

## Important equations

(1.1)–(1.5), (1.7)–(1.9), (1.11)–(1.14), Conjecture 1, (2.6)–(2.8), (2.11)–(2.19), Defs. 1–2, (3.1), (3.5), (3.11)–(3.14), (3.18)–(3.22), (3.29), (3.46)–(3.49), (3.52)–(3.56), (4.1)–(4.5), (5.1)–(5.2), (5.8)–(5.11).

## Numerical methods

Exact diagonalisation of SYK-type models up to $N\sim16$ (sparse), two-flavor $N=9$ and $K=2$, $N=4$ (16 complex fermions); projected-operator spectra, nearest-neighbour spacing vs GUE, Gaussian-filtered spectral form factor; Haar-averaged information/entanglement entropies; continuous-$N$ interpolation (4.2)–(4.3); sampling over sparse couplings.

## Relevant figures/results

Fig. 2 (concentration window); Fig. 3 ($N=9$ two-flavor BPS table); Fig. 5 (mono vs fortuitous lines); Figs. 6–7 (LMRS); Fig. 8–9 (following $N$); Fig. 11–12 (sparse deconcentration); Table 1 (gapped/gapless sectors).

## Limitations

- Conjecture 1 is conjectural; tested only on disordered models. Genericity is not defined for non-disordered theories with extra symmetry.
- Large-$N$ statements locate $R_c$ only approximately; concentration is a microscopic (finite-$N$) statement.
- The Lie-algebra model (5.9) is left unstudied; the paper does not treat cubic purely-creation matrix supercharges like Chen's $\mathrm{Tr}\Psi^3$.
- No $N$-scaling of the Thouless time (footnote 24).

## Relationship to our project

Supplies the **precise definitions** we must use: fortuity vs monotony (Defs. 1–2), irreducible cochain complex, and R-charge concentration (single degree per irreducible complex with macroscopic index). Two consequences matter for us:

1. For a cubic supercharge, $\mathbb Z_3$ ($N_\Psi$ mod 3) refines the complex; **for the matrix models, the $SU(N)$ (and flavor) symmetry commuting with $Q$ refines it further** (irreps are additional labels). A BPS window of *three consecutive* charges is therefore consistent with strict concentration; a window of $\ge4$ consecutive charges is not. This is exactly why Chen's single-matrix model ($N+1$ sectors) fails for $N\ge3$ and why the SYK comparison window is $\frac{N-q}2\le p\le\frac{N+q}2$.
2. The refined index (footnote 11 style, computed per $\mathbb Z_3$ class and per $SU(N)$ irrep) gives *rigorous* lower bounds on BPS counts that our bootstrap results must respect; conversely, the "index saturation" criterion (deconcentration = index no longer saturated, §5.1) is a clean diagnostic we can compute combinatorially for the 3-matrix model.

The chaos diagnostics (LMRS, $S_{\rm info}$, ramp) define what "chaotic" should mean for Chen's proposal; the sparseness discussion (§5.1, (5.8)) provides the natural way to think about *why* a structured cubic matrix supercharge might or might not concentrate. The relation between the single-matrix model and the "phantom saddle"/null-state picture is developed further in Chen 2025.
