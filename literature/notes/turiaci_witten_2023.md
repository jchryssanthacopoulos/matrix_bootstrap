# Turiaci, Witten (2023) — $\mathcal N=2$ JT Supergravity and Matrix Models

**File:** `papers/turiaci_witten_2023.pdf` (126 pp.; only §§1–3 read in full, §§4–7 and appendices skimmed)

## Full citation

Gustavo J. Turiaci, Edward Witten, *$\mathcal N=2$ JT Supergravity and Matrix Models*, arXiv:2305.19438v4 [hep-th] (22 Nov 2023).

## Main question

Which random matrix ensemble is dual to $\mathcal N=2$ JT supergravity on surfaces of arbitrary topology? Answer: supermultiplets of different R-charge are statistically independent, each governed by an Altland–Zirnbauer $(\alpha,\beta)=(1+2\nu,2)$ ensemble; BPS states exist, with fixed (non-fluctuating) numbers, only in a window of R-charges of width $\hat q$.

## Physical system

Global $\mathcal N=2$ algebra $Q^2=\bar Q^2=0$, $\{Q,\bar Q\}=H$ (2.1) with $U(1)_R$ generator $J$, $[J,Q]=\hat qQ$, $\hat q$ odd, $(-1)^F=e^{\pm i\pi J}$ (2.2)–(2.3); possible anomalous shift $J\in\delta+\mathbb Z$ (§2.4). $\mathcal N=2$ JT supergravity = $\mathrm{SU}(1,1|1)$ BF theory; boundary $\mathcal N=2$ Schwarzian. Also time-reversal variants (§6) and $\mathcal N=4$ (§7).

## Important definitions

- Charge-sector decomposition $\mathcal H=\oplus_k\mathcal H_k$, $Q=\sum_kQ_k$, $Q_k:\mathcal H_k\to\mathcal H_{k+\hat q}$; $\mathcal H_k=\mathcal H_k^0\oplus\mathcal H_k^+\oplus\mathcal H_k^-$ (BPS, lower member of a $(k,k+\hat q)$ multiplet, upper member of a $(k-\hat q,k)$ multiplet) (2.19)–(2.20).
- **Index constraint** $\sum_k(-1)^kL^0_k=\sum_k(-1)^kL_k$ (2.22).
- **Genericity**: for a generic solution of $Q^2=0$, "for any $k$, either $L^0_k$ or $L^0_{k+1}$ vanishes" (adjacent in units of $\hat q$), since otherwise a small perturbation pairs them into a long multiplet. Components of the solution space are labelled by the $L^0_k$; JT picks one (no wormhole fluctuations of $N_{\rm BPS}$, §3.2).
- **Ensemble** (2.23), (2.33): measure $\prod dQ_k\prod\delta(Q_{s+1}Q_s)\exp(-\sum\mathrm{Tr}F_r(Q_r^\dagger Q_r))$; after gauge fixing $\prod_kU(L_k)$ the $(r,r+\hat q)$ multiplets decouple and are AZ with $\alpha=\alpha_0+\beta\nu=1+2\nu$, $\beta=2$, $\nu$ = number of BPS states of charge $r$ or $r+\hat q$. Four consecutive charges need BV-type treatment of redundant constraints (2.34)–(2.37) but give the same conclusion.
- Multiplet label: average charge $q=k+\hat q/2$; charges $(q-\hat q/2,q+\hat q/2)$.

## Main assumptions

Large $e^{S_0}$ topological expansion; JT (pure dilaton gravity) with no matter; $\delta$ and $\hat q$ as parameters.

## Main analytical results

- **Disk / $\mathcal N=2$ Schwarzian** (3.7)–(3.9): $Z_{\rm disk}(\beta,\alpha)=\sum_ke^{i\alpha k}N_{\rm BPS}(k)+\sum_q(e^{i\alpha(q-\hat q/2)}+e^{i\alpha(q+\hat q/2)})\int dE\,e^{-\beta E}\rho_q(E)$ with
  $$N_{\rm BPS}(k)=\frac{e^{S_0}}{4\pi^2}\cos\frac{\pi k}{\hat q}\ \ (|k|<\tfrac{\hat q}2),\quad 0\ \text{otherwise}\quad(3.10);\qquad \rho_q(E)=\frac{e^{S_0}}{2\pi}\frac{\sinh\big(2\pi\sqrt{E-E_0(q)}\big)}{4\pi^2E},\ E_0(q)=\frac{q^2}{4\hat q^2}\quad(3.11).$$
  BPS states occur only for $|k|<\hat q/2$ — for $\hat q=3$: $k\in\{-1,0,1\}$ with weights $\cos(\pi k/3)=(\tfrac12,1,\tfrac12)$, i.e. **1:2:1**; never two BPS charges differing by $\hat q$ (they would pair up). Threshold behaviour $\rho\sim\sqrt\varepsilon$ for $E_0(q)>0$, $\rho\sim1/\sqrt\varepsilon$ for the $q=0$ multiplet (AZ with $\alpha=1$); spectral curve $y_q(x)$ (3.12) with a pole at $x=0$ whose residue counts BPS states (3.13).
- Cylinder (§3.2): consistent with statistically independent multiplets and non-fluctuating $N_{\rm BPS}$.
- Table (Fig. 1): $\mathcal N=2$ has BPS states for all $(\hat q,\delta)\ne(1,\tfrac12)$, a gap for $\delta\ne\tfrac12$, vanishing plain index but non-vanishing refined index.
- Higher topology: $\mathcal N=2$ Weil–Petersson-type volumes with a Mirzakhani recursion (5.31), matching the loop equations of the proposed ensemble (§5.4).
- **Check on $\mathcal N=2$ SYK** (§2.5, Figs. 3–4): with $\hat q=3$, $N=10,11$, the ratio statistic $r=(\lambda_{j+2}-\lambda_{j+1})/(\lambda_{j+1}-\lambda_j)$ of the *singular values of $Q_k:\mathcal H_k\to\mathcal H_{k+3}$* follows the $\beta=2$ surmise $P_\beta(r)\propto(r+r^2)^\beta/(1+r+r^2)^{1+3\beta/2}$ for every multiplet, except the $CT$-self-conjugate $(-\tfrac32,\tfrac32)$ multiplet at odd $N$ ($\beta=1$). Depends on $\beta$ only, not $\alpha$; no unfolding needed.

## Important equations

(2.1)–(2.3), (2.19)–(2.23), (2.33), (3.7)–(3.13).

## Numerical methods

ED of $\mathcal N=2$ SYK ($N=10,11$, 200 samples) with Jordan–Wigner qubits; singular values of $Q_k$ per charge sector.

## Relevant figures/results

Fig. 1 (spectral properties table), Fig. 2 (ensembles with/without time reversal), Figs. 3–4 (SYK level statistics).

## Limitations

Pure JT: no matter, leading large-$e^{S_0}$ formulas; the ensemble is *postulated* from genericity and matched to gravity, not derived from a microscopic model; other components of the $Q^2=0$ solution space (non-concentrated BPS distributions) are not described by JT (cf. Chang et al. §6).

## Relationship to our project

This paper defines what "super-Schwarzian / supercharge chaos" predicts *quantitatively* for a concentrating cubic model, and hence what our bootstrap and ED should be compared with:

1. **Concentration and the 1:2:1 pattern.** With $J=N_\Psi-\tfrac{pN^2}{2}$ and $\hat q=3$, (3.10) says BPS states sit at three consecutive charges $k=-1,0,1$ with counts in ratio $1{:}2{:}1$ (up to the $e^{S_0}$ normalisation and subleading corrections). The project's $N=2$ three-matrix ED gives $243{:}486{:}243$ at $N_\Psi=5,6,7$ — exactly this pattern (and identical to FGMS's exact SYK counts (5.7)). The single-matrix model gives binomial weights $\binom N{\cdot}$ over $N+1$ charges (Chen (2.20)), which coincide with $1{:}2{:}1$ only at $N=2$. So a sharp, quantitative large-$N$ test is: do the three-matrix BPS counts stay in the ratio $\cos(\pi k/3)$ as $N$ grows?
2. **Gap prediction.** $E_0(q)=q^2/(4\hat q^2)$ in Schwarzian units: multiplets straddling the BPS window have gaps growing quadratically with distance from the centre — a target for $E_0(N_\Psi)$ bounds in lifted sectors (after fixing the overall scale $\propto1/N^\#$).
3. **Chaos diagnostic.** The right object is the singular-value spectrum of $Q_k$ restricted to a $(k,k+3)$ pair of charge sectors (and, for us, a fixed $SU(N)$ irrep), tested with the $r$-ratio against $\beta=2$ (or $\beta=1$ for $CT$-self-conjugate multiplets when the total fermion number is odd). This avoids the huge $H$-degeneracies that made the project's $N=2$ level-statistics test inconclusive and is feasible with sparse SVD/Lanczos on $Q_k^\dagger Q_k$ in sectors of dimension $\sim10^4$–$10^5$.
4. **Structural constraint.** Genericity forbids BPS states at charges $k$ and $k+\hat q$ simultaneously; the single-matrix model violates this for $N\ge3$ (BPS at $N_\Psi=3$ and $6$ in the same irrep $r_*$), which is another way to say its supercharge is non-generic (Casimir Hamiltonian). Whether the three-matrix model is "generic in this sense" at $N\ge3$ is the concentration question.
