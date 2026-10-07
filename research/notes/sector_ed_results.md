# Sector ground energies up to the BPS window at $N=2,3$ (exact diagonalisation, 2026-10-06)

Purpose: ground truth near the BPS window, where the super-Schwarzian question lives (Q5) and the bootstrap is
weakest. Up to now the only near-window data were the $N=2$ values. This note adds every sector below the window at
$N=3$, with the irrep of each ground state.

## Method

- **Model.** Chen's three-matrix model ($p=q=3$, cyclic $C$, $U(N)$, no singlet projection), $H$ in the D2.13 form
  $H=E_c+H_1+9\,Y^\dagger Y$ (D16.2).
- **Block.** For each sector $k$, the *zero-weight* block (all $N$ Cartan charges zero). Every irrep occurring
  here has weights in the root lattice and hence a zero-weight state, so the lowest eigenvalue of this block is the
  ground energy of the whole sector. A level of irrep $\lambda$ appears $K_{\lambda,0}$ times per copy
  ($K_{\lambda,0}$ = zero-weight multiplicity: 1 for $\mathbf 1,\mathbf{10},\overline{\mathbf{10}}$, 2 for $\mathbf 8$,
  3 for $\mathbf{27}$; 1 for every $SU(2)$ irrep).
- **Code.** `src/sector_ed.py` (vectorised bit-mask builder `build_fast`, Casimir labels, irrep tables),
  `scripts/run_sector_ed.py` (driver). The reference builder is `free_sectors.sector_operators`.
- **Solver.** Dense diagonalisation for blocks $\le3000$; otherwise `scipy` `eigsh` (ARPACK Lanczos, `which='SA'`,
  12 eigenvalues, `ncv=40`, `tol=1e-10`, seeded start vector, seed 1).
- **Labels.** $\langle\hat C_2\rangle$ evaluated on each eigenvector, matched to the occurring irreps.
- **Validation.**
  - Every run checks the D2.13 builder against $\{Q,Q^\dagger\}$ applied directly to a random 40-state vector. The
    largest deviation was $1.4\times10^{-12}$.
  - The two builders give identical spectra at $N=3$, $k=8,9$ ($E_0$ agrees to $5\times10^{-15}$).
  - Known values are reproduced: $N=2$ full ED (98, 22, 5.16536, 1.22706, 0.08796); $N=3$ few-body ED (75.7917 at
    $k=3$, 46.271 at $k=4$); D16 (273, 159).
  - All Lanczos residuals are $\le3\times10^{-10}$.
- **Caveat.** Plain Lanczos can miss copies of exactly degenerate levels, so the multiplicities it reports are lower
  bounds; true degeneracies follow from the irrep label. Irreps sharing a Casimir value ($\mathbf{10}$ and
  $\overline{\mathbf{10}}$) are not separated.
- **Data.** `results/data/sector_ed_N{2,3}_k*_fast.json`, plus `_python.json` reference runs at $N=3$, $k=8,9$.
  The JSON files record parameters, timings, peak memory and the git revision (`0afd4c2+uncommitted`).
- **Cost.** At $N=3$, $k=11$: block 750,699 of the sector's 13,037,895 states; build 10 s; Lanczos 42 min; peak
  3.0 GB.

## Results

$J\equiv N_\Psi-pN^2/2$. Each ground state below the window is the lower member of a $Q$-multiplet $(k,k+3)$. Its
partner cannot sit in sector $k-3$, because that would need an energy below $E_0(k-3)$, which never happens in the
table. Its multiplet charge is therefore $q=k+\tfrac32-\tfrac{pN^2}2$, Turiaci–Witten's label.

**$N=3$** (BPS window $k=12..15$; free sectors $k\le2$ by D16):

| $k$ | $|q|$ | block | $E_0(k;3)$ | ground-state irrep | next levels (irrep) |
|---|---|---|---|---|---|
| 0 | 12 | 1 | 387 | $\mathbf 1$ | |
| 1 | 11 | 9 | 273 | $\mathbf 8$ | 372 ($\mathbf 8$) |
| 2 | 10 | 63 | 159 | $\mathbf{10}\oplus\overline{\mathbf{10}}$ | 247.69 ($\mathbf 8$) |
| 3 | 9 | 381 | 75.79168 | $\mathbf{27}$ | 133.78 ($\mathbf{10}$) |
| 4 | 8 | 1,854 | 46.27054 | $\mathbf{27}$ | |
| 5 | 7 | 7,254 | 22.39999 | $\mathbf{27}$ | |
| 6 | 6 | 23,388 | 12.56019 | $\mathbf{10}\oplus\overline{\mathbf{10}}$ | 16.14 ($C_2=15$), 16.16 ($\mathbf{27}$) |
| 7 | 5 | 63,486 | 6.63510 | $\mathbf 8$ | 8.19 ($C_2=12$), 8.32 ($\mathbf{10}$) |
| 8 | 4 | 146,943 | 2.73751 | $\mathbf 1$ | 3.634 ($\mathbf 8$), 3.739 ($\mathbf{27}$), 4.280 ($C_2=18$) |
| 9 | 3 | 292,179 | 1.14573 | $\mathbf 1$ | 1.497 ($\mathbf 8$), 1.713 ($\mathbf{27}$), 1.771 ($C_2=15$) |
| 10 | 2 | 502,245 | 0.345019 | $\mathbf 8$ | 0.464 ($\mathbf{10}$), 0.471 ($C_2=15$), 0.502 ($C_2=12$) |
| 11 | 1 | 750,699 | 0.061408 | $\mathbf{10}\oplus\overline{\mathbf{10}}$ | 0.0651 ($\mathbf 8$), 0.0696 ($\mathbf{10}$), 0.0710 ($\mathbf{27}$), 0.0839 ($C_2=12$) |

**$N=2$** (window $k=5..7$): $E_0=98,\ 22,\ 5.16536,\ 1.22706,\ 0.087961$ at $k=0..4$. The ground states are
$\mathbf 1$, $\mathbf 3$, $\mathbf 3$, $\mathbf 1$, $\mathbf 1$, at $|q|=4.5,3.5,2.5,1.5,0.5$. At $k=5,6$ all eight
requested eigenvalues are zero (BPS). *Correction (2026-10-06): an earlier version said "8 zero modes each"; 8 was
only the number of eigenvalues requested (`n_eig`). Dense diagonalisation (`src/bps_ritz.py`) gives 63 and 126
zero-weight zero modes at $k=5,6$, matching `cohomology_N2_full.json` ($h=9,27,18,9$ and $18,54,36,18$ in
$C_2=0,2,6,12$, each with zero-weight multiplicity 1).*

## Comparison with the super-Schwarzian gap law

Turiaci–Witten (3.11) give $E_0(q)=q^2/(4\hat q^2)$ in Schwarzian units for the multiplet of charge $q$, so
$E_0/q^2$ should be constant. Figure: `results/figures/near_window_gaps.pdf` (log–log, with the $q^2$ law through
the point nearest the window).

| | $|q|$ | $E_0$ | $E_0/q^2$ | local exponent $d\ln E_0/d\ln|q|$ to the next point out |
|---|---|---|---|---|
| $N=2$ | 0.5 | 0.087961 | 0.352 | 2.40 |
| | 1.5 | 1.22706 | 0.545 | 2.81 |
| | 2.5 | 5.16536 | 0.826 | |
| $N=3$ | 1 | 0.061408 | 0.0614 | 2.49 |
| | 2 | 0.345019 | 0.0863 | 2.96 |
| | 3 | 1.14573 | 0.127 | 3.03 |
| | 4 | 2.73751 | 0.171 | 3.97 |
| | 5 | 6.63510 | 0.265 | |

**What the data show (observations, not fits; two values of $N$ only).**
- **Faster than $q^2$.** At both $N$ the energies rise faster than $q^2$ away from the window. The local exponent is
  closest to 2 next to the window (2.40 at $N=2$, 2.49 at $N=3$) and steepens outward.
- **The near-window scale drops with $N$.** At comparable $|q|$ the $N=3$ energies are an order of magnitude below
  $N=2$ (e.g. $|q|=2.5$: 5.17 at $N=2$, against 0.345 and 1.15 at $|q|=2,3$ for $N=3$), while the vacuum energy
  grows from 98 to 387. In "Schwarzian units" $36E_0/q^2$ at the point nearest the window is 12.7 at $N=2$ and 2.2 at
  $N=3$.
- **Irreps.** The ground states near the window are low-Casimir irreps ($C_2\le6$), and their order changes from
  sector to sector: singlet at $k=8,9$, adjoint at $k=10$, $\mathbf{10}\oplus\overline{\mathbf{10}}$ at $k=11$. At
  $k=11$, the four lowest levels ($C_2=6,3,6,8$) lie within 16% of each other. At $N=2$ the two sectors nearest the window have singlet
  ground states. **Correction to a remark made earlier on 2026-10-06:** the near-window ground states are not all
  gauge singlets at $N=3$.
- **Not the maximal irrep.** There the spectrum is exactly $\{0,1,\dots,N\}$ (D17), so its gap next to the window is
  1 at every $N$. The small near-window energies all come from low-Casimir irreps.

**How to read this.** The $q^2$ law is a large-$N$ prediction for a model whose BPS states concentrate in
$\hat q=3$ charges. $N=2$ concentrates but is tiny (12 modes). $N=3$ has BPS states in four charges, $J=\pm\tfrac12,\pm\tfrac32$,
which a generic $\hat q=3$ ensemble forbids. So neither point is in the regime where the law is supposed to hold,
and agreement was not expected. What the data do establish is that the energy scale next to the window falls
steeply from $N=2$ to $N=3$, and that the near-window spectrum is a dense mix of low-Casimir irreps. Do not infer
an $N$-scaling exponent from these two points (D1 lesson).

## Next

- **Gaps inside the window** ($k=12..15$): the lowest non-BPS energies above the zero modes. A plain Lanczos run is
  swamped by the BPS kernel, so this needs either the cohomology dimensions to deflate or a block where the kernel
  is small.
- **Finer blocks.** Resolve the $\mathbb Z_3$ flavour charge and the cubic Casimir to split $\mathbf{10}$ from
  $\overline{\mathbf{10}}$, and to get level statistics within a single symmetry sector (Q2).
- **Bootstrap against these values.** Run the Casimir-resolved level-3 bootstrap at $k=9..11$ against the table
  above.
