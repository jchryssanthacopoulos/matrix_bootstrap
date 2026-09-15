# Local validation scripts — how to run and what to report back

Four files, keep them in the same folder. Requires only `numpy` and `scipy`
(`pip install numpy scipy`).

| File | Role |
|---|---|
| `mqm_fermion_ed.py` | shared engine (matrix-free, R-charge-resolved $E_0$). Don't run directly. |
| `run_single_matrix.py` | **MODEL: Chen's single matrix** $Q=\mathrm{Tr}[\Psi^3]$ (validation) |
| `run_matrix_syk.py` | **MODEL: the 3-matrix SYK target** $Q=\sum_{i\le j\le k}\mathrm{Tr}[\Psi_i\Psi_j\Psi_k]$ |

Both compute the ground-state energy in each fixed $R$-charge sector $E_0(N_\Psi)$ and
report the **BPS window** = the set of $N_\Psi$ where $E_0=0$. Each writes a JSON file
and prints a `REPORT BACK` block — just paste that block to me.

The engine never builds the $2^{pN^2}$ Hilbert space; it works one charge sector at a
time with matrix-free Lanczos, so it reaches further than dense ED. It is validated
in-session: it reproduces the single matrix at $N=2,3$ and the 3-matrix window
$\{5,6,7\}$ at $N=2$ exactly.

---

## 1. Single matrix (the "does-not-concentrate" baseline)

```
python run_single_matrix.py 4          # N=2,3,4   (seconds–minutes)
python run_single_matrix.py 5 --only 5 # just N=5  (minutes/charge, a few GB RAM)
```

Exact expected answer: BPS window $=[N(N{-}1)/2,\,N(N{+}1)/2]$, i.e. **$N+1$ sectors**.
The interesting run is **$N=5$** (window should be $[10,15]$, 6 sectors): it confirms the
window keeps growing $\sim N$. We have this to $N=4$; $N=5$ is the new data point.
$N=6$ is not feasible ($\binom{36}{18}\approx9\times10^{9}$).

## 2. Matrix SYK target (the concentration question)

```
python run_matrix_syk.py --N 2                  # p=3, N=2: reproduces window {5,6,7}
python run_matrix_syk.py --N 3 --p 2            # 2-matrix cubic model, N=3: FEASIBLE probe
python run_matrix_syk.py --N 3 --center 4       # p=3, N=3: only charges within 4 of center
```

The physics question: **does the BPS window stay $O(1)$ (concentration) or grow like the
single-matrix baseline ($N+1$)?**

- **`--N 2`** (p=3): fast; window width 3 $=$ baseline 3 (they only differ at higher $N$).
- **`--p 2 --N 3`** (2-matrix cubic, $2^{18}$): fully feasible in minutes. A genuine
  higher-$N$ cubic-model data point. Report the window; if its width stays $\sim$3 while
  the single-matrix baseline at $N=3$ is 4 and growing, that is concentration in a cubic
  matrix SYK model at a new $N$.
- **`--N 3` (p=3, $2^{27}$)**: the real target but heavy — each central sector is
  $\binom{27}{13}\approx2\times10^{7}$, needing tens of GB and hours per charge. Only on a
  large workstation, and only scan a few charges around the center
  (`--center 3` or `4`). Even confirming $E_0=0$ at $N_\Psi=13,14$ and $E_0>0$ at, say,
  $N_\Psi=10,11$ would pin the window width.

### What to report back
Paste the `REPORT BACK` block (or the JSON). The key numbers are the **BPS window** and
its **width vs the single-matrix baseline** at that $N$.

---

## Notes / knobs
- Speed: increase the Lanczos tolerance for faster (slightly rougher) runs by editing
  `tol` in `mqm_fermion_ed.E0_charge` (default `1e-6`); `1e-4` is usually plenty to tell
  $E_0=0$ from $E_0>0$.
- If a large sector runs out of RAM, scan fewer charges (`--center k`) so only the
  needed sectors ($k-3,k,k+3$) are built.
- `E_0` is clamped at 0 (the true ground energy is $\ge0$ by SUSY); a value like
  $10^{-5}$ means "BPS", a value like $0.09$ or $1.2$ means "lifted".
