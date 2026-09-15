#!/usr/bin/env python3
"""
run_matrix_syk.py   --   MODEL: the 3-matrix "matrix SYK" TARGET
                         Q = sum_{1<=i<=j<=k<=p} Tr[Psi_i Psi_j Psi_k]  (default p=3).

What it does
------------
Computes E_0(N_Psi) per R-charge sector and reports the BPS window (E_0 = 0) and a
concentration metric.  THIS is the physics question: does the BPS window stay O(1)
(R-charge concentration) or grow with N like the single-matrix baseline?

Known (from in-session ED):
    p=3, N=2 : BPS window = {5,6,7}   (3 sectors, counts 243:486:243)  -- reproduced here.
The prize is p=3, N=3: if the window is still ~O(1) (e.g. {12..15}) while the single
matrix at N=3 already spans 4 sectors and keeps growing, that is R-charge concentration
demonstrated one N beyond what dense ED can reach.

Usage
-----
    python run_matrix_syk.py --N 2                 # p=3, N=2  (fast, reproduces window)
    python run_matrix_syk.py --N 3 --center 4      # p=3, N=3, charges within +-4 of center
    python run_matrix_syk.py --N 3 --p 2           # lighter 2-matrix variant (feasible)

Feasibility (sector dim = C(p*N^2, center)):
    p=3,N=2 : 924            seconds
    p=2,N=3 : ~48,600        minutes            <- good feasible probe of a cubic model
    p=3,N=3 : ~20,000,000    HOURS/charge, tens of GB RAM  <- only on a big workstation;
              run just a few charges around the center with --center 3 or 4.
Needs mqm_fermion_ed.py in the same folder.  Requires numpy, scipy.
"""
import sys, json, time
from math import comb
import mqm_fermion_ed as ed

def parse(flag, default):
    if flag in sys.argv:
        return int(sys.argv[sys.argv.index(flag) + 1])
    return default

if __name__ == "__main__":
    N = parse("--N", 2)
    p = parse("--p", 3)
    M = p * N * N
    center = M // 2
    half = parse("--center", M)          # how many charges each side of center to scan
    lo = max(0, center - half); hi = min(M, center + half)
    charges = list(range(lo, hi + 1))

    print(f"===== MATRIX SYK (p={p}), N={N}, {M} complex fermions, full dim=2^{M} =====")
    print(f"  scanning N_Psi in [{lo},{hi}] (center {center}); sector dims up to C({M},{center})={comb(M,center):,}")
    t0 = time.time()
    prof = ed.profile(N, p=p, charges=charges)
    window = sorted(k for k in charges if prof[k]['E0'] is not None and prof[k]['E0'] < 1e-4)
    single_matrix_window_width = N + 1     # baseline for comparison at this N
    out = dict(kind="matrix_syk", model=f"p{p}_N{N}", p=p, N=N, fermions=M,
               scanned=charges, bps_window=window,
               window_width=len(window),
               single_matrix_baseline_width=single_matrix_window_width,
               profile={int(k): prof[k] for k in charges},
               elapsed_sec=round(time.time() - t0, 1))
    fname = f"results_matrix_syk_p{p}_N{N}.json"
    with open(fname, "w") as f:
        json.dump(out, f, indent=2)

    print("\n==================  REPORT BACK  ==================")
    print(f"MODEL: matrix SYK p={p}, N={N}.  Paste the block below:")
    print(f"  BPS window (E_0=0 sectors) = {window}")
    print(f"  window width = {len(window)}   vs single-matrix baseline width at this N = {N+1}")
    if window:
        print(f"  -> {'CONCENTRATED (narrower than baseline)' if len(window) < N+1 else 'not narrower than baseline'}")
    print(f"Saved -> {fname}   ({out['elapsed_sec']}s)")
