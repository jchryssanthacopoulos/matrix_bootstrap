#!/usr/bin/env python3
"""
run_single_matrix.py   --   MODEL: Chen's SINGLE-matrix model  Q = Tr[Psi^3]
                             (arXiv:2511.00790).  This is the exactly-solvable
                             VALIDATION model, NOT the 3-matrix SYK target.

What it does
------------
Computes the R-charge-resolved ground-state energy  E_0(N_Psi)  by matrix-free
Lanczos, and locates the BPS window (where E_0 = 0).

Known exact answer (to validate against): the BPS window is
        [ N(N-1)/2 ,  N(N+1)/2 ]      ->  exactly  N+1  charge sectors.
The point of running higher N is to confirm the window keeps GROWING as ~N
(this is the "not concentrated" baseline that the 3-matrix model must beat).

Usage
-----
    python run_single_matrix.py 5           # run N=2..5
    python run_single_matrix.py 5 --only 5  # run only N=5

Feasibility (charge-sector dim = C(N^2, N^2/2)):
    N=4  -> 12,870      seconds
    N=5  -> 5,200,300   minutes per charge, a few GB RAM
    N=6  -> ~9e9        NOT feasible
Needs mqm_fermion_ed.py in the same folder.  Requires numpy, scipy.
"""
import sys, json, time
from math import comb
import mqm_fermion_ed as ed

def run(N):
    lo, hi = N * (N - 1) // 2, N * (N + 1) // 2
    print(f"\n===== SINGLE MATRIX (p=1), N={N} =====")
    print(f"  predicted BPS window [N(N-1)/2, N(N+1)/2] = [{lo},{hi}]  ({hi-lo+1} sectors)")
    prof = ed.profile(N, p=1)
    window = sorted(k for k, d in prof.items() if d['E0'] is not None and d['E0'] < 1e-4)
    ok = (window == list(range(lo, hi + 1)))
    print(f"  measured BPS window = {window}   matches prediction: {ok}")
    return dict(model="single_matrix_p1", N=N, dim=2 ** (N * N),
                predicted_window=[lo, hi], measured_window=window,
                matches=bool(ok), profile={int(k): v for k, v in prof.items()})

if __name__ == "__main__":
    Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    only = None
    if "--only" in sys.argv:
        only = int(sys.argv[sys.argv.index("--only") + 1])
    Ns = [only] if only else list(range(2, Nmax + 1))
    results = []
    t0 = time.time()
    for N in Ns:
        results.append(run(N))
    out = dict(kind="single_matrix_baseline", results=results,
               elapsed_sec=round(time.time() - t0, 1))
    with open("results_single_matrix.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\n==================  REPORT BACK  ==================")
    print("MODEL: single matrix (validation).  Paste the block below:")
    for r in results:
        print(f"  N={r['N']}: window={r['measured_window']}  matches={r['matches']}")
    print(f"Saved -> results_single_matrix.json   ({out['elapsed_sec']}s)")
