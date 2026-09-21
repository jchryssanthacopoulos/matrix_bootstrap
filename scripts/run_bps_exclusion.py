"""BPS-exclusion bootstrap: feasibility margin of the sector-k functional with the BPS rows phi(XQ) = phi(QX) = 0
(and the Qbar conjugates), see docs/trace_bootstrap_plan.md and research/notes/bps_exclusion_results.md.

    cd src && PYTHONPATH=. ../.venv/bin/python ../scripts/run_bps_exclusion.py --p 3 --N 3 --k_list 0-27 \
        --L_adj 2 --L_sing 2 --L_eom 3 --out ../results/data/bps_exclusion.jsonl
Margin t* = max t s.t. every cone >= t*1: t* ~ 0 -> the sector cannot be excluded at this level;
t* < 0 (clearly) -> no BPS state in sector k (numerical certificate; rows are exact operator identities).
The same constraint set is evaluated at each N in --N_list (the sector index k is just a number).
"""
import argparse, json, os, sys, time, resource
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from fermion_matrix_model import chen_C
from trace_bootstrap import TraceSDP


def parse_list(spec):
    out = []
    for part in spec.split(','):
        if '-' in part:
            a, b = part.split('-'); out += list(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--p', type=int, default=3); ap.add_argument('--N_list', default='3'); ap.add_argument('--k_list', required=True)
    ap.add_argument('--L_adj', type=int, default=2); ap.add_argument('--L_sing', type=int, default=2); ap.add_argument('--L_eom', type=int, default=3)
    ap.add_argument('--L_bps', type=int, default=None)
    ap.add_argument('--no_bps', action='store_true', help='control run without the BPS rows')
    ap.add_argument('--solver', default='clarabel', choices=['clarabel', 'scs']); ap.add_argument('--eps', type=float, default=1e-8)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    C = chen_C(a.p) if a.p > 1 else np.ones((1, 1, 1))
    Ns = parse_list(a.N_list); ks = parse_list(a.k_list)
    for k in ks:
        t0 = time.time()
        S = TraceSDP(C, a.p, k, a.L_adj, a.L_sing, a.L_eom, bps=not a.no_bps, L_bps=a.L_bps, verbose=False)
        for N in Ns:
            if k > a.p * N * N:
                continue
            r = S.bps_margin(N, solver=a.solver, eps=a.eps)
            rec = dict(p=a.p, N=N, k=k, L_adj=a.L_adj, L_sing=a.L_sing, L_eom=a.L_eom, L_bps=a.L_bps, bps=not a.no_bps,
                       solver=a.solver, margin=r['margin'], status=r['status'], n=r['n'], n_eq=r['n_eq'], n_rows=r['n_rows'],
                       n_bps_rows=S.n_bps_rows, n_monomials=len(S.monos), build_s=round(S.build_time), solve_s=round(r['solve_time'], 1),
                       peak_rss_gb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 2), date=time.strftime('%Y-%m-%d'))
            verdict = 'EXCLUDED' if r['margin'] < -1e-4 else 'not excluded'
            print(f"N={N} k={k:2d} level ({a.L_adj},{a.L_sing},{a.L_eom}) bps={not a.no_bps}: t* = {r['margin']:+.5f}  {verdict}  [{r['status']}, {r['solve_time']:.0f}s]", flush=True)
            if a.out:
                with open(a.out, 'a') as f:
                    f.write(json.dumps(rec) + '\n')


if __name__ == '__main__':
    main()
