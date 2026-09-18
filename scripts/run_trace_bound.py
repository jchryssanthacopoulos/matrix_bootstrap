"""Exact finite-N trace bootstrap: rigorous lower bound on E_0(k; N) (docs/trace_bootstrap_plan.md).

The constraint set is built once (coefficients polynomial in N) and evaluated at every N in --N_list, e.g.
    cd src && PYTHONPATH=. ../.venv/bin/python ../scripts/run_trace_bound.py --p 3 --k 2 --L_adj 3 --L_sing 3 --L_eom 4 \
        --N_list 2,3,4,10,100 --solver scs --out ../results/data/trace_bounds.jsonl
--fixed_N --L_eom_gram 6 --workers 8 selects the parallel fixed-N build used for level 4 (EOM/sector rows only for
Gram entries of total length <= L_eom_gram; the length-8 entries' commutators are almost information-free).
--finiteN_len L adds the finite-N trace relations (antisymmetriser over N+1 indices, total length <= L) for the
integer --N_rel (they are specific to that N; use a single N in --N_list then).  --check compares with the exact
sector ground state when the Fock space is small (N=2, or p=1 with N<=4).  One JSON line per (N) is appended.
"""
import argparse, json, os, resource, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from fermion_matrix_model import build_model, chen_C
from trace_bootstrap import TraceSDP, MemoryBudgetExceeded


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--p', type=int, default=3); ap.add_argument('--k', type=int, required=True)
    ap.add_argument('--N_list', default='2')
    ap.add_argument('--L_adj', type=int, default=2); ap.add_argument('--L_sing', type=int, default=2); ap.add_argument('--L_eom', type=int, default=3)
    ap.add_argument('--finiteN_len', type=int, default=0); ap.add_argument('--N_rel', type=int, default=None)
    ap.add_argument('--plain_adjoint', action='store_true', help='plain adjoint channel instead of the Cho et al. split')
    ap.add_argument('--no_fourier', action='store_true')
    ap.add_argument('--solver', default='clarabel', choices=['clarabel', 'scs'])
    ap.add_argument('--eps', type=float, default=1e-8); ap.add_argument('--max_iters', type=int, default=50000)
    ap.add_argument('--budget_gb', type=float, default=10.0)
    ap.add_argument('--check', action='store_true', help='evaluate all constraints on the exact ground state (small N only)')
    ap.add_argument('--fixed_N', action='store_true', help='parallel build at the single N in --N_list with coefficients evaluated (level 4)')
    ap.add_argument('--L_eom_gram', type=int, default=None, help='EOM/sector rows only for Gram entries of total length <= this (fixed_N path)')
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    C = chen_C(a.p) if a.p > 1 else np.ones((1, 1, 1))
    Ns = [int(x) for x in a.N_list.split(',')]
    t0 = time.time()
    if a.fixed_N:
        assert len(Ns) == 1 and not a.finiteN_len, "--fixed_N takes a single N and no finite-N relations"
        S = TraceSDP.fixed_N(C, a.p, a.k, Ns[0], a.L_adj, a.L_sing, a.L_eom, L_eom_gram=a.L_eom_gram,
                             adjoint_projected=not a.plain_adjoint, fourier=not a.no_fourier, workers=a.workers, verbose=True)
    else:
        S = TraceSDP(C, a.p, a.k, a.L_adj, a.L_sing, a.L_eom, finiteN_len=a.finiteN_len, N_for_relations=a.N_rel,
                     adjoint_projected=not a.plain_adjoint, fourier=not a.no_fourier, verbose=True)
    for N in Ns:
        try:
            res = S.solve(N, solver=a.solver, eps=a.eps, max_iters=a.max_iters, budget_gb=a.budget_gb)
        except MemoryBudgetExceeded as e:
            print("refused:", e); continue
        rec = dict(p=a.p, N=N, k=a.k, L_adj=a.L_adj, L_sing=a.L_sing, L_eom=a.L_eom, finiteN_len=a.finiteN_len, N_rel=a.N_rel,
                   fixed_N=a.fixed_N, L_eom_gram=a.L_eom_gram,
                   adjoint_projected=not a.plain_adjoint, fourier=S.fourier, solver=a.solver, eps=a.eps,
                   n_monomials=len(S.monos), n_vars=res['n'], n_eq=res['n_eq'], n_rows=res['n_rows'], nnz=res['nnz'],
                   cones=[len(c['labels']) for c in S.cones], bound=res['value'], status=res['status'],
                   build_s=round(S.build_time), solve_s=round(res['solve_time'], 1),
                   peak_rss_gb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 2), date=time.strftime('%Y-%m-%d'))
        if a.k == 0:
            rec['exact'] = 16 * N ** 3 - 15 * N if a.p == 3 else None
        elif a.k == 1 and a.p == 3:
            rec['exact'] = 16 * N ** 3 - 53 * N
        if a.check and ((a.p == 3 and N == 2) or (a.p == 1 and N <= 4)):
            model = build_model(N, a.p, C, fourier=S.fourier)
            chk = S.check_exact(model)
            rec.update(exact=chk['E0'], exactGS_max_row_residual=chk['max_row_residual'],
                       exactGS_min_cone_eig=min(chk['cone_min_eig']))
        ex = rec.get('exact')
        print(f"RESULT p={a.p} N={N} k={a.k} level ({a.L_adj},{a.L_sing},{a.L_eom}) finiteN_len={a.finiteN_len}: bound {rec['bound']:.6f}"
              + (f"  exact {ex:.6f}" if ex is not None else '') + f"  [{rec['status']}, {a.solver}, {rec['solve_s']}s] peak {rec['peak_rss_gb']} GB", flush=True)
        if a.out:
            with open(a.out, 'a') as f:
                f.write(json.dumps(rec) + '\n')


if __name__ == '__main__':
    main()
