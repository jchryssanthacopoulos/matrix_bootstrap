"""Rigorous lower bound on E_0(k) from the sector bootstrap (docs/derivations.md D3/D4), one sector per call.

Example (three-matrix Chen model, N=2, sector k=3, level 3, symmetry-reduced, 10 GB budget):
    cd src && PYTHONPATH=. ../.venv/bin/python ../scripts/run_sector_bound.py --N 2 --p 3 --k 3 \
        --L_adj 3 --L_sing 3 --L_eom 4 --symmetry --solver SCS --budget_gb 10 --out ../results/data/bounds.jsonl
Each run appends one JSON line (parameters, bound, exact value, sizes, timings, peak RSS) to --out.
--dry_run prints the plan (cone sizes, memory estimate) and exits without building anything.
"""
import argparse, json, os, resource, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from fermion_matrix_model import build_model, chen_C
from sector_bootstrap import SectorSDP, plan_sector, MemoryBudgetExceeded, solver_options


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--N', type=int, default=2); ap.add_argument('--p', type=int, default=3)
    ap.add_argument('--k', type=int, required=True)
    ap.add_argument('--L_adj', type=int, default=3); ap.add_argument('--L_sing', type=int, default=3)
    ap.add_argument('--L_eom', type=int, default=4)
    ap.add_argument('--gs', action='store_true', help='add ground-state positivity blocks')
    ap.add_argument('--symmetry', action='store_true', help='SU(2)_gauge x Z_p invariant functional (N=2 only)')
    ap.add_argument('--no_fourier', action='store_true', help='keep the flavor basis (no Z_p grading)')
    ap.add_argument('--adjoint_projected', action='store_true', help='Cho et al. adjoint tensor structure (needs --symmetry)')
    ap.add_argument('--solver', default='SCS'); ap.add_argument('--eps', type=float, default=1e-7)
    ap.add_argument('--max_iters', type=int, default=100000)
    ap.add_argument('--budget_gb', type=float, default=10.0)
    ap.add_argument('--dry_run', action='store_true')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    model = build_model(a.N, a.p, chen_C(a.p) if a.p > 1 else np.ones((1, 1, 1)), fourier=(a.p > 1 and not a.no_fourier))
    t0 = time.time()
    plan = plan_sector(model, a.k, a.L_adj, a.L_sing, a.L_eom, True, a.gs, a.symmetry)
    print(f"plan: d={plan['d']} D={plan['D']} cones(max {max(plan['cones'])}, n={len(plan['cones'])}) sum m^2={plan['sum_m2']} "
          f"Gram entries={plan['n_entries']} est {plan['est_gb']:.2f} GB (budget {a.budget_gb})", flush=True)
    if a.dry_run:
        print(sorted(plan['cones'], reverse=True)); return
    try:
        S = SectorSDP(model, a.k, a.L_adj, a.L_sing, a.L_eom, add_Q=True, gs=a.gs, verbose=True, budget_gb=a.budget_gb,
                      symmetry=a.symmetry, adjoint_projected=a.adjoint_projected)
    except MemoryBudgetExceeded as e:
        print("refused:", e); return
    t1 = time.time()
    prob, y = S._problem(S.cH, 'min')
    solver = a.solver
    try:
        prob.solve(solver=solver, verbose=False, **solver_options(solver, a.eps, a.max_iters))
    except Exception as e:                       # CLARABEL fails on some (not strictly feasible) instances
        print(f"{solver} failed ({str(e)[:80]}); retrying with SCS", flush=True)
        solver = 'SCS'
        prob.solve(solver='SCS', verbose=False, eps=a.eps, max_iters=a.max_iters)
    t2 = time.time()
    peak_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9
    rec = dict(model=f"N={a.N} p={a.p}", k=a.k, L_adj=a.L_adj, L_sing=a.L_sing, L_eom=a.L_eom, gs=a.gs, symmetry=a.symmetry,
               fourier=bool(model.get('fourier')), adjoint_projected=a.adjoint_projected, solver=solver,
               d=S.d, D=(S.reducer.D if S.reducer else None), r=S.r, n_eom_raw=S.n_eom_raw, n_eom=int(S.eom_rows.shape[0]),
               cones=[m for m, _ in S.block_coeffs], bound=float(prob.value), exact=S.exact, status=prob.status,
               build_s=round(t1 - t0), solve_s=round(t2 - t1), peak_rss_gb=round(peak_gb, 2), date=time.strftime('%Y-%m-%d'))
    print(f"RESULT k={a.k} level ({a.L_adj},{a.L_sing},{a.L_eom}) gs={a.gs} sym={a.symmetry}: bound {rec['bound']:.5f}  exact {S.exact:.5f}  "
          f"[{prob.status}] r={S.r} EOM {S.n_eom_raw}->{rec['n_eom']} build {rec['build_s']}s solve {rec['solve_s']}s peak {peak_gb:.2f} GB", flush=True)
    if a.out:
        with open(a.out, 'a') as f:
            f.write(json.dumps(rec) + '\n')


if __name__ == '__main__':
    main()
