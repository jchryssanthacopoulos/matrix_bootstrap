#!/usr/bin/env python3
"""
Phase 2a of docs/nearbps_bootstrap_plan.md: the supermultiplet eigen-bootstrap for one R-charge sector at small N
(Hilbert-space engine, src/sector_bootstrap.py, symmetry-reduced, Fourier flavour letters).

Constraints on the functional phi (multiplet-averaged eigenprojector of the target level):
  * the usual sector bootstrap (adjoint/singlet Gram cones, EOM, sector), optional
  * lower-member rows  phi(Q X) = 0  for every component X of every charge -3 word   (Qbar psi = 0),
  * eigen rows         phi(X H) = E phi(X)  for every registered neutral X, with phi(H) = E fixed,
and E is scanned over E > 0.  E = 0 is excluded by the exact cohomology input (no BPS state in the sector).  A sector
ground state below the window is a lower member, and E_0(k) = min(L(k), L(k-3)) with L the lowest lower-member
energy, so a scan with lower-member rows is a rigorous probe of E_0(k) once L(k-3) is known to be larger.

    python scripts/run_supermultiplet_scan.py --N 2 --k 4 --level 2 --lower_member --eigenstate \
        --out results/data/supermultiplet_scan_N2_k4_L2.json

Prints the convex bound min phi(H) (with the chosen rows, no fixed E), then feasibility at each E of the grid.
"""
import argparse, json, os, resource, sys, time, warnings
import numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
from fermion_matrix_model import build_model, chen_C                          # noqa: E402
from sector_bootstrap import SectorSDP, plan_sector                          # noqa: E402

LEVELS = {2: (2, 2, 3), 3: (3, 3, 4)}


def peak_gb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1e9 if sys.platform == 'darwin' else r / 1e6


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--N', type=int, default=2)
    ap.add_argument('--k', type=int, required=True)
    ap.add_argument('--level', type=int, default=2, choices=sorted(LEVELS))
    ap.add_argument('--lower_member', action='store_true')
    ap.add_argument('--eigenstate', action='store_true')
    ap.add_argument('--E', default=None, help='comma-separated energies; default: a grid in (0, 0.45] plus exact levels')
    ap.add_argument('--solver', default='SCS')
    ap.add_argument('--budget_gb', type=float, default=8.0)
    ap.add_argument('--skip_convex', action='store_true', help='skip the convex min phi(H) solve (scan only)')
    ap.add_argument('--dry_run', action='store_true')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    L_adj, L_sing, L_eom = LEVELS[a.level]
    model = build_model(a.N, 3, chen_C(3), fourier=True)
    plan = plan_sector(model, a.k, L_adj, L_sing, L_eom, True, False, True)
    print(f"plan: d={plan['d']}, D={plan['D']}, cones={plan['cones']}, est {plan['est_gb']:.2f} GB", flush=True)
    if a.dry_run:
        return
    t0 = time.time()
    S = SectorSDP(model, a.k, L_adj, L_sing, L_eom, add_Q=True, gs=False, eigenstate=a.eigenstate, verbose=True,
                  budget_gb=a.budget_gb, symmetry=True, adjoint_projected=True, lower_member=a.lower_member)
    lam = np.linalg.eigvalsh(S.Hk.toarray())
    levels = sorted(set(np.round(lam, 8)))[:6]
    conv = dict(value=None, status='skipped') if a.skip_convex else S.solve('H', 'min', solver=a.solver)
    cv = 'skipped' if conv['value'] is None else f"{conv['value']:.6f}"
    print(f"N={a.N} k={a.k} level {a.level} lower_member={a.lower_member} eigenstate={a.eigenstate}: "
          f"convex bound min phi(H) = {cv} [{conv['status']}]; exact levels {[float(x) for x in levels]}  "
          f"[build {S.build_time:.0f}s, peak {peak_gb():.2f} GB]", flush=True)
    if a.E:
        grid = [float(x) for x in a.E.split(',')]
    else:
        grid = sorted(set([1e-4, 1e-3, 5e-3, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.085, 0.0875, 0.0885, 0.09,
                           0.1, 0.12, 0.15, 0.2, 0.25, 0.3, 0.4, 0.45] + [float(x) for x in levels if 0 < x < 0.45]))
    rows = []
    for E in grid:
        f = S.solve(None, 'min', fixE=E, solver=a.solver)
        feas = f['status'] in ('optimal', 'optimal_inaccurate')
        rows.append(dict(E=E, feasible=feas, status=f['status']))
        tag = '  <- exact level' if any(abs(E - x) < 1e-7 for x in levels) else ''
        print(f"  E={E:11.8f}: {'FEASIBLE  ' if feas else 'infeasible'} [{f['status']}]{tag}", flush=True)
    rec = dict(params=vars(a), level=(L_adj, L_sing, L_eom), d=S.d, r=S.r, n_eq=int(S.eom_rows.shape[0]),
               convex_bound=conv['value'], convex_status=conv['status'], exact_levels=[float(x) for x in levels],
               scan=rows, build_s=round(S.build_time, 1), total_s=round(time.time() - t0, 1),
               peak_rss_gb=round(peak_gb(), 2))
    if a.out:
        os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
        json.dump(rec, open(a.out, 'w'), indent=1)
        print('saved', a.out)


if __name__ == '__main__':
    main()
