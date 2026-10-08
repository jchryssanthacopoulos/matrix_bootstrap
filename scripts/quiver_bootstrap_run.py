#!/usr/bin/env python3
"""
Quiver singlet bootstrap runs (src/quiver_bootstrap.py; research/notes/quiver_project.md section 12).

For each edge occupation m (degree k = 3m) builds the singlet-sector SDP and runs the requested modes:
  energy       min phi(H) without BPS rows: a lower bound on the lowest singlet energy E_0(k);
  margin       BPS rows on, t* = max t with every cone >= t 1 (t* clearly < 0: no singlet BPS state at k);
  feasibility  BPS rows on, pure feasibility with a verified Farkas certificate.
One JSON record per (m, mode) is appended to --out.  Couplings: integer couplings in 1..5 from
np.random.default_rng(seed) (seed 3 reproduces the couplings of the stored exact singlet spectra).

    python scripts/quiver_bootstrap_run.py --n 2 --p 2 --m 1-3 --modes energy margin --out results/data/quiver_bootstrap.jsonl
"""
import argparse, json, os, resource, subprocess, sys, time
import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
from quiver_bootstrap import QuiverSDP                                         # noqa: E402


def parse_list(spec):
    out = []
    for part in spec.split(','):
        if '-' in part:
            a, b = part.split('-'); out += list(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


def git_rev():
    try:
        return subprocess.check_output(['git', '-C', ROOT, 'rev-parse', '--short', 'HEAD'], text=True).strip()
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--p', type=int, required=True)
    ap.add_argument('--m', required=True, help='edge occupations, e.g. 1-5 or 2,4')
    ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('--modes', nargs='+', default=['energy'], choices=['energy', 'margin', 'feasibility'])
    ap.add_argument('--L-adj', type=int, default=2)
    ap.add_argument('--L-sing', type=int, default=3)
    ap.add_argument('--L-eom-gram', default='auto')
    ap.add_argument('--L-cas', default='auto')
    ap.add_argument('--L-bps', type=int, default=None)
    ap.add_argument('--L-bpsH', default='auto')
    ap.add_argument('--no-casimir', action='store_true')
    ap.add_argument('--gauss', type=int, nargs=2, default=None, metavar=('L_Y', 'L_W'))
    ap.add_argument('--finite-n-len', type=int, default=0)
    ap.add_argument('--cone-casimir', action='store_true', help='sandwiched-Casimir (cone irrep) rows, D21.4b')
    ap.add_argument('--L-ccas', default='auto')
    ap.add_argument('--z3', action='store_true', help='impose the quiver-rotation symmetry (needs cyclic couplings)')
    ap.add_argument('--solver', default='clarabel', choices=['clarabel', 'scs'])
    ap.add_argument('--eps', type=float, default=1e-8)
    ap.add_argument('--max-iters', type=int, default=100000)
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--exact', default=None, help='JSON dict {k: exact singlet E_0} for the record (optional)')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    cap = lambda s: None if s in (None, 'None', 'none') else ('auto' if s == 'auto' else int(s))
    C = np.random.default_rng(a.seed).integers(1, 6, size=(a.p, a.p, a.p)).astype(float)
    exact = {int(k): v for k, v in json.loads(a.exact).items()} if a.exact else {}
    for m in parse_list(a.m):
        k = 3 * m
        builds = {}
        for mode in a.modes:
            bps = mode != 'energy'
            if bps not in builds:
                S = QuiverSDP(C, a.n, m, L_adj=a.L_adj, L_sing=a.L_sing, L_eom_gram=cap(a.L_eom_gram), casimir=not a.no_casimir,
                              L_cas=cap(a.L_cas), gauss=tuple(a.gauss) if a.gauss else None, bps=bps, L_bps=a.L_bps,
                              L_bpsH=cap(a.L_bpsH), finite_n_len=a.finite_n_len, cone_casimir=a.cone_casimir,
                              L_ccas=cap(a.L_ccas), z3=a.z3, workers=a.workers, verbose=True)
                builds[bps] = S
            S = builds[bps]
            fn = dict(energy=S.energy, margin=S.margin, feasibility=S.feasibility)[mode]
            r = fn(solver=a.solver, eps=a.eps, max_iters=a.max_iters)
            rec = dict(n=a.n, p=a.p, m=m, k=k, seed=a.seed, couplings=C.astype(int).tolist(), mode=mode, solver=a.solver,
                       eps=a.eps, level=S.level, n_monomials=len(S.monos), n_rows=len(S.rows),
                       cones=sorted((c['size'] for c in S.cones), reverse=True), n_vars=r['n'], n_eq=r['n_eq'],
                       status=r['status'], solve_s=round(r['solve_time'], 1), build_s=round(S.build_time, 1),
                       peak_rss_gb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 2),
                       exact_E0=exact.get(k), git=git_rev(), date=time.strftime('%Y-%m-%d %H:%M'))
            if mode == 'energy':
                rec['bound'] = r['obj']
                msg = f"E_0 >= {r['obj']:.6f}" + (f"  (exact {exact[k]:.6f})" if k in exact else '')
            elif mode == 'margin':
                rec['margin'] = r['margin']
                msg = f"t* = {r['margin']:+.3e}  " + ('EXCLUDED' if r['margin'] < -1e-5 else 'not excluded')
            else:
                for key in ('cert_ATy', 'cert_by', 'cert_min_eig'):
                    rec[key] = r.get(key)
                infeas = 'infeasible' in r['status'].lower()
                msg = f"{r['status']}  cert: |A^T y| {r.get('cert_ATy', float('nan')):.1e}, b.y {r.get('cert_by', float('nan')):+.2e}, " \
                      f"min eig {r.get('cert_min_eig', float('nan')):.1e}" + ('  EXCLUDED' if infeas else '')
            print(f"(n,p)=({a.n},{a.p}) k={k:2d} {mode:11s}: {msg}  [{r['status']}, {r['solve_time']:.0f}s, "
                  f"{r['n']} vars, {r['n_eq']} eq, cones {rec['cones'][:4]}]", flush=True)
            if a.out:
                with open(a.out, 'a') as f:
                    f.write(json.dumps(rec) + '\n')


if __name__ == '__main__':
    main()
