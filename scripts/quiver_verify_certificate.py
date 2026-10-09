#!/usr/bin/env python3
"""
Solve the quiver singlet bootstrap's pure feasibility problem (BPS rows on) at (n, p, m) and verify the resulting
Farkas certificate in exact integer arithmetic (src/quiver_certificate.py; research/notes/quiver_project.md section
12).  A verified certificate proves that the singlet sector of degree k = 3m has no BPS state, for these couplings
and hence (by upper semicontinuity of h^k in C) for generic couplings -- given the correctness of the constraint
generation, which is tested separately against explicit Fock-space operators.

    python scripts/quiver_verify_certificate.py --n 3 --p 3 --m 6 --solver clarabel --out results/data/quiver_certificate_checks.jsonl
"""
import argparse, json, os, resource, subprocess, sys, time
import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
from quiver_bootstrap import QuiverSDP                                          # noqa: E402
from quiver_certificate import verify_farkas_exact                              # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--p', type=int, required=True)
    ap.add_argument('--m', type=int, required=True)
    ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('--solver', default='clarabel', choices=['clarabel', 'scs'])
    ap.add_argument('--eps', type=float, default=1e-8)
    ap.add_argument('--max-iters', type=int, default=200000)
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--scale', action='store_true', help='column scaling x = n^w y (needed for conditioning at n >= 4 or p >= 4)')
    ap.add_argument('--time-limit', type=float, default=None)
    ap.add_argument('--no-exact', action='store_true', help='skip the exact verification (fast feasibility verdict only)')
    ap.add_argument('--save-cert', default=None, help='npz file for the certificate vector and its provenance maps')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    C = np.random.default_rng(a.seed).integers(1, 6, size=(a.p, a.p, a.p)).astype(float)
    t0 = time.time()
    S = QuiverSDP(C, a.n, a.m, bps=True, workers=a.workers, verbose=True)
    data = S.assemble('feasibility', a.solver, prune=(a.solver == 'clarabel'), scale=a.scale)
    out = S._solve(data, a.solver, a.eps, a.max_iters, False, time_limit=a.time_limit)
    cert = S.check_certificate(data, out['y'], a.solver)
    print(f"(n,p)=({a.n},{a.p}) k={3*a.m} {a.solver}: {out['status']} in {out['solve_time']:.0f}s; float certificate "
          f"|A^T y|/|y| {cert['cert_ATy']:.1e}, b.y/|y| {cert['cert_by']:+.2e}, min dual eig {cert['cert_min_eig']:.1e}", flush=True)
    excluded = ('Infeasible' in out['status'] and cert['cert_by'] < -1e-6 and cert['cert_ATy'] < 1e-6
                and cert['cert_min_eig'] > -1e-8)
    rep = verify_farkas_exact(S, data, out['y']) if (excluded and not a.no_exact) else None
    try:
        rev = subprocess.check_output(['git', '-C', ROOT, 'rev-parse', '--short', 'HEAD'], text=True).strip()
    except Exception:
        rev = None
    rec = dict(n=a.n, p=a.p, m=a.m, k=3 * a.m, seed=a.seed, couplings=C.astype(int).tolist(), solver=a.solver, eps=a.eps, scaled=a.scale,
               status=out['status'], solve_s=round(out['solve_time'], 1), n_vars=out['n'], n_eq=out['n_eq'],
               float_certificate=cert, excluded=bool(excluded), exact=rep, total_s=round(time.time() - t0, 1), f=a.m / (a.p * a.n * a.n),
               peak_rss_gb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 2), git=rev,
               date=time.strftime('%Y-%m-%d %H:%M'))
    if a.save_cert:
        np.savez_compressed(a.save_cert, y=out['y'], row_map=data['row_map'], col_map=data['col_map'], rn=data['rn'],
                            dims=np.array(data['dims']))
    if a.out:
        with open(a.out, 'a') as f:
            f.write(json.dumps(rec) + '\n')


if __name__ == '__main__':
    main()
