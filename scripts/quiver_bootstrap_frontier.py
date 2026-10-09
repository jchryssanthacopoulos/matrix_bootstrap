#!/usr/bin/env python3
"""
Locate the level-2 BPS-exclusion frontier of the quiver singlet bootstrap at (n, p) (research/notes/quiver_project.md
section 12; docs/derivations.md D21).

  1. coarse scan with pure feasibility solves at edge occupations m chosen from the filling fractions --coarse-f
     (f = m / p n^2): excluded = 'PrimalInfeasible' with a Farkas vector of b.y/|y| < -1e-6 and |A^T y|/|y| < 1e-6;
     optionally (--margins) also the margin t* at the coarse points;
  2. bisection of the crossing with feasibility solves;
  3. exact verification (src/quiver_certificate.py) of the certificate at the last excluded m and at every --targets m.
Every solve is one JSON record in --out.  The SDP has the same size at every n (only coefficients change).  All solves
use the column scaling x_j = n^{w_j} y_j (QuiverSDP.assemble(scale=True)): without it the interior-point method
needs hours per solve at n >= 4 (monomial values span ~ n^8); with it, ~9 iterations and ~1 min.

    python scripts/quiver_bootstrap_frontier.py --n 4 --p 3 --coarse-f 0.125,0.21,0.27,0.31,0.375 --targets 13,14 \
        --out results/data/quiver_bootstrap_frontier.jsonl
"""
import argparse, json, os, resource, subprocess, sys, time
import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
from quiver_bootstrap import QuiverSDP                                          # noqa: E402
from quiver_certificate import verify_farkas_exact                              # noqa: E402


def git_rev():
    try:
        return subprocess.check_output(['git', '-C', ROOT, 'rev-parse', '--short', 'HEAD'], text=True).strip()
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--p', type=int, required=True)
    ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('--coarse-f', default='0.125,0.21,0.27,0.31,0.375')
    ap.add_argument('--targets', default='', help='edge occupations m to certify in any case (e.g. fortuity degrees)')
    ap.add_argument('--no-exact', action='store_true', help='skip the exact certificate verification')
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--eps', type=float, default=1e-8)
    ap.add_argument('--margins', action='store_true', help='also compute the margin t* at the coarse points (slow)')
    ap.add_argument('--time-limit', type=float, default=1800.0, help='seconds per solve')
    ap.add_argument('--subprocess', action='store_true',
                    help='run every point in a fresh process (scripts/quiver_verify_certificate.py); needed at p >= 4, '
                         'where a second in-process build forks workers from a large parent and the summed RSS exceeds the guard')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    n, p = a.n, a.p
    pn2 = p * n * n
    C = np.random.default_rng(a.seed).integers(1, 6, size=(p, p, p)).astype(float)
    ms = sorted({max(1, min(int(round(float(f) * pn2)), (pn2 - 1) // 2)) for f in a.coarse_f.split(',')})
    targets = [int(x) for x in a.targets.split(',') if x.strip()]

    def record(m, mode, r, S, extra=None):
        rec = dict(n=n, p=p, m=m, k=3 * m, f=m / pn2, seed=a.seed, mode=mode, status=r['status'], solve_s=round(r['solve_time'], 1),
                   n_vars=r['n'], n_eq=r['n_eq'], n_monomials=len(S.monos), level=S.level, git=git_rev(),
                   peak_rss_gb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 2), date=time.strftime('%Y-%m-%d %H:%M'))
        rec.update(extra or {})
        with open(a.out, 'a') as fh:
            fh.write(json.dumps(rec) + '\n')
        return rec

    def margin(m):
        S = QuiverSDP(C, n, m, bps=True, workers=a.workers)
        r = S.margin(eps=a.eps, scale=True, time_limit=a.time_limit)
        rec = record(m, 'margin', r, S, dict(margin=r['margin']))
        print(f"(n,p)=({n},{p}) k={3*m:3d} f={m/pn2:.3f} margin   t* = {r['margin']:+.4e}  [{r['status']}, {r['solve_time']:.0f}s]", flush=True)
        return r['margin'] < -1e-4 and 'olved' in r['status'], rec

    def feasibility_subprocess(m, exact):
        cmd = [sys.executable, os.path.join(ROOT, 'scripts', 'quiver_verify_certificate.py'), '--n', str(n), '--p', str(p),
               '--m', str(m), '--seed', str(a.seed), '--scale', '--workers', str(a.workers), '--eps', str(a.eps),
               '--time-limit', str(a.time_limit), '--out', a.out] + ([] if exact else ['--no-exact'])
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        rec = [json.loads(l) for l in open(a.out)][-1]
        assert rec['m'] == m and rec['n'] == n and rec['p'] == p
        ex, cert = rec['excluded'], rec['float_certificate']
        exr = rec.get('exact')
        print(f"(n,p)=({n},{p}) k={3*m:3d} f={m/pn2:.3f} feasib.  {rec['status']}: b.y/|y| {cert['cert_by']:+.2e}, |A^Ty|/|y| {cert['cert_ATy']:.1e}"
              + (f", exact {'VERIFIED' if exr['verified'] else 'NOT verified'} (ratio {exr['ratio_bound_over_b']:.2e})" if exr else '')
              + f"  [{rec['solve_s']:.0f}s solve, {rec['total_s']:.0f}s total, {rec['peak_rss_gb']} GB]", flush=True)
        return ex, (exr['verified'] if exr else None)

    def feasibility(m, exact):
        if a.subprocess:
            return feasibility_subprocess(m, exact)
        S = QuiverSDP(C, n, m, bps=True, workers=a.workers)
        data = S.assemble('feasibility', 'clarabel', prune=True, scale=True)
        r = S._solve(data, 'clarabel', a.eps, 100000, False, time_limit=a.time_limit)
        cert = S.check_certificate(data, r['y'], 'clarabel')
        excluded = 'Infeasible' in r['status'] and cert['cert_by'] < -1e-6 and cert['cert_ATy'] < 1e-6 and cert['cert_min_eig'] > -1e-8
        extra = dict(float_certificate=cert, excluded=bool(excluded), iterations=r.get('iterations'), scaled=True)
        if excluded and exact:
            extra['exact'] = verify_farkas_exact(S, data, r['y'], verbose=False)
        record(m, 'feasibility', r, S, extra)
        ex = extra.get('exact')
        print(f"(n,p)=({n},{p}) k={3*m:3d} f={m/pn2:.3f} feasib.  {r['status']}: b.y/|y| {cert['cert_by']:+.2e}, |A^Ty|/|y| {cert['cert_ATy']:.1e}"
              + (f", exact {'VERIFIED' if ex['verified'] else 'NOT verified'} (ratio {ex['ratio_bound_over_b']:.2e})" if ex else '')
              + f"  [{r['solve_time']:.0f}s]", flush=True)
        return excluded, (ex['verified'] if ex else None)

    t0 = time.time()
    status = {}
    for m in ms:
        status[m] = feasibility(m, exact=False)[0]
        if a.margins:
            margin(m)
    excl = [m for m in ms if status[m]]
    lo = max(excl) if excl else 0
    hi_c = [m for m in ms if m > lo and not status[m]]
    hi = min(hi_c) if hi_c else (pn2 + 1) // 2
    if any(not status[m] for m in ms if m < lo):
        print('WARNING: non-monotone coarse scan', {m: status[m] for m in ms}, flush=True)
    certified = {}
    while hi - lo > 1:
        mid = (lo + hi) // 2
        ex, ver = feasibility(mid, exact=False)
        if ex:
            lo = mid
        else:
            hi = mid
    if lo > 0 and not a.no_exact:
        certified[lo] = feasibility(lo, exact=True)
    for m in targets:
        if m not in certified:
            certified[m] = feasibility(m, exact=not a.no_exact)
    print(f"(n,p)=({n},{p}): frontier: excluded through k={3*lo} (f={lo/pn2:.3f}), not excluded at k={3*hi} (f={hi/pn2:.3f}); "
          f"certified {certified}  [{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    main()
