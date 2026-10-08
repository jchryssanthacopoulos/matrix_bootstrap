#!/usr/bin/env python3
"""
Test of the observed pattern |I_0(n, p)| = |Dixon(n p)| = (3np/2)! / ((np/2)!)^3 for the singlet index of the U(n)^3
fermionic quiver (research/notes/quiver_project.md section 7).  The index is computed exactly as
I_0 = n(t=-1) = Tr[R(-1)^3] with the edge transfer matrix over U(n) irreps (src/quiver_index.singlet_index_exact).

    scripts/run_guarded.sh 4 log python3 scripts/quiver_index_dixon.py --cases 2,2 2,3 3,2 2,4 4,2 3,3
"""
import argparse, json, os, resource, sys, time
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'scripts')); sys.path.insert(0, os.path.join(ROOT, 'src'))
from quiver_index import singlet_index_exact, dixon                         # noqa: E402


def peak_gb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1e9 if sys.platform == 'darwin' else r / 1e6


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cases', nargs='+', required=True, help='n,p pairs')
    ap.add_argument('--out', default=os.path.join(ROOT, 'results', 'data', 'quiver_singlet_index_dixon.json'))
    a = ap.parse_args()
    rec = json.load(open(a.out)) if os.path.exists(a.out) else {}
    for cs in a.cases:
        n, p = (int(x) for x in cs.split(','))
        t0 = time.time()
        r = singlet_index_exact(n, p, verbose=False)
        d = dixon(n * p)
        r.update(n=n, p=p, dixon_np=d, matches=abs(r['index']) == abs(d), time_s=round(time.time() - t0, 1),
                 peak_rss_gb=round(peak_gb(), 2))
        rec[f'{n},{p}'] = r
        print(f"(n,p)=({n},{p}): I_0 = {r['index']}, Dixon(np={n*p}) = {d}  -> {'MATCH' if r['matches'] else 'differs'}"
              f"  [{r['n_tuples']} tuples, {r['n_nu']} irreps, LR rounding {r['lr_rounding_error']:.1e}, "
              f"max|R| {r['max_abs_R']}, {r['time_s']}s, {r['peak_rss_gb']} GB]", flush=True)
        json.dump(rec, open(a.out, 'w'), indent=1)
    print('saved', os.path.relpath(a.out, ROOT))


if __name__ == '__main__':
    main()
