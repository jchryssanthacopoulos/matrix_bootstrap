#!/usr/bin/env python3
"""
Degree-resolved gauge-singlet counts n(k) of the U(n)^3 quiver (research/notes/quiver_project.md section 9), from the
flavour transfer matrix with a fugacity per fermion (src/quiver_index.singlet_series_recursive), with cross-checks
against the stored totals and indices (results/data/quiver_singlet_index_dixon.json) and palindromy.

Used to size the next exact computations: a singlet-adapted basis needs dim = n(k) per degree.

    python scripts/quiver_singlet_series.py --pairs 2,2 2,3 2,4 2,5 3,2 3,3 4,2
"""
import argparse, json, os, sys, time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
from quiver_index import singlet_series_recursive          # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--pairs', nargs='+', default=['2,2', '2,3', '2,4', '2,5', '3,2', '3,3', '4,2'])
    ap.add_argument('--out', default=os.path.join(ROOT, 'results', 'data', 'quiver_singlet_series.json'))
    a = ap.parse_args()
    stored = json.load(open(os.path.join(ROOT, 'results', 'data', 'quiver_singlet_index_dixon.json')))
    rec = {}
    for pr in a.pairs:
        n, p = (int(x) for x in pr.split(','))
        t0 = time.time()
        c, resid = singlet_series_recursive(n, p)
        total = sum(c)
        index = sum((-1) ** (k // 3) * c[k] for k in range(0, len(c), 3))
        st = stored.get(f'{n},{p}', {})
        checks = dict(palindromic=c == c[::-1], only_k_multiple_of_3=all(x == 0 for k, x in enumerate(c) if k % 3),
                      total_matches=(st.get('singlet_total') is None or abs(total - st['singlet_total']) <= 1e-9 * total),
                      index_matches=(st.get('index') is None or index == st['index']))
        rec[pr] = dict(n=n, p=p, counts=c, total=total, index=index, fft_rounding_residual=resid, checks=checks,
                       half_filling=3 * p * n * n / 2, time_s=round(time.time() - t0, 2))
        print(f"({n},{p}): total {total:.4g}, index {index}, checks {checks}, residual {resid:.1e}; "
              f"n(k) at half filling {c[3 * p * n * n // 2] if (3 * p * n * n) % 2 == 0 else 'n/a'}", flush=True)
    json.dump(rec, open(a.out, 'w'), indent=1)
    print('saved', os.path.relpath(a.out, ROOT))


if __name__ == '__main__':
    main()
