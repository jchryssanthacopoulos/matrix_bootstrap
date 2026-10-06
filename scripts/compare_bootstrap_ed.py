#!/usr/bin/env python3
"""
Grade trace-bootstrap lower bounds on sector ground energies against exact diagonalisation.

Reads bound records (JSON lines from scripts/run_trace_bound.py) and the exact values in
results/data/sector_ed_N{N}_k{k}_fast.json (scripts/run_sector_ed.py).  For a plain-sector bound the target is
the sector ground energy E_0(k); for a Casimir-resolved bound (rows phi((C2-c)X)=0) the target is the lowest
exact level with that Casimir among the computed levels (marked n/a if none was computed).  Prints a table and
optionally plots bound/exact against the sector.

    python scripts/compare_bootstrap_ed.py --N 3 --bounds results/data/trace_bounds_N3_level3.jsonl \
        --fig results/figures/bootstrap_vs_ed_N3
"""
import argparse
import json
import os

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')


def exact_targets(N, k):
    f = os.path.join(ROOT, 'results', 'data', f'sector_ed_N{N}_k{k}_fast.json')
    if not os.path.exists(f):
        return None
    rec = json.load(open(f))
    return dict(E0=rec['levels'][0]['E'], levels=[(g['E'], g['C2_mean']) for g in rec['levels']])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--N', type=int, required=True)
    ap.add_argument('--bounds', required=True)
    ap.add_argument('--fig', default=None)
    a = ap.parse_args()
    N, p = a.N, 3
    rows = []
    for line in open(a.bounds):
        r = json.loads(line)
        if r['N'] != N:
            continue
        ex = exact_targets(N, r['k'])
        if ex is None:
            continue
        c = r.get('casimir')
        if c is None:
            target = ex['E0']
        else:
            cand = [E for E, cc in ex['levels'] if abs(cc - c) < 1e-4]
            target = min(cand) if cand else None
        rows.append(dict(k=r['k'], c=c, level=(r['L_adj'], r['L_sing'], r['L_eom']), bound=r['bound'],
                         status=r['status'], target=target, solve_s=r['solve_s'], peak=r['peak_rss_gb']))
    rows.sort(key=lambda r: (r['k'], -1 if r['c'] is None else r['c']))
    print(f"N={N}: trace-bootstrap lower bounds against exact diagonalisation")
    print(f"{'k':>3} {'x=k/N^2':>8} {'|q|':>5} {'C2':>5} {'level':>9} {'bound':>12} {'exact':>12} {'bound/exact':>11}  status")
    for r in rows:
        q = p * N * N / 2 - r['k'] - 1.5
        ratio = r['bound'] / r['target'] if r['target'] else float('nan')
        cs = 'all' if r['c'] is None else f"{r['c']:g}"
        tg = f"{r['target']:12.6f}" if r['target'] is not None else f"{'n/a':>12}"
        print(f"{r['k']:>3} {r['k'] / N ** 2:8.3f} {q:5.1f} {cs:>5} {str(r['level']):>9} {r['bound']:12.6f} {tg} "
              f"{ratio:11.4f}  {r['status']} [{r['solve_s']}s, {r['peak']} GB]")
    if a.fig:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.rcParams.update({'font.size': 11})
        fig, ax = plt.subplots(figsize=(6.2, 3.9))
        for sel, mk, lab in ((lambda r: r['c'] is None, 'o', 'whole sector'),
                             (lambda r: r['c'] is not None, 's', 'ground-state irrep (Casimir rows)')):
            pts = [r for r in rows if sel(r) and r['target']]
            if pts:
                ax.plot([r['k'] for r in pts], [r['bound'] / r['target'] for r in pts], mk + '-', ms=6, label=lab)
        ax.axhline(1, color='k', lw=0.8, ls=':')
        ax.set_xlabel(r'sector $k=N_\Psi$')
        ax.set_ylabel('bootstrap lower bound / exact $E_0$')
        ax.set_ylim(-0.05, 1.05)
        ax.legend(frameon=False, fontsize=9)
        ax.grid(alpha=0.2)
        fig.tight_layout()
        for ext in ('pdf', 'png'):
            fig.savefig(f'{a.fig}.{ext}', dpi=200)
        print(f'saved {a.fig}.pdf/.png')


if __name__ == '__main__':
    main()
