#!/usr/bin/env python3
"""
Sector ground energies below the BPS window against the multiplet charge |q|, at N = 2 and 3, with the
super-Schwarzian law E_0(q) = q^2/(4 qhat^2) (Turiaci-Witten (3.11)) drawn through the point nearest the window.

Each sector ground state below the window is the lower member of a Q-multiplet (k, k+3) (its partner in sector
k-3 would need an energy below E_0(k-3), which never happens in the data), so its charge relative to half filling
is q = k + 3/2 - p N^2/2.  Reads results/data/sector_ed_N{N}_k{k}_fast.json (scripts/run_sector_ed.py).

    python scripts/plot_near_window.py --out results/figures/near_window_gaps
    python scripts/plot_near_window.py --no-title --size 6.2 4.0 --fontsize 10.5 --free-label 'free sectors (exact)' --out research/tex/report_figs/near_window_gaps
"""
import argparse
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                               # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
NAMES = {(0, 0): r'$\mathbf{1}$', (1, -1): r'$\mathbf{3}$', (2, -2): r'$\mathbf{5}$', (3, -3): r'$\mathbf{7}$',
         (0, 0, 0): r'$\mathbf{1}$', (1, 0, -1): r'$\mathbf{8}$', (2, -1, -1): r'$\mathbf{10}\oplus\overline{\mathbf{10}}$',
         (1, 1, -2): r'$\mathbf{10}\oplus\overline{\mathbf{10}}$', (2, 0, -2): r'$\mathbf{27}$'}


def load(N, kmax, p=3):
    rows = []
    for k in range(1, kmax + 1):
        f = os.path.join(ROOT, 'results', 'data', f'sector_ed_N{N}_k{k}_fast.json')
        if not os.path.exists(f):
            continue
        rec = json.load(open(f))
        g = rec['levels'][0]
        lam = tuple(g['candidates'][0]['lam']) if g['candidates'] else None
        q = k + 1.5 - p * N * N / 2
        rows.append(dict(k=k, q=q, E=g['E'], lam=lam, free=k <= N * N // 4))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default=os.path.join(ROOT, 'results', 'figures', 'near_window_gaps'))
    ap.add_argument('--no-title', action='store_true', help='omit the title (report version)')
    ap.add_argument('--size', type=float, nargs=2, default=(6.4, 4.6), help='width height in inches')
    ap.add_argument('--fontsize', type=float, default=12)
    ap.add_argument('--free-label', default='free sectors (D16)')
    a = ap.parse_args()
    data = {2: load(2, 4), 3: load(3, 11)}
    colors = {2: '#1f77b4', 3: '#d62728'}
    plt.rcParams.update({'font.size': a.fontsize})
    fig, ax = plt.subplots(figsize=tuple(a.size))
    for N, rows in data.items():
        aq = [abs(r['q']) for r in rows]
        E = [r['E'] for r in rows]
        inter = [r for r in rows if not r['free']]
        free = [r for r in rows if r['free']]
        ax.plot(aq, E, '-', color=colors[N], lw=1, alpha=0.5)
        ax.plot([abs(r['q']) for r in inter], [r['E'] for r in inter], 'o', color=colors[N], ms=7,
                label=f'$N={N}$, exact diagonalisation')
        if free:
            ax.plot([abs(r['q']) for r in free], [r['E'] for r in free], 'o', mfc='white', color=colors[N], ms=7,
                    label=f'$N={N}$, {a.free_label}')
        near = min(rows, key=lambda r: abs(r['q']))
        qq = [0.4 + 0.05 * i for i in range(240)]
        qq = [x for x in qq if x <= max(aq) + 0.5]
        ax.plot(qq, [near['E'] * (x / abs(near['q'])) ** 2 for x in qq], '--', color=colors[N], lw=1.2,
                label=rf'$\propto q^2$ through $|q|={abs(near["q"]):g}$')
        for r in rows:
            if abs(r['q']) <= (2.5 if N == 2 else 5) and r['lam'] in NAMES:
                ax.annotate(NAMES[r['lam']], (abs(r['q']), r['E']), textcoords='offset points',
                            xytext=(7, -11 if N == 3 else 5), fontsize=9, color=colors[N])
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(r'multiplet charge $|q|=\frac{pN^2}{2}-k-\frac{3}{2}$')
    ax.set_ylabel(r'sector ground energy $E_0(k;N)$')
    if not a.no_title:
        ax.set_title(r'Approach to the BPS window, $p=q=3$, $U(N)$', fontsize=a.fontsize)
    ax.legend(fontsize=8.5, loc='upper left', frameon=False)
    ax.grid(True, which='both', alpha=0.2)
    fig.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(f'{a.out}.{ext}', dpi=200)
    print(f'saved {a.out}.pdf/.png')


if __name__ == '__main__':
    main()
