#!/usr/bin/env python3
"""
Concentration figure for the quiver progress report (research/notes/quiver_project.md sections 4-6;
research/tex/quiver_progress_report.tex).

  (a) dimension n(k) of the gauge-singlet sector in each degree k (dual-Cauchy series, src/quiver_index.py) for
      (n,p) = (2,2) and (2,3), with the singlet BPS states: 90 at k = 12 (exact ranks) and 1680 at k = 18
      (Lanczos/dense diagonalisation of H + mu * sum_v J_v^2; none in any other degree).
  (b) lowest non-BPS singlet energy in each degree, E_min(k) = min(E_0 of the Q-pair (k, k+3), E_0 of (k-3, k)),
      from the pair edges assembled by scripts/plot_quiver_step1.edges (dense pairs, Lanczos edge files, and the
      10,024-level pair (9,12) at (2,3)); particle-hole symmetry k -> 3pn^2 - k is used for the upper half.

    python scripts/plot_quiver_concentration.py --out results/figures/quiver_concentration
"""
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
sys.path.insert(0, os.path.join(ROOT, 'src'))
from quiver_index import singlet_series          # noqa: E402
from plot_quiver_step1 import edges               # noqa: E402

COL = {2: '#0072B2', 3: '#D55E00'}
MK = {2: 's', 3: 'o'}
BPS = {2: (12, 90), 3: (18, 1680)}               # (degree, number of singlet BPS states) at n = 2


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default=os.path.join(ROOT, 'results', 'figures', 'quiver_concentration'))
    ap.add_argument('--fontsize', type=float, default=10.5)
    a = ap.parse_args()
    plt.rcParams.update({'font.size': a.fontsize})
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.6))
    n = 2
    for p in (2, 3):
        nk, _, _ = singlet_series(n, p)
        top = 3 * p * n * n
        ks = np.arange(0, top + 1, 3)
        dims = np.array([int(nk[k]) for k in ks])
        kb, hb = BPS[p]
        ax = axs[0]
        ax.plot(ks, dims, MK[p] + '-', color=COL[p], ms=5, lw=1.2, label=f'$(n,p)=(2,{p})$: singlet states $n(k)$')
        ax.plot([kb], [hb], '*', color=COL[p], ms=16, mec='k', mew=0.7, zorder=5,
                label=f'$(2,{p})$: all {hb} singlet BPS states, $k={kb}$')
        # lowest non-BPS singlet energy per degree
        E = edges(p)                                            # |q| -> E_0 of the pair with that |q|
        half = 3 * p * n * n / 2

        def pair_E0(k0):
            q = abs(k0 + 1.5 - half)
            return E.get(q)
        emin = []
        for k in ks:
            cands = [pair_E0(k0) for k0 in (k, k - 3) if 0 <= k0 and k0 + 3 <= top]
            cands = [c for c in cands if c is not None]
            emin.append(min(cands) if cands else np.nan)
        ax = axs[1]
        ax.plot(ks, emin, MK[p] + '-', color=COL[p], ms=5, lw=1.2, label=f'$(n,p)=(2,{p})$')
        ax.axvline(kb, color=COL[p], ls=':', lw=1)
        ax.text(kb + 0.4, 2.5e3 if p == 2 else 2.5e1, f'BPS degree\n$k={kb}$', color=COL[p], fontsize=a.fontsize - 3,
                va='center')
        print(f"(2,{p}): n(k) = {dict(zip(ks.tolist(), dims.tolist()))}")
        print(f"(2,{p}): E_min(k) = {dict(zip(ks.tolist(), [round(x, 5) for x in emin]))}")
    ax = axs[0]
    ax.set_yscale('log')
    ax.set_xlabel('fermion number (R-charge) $k$')
    ax.set_ylabel('number of states')
    ax.set_title('(a) gauge-singlet sector and its BPS states')
    ax.set_xticks(range(0, 37, 6))
    ax.legend(fontsize=a.fontsize - 2.5, loc='upper center', ncol=2)
    ax.set_ylim(0.5, 5e8)
    ax = axs[1]
    ax.set_yscale('log')
    ax.set_xlabel('fermion number (R-charge) $k$')
    ax.set_ylabel('lowest non-BPS singlet energy $E_{\\min}(k)$')
    ax.set_title('(b) singlet gap per degree')
    ax.set_xticks(range(0, 37, 6))
    ax.set_ylim(1e-2, 5e3)
    ax.legend(fontsize=a.fontsize - 2.5, loc='lower left')
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(f'{a.out}.{ext}', dpi=160)
    print('saved', a.out + '.{pdf,png}')


if __name__ == '__main__':
    main()
