#!/usr/bin/env python3
"""
Figure for the large-n saddle of the quiver singlet index and singlet count (research/notes/quiver_project.md section 8;
docs/derivations.md D18).  Input: results/data/quiver_index_saddle.json (scripts/quiver_index_saddle.py).

  (a) ln|I_0(n,p)| / n^2 against n: exact indices (filled), the corrected Laplace estimate about the discrete saddle
      (Gaussian fluctuations minus the universal beta = 2 crystal excess 3n ln(e/sqrt(2pi))), and the continuum
      asymptote F*(p) (dashed).  p = 2: the exact closed form (3n)!/(n!)^3 (F* = 0).
  (b) the same for the total singlet count, asymptote G*(p).
  (c) leading coefficients per fermion mode, F*(p)/(3p) and G*(p)/(3p), from the continuum solution, with the
      threshold forms (p/2 - 1) 39 zeta(3)/(2pi^2) and (p - 1) 21 zeta(3)/(2pi^2) and the large-p (semicircle) forms;
      reference lines ln 2 (all states) and (1/2) ln 3.
  (d) equilibrium densities at p = 3: the three index arcs (one per node) and the coincident count density, with
      histograms of the discrete maximum at n = 512.

    python scripts/plot_quiver_index_saddle.py --out results/figures/quiver_index_saddle
"""
import argparse, json, math, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import zeta

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
from quiver_saddle import Continuum, newton_reduced, expand          # noqa: E402

COL = {2: '#0072B2', 3: '#D55E00', 4: '#009E73', 5: '#CC79A7', 6: '#E69F00'}
MK = {2: 's', 3: 'o', 4: '^', 5: 'D', 6: 'v'}


def per_n2_panel(ax, rec, kind, ps, fs):
    for p in ps:
        rows = rec['discrete'].get(f'{kind},{p}', [])
        if not rows:
            continue
        cont = rec['continuum'][kind].get(str(float(p)))
        Fs = cont['F'] if cont else (0.0 if kind == 'index' and p == 2 else None)
        lap = [(r['n'], r['laplace_corr'] / r['n'] ** 2) for r in rows
               if r.get('laplace_corr') is not None and np.isfinite(r['laplace_corr'])]
        if kind == 'index' and p == 2:
            nn = np.unique(np.round(np.logspace(0, np.log10(256), 60)).astype(int))
            ax.plot(nn, [(math.lgamma(3 * n + 1) - 3 * math.lgamma(n + 1)) / n ** 2 for n in nn], '-', color=COL[p], lw=1.4,
                    label=f'$p={p}$: exact $(3n)!/(n!)^3$')
        elif lap:
            x, y = zip(*lap)
            ax.plot(x, y, '-', color=COL[p], lw=1.4, label=f'$p={p}$: saddle + fluctuations')
        ex = [(r['n'], r['exact'] / r['n'] ** 2) for r in rows if r.get('exact') is not None]
        if ex:
            x, y = zip(*ex)
            ax.plot(x, y, MK[p], color=COL[p], ms=7, mec='k', mew=0.6, zorder=4,
                    label=f'$p={p}$: exact' if not (kind == 'index' and p == 2) else None)
        if Fs is not None:
            ax.axhline(Fs, color=COL[p], ls='--', lw=1, alpha=0.8)
            ax.text(330, Fs + 0.06, f'{"$F^*$" if kind == "index" else "$G^*$"}$={Fs:.4f}$', color=COL[p], va='bottom',
                    fontsize=fs - 3, bbox=dict(fc='white', ec='none', pad=0.3))
    ax.set_xscale('log')
    ax.set_xlim(0.9, 1100)
    ax.set_xlabel('rank $n$ of each gauge group')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--data', default=os.path.join(ROOT, 'results', 'data', 'quiver_index_saddle.json'))
    ap.add_argument('--out', default=os.path.join(ROOT, 'results', 'figures', 'quiver_index_saddle'))
    ap.add_argument('--fontsize', type=float, default=10.5)
    ap.add_argument('--n-hist', type=int, default=512)
    a = ap.parse_args()
    rec = json.load(open(a.data))
    fs = a.fontsize
    plt.rcParams.update({'font.size': fs})
    fig, axs2 = plt.subplots(2, 2, figsize=(13, 10))
    axs = axs2.ravel()

    # (a) index
    ax = axs[0]
    per_n2_panel(ax, rec, 'index', (2, 3, 4, 5, 6), fs)
    ax.set_ylabel(r'$\ln|I_0(n,p)|\,/\,n^2$')
    ax.set_ylim(0, 8.5)
    ax.set_title('(a) singlet index: exact vs large-$n$ saddle')
    ax.legend(fontsize=fs - 3, loc='upper right', ncol=2)

    # (b) count
    ax = axs[1]
    per_n2_panel(ax, rec, 'count', (2, 3, 4, 5, 6), fs)
    ax.set_ylabel(r'$\ln(\#\,\mathrm{singlets})\,/\,n^2$')
    ax.set_ylim(0, 12)
    ax.set_title('(b) total singlet count: exact vs large-$n$ saddle')
    ax.legend(fontsize=fs - 3, loc='upper right', ncol=2)

    # (c) coefficients per mode
    ax = axs[2]
    C_arc, C_half = rec['C_arc'], rec['C_half']
    for kind, col, lab, thr_slope, p0 in (('index', '#D55E00', r'index $F^*(p)/3p$', C_arc / 2, 2.0),
                                          ('count', '#0072B2', r'singlets $G^*(p)/3p$', C_half, 1.0)):
        cs = rec['continuum'][kind]
        ps = np.array(sorted(float(k) for k in cs))
        F = np.array([cs[str(p)]['F'] for p in ps])
        ax.plot(ps, F / (3 * ps), '-', color=col, lw=1.8, label=lab + ' (continuum)')
        integ = [p for p in ps if abs(p - round(p)) < 1e-9]
        ax.plot(integ, [cs[str(p)]['F'] / (3 * p) for p in integ], 'o', color=col, ms=5)
        pp = np.linspace(p0, p0 + 1.2, 40)
        if kind == 'index':
            thr = (pp / 2 - 1) * C_arc
            pl = np.logspace(np.log10(4), np.log10(128), 60)
            large = 1.5 * pl * np.log(3) - 1.5 * np.log(2 * pl / 3) - 2.25
        else:
            thr = (pp - 1) * C_half
            pl = np.logspace(np.log10(3), np.log10(128), 60)
            large = 3 * pl * np.log(2) - 1.5 * np.log(pl / 2) - 2.25
        ax.plot(pp, thr / (3 * pp), ':', color=col, lw=1.4)
        ax.plot(pl, large / (3 * pl), '--', color=col, lw=1, alpha=0.8)
    ax.plot([], [], ':', color='0.3', label='threshold form (lower bound)')
    ax.plot([], [], '--', color='0.3', label='large-$p$ semicircle form')
    ax.axhline(np.log(2), color='0.4', lw=0.8, ls='-.')
    ax.text(1.05, np.log(2) + 0.012, r'$\ln 2$ (all $2^{3pn^2}$ states)', fontsize=fs - 3, color='0.3')
    ax.axhline(0.5 * np.log(3), color='0.6', lw=0.8, ls='-.')
    ax.text(1.05, 0.5 * np.log(3) + 0.012, r'$\frac{1}{2}\ln 3$ (index, $p\to\infty$)', fontsize=fs - 3, color='0.4')
    ax.axvline(2, color='0.8', lw=0.8)
    ax.set_xscale('log'); ax.set_xlim(1, 130); ax.set_ylim(0, 0.75)
    ax.set_xticks([1, 2, 3, 4, 6, 10, 16, 32, 64, 128]); ax.set_xticklabels(['1', '2', '3', '4', '6', '10', '16', '32', '64', '128'])
    ax.minorticks_off()
    ax.set_xlabel('flavours per edge $p$')
    ax.set_ylabel('coefficient of $n^2$ per fermion mode')
    ax.set_title('(c) leading large-$n$ coefficients')
    ax.legend(fontsize=fs - 3, loc='lower right')

    # (d) densities at p = 3 with the discrete maximum
    ax = axs[3]
    th = np.linspace(-np.pi, np.pi, 2000)
    ci = Continuum(3, 'index', 60, 700); ci.optimise()
    cc = Continuum(3, 'count', 60, 700); cc.optimise()
    n = a.n_hist
    Li, thi, _, _ = newton_reduced(ci.quantiles(n), 3, 'index')
    Lc, thc, _, _ = newton_reduced(cc.quantiles(n), 3, 'count')
    full_i = expand(thi, 'index')
    bins = np.linspace(-np.pi, np.pi, 73)
    node_cols = ['#D55E00', '#E69F00', '#CC79A7']
    for v in range(3):
        ang = np.mod(full_i[v * n:(v + 1) * n] + np.pi, 2 * np.pi) - np.pi
        ax.hist(ang, bins=bins, density=False, weights=np.full(n, 1 / (n * (bins[1] - bins[0]))), color=node_cols[v],
                alpha=0.35)
        shifted = np.mod(th - 2 * np.pi * v / 3 + np.pi, 2 * np.pi) - np.pi
        ax.plot(th, ci.density(shifted), '-', color=node_cols[v], lw=1.6,
                label=f'index saddle, node {v}' + (f' (arc half-width {ci.a:.4f})' if v == 0 else ''))
    ax.hist(thc, bins=bins, weights=np.full(n, 1 / (n * (bins[1] - bins[0]))), color='#0072B2', alpha=0.2)
    ax.plot(th, cc.density(th), '-', color='#0072B2', lw=1.6, label=f'count saddle, all nodes (half-width {cc.a:.4f})')
    ax.set_xlim(-np.pi, np.pi)
    ax.set_xticks([-np.pi, -2 * np.pi / 3, -np.pi / 3, 0, np.pi / 3, 2 * np.pi / 3, np.pi])
    ax.set_xticklabels([r'$-\pi$', r'$-\frac{2\pi}{3}$', r'$-\frac{\pi}{3}$', '0', r'$\frac{\pi}{3}$', r'$\frac{2\pi}{3}$', r'$\pi$'])
    ax.set_xlabel(r'eigenvalue phase $\theta$ of $U_v$')
    ax.set_ylabel(r'eigenvalue density $\rho_v(\theta)$')
    ax.set_title(f'(d) $p=3$ saddles: continuum (lines), discrete $n={n}$ (bars)')
    ax.legend(fontsize=fs - 3, loc='upper left')
    ax.set_ylim(0, max(ci.density(np.linspace(-ci.a, ci.a, 400)).max(), cc.density(np.linspace(-cc.a, cc.a, 400)).max()) * 1.35)

    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(f'{a.out}.{ext}', dpi=160)
    print('saved', a.out + '.{pdf,png}')


if __name__ == '__main__':
    main()
