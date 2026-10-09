#!/usr/bin/env python3
"""
Figure for the large-p Schwinger-Dyson analysis (docs/derivations.md D22; research/notes/quiver_project.md section 14):
the n = 1 quiver (tripartite N = 2 SYK, N = 3p complex fermions, J = p^2 <C^2>) against the large-N saddle.

  (a) entropy S/N vs beta J: full ED (p = 3..6, mean over seeds) and the SD solution (2 x per-Majorana N = 1 values);
  (b) energy E/(N J) vs beta J, with the Schwarzian asymptote 2 pi^2 alpha_s / (beta J)^2;
  (c) Var(H)/(N J^2) at infinite temperature vs 1/p: ED seeds, the exact disorder average (1/16)(1 + 1/p + 1/p^2),
      and the SD value 1/16;
  (d) number of BPS states / 3^(N/2): ED and sum over flavour sectors of |index|, with the index asymptotes 4/3, 2/sqrt3;
  (e) R-charge distribution of BPS states: fraction at J_R = 0 (even N; Turiaci-Witten 1/2) and fraction outside the
      window |J_R| < 3/2, from the sector indices (lines) and ED (markers);
  (f) singlet-sector edges E_0(q) N / J vs q^2 (seed 1) against E_0(q) = J q^2 / (8 qh^2 alpha_s N).

    python scripts/plot_n2syk_large_p.py --out results/figures/n2syk_large_p
    python scripts/plot_n2syk_large_p.py --layout 3x2 --figsize 10 11.6 --out research/tex/report_figs/n2syk_large_p
"""
import argparse, glob, json, math, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
from n2syk_sd import thermodynamics                                             # noqa: E402

D = os.path.join(ROOT, 'results', 'data')


def sector_table(p):
    """Exact Euler characteristics chi(s,t) of the flavour sectors of the full n = 1 model, with their window charge."""
    out = []
    for s in range(-p, p + 1):
        for t in range(-p, p + 1):
            chi, Js = 0, []
            for z in range(p + 1):
                x, y = z + s + t, z + t
                if 0 <= x <= p and 0 <= y <= p:
                    chi += (-1) ** (x + y + z) * math.comb(p, x) * math.comb(p, y) * math.comb(p, z)
                    Js.append(x + y + z - 1.5 * p)
            w = [j for j in Js if abs(j) < 1.5]
            out.append((abs(chi), w[0] if w else None))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default=os.path.join(ROOT, 'results', 'figures', 'n2syk_large_p'))
    ap.add_argument('--figsize', type=float, nargs=2, default=[14.0, 8.4])
    ap.add_argument('--fontsize', type=float, default=11.5)
    ap.add_argument('--layout', choices=['2x3', '3x2'], default='2x3', help='panel grid (3x2 for the report page width)')
    a = ap.parse_args()
    plt.rcParams.update({'font.size': a.fontsize, 'axes.labelsize': a.fontsize, 'legend.fontsize': a.fontsize - 1.5,
                         'xtick.labelsize': a.fontsize - 1, 'ytick.labelsize': a.fontsize - 1})
    sd = json.load(open(os.path.join(D, 'n2syk_thermo_qh3.json')))
    aS2 = sd['alpha_s_N2']
    f_sd = thermodynamics([r['betaJ'] for r in sd['rows']], [r['E'] for r in sd['rows']], 3, sd['alpha_s_N1'])
    ed = [json.loads(l) for l in open(os.path.join(D, 'quiver_n1_full_thermo.jsonl'))]
    mom = [json.loads(l) for l in open(os.path.join(D, 'quiver_n1_full_moments.jsonl'))]
    ps = sorted({r['p'] for r in ed if r['p'] >= 3})
    cmap = plt.get_cmap('viridis')
    col = {p: cmap(0.08 + 0.8 * (p - 3) / 7) for p in range(3, 11)}
    col[2] = cmap(0.0)
    rows, cols = (2, 3) if a.layout == '2x3' else (3, 2)
    fig, ax = plt.subplots(rows, cols, figsize=a.figsize)
    ax = ax.ravel()

    # (a), (b): thermodynamics
    bJ = np.geomspace(0.05, 200, 400)
    lz, E1, S1 = f_sd(bJ)
    for p in ps:
        rr = [r for r in ed if r['p'] == p]
        b = np.array([x['betaJ'] for x in rr[0]['rows']])
        S = np.mean([[x['S_per_N'] for x in r['rows']] for r in rr], axis=0)
        E = np.mean([[x['E_per_N'] / r['J'] for x in r['rows']] for r in rr], axis=0)
        lab = f'ED $p={p}$, $N={3*p}$'
        ax[0].plot(b, S, color=col[p], lw=1.8, label=lab)
        ax[1].plot(b, E, color=col[p], lw=1.8, label=lab)
    ax[0].plot(bJ, 2 * S1, 'k-', lw=2.2, label='SD, $N\\to\\infty$')
    ax[0].axhline(np.log(np.sqrt(3)), color='0.45', ls='--', lw=1.2)
    ax[0].text(0.06, np.log(np.sqrt(3)) - 0.012, '$S_0/N=\\ln\\sqrt{3}$', color='0.3', va='top')
    ax[0].axhline(np.log(2), color='0.45', ls=':', lw=1.2)
    ax[0].text(0.06, np.log(2) - 0.004, '$\\ln 2$', color='0.3', va='top')
    ax[0].set_xscale('log'); ax[0].set_xlabel('$\\beta J$'); ax[0].set_ylabel('entropy $S/N$')
    ax[0].set_ylim(0.52, 0.70); ax[0].legend(loc='lower left', bbox_to_anchor=(0.0, 0.14), frameon=False)
    ax[0].set_title('(a) entropy per complex fermion', loc='left')
    ax[1].plot(bJ, 2 * E1, 'k-', lw=2.2, label='SD, $N\\to\\infty$')
    ax[1].plot(bJ[bJ > 5], 2 * np.pi ** 2 * aS2 / bJ[bJ > 5] ** 2, color='0.45', ls='--', lw=1.2,
               label=f'Schwarzian $2\\pi^2\\alpha_s/(\\beta J)^2$, $\\alpha_s={aS2:.5f}$')
    ax[1].set_xscale('log'); ax[1].set_yscale('log'); ax[1].set_ylim(1e-6, 0.12)
    ax[1].set_xlabel('$\\beta J$'); ax[1].set_ylabel('energy $E/(NJ)$')
    ax[1].legend(loc='lower left', frameon=False)
    ax[1].set_title('(b) energy per complex fermion', loc='left')

    # (c): variance at infinite temperature
    x = np.linspace(0, 0.55, 100)
    ax[2].plot(x, (1 + x + x ** 2) / 16, color='0.25', lw=1.6, label='exact disorder average $\\frac{1}{16}(1+\\frac{1}{p}+\\frac{1}{p^2})$')
    for p in sorted({r['p'] for r in mom} | {r['p'] for r in ed}):
        # p <= 5: 20 seeds in nominal units (sigma = 1/p, J = 1), whose mean estimates the exact average directly;
        # p = 6: the few full-thermo seeds, in units of their realised J (bias ~ 2/p^3 < 1%)
        v = [r['var_H_per_N'] for r in mom if r['p'] == p] or [r['var_H_per_N'] / r['J'] ** 2 for r in ed if r['p'] == p]
        jit = 0.004 * np.linspace(-1, 1, len(v))
        ax[2].plot(1 / p + jit, v, 'o', ms=3.5, color=col.get(p, 'C0'), alpha=0.45, mec='none')
        m, e = np.mean(v), (np.std(v, ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0)
        ax[2].errorbar(1 / p, m, yerr=e, fmt='s', ms=7, color=col.get(p, 'C0'), mec='white', mew=1.2, capsize=3, zorder=5)
    ax[2].plot(0, 1 / 16, '*', ms=15, color='k', mec='white', zorder=6, label='SD, $N\\to\\infty$: $1/16$')
    ax[2].set_xlabel('$1/p$'); ax[2].set_ylabel('$\\mathrm{Var}(H)/(N J^2)$ at $T=\\infty$')
    ax[2].set_xlim(-0.02, 0.55); ax[2].set_ylim(0, 0.215)
    ax[2].legend(loc='upper left', frameon=False)
    ax[2].set_title('(c) second moment: $1/p$ approach', loc='left')

    # (d), (e): BPS counts and R-charge distribution from the sector indices
    pp = list(range(2, 25))
    tot, f0, fout = {}, {}, {}
    for p in pp:
        tab = sector_table(p)
        T = sum(c for c, _ in tab)
        tot[p] = T
        f0[p] = sum(c for c, w in tab if w == 0) / T
        fout[p] = sum(c for c, w in tab if w is None) / T
    for par, mk, lab, asym in ((0, 'o', 'even $N$', 4 / 3), (1, 'D', 'odd $N$', 2 / np.sqrt(3))):
        sel = [p for p in pp if (3 * p) % 2 == par]
        ax[3].plot([3 * p for p in sel], [tot[p] / 3 ** (1.5 * p) for p in sel], mk, mfc='none', color='0.35', ms=6,
                   label=f'$\\Sigma|\\chi|$ over sectors, {lab}')
        ax[3].axhline(asym, color='0.55', ls='--', lw=1.1)
    ax[3].text(73, 4 / 3 - 0.006, '$4/3$', ha='right', va='top', color='0.3')
    ax[3].text(73, 2 / np.sqrt(3) - 0.008, '$2/\\sqrt{3}$', ha='right', va='top', color='0.3')
    nb = {}
    for r in ed + mom:
        nb.setdefault(r['p'], set()).add(r['n_zero'])
    for p, v in sorted(nb.items()):
        for n in v:
            ax[3].plot(3 * p, n / 3 ** (1.5 * p), 'o' if p % 2 == 0 else 'D', ms=6.5, color=col.get(p, col[2]),
                       mec='k', mew=0.6, zorder=5)
    ax[3].plot([], [], 'o', color=col[4], mec='k', mew=0.6, label='ED (every seed)')
    ax[3].annotate('$p=3$: 2+2 singlet states\nbeyond $|\\chi|$ (index 0)', xy=(9, 172 / 3 ** 4.5), xytext=(20, 1.22),
                   arrowprops=dict(arrowstyle='->', color='0.3', lw=0.9), color='0.25', fontsize=a.fontsize - 2,
                   bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none'))
    ax[3].set_xlabel('$N=3p$'); ax[3].set_ylabel('#BPS $/\\,3^{N/2}$')
    ax[3].set_ylim(1.1, 1.47); ax[3].legend(loc='upper right', frameon=False)
    ax[3].set_title('(d) BPS count vs index asymptotics', loc='left')

    ev = [p for p in pp if p % 2 == 0]; od = [p for p in pp if p % 2 == 1]
    c0, c1 = '0.15', 'tab:orange'
    nz = lambda v: v if v > 0 else np.nan                 # zero fractions (p = 2) cannot be drawn on the log axis
    ax[4].plot([3 * p for p in ev], [f0[p] for p in ev], '-', color=c0, lw=1.6, label='fraction at $J_R=0$ (even $N$), sector indices')
    ax[4].plot([3 * p for p in ev], [nz(fout[p]) for p in ev], '-', color=c1, lw=1.6, label='fraction outside $|J_R|<3/2$, even $N$')
    ax[4].plot([3 * p for p in od], [fout[p] for p in od], '--', color=c1, lw=1.6, label='fraction outside $|J_R|<3/2$, odd $N$')
    ax[4].axhline(0.5, color='0.55', ls=':', lw=1.1)
    ax[4].text(70, 0.53, 'Turiaci–Witten $\\cos(\\pi k/3)$: $1/2$', ha='right', va='bottom', color='0.3', fontsize=a.fontsize - 1.5)
    for r in ed:
        byJ = r['bps_by_JR']; n = r['n_zero']
        out = sum(v for k, v in byJ.items() if abs(float(k)) >= 1.5) / n
        if r['p'] % 2 == 0:
            ax[4].plot(r['N'], byJ.get('+0.0', 0) / n, 'o', ms=6, color=c0, mec='white', mew=0.8, zorder=5)
        if out > 0:
            ax[4].plot(r['N'], out, 'o' if r['p'] % 2 == 0 else 'D', ms=6, color=c1, mec='k', mew=0.6, zorder=5)
    ax[4].plot([], [], 'o', color='0.6', mec='k', mew=0.6, label='ED')
    ax[4].set_yscale('log'); ax[4].set_ylim(1e-5, 1.5)
    ax[4].set_xlabel('$N=3p$'); ax[4].set_ylabel('fraction of BPS states')
    ax[4].legend(loc='lower left', frameon=False, fontsize=a.fontsize - 2.5)
    ax[4].set_title('(e) R-charge window of BPS states', loc='left')

    # (f): singlet-sector edges vs the Schwarzian threshold
    edges = {}
    for fpath in sorted(glob.glob(os.path.join(D, 'quiver_n1_spectrum_gauss_seed1_*.json'))):
        d = json.load(open(fpath))
        res = d['results'] if 'results' in d else {str(d['p']): d}
        for p, r in res.items():
            for pr in r.get('pairs', {}).values():
                if abs(pr['q']) > 0:
                    edges.setdefault(int(p), {})[abs(pr['q'])] = pr['E0']
    qq = np.linspace(1.2, 11.5, 50) ** 2
    slope = 1 / (8 * 9 * aS2)
    ax[5].plot(qq, slope * qq, 'k-', lw=2.0, label=f'saddle: $E_0={slope:.3f}\\,Jq^2/N$')
    for p in sorted(edges):
        if p < 4:
            continue
        C = np.random.default_rng(1).standard_normal((p, p, p)) / p
        J = p ** 2 * (C ** 2).mean(); N = 3 * p
        q = np.array(sorted(edges[p])); e = np.array([edges[p][x] for x in q])
        ax[5].plot(q ** 2, e * N / J, 'o-' if p % 2 == 0 else 'D--', color=col[p], ms=6, mec='k', mew=0.5, lw=1.2,
                   label=f'$p={p}$')
    ax[5].set_xscale('log'); ax[5].set_yscale('log')
    ax[5].set_xlabel('$q^2$ (multiplet R-charge)'); ax[5].set_ylabel('$E_0\\,N/J$ (singlet sector)')
    ax[5].legend(loc='upper left', frameon=False, ncol=1, fontsize=a.fontsize - 2.5)
    ax[5].set_title('(f) near-BPS edges vs threshold $E_0(q)$', loc='left')
    for x in ax:
        x.grid(alpha=0.25, lw=0.6)
    fig.tight_layout()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(f'{a.out}.{ext}', dpi=170)
    print('saved', a.out + '.{pdf,png}')


if __name__ == '__main__':
    main()
