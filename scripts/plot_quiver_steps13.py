#!/usr/bin/env python3
"""
Figure for steps 1 and 3 of the quiver plan (research/notes/quiver_project.md sections 10-11).

  (a) n = 1 (three-species N = 2 SYK, the large-p test): singlet BPS states per block against |I_0|, per p.
  (b) n = 1: near-BPS edges E_0(|q|) of the innermost pairs against p (log-log), with the Turiaci-Witten ratio
      E_0(4.5)/E_0(1.5) = 9 indicated.
  (c) n = 1: <r> of the largest pair spectra against p, with orthogonal-class and Poisson values.
  (d) next rank in the p = 2 corner: lowest non-BPS singlet energy per degree at (n,p) = (2,2) and (3,2), in units
      of the vacuum energy E_vac = n^3 sum |C|^2, against k / (3pn^2); BPS states only at half filling.

Inputs: results/data/quiver_n1_spectrum_*.json (merged), quiver_singlet_sw_n3_p2_seed3.json, and the (2,2) pair
data quiver_singlet_pairs_n2_p2.json.

    python scripts/plot_quiver_steps13.py --n1 results/data/quiver_n1_spectrum_gauss_seed1.json --out results/figures/quiver_steps13
"""
import argparse, glob, json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
COLS = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#000000']


def load_n1(paths):
    res = {}
    for pth in paths:
        d = json.load(open(pth))
        for p, r in d['results'].items():
            res.setdefault(int(p), {}).update({k: v for k, v in r.items()})
            for key in ('dims', 'bps', 'lowest', 'pairs', 'pair_spectra'):
                if key in r:
                    res[int(p)].setdefault(key, {}).update(r[key])
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--n1', nargs='+', default=sorted(glob.glob(os.path.join(ROOT, 'results', 'data', 'quiver_n1_spectrum_*.json'))))
    ap.add_argument('--sw32', default=os.path.join(ROOT, 'results', 'data', 'quiver_singlet_sw_n3_p2_seed3.json'))
    ap.add_argument('--out', default=os.path.join(ROOT, 'results', 'figures', 'quiver_steps13'))
    ap.add_argument('--fontsize', type=float, default=10.5)
    a = ap.parse_args()
    plt.rcParams.update({'font.size': a.fontsize})
    fig, axs2 = plt.subplots(2, 2, figsize=(13, 9.5))
    axs = axs2.ravel()
    n1 = load_n1(a.n1)
    ps = sorted(n1)

    # (a) BPS per block
    ax = axs[0]
    for i, p in enumerate(ps):
        r = n1[p]
        ms = sorted(int(m) for m in r['bps'])
        vals = [r['bps'][str(m)] if str(m) in r['bps'] else r['bps'][m] for m in ms]
        vals = [v if isinstance(v, (int, float)) else np.nan for v in vals]
        x = np.array(ms) - p / 2
        ax.plot(x, np.maximum(vals, 0.5), 'o-', color=COLS[i % len(COLS)], ms=5, lw=1,
                label=f'$p={p}$ ($|I_0|={abs(r["index"])}$)')
    ax.set_yscale('log')
    ax.set_ylim(0.4, 5e4)
    ax.set_xlabel('block $m-p/2$ (charge-balanced, $k=3m$)')
    ax.set_ylabel('singlet BPS states (0 plotted at 0.5)')
    ax.set_title('(a) $n=1$: BPS states per block')
    ax.legend(fontsize=a.fontsize - 3, ncol=2)

    # (b) edges vs p
    ax = axs[1]
    for qabs, mk in ((1.5, 'o'), (4.5, 's'), (7.5, '^')):
        xs, ys = [], []
        for p in ps:
            for key, v in n1[p].get('pairs', {}).items():
                if abs(abs(v['q']) - qabs) < 1e-9 and p % 2 == 0:
                    xs.append(p); ys.append(v['E0'])
                    break
        if xs:
            ax.plot(xs, ys, mk + '-', ms=6, label=f'$E_0$, $|q|={qabs}$ (even $p$)')
    pe = [p for p in ps if p % 2 == 0]
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xticks(pe); ax.set_xticklabels([str(p) for p in pe]); ax.minorticks_off()
    xx = np.array([min(pe), max(pe)], float)
    ax.plot(xx, 0.25 * (xx / 4) ** -1.0, ':', color='0.5', label='$\\propto p^{-1}$ (Schwarzian, $N_{\\rm eff}=3p$)')
    ax.set_xlabel('flavours $p$ (at $n=1$, $\\langle C^2\\rangle=1/p^2$)')
    ax.set_ylabel('lowest multiplet energy $E_0$')
    ratios = []
    for p in pe:
        e = {abs(v['q']): v['E0'] for v in n1[p].get('pairs', {}).values()}
        if 1.5 in e and 4.5 in e:
            ratios.append(f'{p}: {e[4.5] / e[1.5]:.1f}')
    ax.set_title('(b) $n=1$: near-BPS edges; $E_0(4.5)/E_0(1.5)$ = ' + ', '.join(ratios) + ' (TW: 9)',
                 fontsize=a.fontsize - 1)
    ax.legend(fontsize=a.fontsize - 3)

    # (c) <r>
    ax = axs[2]
    for p in ps:
        for key, v in n1[p].get('pair_spectra', {}).items():
            if v['n_levels'] >= 100:
                st = v['stats']
                ax.errorbar(p + 0.08 * (float(key) - p / 2), st['mean_r'], yerr=st['sem_r'], fmt='o',
                            color=COLS[ps.index(p) % len(COLS)], ms=4 + np.log10(v['n_levels']))
    ax.axhline(0.5307, color='0.3', ls='--', lw=1, label='orthogonal class (0.531)')
    ax.axhline(0.3863, color='0.6', ls=':', lw=1, label='Poisson (0.386)')
    ax.set_xlabel('flavours $p$')
    ax.set_ylabel('$\\langle r\\rangle$ of pair spectra ($\\geq$100 levels)')
    ax.set_title('(c) $n=1$: level statistics per multiplet sector')
    ax.legend(fontsize=a.fontsize - 3)

    # (d) next rank at p = 2
    ax = axs[3]
    pr = json.load(open(os.path.join(ROOT, 'results', 'data', 'quiver_singlet_pairs_n2_p2.json')))['pairs']
    C2 = np.array([[[5, 1], [1, 2]], [[1, 5], [5, 3]]], float)
    for n, col, mk, edges in ((2, COLS[0], 's', {int(k.split(',')[0]): v['E0'] for k, v in pr.items() if v['E0'] is not None}),
                              (3, COLS[1], 'o', None)):
        if n == 3:
            if not os.path.exists(a.sw32):
                continue
            d = json.load(open(a.sw32))
            edges = {int(k): v['E0'] for k, v in d['pairs'].items()}
        top = 3 * 2 * n * n
        Evac = n ** 3 * (C2 ** 2).sum()
        ks, em = [], []
        for k in range(0, top // 2 + 1, 3):
            c = [edges.get(k0) for k0 in (k, k - 3) if edges.get(k0) is not None]
            if c:
                ks.append(k); em.append(min(c) / Evac)
        ks = np.array(ks)
        ax.plot(ks / top, em, mk + '-', color=col, ms=6, label=f'$(n,p)=({n},2)$, $E_{{\\rm vac}}={Evac:.0f}$')
        ax.plot(1 - ks / top, em, mk + '-', color=col, ms=6, alpha=0.35)
    ax.axvline(0.5, color='0.5', ls=':', lw=1)
    ax.text(0.505, 0.9, 'half filling: all BPS singlets', fontsize=a.fontsize - 3, color='0.4', transform=ax.get_xaxis_transform())
    ax.set_yscale('log')
    ax.set_xlabel('$k/(3pn^2)$')
    ax.set_ylabel('lowest non-BPS singlet energy / $E_{\\rm vac}$')
    ax.set_title('(d) next rank, $p=2$: singlet gap per degree')
    ax.legend(fontsize=a.fontsize - 3)
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(f'{a.out}.{ext}', dpi=160)
    print('saved', a.out + '.{pdf,png}')


if __name__ == '__main__':
    main()
