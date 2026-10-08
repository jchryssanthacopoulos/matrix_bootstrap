#!/usr/bin/env python3
"""
Figure for the quiver singlet bootstrap (research/notes/quiver_project.md section 12; docs/derivations.md D21).

  (a) BPS-exclusion margin t* (level 2, BPS rows on) against the filling fraction k/(3pn^2) at (n,p) = (2,2), (2,3),
      (3,3).  t* < 0: no singlet BPS state at that degree.  Hollow markers: degrees where BPS singlets exist (exact).
  (b) Lowest singlet energy per degree, in units of E_vac = n^3 sum C^2: level-2 lower bounds (open markers, dashed;
      zero bounds omitted) against the exact values (filled markers, solid) where known; BPS degrees marked at the
      bottom edge.

Inputs: results/data/quiver_bootstrap_n{n}_p{p}_L2.jsonl, quiver_singlet_pairs_n2_p2.json,
quiver_singlet_pairs_n2_p3.json, quiver_singlet_edge_n2_p3_k{9,12,15}.json, quiver_singlet_sw_n3_p3_seed3.json.

    python scripts/plot_quiver_bootstrap.py --out results/figures/quiver_bootstrap
    python scripts/plot_quiver_bootstrap.py --figsize 10 4.1 --fontsize 11.5 --out research/tex/report_figs/quiver_bootstrap
"""
import argparse, json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
DATA = os.path.join(ROOT, 'results', 'data')
COLS = {(2, 2): '#0072B2', (2, 3): '#D55E00', (3, 3): '#009E73'}
MARK = {(2, 2): 's', (2, 3): 'o', (3, 3): '^'}


def load_runs(n, p):
    pth = os.path.join(DATA, f'quiver_bootstrap_n{n}_p{p}_L2.jsonl')
    recs = [json.loads(l) for l in open(pth)] if os.path.exists(pth) else []
    margin, bound = {}, {}
    for r in recs:
        if r['mode'] == 'margin' and not r['level'].get('cone_casimir'):
            margin[r['k']] = r['margin']
        elif r['mode'] == 'energy':
            bound[r['k']] = r['bound']
    return margin, bound, (np.array(recs[0]['couplings'], float) if recs else None)


def exact_E0(n, p):
    """Lowest singlet energy per degree (exact), from the stored spectra."""
    out = {}
    if (n, p) == (2, 2):
        pr = json.load(open(os.path.join(DATA, 'quiver_singlet_pairs_n2_p2.json')))['pairs']
        e = {tuple(int(x) for x in key.split(',')): v['E0'] for key, v in pr.items() if v['E0'] is not None}
        for k in range(0, 13, 3):
            c = [e[(a, b)] for (a, b) in e if k in (a, b)]
            out[k] = min(c) if c else None
        out[12] = 0.0
    elif (n, p) == (2, 3):
        pr = json.load(open(os.path.join(DATA, 'quiver_singlet_pairs_n2_p3.json')))['pairs']
        e = {tuple(int(x) for x in key.split(',')): v['E0'] for key, v in pr.items() if v['E0'] is not None}
        for k in (0, 3, 6):
            out[k] = min(e[(a, b)] for (a, b) in e if k in (a, b))
        for k in (9, 12, 15):
            d = json.load(open(os.path.join(DATA, f'quiver_singlet_edge_n2_p3_k{k}.json')))
            out[k] = min(l['E'] for l in d['levels'] if l['kind'] == 'singlet')
        out[18] = 0.0
    elif (n, p) == (3, 3):
        d = json.load(open(os.path.join(DATA, 'quiver_singlet_sw_n3_p3_seed3.json')))
        for k, v in d['degrees'].items():
            out[int(k)] = v['lowest_eigenvalue']
    return {k: v for k, v in out.items() if v is not None}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default=os.path.join(ROOT, 'results', 'figures', 'quiver_bootstrap'))
    ap.add_argument('--fontsize', type=float, default=10.5)
    ap.add_argument('--figsize', type=float, nargs=2, default=(12.5, 4.8), metavar=('W', 'H'),
                    help='inches; the report version uses 10 4.1 with --fontsize 11.5')
    a = ap.parse_args()
    plt.rcParams.update({'font.size': a.fontsize})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=tuple(a.figsize))
    for (n, p) in ((2, 2), (2, 3), (3, 3)):
        margin, bound, C = load_runs(n, p)
        if not margin:
            continue
        top = 3 * p * n * n
        col, mk = COLS[(n, p)], MARK[(n, p)]
        ks = sorted(margin)
        x = np.array(ks) / top
        t = np.array([margin[k] for k in ks])
        ax1.plot(x, t, '-', color=col, lw=1.2)
        exc = t < -1e-4
        ax1.plot(x[exc], t[exc], mk, color=col, ms=6, label=f'$(n,p)=({n},{p})$: excluded through $k={max(np.array(ks)[exc])}$')
        ax1.plot(x[~exc], t[~exc], mk, color=col, ms=6, mfc='white')
        # (b) energies
        Evac = n ** 3 * (C ** 2).sum()
        ex = exact_E0(n, p)
        kb = sorted(k for k in bound if bound[k] > 1e-6 * Evac)          # a zero bound carries no information
        ax2.plot(np.array(kb) / top, [bound[k] / Evac for k in kb], mk + '--', color=col, mfc='white', ms=6, lw=1,
                 label=f'$({n},{p})$ level-2 lower bound')
        ke = sorted(k for k in ex if k <= top / 2 and ex[k] > 0)
        ax2.plot(np.array(ke) / top, [ex[k] / Evac for k in ke], mk + '-', color=col, ms=6, lw=1.2,
                 label=f'$({n},{p})$ exact ($E_{{\\rm vac}}={Evac:.0f}$)')
        kz = [k for k in ex if k <= top / 2 and ex[k] == 0]                  # BPS degrees (exact zero)
        if kz:
            ax2.plot(np.array(kz) / top, [1.6e-5] * len(kz), 'v', color=col, ms=8, mfc=col if n == 2 and p == 2 else 'none')
    ax1.axhline(0, color='0.4', lw=0.8)
    ax1.axvline(0.5, color='0.5', ls=':', lw=1)
    ax1.set_xlabel('filling fraction $k/(3pn^2)$')
    ax1.set_ylabel('BPS-exclusion margin $t^*$ (level 2)')
    ax1.set_title('(a) singlet BPS exclusion: $t^*<0$ means no BPS singlet at degree $k$', fontsize=a.fontsize - 0.5)
    ax1.text(0.49, -2.0, 'half filling', rotation=90, ha='right', va='center', fontsize=a.fontsize - 2, color='0.4')
    ax1.legend(fontsize=a.fontsize - 2, loc='lower right')
    ax1.set_xlim(0, 0.55)
    ax2.set_yscale('log')
    ax2.set_xlabel('filling fraction $k/(3pn^2)$')
    ax2.set_ylabel('lowest singlet energy / $E_{\\rm vac}$')
    ax2.set_title('(b) singlet gap per degree: bootstrap bound vs exact', fontsize=a.fontsize - 0.5)
    ax2.set_xlim(0, 0.55)
    ax2.set_ylim(1e-5, 1.5)
    ax2.axvline(0.5, color='0.5', ls=':', lw=1)
    ax2.annotate('$E_0=0$: BPS singlets at\nhalf filling, (2,2) and (2,3)', (0.5, 1.6e-5), xytext=(0.545, 1.2e-2),
                 textcoords='data', ha='right', va='center', fontsize=a.fontsize - 2.5, color='0.35',
                 bbox=dict(fc='white', ec='none', pad=1.0),
                 arrowprops=dict(arrowstyle='->', color='0.5', lw=0.8, shrinkB=6))
    ax2.legend(fontsize=a.fontsize - 2.5, ncol=1, loc='lower left')
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(f'{a.out}.{ext}', dpi=160)
    print('saved', a.out + '.{pdf,png}')


if __name__ == '__main__':
    main()
