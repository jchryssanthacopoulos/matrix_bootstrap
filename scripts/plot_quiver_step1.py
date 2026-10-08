#!/usr/bin/env python3
"""
Figure for step 1 of the quiver programme (research/notes/quiver_project.md section 6): near-BPS singlet multiplet
edges and level statistics of the n = 2 quiver.

  (a) lowest singlet multiplet energy E_0 of each Q-pair (k, k+3) against the multiplet charge |q| = |k + 3/2 - 3pn^2/2|,
      for (n,p) = (2,2) (hidden N = 4 on singlets) and (2,3); the Turiaci-Witten law E_0 = c q^2 is drawn through the
      |q| = 1.5 point of each model.
  (b) <r> of every symmetry-resolved ensemble with its size-matched references.
  (c) unfolded spacing distribution of the largest ensemble against the GOE surmise and Poisson.
  (d) near-edge counting functions N_q(E) of the two innermost (2,3) multiplet sectors (|q| = 1.5, 4.5) against the
      Turiaci-Witten density rho_q(E) ~ sinh(2 pi sqrt((E - E_0(q))/E_s)) / E with E_0(q) = q^2 E_s / (4 qhat^2),
      qhat = 3: E_s fixed by E_0(1.5), the normalisation by the |q| = 1.5 count; the |q| = 4.5 curve is a prediction.

Inputs: results/data/quiver_singlet_pairs_n2_p{2,3}.json, quiver_singlet_edge_n2_p3_k{9,12,15}.json,
quiver_p2_n4_*.json, quiver_singlet_pair_n2_p3_k9.json (optional).

    python scripts/plot_quiver_step1.py --out results/figures/quiver_step1
"""
import argparse, glob, json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
D = lambda f: os.path.join(ROOT, 'results', 'data', f)
COL = {2: '#1f77b4', 3: '#d62728'}
MK = {2: 's', 3: 'o'}


def load(f):
    p = D(f)
    return json.load(open(p)) if os.path.exists(p) else None


def edges(p):
    """|q| -> E_0 from the dense pairs file, overridden/extended by Lanczos edge files (p = 3)."""
    out = {}
    pr = load(f'quiver_singlet_pairs_n2_p{p}.json')
    for key, v in pr['pairs'].items():
        if v['E0'] is not None and v['q'] < 0:
            out[abs(v['q'])] = v['E0']
    for k in (9, 12, 15):
        e = load(f'quiver_singlet_edge_n2_p{p}_k{k}.json')
        if not e:
            continue
        low = [r for r in e['levels'] if r['kind'] == 'singlet' and r['member'] == 'lower']
        if low:
            out[abs(low[0]['q'])] = min(r['E'] for r in low)
    big = load(f'quiver_singlet_pair_n2_p{p}_k9_f64.json') or load(f'quiver_singlet_pair_n2_p{p}_k9.json')
    if big:
        out[abs(big['q'])] = big['E0']
    return dict(sorted(out.items()))


def local_unfold(levels, trim=0.1, w=8):
    """Spacings normalised by the local mean spacing over +-w neighbours (central part of the spectrum)."""
    E = np.sort(np.asarray(levels, float))
    lo, hi = int(trim * len(E)), int((1 - trim) * len(E))
    out = []
    for i in range(max(lo, w), min(hi, len(E) - w - 1)):
        rho = 2 * w / (E[i + w] - E[i - w])
        out.append((E[i + 1] - E[i]) * rho)
    out = np.array(out)
    return out / out.mean()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default=os.path.join(ROOT, 'results', 'figures', 'quiver_step1'))
    ap.add_argument('--fontsize', type=float, default=10.5)
    a = ap.parse_args()
    plt.rcParams.update({'font.size': a.fontsize})
    fig, axs2 = plt.subplots(2, 2, figsize=(13, 9.5))
    axs = axs2.ravel()

    # (a) multiplet edges
    ax = axs[0]
    for p in (2, 3):
        E = edges(p)
        qs, es = np.array(list(E.keys())), np.array(list(E.values()))
        lab = '$(n,p)=(2,2)$, hidden $\\mathcal{N}=4$' if p == 2 else '$(n,p)=(2,3)$, $\\mathcal{N}=2$'
        short = '$(2,2)$' if p == 2 else '$(2,3)$'
        ax.plot(qs, es, MK[p], color=COL[p], ms=7, label=lab, zorder=3)
        ax.plot(qs, es, '-', color=COL[p], lw=1, alpha=0.4)
        c = es[qs == 1.5][0] / 1.5 ** 2
        qq = np.linspace(1.2, 17, 100)
        ax.plot(qq, c * qq ** 2, '--', color=COL[p], lw=1, alpha=0.8,
                label=f'Turiaci–Witten $E_0\\propto q^2$ through $|q|=1.5$, {short}')
    ax.set_yscale('log'); ax.set_xscale('log')
    ticks = [1.5, 4.5, 7.5, 10.5, 13.5, 16.5]
    ax.set_xticks(ticks); ax.set_xticklabels([str(t) for t in ticks]); ax.minorticks_off()
    ax.set_xlabel('multiplet charge $|q|=|k+3/2-3pn^2/2|$')
    ax.set_ylabel('lowest singlet multiplet energy $E_0$')
    ax.set_title('(a) near-BPS edges (gauge singlets)')
    ax.legend(fontsize=a.fontsize - 3, loc='lower right')

    # (b) <r> summary
    ax = axs[1]
    rows = []
    pr3 = load('quiver_singlet_pairs_n2_p3.json')
    v = pr3['pairs']['6,9']; ref = pr3['references']['6,9']
    rows.append(('(2,3) pair (6,9), $|q|$=10.5', v['stats'], ref['loe'], ref['poisson'], 3))
    big = load('quiver_singlet_pair_n2_p3_k9_f64.json') or load('quiver_singlet_pair_n2_p3_k9.json')
    if big:
        rows.append(('(2,3) pair (9,12), $|q|$=7.5', big['stats'], big['references']['loe'], big['references']['poisson'], 3))
    for f in sorted(glob.glob(D('quiver_p2_n4_*_seed*.json'))):
        r = json.load(open(f))
        tag = os.path.basename(f).replace('quiver_p2_n4_', '').replace('.json', '')
        for k in ('9', '6'):
            b = r['bottoms'].get(k, {})
            st = b.get('stats_+')
            key = f"{k}+"
            if st and key in r.get('references', {}):
                rr = r['references'][key]
                rm = [kk for kk in rr if kk != 'poisson'][0]
                rows.append((f'(2,2) bottom k={k}, {tag}', st, rr[rm], rr['poisson'], 2))
    y = np.arange(len(rows))[::-1]
    for yi, (lab, st, rm, rp, p) in zip(y, rows):
        ax.errorbar(st['mean_r'], yi, xerr=st['sem_r'], fmt=MK[p], color=COL[p], ms=6, capsize=3, zorder=3)
        ax.plot([rm['mean'] - rm['std'], rm['mean'] + rm['std']], [yi, yi], color='0.25', lw=5, alpha=0.25, solid_capstyle='butt')
        ax.plot([rp['mean'] - rp['std'], rp['mean'] + rp['std']], [yi, yi], color='0.6', lw=5, alpha=0.35, solid_capstyle='butt')
        ax.text(rm['mean'], yi + 0.3, rm['kind'].upper().replace('LOE', 'GOE'), ha='center', fontsize=a.fontsize - 4, color='0.25')
    ax.axvline(0.386, color='0.6', ls=':', lw=1); ax.text(0.388, -0.9, 'Poisson', ha='left', fontsize=a.fontsize - 3, color='0.4')
    ax.set_ylim(-1.2, len(rows) - 0.3)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=a.fontsize - 3)
    ax.set_xlabel('mean ratio of consecutive spacings $\\langle r\\rangle$')
    ax.set_title('(b) level statistics (bars: size-matched references)')
    ax.set_xlim(0.25, 0.75)

    # (c) spacing distribution of the largest ensemble
    ax = axs[2]
    src = big if big else pr3['pairs']['6,9']
    s = local_unfold(src['levels'])
    lab = '(2,3) pair (9,12)' if big else '(2,3) pair (6,9)'
    ax.hist(s, bins=np.linspace(0, 4, 33), density=True, color=COL[3], alpha=0.45, label=f'{lab}, {len(s)} spacings')
    x = np.linspace(0, 4, 400)
    ax.plot(x, np.pi * x / 2 * np.exp(-np.pi * x ** 2 / 4), '-', color='0.15', lw=1.5, label='GOE (Wigner surmise)')
    ax.plot(x, np.exp(-x), '--', color='0.45', lw=1.5, label='Poisson')
    ax.set_xlabel('unfolded spacing $s$'); ax.set_ylabel('$P(s)$')
    ax.set_title(f'(c) spacing distribution (local unfolding, $\\pm$8 levels)')
    ax.legend(fontsize=a.fontsize - 2)
    # (d) Turiaci-Witten near-edge test at (2,3)
    ax = axs[3]
    lev = {}
    for k, qabs in ((15, 1.5), (12, 4.5)):
        e = load(f'quiver_singlet_edge_n2_p3_k{k}.json')
        lev[qabs] = np.sort([r['E'] for r in e['levels'] if r['kind'] == 'singlet' and r['member'] == 'lower'])
    E0a = lev[1.5][0]
    Es = E0a * 36 / 1.5 ** 2                                   # E_0(q) = q^2 E_s / 36
    def NTW(E, q, A):
        e0 = q ** 2 * Es / 36
        x = np.linspace(e0, E, 4000)
        f = np.sinh(2 * np.pi * np.sqrt(np.maximum(x - e0, 0) / Es)) / x
        return A * np.trapezoid(f, x)
    A = (len(lev[1.5]) - 0.5) / NTW(lev[1.5][-1], 1.5, 1.0)   # midpoint of the last step
    for qabs, mk, ls in ((1.5, 'o', '-'), (4.5, '^', '--')):
        L = lev[qabs]
        ax.step(np.concatenate([[L[0]], L]), np.concatenate([[0], np.arange(1, len(L) + 1)]), where='post',
                color=COL[3], lw=1.5, ls=ls, label=f'measured, $|q|={qabs}$ ({len(L)} lowest levels)')
        grid = np.linspace(qabs ** 2 * Es / 36 * 1.0001, max(L[-1], qabs ** 2 * Es / 36) * 1.05, 120)
        ax.plot(grid, [NTW(g, qabs, A) for g in grid], ls, color='0.2', lw=1,
                label=f'Turiaci–Witten, $|q|={qabs}$' + (' (prediction)' if qabs == 4.5 else ' (fit)'))
    ax.axvline(4.5 ** 2 * Es / 36, color='0.5', ls=':', lw=1)
    ax.text(4.5 ** 2 * Es / 36, 1, ' predicted $E_0(4.5)$', fontsize=a.fontsize - 3, color='0.4', rotation=90, va='bottom')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_ylim(0.7, 3e3); ax.set_xlim(0.04, 8)
    ax.set_xlabel('singlet multiplet energy $E$ at $(n,p)=(2,3)$'); ax.set_ylabel('number of multiplets below $E$')
    ax.set_title(f'(d) near-edge counting vs Turiaci–Witten ($E_s$={Es:.3f})')
    ax.legend(fontsize=a.fontsize - 3, loc='upper left')
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(f'{a.out}.{ext}', dpi=160)
    print('saved', a.out + '.{pdf,png}')


if __name__ == '__main__':
    main()
