#!/usr/bin/env python3
"""Figures from the sector scans produced by scripts/scan_sector.py (results/data/scan_*.json).

  python scripts/plot_scans.py            -> results/figures/sector_scan_single_matrix_N3.png
                                             results/figures/sector_scan_three_matrix_N2_k2.png
"""
import json, os, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
DATA = os.path.join(ROOT, 'results', 'data'); FIG = os.path.join(ROOT, 'results', 'figures')
C = {'blue': '#2a78d6', 'orange': '#eb6834', 'aqua': '#1baf7a', 'violet': '#4a3aa7', 'ink': '#0b0b0b', 'muted': '#52514e'}
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True,
                     'grid.alpha': 0.25, 'legend.frameon': False})


def load(name):
    p = os.path.join(DATA, name)
    return json.load(open(p)) if os.path.exists(p) else None


def feasible_intervals(d, step=None):
    Es = [r['E'] for r in d['scan'] if r['feasible']]
    if not Es:
        return []
    grid = np.diff([r['E'] for r in d['scan']])[0]
    iv = [[Es[0], Es[0]]]
    for E in Es[1:]:
        if abs(E - iv[-1][1] - grid) < 1e-9:
            iv[-1][1] = E
        else:
            iv.append([E, E])
    return iv, grid


def region(ax, d, name, color, label):
    E = np.array([r['E'] for r in d['scan'] if r['feasible'] and name in r])
    lo = np.array([r[name][0] for r in d['scan'] if r['feasible'] and name in r])
    hi = np.array([r[name][1] for r in d['scan'] if r['feasible'] and name in r])
    ax.fill_between(E, lo, hi, color=color, alpha=0.25, lw=0, label=label)
    ax.plot(E, lo, color=color, lw=1.5); ax.plot(E, hi, color=color, lw=1.5)


def archipelago_row(ax, d, y, color, label, exact=None):
    iv, grid = feasible_intervals(d)
    for a, b in iv:
        if b - a < grid / 2:
            ax.plot([a], [y], 'o', color=color, ms=7, mec='white', mew=1)
        else:
            ax.plot([a, b], [y, y], color=color, lw=6, solid_capstyle='butt', alpha=0.85)
    if exact is not None:
        ax.plot(sorted(set(np.round(exact, 4))), [y] * len(set(np.round(exact, 4))), 'x', color=C['ink'], ms=8, mew=1.5, zorder=5)
    ax.text(-3.5, y, label, ha='right', va='center', fontsize=9, color=C['muted'])


def fig_single_matrix():
    k2 = load('scan_p1_N3_k2_L2.json'); k3 = load('scan_p1_N3_k3_L2.json')
    e_k2 = load('scan_p1_N3_k2_L2_eig.json'); e_k3_2 = load('scan_p1_N3_k3_L2_eig.json'); e_k3_3 = load('scan_p1_N3_k3_L3_eig.json')
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axs[0]
    rows = [(e_k2, 3, C['blue'], r'$N_\Psi=2$, level 2'), (e_k3_2, 2, C['orange'], r'$N_\Psi=3$, level 2'), (e_k3_3, 1, C['aqua'], r'$N_\Psi=3$, level 3')]
    for d, y, col, lab in rows:
        if d: archipelago_row(ax, d, y, col, lab, exact=d['exact']['E'])
    ax.set_yticks([]); ax.set_ylim(0.3, 3.7); ax.set_xlim(-4, 76); ax.set_xlabel(r'$E$'); ax.grid(axis='y', alpha=0)
    ax.set_title(r'(a) Feasible energies with eigenstate constraints $\phi(OH)=E\,\phi(O)$', fontsize=10, loc='left')
    ax.plot([], [], 'x', color=C['ink'], label='exact eigenvalues (ED)'); ax.plot([], [], color=C['muted'], lw=6, alpha=0.85, label='feasible'); ax.legend(loc='upper right', fontsize=8)
    ax = axs[1]
    if k2: region(ax, k2, 'C2gauge', C['blue'], r'$N_\Psi=2$ allowed region (level 2)')
    if k3: region(ax, k3, 'C2gauge', C['orange'], r'$N_\Psi=3$ allowed region (level 2)')
    Eg = np.linspace(0, 72, 2); ax.plot(Eg, (72 - Eg) / 9, '--', color=C['ink'], lw=1, label=r'exact: $E=3N(N^2-1)-9C_2$')
    for d, col in ((k2, C['blue']), (k3, C['orange'])):
        if d: ax.plot(d['exact']['E'], d['exact']['C2gauge'], 'x', color=col, ms=8, mew=1.5)
    ax.set_xlabel(r'$E$'); ax.set_ylabel(r'$\langle \hat C_2^{\rm gauge}\rangle$'); ax.set_title(r'(b) Allowed region in the $(E,\langle\hat C_2\rangle)$ plane', fontsize=10, loc='left'); ax.legend(fontsize=8)
    fig.suptitle(r'Single-matrix model $Q=\mathrm{Tr}\,\Psi^3$, $N=3$: sector bootstrap scans vs exact solution', fontsize=11)
    fig.tight_layout(); out = os.path.join(FIG, 'sector_scan_single_matrix_N3.png'); fig.savefig(out, dpi=170); print('wrote', out)


def fig_three_matrix():
    L2 = load('scan_p3_N2_k2_L2.json')
    if not L2:
        print('no three-matrix scan data yet'); return
    ex = L2['exact']
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.4))
    for ax, name, ylab in zip(axs, ['C2gauge', 'Ns'], [r'$\langle\hat C_2^{\rm gauge}\rangle$', r'$\langle N_s\rangle$  (flavor-symmetric occupation)']):
        region(ax, L2, name, C['blue'], 'allowed region, level 2')
        note = '' if ex.get(name + '_commutes', True) else ' (range within each level)'
        ax.plot(ex['E'], ex[name], 'x', color=C['ink'], ms=6, mew=1.2, label='exact eigenstates (ED)' + note)
        ax.axvline(ex['E'][0], color=C['orange'], lw=1, ls='--'); ax.text(ex['E'][0] + 1, ax.get_ylim()[1] * 0.02 + ax.get_ylim()[0], r'$E_0^{\rm exact}$', color=C['orange'], fontsize=8)
        ax.set_xlabel(r'$E$'); ax.set_ylabel(ylab)
    h, l = axs[1].get_legend_handles_labels(); fig.legend(h, l, loc='lower center', ncol=2, fontsize=8, bbox_to_anchor=(0.5, -0.01))
    axs[0].set_title(r'(a) gauge Casimir vs energy', loc='left', fontsize=10); axs[1].set_title(r'(b) flavor-symmetric occupation vs energy', loc='left', fontsize=10)
    axs[1].set_ylim(-0.05, 2.1)
    fig.suptitle(r'Three-matrix model (Chen), $N=2$, sector $N_\Psi=2$ ($\dim=66$): level-2 allowed regions vs exact spectrum', fontsize=11)
    fig.tight_layout(rect=(0, 0.05, 1, 1)); out = os.path.join(FIG, 'sector_scan_three_matrix_N2_k2.png'); fig.savefig(out, dpi=170); print('wrote', out)


if __name__ == '__main__':
    fig_single_matrix(); fig_three_matrix()
