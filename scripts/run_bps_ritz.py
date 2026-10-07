#!/usr/bin/env python3
"""
Phase 1 of docs/nearbps_bootstrap_plan.md (small blocks, dense references): Rayleigh-Ritz upper bounds on the sector
ground energy E_0(k) from BPS states of window sectors k_B dressed by adjoint-valued words of charge k - k_B.

    python scripts/run_bps_ritz.py --N 2 --k 4 --kB 5,6,7 --max-len 3 --refs all
    python scripts/run_bps_ritz.py --N 2 --k 4 --kB 5 --max-len 5 --refs random --n-random 1 --seed 1

Trial space (src/bps_ritz.py): span{O_ii B, O_ij E_ji B} on the zero-weight block of sector k, for every selected
reference B (zero-weight zero mode of H in sector k_B, Casimir-resolved) and every word of length <= L with the right
charge, for L = |q|, |q|+2, ..., max-len.  Reported per L: number of trial vectors, rank, lowest Ritz values, the
Casimir <C_2> of the lowest Ritz vector, and (blocks <= 3000 states) the exact E_0 and the weight of the exact
ground space inside the trial space.  Every Ritz value is a rigorous upper bound on E_0(k) (up to round-off).
"""
import argparse
import json
import os
import resource
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
from fermion_matrix_model import chen_C                                       # noqa: E402
import bps_ritz as br                                                         # noqa: E402


def peak_gb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1e9 if sys.platform == 'darwin' else r / 1e6


def git_rev():
    try:
        rev = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'], cwd=ROOT, text=True).strip()
        dirty = subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True).strip()
        return rev + ('+uncommitted' if dirty else '')
    except Exception:
        return 'unknown'


def select_refs(B, c2, how, n_random, rng):
    """Columns of the Casimir-resolved BPS basis to use as references."""
    if how == 'all':
        return [B[:, n] for n in range(B.shape[1])], 'all'
    if how.startswith('casimir='):
        c = float(how.split('=')[1])
        sel = np.abs(c2 - c) < 1e-6
        return [B[:, n] for n in np.nonzero(sel)[0]], how
    if how == 'singlet':
        sel = np.abs(c2) < 1e-6
        return [B[:, n] for n in np.nonzero(sel)[0]], how
    if how in ('random', 'random-singlet'):
        cols = B if how == 'random' else B[:, np.abs(c2) < 1e-6]
        out = []
        for _ in range(n_random):
            v = cols @ rng.standard_normal(cols.shape[1])
            out.append(v / np.linalg.norm(v))
        return out, f'{how}x{n_random}'
    raise ValueError(how)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--N', type=int, default=2)
    ap.add_argument('--k', type=int, required=True, help='target sector')
    ap.add_argument('--kB', default=None, help='comma-separated reference (window) sectors')
    ap.add_argument('--max-len', type=int, default=3)
    ap.add_argument('--refs', default='all', help="all | singlet | casimir=<c> | random | random-singlet")
    ap.add_argument('--n-random', type=int, default=1)
    ap.add_argument('--no-offdiag', action='store_true', help='only O_ii B (drop O_ij E_ji B)')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--budget-gb', type=float, default=4.0,
                    help='refuse to build a trial set whose dense storage + Gram work exceeds this estimate')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    N, p, k = a.N, 3, a.k
    C = chen_C()
    rng = np.random.default_rng(a.seed)
    kBs = [int(x) for x in a.kB.split(',')]
    t0 = time.time()
    target = br.Block(N, p, k)
    Hk = br.SectorH(N, C, target)
    rec = dict(params=vars(a), git=git_rev(), date=time.strftime('%Y-%m-%d %H:%M'), target_dim=target.dim)
    exact = None
    if target.dim <= 3000:
        lam, U = np.linalg.eigh(Hk.dense())
        E0 = lam[0]
        G = U[:, lam < E0 + 1e-8]
        exact = dict(E0=float(E0), degeneracy=int(G.shape[1]), levels=[float(x) for x in lam[:8]])
        print(f"N={N} k={k}: zero-weight block {target.dim}, exact E_0 = {E0:.8f} (x{G.shape[1]} in block)")
    rec['exact'] = exact
    refs = {}
    for kB in kBs:
        blk = br.Block(N, p, kB)
        B, c2, _ = br.bps_basis_dense(N, C, blk)
        vals, cnt = np.unique(np.round(c2, 6), return_counts=True)
        chosen, tag = select_refs(B, c2, a.refs, a.n_random, rng)
        refs[kB] = [blk.sparse(v) for v in chosen]
        print(f"  k_B={kB}: block {blk.dim}, {B.shape[1]} zero-weight BPS states, Casimir content "
              f"{dict(zip(vals.tolist(), cnt.tolist()))}; using {len(chosen)} ({tag})")
    rows = []
    max_len = a.max_len
    rec['budget_gb'] = a.budget_gb
    Wacc = np.zeros((target.dim, 0))
    lengths_done = {kB: 0 for kB in kBs}
    stop = False
    for L in range(1, max_len + 1):
        new = []
        for kB in kBs:
            q = k - kB
            if abs(q) > L or (L - q) % 2:
                continue
            words = br.words_of_charge(q, L, p)
            m_est = Wacc.shape[1] + sum(len(r) for r in [refs[kB]]) * len(words) * N * N
            gb = 8 * (2 * m_est * target.dim + 3 * min(m_est, target.dim) ** 2) / 1e9
            if gb > a.budget_gb:
                print(f"  L<={L}: up to {m_est} trial vectors on {target.dim} states, est {gb:.1f} GB > budget "
                      f"{a.budget_gb} GB; stopping", flush=True)
                new, stop = [], True
                break
            W, _ = br.trial_vectors(N, p, refs[kB], words, target, include_offdiag=not a.no_offdiag)
            new.append(W)
            lengths_done[kB] = L
        if stop:
            break
        if not new:
            continue
        Wacc = np.column_stack([Wacc] + new)
        rr = br.rayleigh_ritz(Hk, Wacc)
        v0 = rr['vectors'][:, 0]
        c2v = float(v0 @ br.casimir_apply(br.casimir_W(N, p, target), v0))
        row = dict(L=L, n_trial=int(Wacc.shape[1]), rank=rr['rank'], ritz=[float(x) for x in rr['values']],
                   residual0=float(rr['residuals'][0]), C2_lowest=c2v)
        msg = (f"  L<={L}: {Wacc.shape[1]} trial vectors, rank {rr['rank']} / {target.dim}; Ritz "
               f"{', '.join(f'{x:.6f}' for x in rr['values'][:4])}; <C2> of lowest {c2v:.4f}")
        if exact is not None:
            # weight of the exact ground space inside the trial space, and of the lowest Ritz vector on it
            Us, sv, _ = np.linalg.svd(Wacc, full_matrices=False)
            basis = Us[:, sv > 1e-10 * sv[0]]
            cap = float(np.linalg.norm(basis.T @ G) ** 2 / G.shape[1])
            ov = float(np.linalg.norm(G.T @ v0) ** 2)
            row.update(ground_weight_in_V=cap, lowest_ritz_overlap=ov, ratio=float(rr['values'][0] / exact['E0']))
            msg += f"; ratio to exact {rr['values'][0] / exact['E0']:.4f}; |P_V psi_0|^2 = {cap:.4f}"
        print(msg, flush=True)
        rows.append(row)
    rec['rows'] = rows
    rec['time_s'] = round(time.time() - t0, 1)
    rec['peak_rss_gb'] = round(peak_gb(), 3)
    print(f"  [{rec['time_s']} s, peak {rec['peak_rss_gb']} GB]")
    if a.out:
        json.dump(rec, open(a.out, 'w'), indent=1)
        print('  saved', os.path.relpath(a.out, ROOT))


if __name__ == '__main__':
    main()
