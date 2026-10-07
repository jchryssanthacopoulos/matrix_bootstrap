#!/usr/bin/env python3
"""
Phase 1 at N=3 (docs/nearbps_bootstrap_plan.md): Rayleigh-Ritz upper bounds on E_0(k) for k = 9, 10, 11 from a
singlet BPS state of the window centre k_B = 13, dressed by adjoint-valued words of charge k - 13.

Stage 'bps': a singlet BPS state B in the zero-weight block of sector 13 as the lowest eigenvector of H + mu C_2
(ground space = singlet BPS states, 243-dimensional: refined index |I_{1,1,w}| = 81 per Z_3 charge w, and class
k = 1 mod 3 has no other degree carrying cohomology since E_0(k) > 0 for k <= 11 (ED) and by particle-hole symmetry).
Lanczos from a random start returns the (normalised) projection of the start onto that space.  Saved as .npy.

Stage 'ritz': trial vectors O_ii B (B is a singlet, so E_ji B = 0) for all words of length L of charge q = k - 13,
on the zero-weight block of sector k; lowest Ritz values and <C_2> of the lowest Ritz vector.

    scripts/run_guarded.sh 9 log python3 scripts/run_bps_ritz_n3.py bps --seed 1
    scripts/run_guarded.sh 9 log python3 scripts/run_bps_ritz_n3.py ritz --k 11 --lengths 2

Memory (estimated from the k=11 ED, 3.0 GB peak at 750 699 states): stage bps ~5-7 GB (block ~1.15M states,
H ~2.5 GB, Casimir maps ~0.8 GB, Lanczos basis ~0.4 GB); stage ritz ~ 8 * dim * (2 * n_trial + 3) bytes + H.
"""
import argparse
import json
import os
import resource
import sys
import time

import numpy as np
import scipy.sparse.linalg as sla

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
from fermion_matrix_model import chen_C                                       # noqa: E402
import bps_ritz as br                                                         # noqa: E402
from run_bps_ritz import git_rev                                              # noqa: E402

N, P, KB = 3, 3, 13
EXACT = {9: (1.14573, '1'), 10: (0.345019, '8'), 11: (0.061408, '10+10bar')}   # sector_ed_results.md


def peak_gb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1e9 if sys.platform == 'darwin' else r / 1e6


def stage_bps(a):
    t0 = time.time()
    blk = br.Block(N, P, KB)
    print(f"k_B={KB}: zero-weight block {blk.dim} [{time.time() - t0:.0f}s, {peak_gb():.2f} GB]", flush=True)
    H = br.SectorH(N, chen_C(), blk)
    print(f"  H built: H1 nnz {H.H1.nnz}, Ydag nnz {H.Yd.nnz} [{time.time() - t0:.0f}s, {peak_gb():.2f} GB]", flush=True)
    Ws = br.casimir_W(N, P, blk)
    print(f"  Casimir maps nnz {sum(W.nnz for W in Ws)} [{time.time() - t0:.0f}s, {peak_gb():.2f} GB]", flush=True)
    op = sla.LinearOperator((blk.dim, blk.dim), dtype=float, matvec=lambda x: H(x) + a.mu * br.casimir_apply(Ws, x))
    rng = np.random.default_rng(a.seed)
    vals, vecs = sla.eigsh(op, k=a.n_eig, which='SA', tol=a.tol, ncv=a.ncv, v0=rng.standard_normal(blk.dim))
    order = np.argsort(vals); vals, vecs = vals[order], vecs[:, order]
    B = vecs[:, 0]
    e_H = float(B @ H(B)); c2 = float(B @ br.casimir_apply(Ws, B)); res = float(np.linalg.norm(H(B)))
    print(f"  eigsh: lowest of H + {a.mu} C2: {', '.join(f'{x:.3e}' for x in vals)}; B: <H> = {e_H:.2e}, "
          f"|H B| = {res:.2e}, <C2> = {c2:.2e} [{time.time() - t0:.0f}s, {peak_gb():.2f} GB]", flush=True)
    out = os.path.join(ROOT, 'results', 'data', f'bps_N3_k13_singlet_seed{a.seed}')
    np.save(out + '.npy', B)
    json.dump(dict(N=N, k=KB, seed=a.seed, mu=a.mu, tol=a.tol, ncv=a.ncv, dim=blk.dim, lowest=[float(x) for x in vals],
                   H_expect=e_H, H_residual=res, C2=c2, git=git_rev(), time_s=round(time.time() - t0, 1),
                   peak_rss_gb=round(peak_gb(), 2)), open(out + '.json', 'w'), indent=1)
    print(f"  saved {os.path.relpath(out, ROOT)}.npy/.json")


def stage_ritz(a):
    t0 = time.time()
    src = br.Block(N, P, KB)
    B = np.load(os.path.join(ROOT, 'results', 'data', f'bps_N3_k13_singlet_seed{a.seed}.npy'))
    assert len(B) == src.dim
    Bs = src.sparse(B)
    del src
    tgt = br.Block(N, P, a.k)
    q = a.k - KB
    lengths = [int(x) for x in a.lengths.split(',')]
    n_words = sum(len(br.words_of_charge(q, L, P)) for L in lengths)
    m = n_words * N
    gb = 8 * tgt.dim * (2 * m + 3) / 1e9
    print(f"k={a.k} (charge {q}): block {tgt.dim}, {n_words} words x {N} diagonal components = {m} trial vectors, "
          f"dense storage est {gb:.2f} GB (+ H)", flush=True)
    if gb > a.budget_gb:
        raise SystemExit(f"refusing: {gb:.1f} GB > budget {a.budget_gb} GB")
    H = br.SectorH(N, chen_C(), tgt)
    print(f"  H built [{time.time() - t0:.0f}s, {peak_gb():.2f} GB]", flush=True)
    Ws = br.casimir_W(N, P, tgt)
    rows, Wacc = [], np.zeros((tgt.dim, 0))
    for L in lengths:
        W, _ = br.trial_vectors(N, P, [Bs], br.words_of_charge(q, L, P), tgt, include_offdiag=False)
        Wacc = np.column_stack([Wacc, W]); del W
        rr = br.rayleigh_ritz(H, Wacc, n_ritz=8)
        c2 = [float(rr['vectors'][:, n] @ br.casimir_apply(Ws, rr['vectors'][:, n])) for n in range(len(rr['values']))]
        E0, irr = EXACT[a.k]
        print(f"  L in {lengths[:lengths.index(L) + 1]}: {Wacc.shape[1]} vectors, rank {rr['rank']}; Ritz "
              f"{', '.join(f'{x:.5f}' for x in rr['values'])}; <C2> {', '.join(f'{x:.2f}' for x in c2)}; "
              f"exact E_0 = {E0} ({irr}), ratio {rr['values'][0] / E0:.2f} [{time.time() - t0:.0f}s, {peak_gb():.2f} GB]",
              flush=True)
        rows.append(dict(lengths=lengths[:lengths.index(L) + 1], n_trial=int(Wacc.shape[1]), rank=rr['rank'],
                         ritz=[float(x) for x in rr['values']], C2=c2, residuals=[float(x) for x in rr['residuals']]))
    out = os.path.join(ROOT, 'results', 'data', f'bps_ritz_N3_k{a.k}_seed{a.seed}.json')
    json.dump(dict(N=N, k=a.k, kB=KB, ref=f'singlet BPS, seed {a.seed}', exact_E0=EXACT[a.k], rows=rows,
                   git=git_rev(), time_s=round(time.time() - t0, 1), peak_rss_gb=round(peak_gb(), 2)),
              open(out, 'w'), indent=1)
    print(f"  saved {os.path.relpath(out, ROOT)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('stage', choices=['bps', 'ritz'])
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--mu', type=float, default=1.0)
    ap.add_argument('--n-eig', type=int, default=3)
    ap.add_argument('--ncv', type=int, default=40)
    ap.add_argument('--tol', type=float, default=1e-10)
    ap.add_argument('--k', type=int, default=11)
    ap.add_argument('--lengths', default='2')
    ap.add_argument('--budget-gb', type=float, default=6.0)
    a = ap.parse_args()
    stage_bps(a) if a.stage == 'bps' else stage_ritz(a)


if __name__ == '__main__':
    main()
