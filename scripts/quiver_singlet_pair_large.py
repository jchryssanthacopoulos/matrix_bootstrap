#!/usr/bin/env python3
"""
Full singlet multiplet spectrum of one Q-pair (k, k+3) of the n = 2 quiver when the singlet sector at k is too large
for a dense singlet basis (research/notes/quiver_project.md section 6; target: (n,p) = (2,3), pair (9,12), about
10^4 multiplets).

Method.
  * Singlet projector on the W-symmetric zero-weight sector at degree k:
        Pi = prod_{v=0,1,2} prod_{j=1}^{jmax} (1 - J_v^2 / (j(j+1))),   J_v^2 = R_v^T R_v  (R_v = E^{(v)}_12 B),
    exact because the three node Casimirs commute and have spectra {j(j+1)}; jmax = 2 p n / 2 = maximal node spin
    (2pn doublets of fermion modes meet each node).
  * Y = Pi X for m = n(k) + oversample random columns X (chunked; stored float32), so span(Y) = singlet sector.
  * Gram B = Y^T Y and A = Y^T Q^T Q Y (chunked);  orthonormalise through B's eigenvectors (cut at tol) and
    diagonalise the compressed A.  Nonzero eigenvalues = squared singular values of Q on singlets = multiplet energies.
  * Checks: rank(B) = n(k) from the dual-Cauchy count; residual non-singlet weight of Y; stats as in
    quiver_singlet_spectrum.py (size-matched LOE / Poisson references).

    scripts/run_guarded.sh 11 log python3 scripts/quiver_singlet_pair_large.py --n 2 --p 3 --k 9 --dry-run
"""
import argparse, json, os, resource, sys, time
import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'scripts')); sys.path.insert(0, os.path.join(ROOT, 'src'))
from quiver_n2 import QuiverN2, sector, Q_sym, raising_sym                        # noqa: E402
from quiver_index import singlet_series                                            # noqa: E402
from quiver_singlet_spectrum import r_stats, reference_r, unfolded_spacings, charge   # noqa: E402


def peak_gb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1e9 if sys.platform == 'darwin' else r / 1e6


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--n', type=int, default=2); ap.add_argument('--p', type=int, default=3)
    ap.add_argument('--k', type=int, default=9); ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('--oversample', type=int, default=300)
    ap.add_argument('--chunk', type=int, default=500, help='columns per projector chunk')
    ap.add_argument('--qchunk', type=int, default=200, help='columns per Q^T Q chunk')
    ap.add_argument('--trim', type=float, default=0.1)
    ap.add_argument('--thr-rel', type=float, default=1e-4, help='zero-eigenvalue cut relative to the largest (float32 storage leaves ~1e-6 relative noise on kernel vectors)')
    ap.add_argument('--dtype', choices=['float32', 'float64'], default='float32',
                    help='storage of the projected vectors Y; float32 accumulates ~1e-5 relative error in the Gram matrices at dim ~5e4')
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    t0 = time.time()
    q = QuiverN2(a.n, a.p, a.seed)
    k = a.k
    nk, _, _ = singlet_series(a.n, a.p)
    S = {kk: sector(q, kk) for kk in (k, k + 3)}
    dim = len(S[k]['sb']['reps']); m = nk[k] + a.oversample
    Rs = [raising_sym(q, S[k]['masks'], S[k]['sb'], v, k) for v in range(3)]
    Ps = [(R.T @ R).tocsr() for R in Rs]
    Qk = Q_sym(q, S[k]['masks'], S[k]['sb'], S[k + 3]['sb'])
    jmax = a.p * a.n                                   # 2pn doublets per node -> max spin pn
    est = (dim * m * (8 if a.dtype == 'float64' else 4) + 2 * m * m * 8 + 3 * a.chunk * dim * 8 + a.qchunk * Qk.shape[0] * 8
           + sum(P.nnz for P in Ps) * 12 + Qk.nnz * 12) / 1e9
    print(f"(n,p)=({a.n},{a.p}) pair ({k},{k+3}): W-sym dim {dim}, n(k) = {nk[k]}, columns m = {m}; nnz J^2_v "
          f"{[P.nnz for P in Ps]}, Q nnz {Qk.nnz}; jmax {jmax}; memory estimate {est:.1f} GB  "
          f"[{time.time()-t0:.0f}s, {peak_gb():.2f} GB]", flush=True)
    if a.dry_run:
        return
    rng = np.random.default_rng(1)
    Y = np.empty((dim, m), dtype=np.float64 if a.dtype == 'float64' else np.float32)
    leak_max = 0.0
    for c0 in range(0, m, a.chunk):
        X = rng.standard_normal((dim, min(a.chunk, m - c0)))
        for P in Ps:
            for j in range(1, jmax + 1):
                X = X - (P @ X) / (j * (j + 1))
        # residual non-singlet weight (should be round-off): sum_v |J_v^2 X| / |X|
        leak = max(float(np.linalg.norm(P @ X) / max(np.linalg.norm(X), 1e-300)) for P in Ps)
        leak_max = max(leak_max, leak)
        Y[:, c0:c0 + X.shape[1]] = X.astype(Y.dtype)
        print(f"  projected columns {c0 + X.shape[1]}/{m}; non-singlet residual {leak:.1e}  [{time.time()-t0:.0f}s, "
              f"{peak_gb():.2f} GB]", flush=True)
    del Ps, Rs
    B = np.zeros((m, m)); A = np.zeros((m, m))
    for c0 in range(0, m, a.qchunk):
        Yc = Y[:, c0:c0 + a.qchunk]
        B[:, c0:c0 + Yc.shape[1]] = (Y.T @ Yc).astype(np.float64)
        Z = Qk.T @ (Qk @ Yc.astype(np.float64))
        A[:, c0:c0 + Yc.shape[1]] = (Y.T @ Z.astype(Y.dtype)).astype(np.float64)
        if (c0 // a.qchunk) % 10 == 0:
            print(f"  Gram/compression columns {c0 + Yc.shape[1]}/{m}  [{time.time()-t0:.0f}s, {peak_gb():.2f} GB]", flush=True)
    del Y
    B = (B + B.T) / 2; A = (A + A.T) / 2
    s, V = np.linalg.eigh(B)
    keep = s > 1e-6 * s.max()
    rank = int(keep.sum())
    T = V[:, keep] / np.sqrt(s[keep])
    del V, B
    Ac = T.T @ A @ T
    sig2 = np.linalg.eigvalsh((Ac + Ac.T) / 2)
    thr = a.thr_rel * sig2.max()
    lev = np.sort(sig2[sig2 > thr])
    print(f"  rank(B) = {rank} (n(k) = {nk[k]}: {'ok' if rank == nk[k] else 'MISMATCH'}); B spectrum gap "
          f"{s[~keep].max() if (~keep).any() else float('nan'):.1e} / {s[keep].min():.3e}; kernel of Q on singlets "
          f"{int((sig2 <= thr).sum())}; multiplets {len(lev)}, E0 = {lev[0]:.6g}", flush=True)
    st = r_stats(lev, a.trim)
    refs = {kd: reference_r(len(lev), kd, a.trim, samples=40 if len(lev) < 3000 else 6) for kd in ('loe', 'poisson')}
    print(f"  <r> = {st['mean_r']:.4f} +- {st['sem_r']:.4f} ({st['n_ratios']} ratios, {st['near_degenerate']} near-degenerate); "
          + ", ".join(f"{kd} {v['mean']:.4f}+-{v['std']:.4f}" for kd, v in refs.items()), flush=True)
    rec = dict(params=vars(a), pair=[k, k + 3], q=charge(k, a.n, a.p), w_sym=dim, n_k=int(nk[k]), rank=rank,
               projector_residual=leak_max, n_levels=int(len(lev)), E0=float(lev[0]), levels=lev.tolist(), stats=st,
               references=refs, unfolded_spacings=unfolded_spacings(lev, a.trim).tolist(),
               time_s=round(time.time() - t0, 1), peak_rss_gb=round(peak_gb(), 2))
    out = a.out or os.path.join(ROOT, 'results', 'data', f'quiver_singlet_pair_n{a.n}_p{a.p}_k{k}' + ('_f64' if a.dtype == 'float64' else '') + '.json')
    json.dump(rec, open(out, 'w'), indent=1)
    print(f"saved {os.path.relpath(out, ROOT)}  [{rec['time_s']}s, peak {rec['peak_rss_gb']} GB]")


if __name__ == '__main__':
    main()
