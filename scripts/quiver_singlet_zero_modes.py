#!/usr/bin/env python3
"""
Singlet BPS zero modes of the U(n)^3 fermionic quiver in one degree k (n = 2), by Lanczos on the Weyl-symmetric
zero-weight sector.  Purpose: test whether the singlet BPS states concentrate in a single degree (research/notes/
quiver_project.md section 4): a U(n)^3 singlet has N_A = N_B = N_C, so all singlets sit at k = 3m, the singlet index
I_0 = sum_m (-1)^m n(3m) lower-bounds the singlet BPS count, with equality iff the singlet cohomology lives in one degree.

Reduction.
  * Zero U(2)^3 weight (contains every singlet).  Masks enumerated through occupation patterns of the 4 cells per
    edge (flavour subsets expanded vectorially).
  * Weyl group W = (S_2)^3 (index swap 0 <-> 1 at each node) acts by signed mode permutations U_g, a genuine gauge
    action, so U_g Q U_g^dag = Q and every singlet is W-invariant.  Basis: normalised orbit sums b_o (orbits whose
    stabiliser acts with a sign -1 carry no invariant vector and are dropped).  Q_sym = B^T Q B, H_sym =
    Q_sym^T Q_sym + Q_sym' Q_sym'^T.
  * Singlet penalty: for a W-invariant zero-weight vector, psi is a singlet iff E^{(v)}_{12} psi = 0 for all three
    nodes v (J_z = 0 and J_+ psi = 0 => spin 0; E_21 psi = 0 then follows by W-invariance).  So
        K = H_sym + mu sum_v (E^{(v)}_{12} B)^T (E^{(v)}_{12} B) >= 0,  and  ker K = singlet BPS states in degree k.
Output: the lowest eigenvalues of K with residuals.  lambda_min clearly > 0 => no singlet BPS state in degree k
(numerical, Lanczos); a zero eigenvalue => singlet BPS states exist there.

    scripts/run_guarded.sh 9 log python3 scripts/quiver_singlet_zero_modes.py --n 2 --p 3 --k 15
"""
import argparse, itertools, json, os, resource, sys, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
U64 = np.uint64


sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
from quiver_n2 import QuiverN2, sym_basis, Q_sym, raising_sym, U64      # noqa: E402  (moved to src/, 2026-10-07)


def peak_gb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1e9 if sys.platform == 'darwin' else r / 1e6


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--n', type=int, default=2); ap.add_argument('--p', type=int, default=3)
    ap.add_argument('--k', type=int, required=True); ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('--mu', type=float, default=1.0); ap.add_argument('--nev', type=int, default=4)
    ap.add_argument('--dense-max', type=int, default=4000)
    ap.add_argument('--v0-seed', type=int, default=1, help='Lanczos start vector seed (robustness checks)')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    t0 = time.time()
    q = QuiverN2(a.n, a.p, a.seed)
    k, zero = a.k, [0] * 6
    blk = {kk: q.block(kk, zero) for kk in (k - 3, k, k + 3) if 0 <= kk <= q.nm}
    sbs = {kk: sym_basis(q, m) for kk, m in blk.items()}
    print(f"(n,p)=({a.n},{a.p}) k={k} seed {a.seed}: zero-weight blocks " +
          ", ".join(f"k={kk}: {len(m)} -> W-sym {len(sbs[kk]['reps'])}" for kk, m in blk.items()) +
          f"  [{time.time()-t0:.0f}s, {peak_gb():.2f} GB]", flush=True)
    ops = []
    if k + 3 in blk:
        Qu = Q_sym(q, blk[k], sbs[k], sbs[k + 3]); ops.append(('up', Qu))
    if k - 3 in blk:
        Qd = Q_sym(q, blk[k - 3], sbs[k - 3], sbs[k]); ops.append(('down', Qd))
    Ws = [raising_sym(q, blk[k], sbs[k], v, k) for v in range(3)]
    print(f"  operators built: " + ", ".join(f"{nm} nnz {M.nnz}" for nm, M in ops) +
          f", raising nnz {[W.nnz for W in Ws]}  [{time.time()-t0:.0f}s, {peak_gb():.2f} GB]", flush=True)
    dim = len(sbs[k]['reps'])
    def Kmv(x):
        y = np.zeros_like(x)
        for nm, M in ops:
            y += (M.T @ (M @ x)) if nm == 'up' else (M @ (M.T @ x))
        for W in Ws:
            y += a.mu * (W.T @ (W @ x))
        return y
    rec = dict(params=vars(a), dims={str(kk): [int(len(m)), int(len(sbs[kk]['reps']))] for kk, m in blk.items()})
    if dim <= a.dense_max:
        Kd = np.column_stack([Kmv(e) for e in np.eye(dim)])
        vals = np.linalg.eigvalsh((Kd + Kd.T) / 2)
        nz = int((vals < 1e-9).sum())
        print(f"  dense: {nz} zero modes (singlet BPS); lowest eigenvalues {np.round(vals[:max(a.nev, nz + 2)], 9).tolist()}")
        rec.update(method='dense', zero_modes=nz, lowest=[float(v) for v in vals[:a.nev + nz]])
    else:
        op = sla.LinearOperator((dim, dim), matvec=Kmv, dtype=float)
        vals, vecs = sla.eigsh(op, k=a.nev, which='SA', tol=1e-10, ncv=max(4 * a.nev, 40),
                               v0=np.random.default_rng(a.v0_seed).standard_normal(dim))
        order = np.argsort(vals); vals, vecs = vals[order], vecs[:, order]
        res = [float(np.linalg.norm(Kmv(vecs[:, i]) - vals[i] * vecs[:, i])) for i in range(len(vals))]
        print(f"  Lanczos: lowest eigenvalues {np.round(vals, 9).tolist()}, residuals {[f'{r:.1e}' for r in res]}")
        rec.update(method='eigsh', lowest=[float(v) for v in vals], residuals=res)
    rec.update(time_s=round(time.time() - t0, 1), peak_rss_gb=round(peak_gb(), 2))
    print(f"  [{rec['time_s']}s, peak {rec['peak_rss_gb']} GB]")
    if a.out:
        json.dump(rec, open(a.out, 'w'), indent=1); print('  saved', a.out)


if __name__ == '__main__':
    main()
