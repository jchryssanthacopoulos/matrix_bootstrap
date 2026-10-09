#!/usr/bin/env python3
"""
Fortuity in the ARROW (flavour) direction at fixed rank n = 1 (research/notes/quiver_project.md section 13): the
fermionic counterpart of the Abelian quiver example of Chryssanthacopoulos-Vegh (2026), Discussion and Table 2, and
of their SYK uplift channels.

Theories: the n = 1 quiver with p flavours, Q = sum C_fgh a_f^+ b_g^+ c_h^+, with NESTED couplings C^(p) =
C^(pmax+1)[:p,:p,:p] drawn once.  Singlet blocks B_m (k = 3m), H^k(p) = BPS states (harmonic representatives).
Adding flavour p (0-based index p) gives two cochain maps:
  horizontal  pi : Lambda(V_{p+1}) -> Lambda(V_p), delete states containing a flavour-p mode (a quotient map);
              pi_* : H^k(p+1) -> H^k(p).  A class lifts horizontally iff it lies in im pi_*.
  diagonal    iota : B_m(p) -> B_{m+1}(p+1), psi -> a_p^+ b_p^+ c_p^+ psi (inclusion of the sub-complex of states
              with flavour p occupied on all three edges);  iota_* : H^k(p) -> H^{k+3}(p+1).
Both are verified to be cochain maps numerically.  The SYK statement "every class uplifts through at least one
channel" becomes  ker iota_*  subset of  im pi_*; we report dim H, rank pi_*, rank iota_*, and the dimension of the
classes lifting through neither channel, dim ker iota_* - dim(ker iota_* cap im pi_*).

    python scripts/quiver_n1_arrow_tower.py --pmax 6 --out results/data/quiver_n1_arrow_tower.json
"""
import argparse, json, os, resource, sys, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
from quiver_n1 import QuiverN1, subsets, dixon_index                          # noqa: E402


def bps_basis(qq, m, dense_max, nev_sparse=40):
    """Orthonormal basis (columns) of the singlet BPS space in block m, and the lowest nonzero level."""
    D = qq.dim(m)
    Qm = qq.Q_sparse(m) if m < qq.p else None
    Qp = qq.Q_sparse(m - 1) if m > 0 else None
    if D <= dense_max:
        H = np.zeros((D, D))
        if Qm is not None:
            H += (Qm.T @ Qm).toarray()            # sparse product first: Q itself can be much larger than D x D
        if Qp is not None:
            H += (Qp @ Qp.T).toarray()
        w, V = np.linalg.eigh(H)
        tol = 1e-9 * max(1.0, w.max())
        nb = int((w < tol).sum())
        return V[:, :nb], (float(w[nb]) if nb < D else None), 'dense'
    ops = []
    if Qm is not None:
        ops.append(lambda v, Q=Qm: Q.T @ (Q @ v))
    if Qp is not None:
        ops.append(lambda v, Q=Qp: Q @ (Q.T @ v))
    Hop = sla.LinearOperator((D, D), matvec=lambda v: sum(op(v) for op in ops), dtype=float)
    k = min(nev_sparse, D - 2)
    w, V = sla.eigsh(Hop, k=k, which='SA', tol=1e-12, ncv=max(3 * k, 120), maxiter=200000)
    o = np.argsort(w); w, V = w[o], V[:, o]
    scale = max(1.0, abs(w).max())
    nb = int((w < 1e-9 * scale).sum())
    if nb == k:
        raise RuntimeError(f'block p={qq.p} m={m}: all {k} computed eigenvalues are zero; raise nev_sparse')
    return V[:, :nb], float(w[nb]), 'lanczos'


def pi_matrix(p, m):
    """pi: B_m(p+1) -> B_m(p) (0/1 selection; states with flavour p present map to 0)."""
    big, small = subsets(p + 1, m), subsets(p, m)
    sidx = {s: i for i, s in enumerate(small)}
    Db, Ds = len(big), len(small)
    rows, cols = [], []
    for ia, Sa in enumerate(big):
        if p in Sa:
            continue
        for ib, Sb in enumerate(big):
            if p in Sb:
                continue
            for ic, Sc in enumerate(big):
                if p in Sc:
                    continue
                rows.append((sidx[Sa] * Ds + sidx[Sb]) * Ds + sidx[Sc]); cols.append((ia * Db + ib) * Db + ic)
    return sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(Ds ** 3, Db ** 3))


def iota_matrix(p, m):
    """iota: B_m(p) -> B_{m+1}(p+1), |S_a,S_b,S_c> -> a_p^+ b_p^+ c_p^+ |...> = +|S_a+p, S_b+p, S_c+p> (p is the largest
    flavour, so the reordering sign is (-1)^{3m} (-1)^{2m} (-1)^m = +1)."""
    small, big = subsets(p, m), subsets(p + 1, m + 1)
    bidx = {s: i for i, s in enumerate(big)}
    Ds, Db = len(small), len(big)
    rows, cols = [], []
    for ia, Sa in enumerate(small):
        ja = bidx[Sa + (p,)]
        for ib, Sb in enumerate(small):
            jb = bidx[Sb + (p,)]
            for ic, Sc in enumerate(small):
                rows.append((ja * Db + jb) * Db + bidx[Sc + (p,)]); cols.append((ia * Ds + ib) * Ds + ic)
    return sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(Db ** 3, Ds ** 3))


def rank(M, tol=1e-8):
    if M.size == 0:
        return 0
    s = np.linalg.svd(M, compute_uv=False)
    return int((s > tol * max(1.0, s.max())).sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--pmax', type=int, default=6, help='largest p whose classes are lifted (needs p+1 too)')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--dense-max', type=int, default=10000)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    t0 = time.time()
    P = a.pmax + 1
    Cbig = np.random.default_rng(a.seed).standard_normal((P, P, P))
    theories = {p: QuiverN1(p, Cbig[:p, :p, :p]) for p in range(1, P + 1)}
    bps, info = {}, {}
    for p in range(1, P + 1):
        for m in range(p + 1):
            B, gap, how = bps_basis(theories[p], m, a.dense_max)
            bps[(p, m)] = B
            info[f'{p},{m}'] = dict(dim=theories[p].dim(m), bps=int(B.shape[1]), lowest_nonzero=gap, method=how)
        counts = {m: bps[(p, m)].shape[1] for m in range(p + 1)}
        print(f'p={p}: BPS by block {counts}  (index {dixon_index(p)})  [{time.time()-t0:.0f}s]', flush=True)
    # cochain checks
    checks = {}
    for p in range(1, P):
        for m in range(p):
            Qs, Qb = theories[p].Q_sparse(m), theories[p + 1].Q_sparse(m)
            d_pi = abs(pi_matrix(p, m + 1) @ Qb - Qs @ pi_matrix(p, m)).max()
            iom, iom1 = iota_matrix(p, m), iota_matrix(p, m + 1)
            Qb1 = theories[p + 1].Q_sparse(m + 1)
            plus, minus = abs(Qb1 @ iom - iom1 @ Qs).max(), abs(Qb1 @ iom + iom1 @ Qs).max()
            checks[f'{p},{m}'] = dict(pi=float(d_pi), iota_commutes=float(plus), iota_anticommutes=float(minus))
    worst_pi = max(c['pi'] for c in checks.values())
    worst_iota = max(min(c['iota_commutes'], c['iota_anticommutes']) for c in checks.values())
    print(f'cochain checks: |pi Q - Q pi| <= {worst_pi:.1e};  iota (anti)commutes with Q to {worst_iota:.1e}', flush=True)
    # the two channels
    rows = []
    for p in range(1, a.pmax + 1):
        for m in range(p + 1):
            Bp = bps[(p, m)]
            h = Bp.shape[1]
            if h == 0:
                continue
            Bh = bps[(p + 1, m)]
            Mh = Bp.T @ (pi_matrix(p, m) @ Bh) if Bh.shape[1] else np.zeros((h, 0))      # im pi_* in the H(p) basis
            Bd = bps[(p + 1, m + 1)] if m + 1 <= p + 1 else np.zeros((0, 0))
            Md = Bd.T @ (iota_matrix(p, m) @ Bp) if Bd.shape[1] else np.zeros((0, h))   # iota_* : H(p) -> H(p+1)
            r_h, r_d = rank(Mh), rank(Md)
            # ker iota_* and im pi_* as subspaces of H^k(p) (orthonormal bases)
            if Md.shape[0]:
                _, s, Vt = np.linalg.svd(Md)
                rk = int((s > 1e-8 * max(1.0, s.max())).sum()) if len(s) else 0
                K = Vt[rk:].T
            else:
                K = np.eye(h)
            Ih = np.linalg.svd(Mh, full_matrices=False)[0][:, :r_h] if r_h else np.zeros((h, 0))
            inter = K.shape[1] + Ih.shape[1] - rank(np.hstack([K, Ih])) if (K.shape[1] and Ih.shape[1]) else 0
            neither = K.shape[1] - inter
            rows.append(dict(p=p, m=m, k=3 * m, h=h, h_next_same=int(Bh.shape[1]), h_next_shifted=int(Bd.shape[1]),
                             rank_pi=r_h, rank_iota=r_d, dim_ker_iota=int(K.shape[1]), lift_neither=int(neither)))
            print(f'p={p} k={3*m:2d}: dim H={h:5d} | horizontal: dim H^k(p+1)={Bh.shape[1]:5d}, rank pi_*={r_h:5d} | '
                  f'diagonal: dim H^(k+3)(p+1)={Bd.shape[1]:5d}, rank iota_*={r_d:5d} | lift through neither: {neither}', flush=True)
    rec = dict(params=vars(a), couplings='nested Gaussian, C^(p) = C[:p,:p,:p], C ~ N(0,1) iid, shape (pmax+1)^3',
               blocks=info, cochain_checks=checks, channels=rows, time_s=round(time.time() - t0, 1),
               peak_rss_gb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 2), date=time.strftime('%Y-%m-%d %H:%M'))
    if a.out:
        json.dump(rec, open(a.out, 'w'), indent=1)
        print('saved', a.out)


if __name__ == '__main__':
    main()
