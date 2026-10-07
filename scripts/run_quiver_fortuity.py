#!/usr/bin/env python3
"""
Fortuity test for the singlet BPS states of the U(n)^3 quiver Q = sum_{abc} C_abc Tr(A^a B^b C^c) at (n, p) = (2, 2)
(the 90 concentrated singlets at k = 12 of research/notes/cohomology_results.md section 4x; same couplings, seed 3).

Definition (Tierz 2026, adopted in literature/notes/tierz_2026.md): the projection pi_{3->2}, deleting every
monomial that contains a third index on any node, is a cochain map, pi Q^(3) = Q^(2) pi; a rank-2 class is monotone
only if it lies in Im Pi_{3->2} (and in the images from every larger rank), fortuitous otherwise.  Dually, the
inclusion iota: H_2 -> H_3 is a cochain map for Q^dag, and iota_*[z] != 0 in Q^dag-cohomology iff [z] is not
orthogonal to Im Pi (harmonic representatives), so one computation covers the Q- and Q^dag-versions of CCSY's
definition.

Leading-order obstruction (derived in the accompanying note).  Write Q^(3) = Q^(2) + Q_new and z = z_0 + z_+ with z_+
the part carrying index-3 modes.  A closed lift exists iff Q_new z_0 is Q^(3)-exact inside the index-3 subcomplex.
Grade that subcomplex by the number m of index-3 modes: Q^(2) preserves m and Q_new raises it by 2 (one index equal
to 3) or 3 (two or three indices equal to 3), never by 1.  At m = 2 the condition is therefore Q^(2) z_2 = -Q_new^(1) z_0,
whose components are single rank-2 letters acting on z_0:
    X z_0 must be Q^(2)-exact for every letter X in {A^a_ij, B^b_jk, C^c_ki}   (generic C_abc),
i.e. the "BPS three-point functions" O_{beta, z} = <beta| X |z> must vanish for every BPS state beta at k = 13 in the
weight block of X.  Any class z with O z != 0 cannot come from rank 3: it is fortuitous.  (Vanishing of O is only
necessary for liftability; higher orders m >= 3 are not tested here.)  By gauge covariance it suffices to test
X_{11} for each edge and flavour.

    scripts/run_guarded.sh 6 log python3 scripts/run_quiver_fortuity.py --dry-run
    scripts/run_guarded.sh 8 log python3 scripts/run_quiver_fortuity.py
"""
import argparse, itertools, json, os, resource, sys, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')


def peak_gb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1e9 if sys.platform == 'darwin' else r / 1e6


def parity_below(masks, m):
    return 1.0 - 2.0 * (np.bitwise_count(masks & np.uint64((1 << m) - 1)) & 1)


class Quiver:
    def __init__(self, n, P, seed):
        self.n, self.P = n, P
        self.nm = 3 * P * n * n
        rng = np.random.default_rng(seed)
        self.C = rng.integers(1, 6, size=(P, P, P))     # identical to run_quiver_irrep_cohomology.py
        W = np.zeros((self.nm, 3 * n), dtype=np.int8)
        for e, f, x, y in itertools.product(range(3), range(P), range(n), range(n)):
            W[self.mode(e, f, x, y), e * n + x] += 1
            W[self.mode(e, f, x, y), ((e + 1) % 3) * n + y] -= 1
        self.W = W
        # all masks with their (k, weight) -- 2^24 at n = P = 2
        masks = np.arange(1 << self.nm, dtype=np.uint64)
        wv = np.zeros((1 << self.nm, 3 * n), dtype=np.int8)
        for i in range(self.nm):
            bit = ((masks >> np.uint64(i)) & np.uint64(1)).astype(np.int8)
            wv += bit[:, None] * W[i][None, :]
        self.masks, self.wv = masks, wv
        self.kk = np.bitwise_count(masks).astype(np.int8)

    def mode(self, e, f, x, y):
        return ((e * self.P + f) * self.n + x) * self.n + y

    def block(self, k, w):
        sel = (self.kk == k) & np.all(self.wv == np.asarray(w, dtype=np.int8)[None, :], axis=1)
        return self.masks[sel]                          # sorted

    def creation_chain(self, src, tgt, chain):
        """Sparse matrix of the product of creation operators c^dag_{chain[0]} ... c^dag_{chain[-1]} (last acts first)."""
        m = src.copy(); val = np.ones(len(src)); alive = np.ones(len(src), bool)
        for mo in reversed(chain):
            bit = np.uint64(1 << mo)
            alive &= (m & bit) == 0
            val = val * parity_below(m, mo)
            m = m | bit
        cols = np.nonzero(alive)[0]
        rows = np.searchsorted(tgt, m[cols])
        assert np.all(tgt[np.minimum(rows, len(tgt) - 1)] == m[cols])
        return sp.csr_matrix((val[cols], (rows, cols)), shape=(len(tgt), len(src)))

    def Q(self, src, tgt):
        n, P = self.n, self.P
        out = sp.csr_matrix((len(tgt), len(src)))
        for a, b, c in itertools.product(range(P), repeat=3):
            for i, j, k in itertools.product(range(n), repeat=3):
                out = out + self.C[a, b, c] * self.creation_chain(
                    src, tgt, [self.mode(0, a, i, j), self.mode(1, b, j, k), self.mode(2, c, k, i)])
        return out.tocsr()

    def lowering(self, src, tgt, node, l, k_):
        """Off-diagonal gauge generator E^{(node)}_{l k} (l != k) restricted to src -> tgt (one-body hops)."""
        n, P = self.n, self.P
        rows, cols, vals = [], [], []
        def hop(m_to, m_from, s):
            bf, bt = np.uint64(1 << m_from), np.uint64(1 << m_to)
            ok = ((src & bf) != 0)
            m = src[ok]; c = np.nonzero(ok)[0]
            sg = parity_below(m, m_from); m = m ^ bf
            ok2 = (m & bt) == 0
            m, c, sg = m[ok2], c[ok2], sg[ok2]
            sg = sg * parity_below(m, m_to); m = m | bt
            r = np.searchsorted(tgt, m)
            assert np.all(tgt[np.minimum(r, len(tgt) - 1)] == m)
            rows.append(r); cols.append(c); vals.append(s * sg)
        for e in range(3):
            for f in range(P):
                for z in range(n):
                    if e == node:                         # fundamental index of edge e lives on node e
                        hop(self.mode(e, f, l, z), self.mode(e, f, k_, z), 1.0)
                    if (e + 1) % 3 == node:               # antifundamental index of edge e lives on node e+1
                        hop(self.mode(e, f, z, k_), self.mode(e, f, z, l), -1.0)
        return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(len(tgt), len(src)))


def null_space(S, dim, expect=None, tol=1e-8):
    """Orthonormal basis of ker S (S sparse, dim columns): dense eigh of S^T S if expect is None (dim <= ~12000,
    8 dim^2 bytes per copy), else shift-invert Lanczos for an expected small kernel."""
    if expect is None:
        G = (S.T @ S).toarray()
        vals, vecs = np.linalg.eigh(G)
        del G
        sel = vals < tol * max(1.0, vals[-1])
        nz = int(sel.sum())
        return vecs[:, sel].copy(), (float(vals[nz - 1]) if nz else float('nan'), float(vals[nz]))
    G = (S.T @ S).tocsc()
    nev = min(dim - 2, (expect or 50) + 20)
    vals, vecs = sla.eigsh(G + 0 * sp.identity(dim), k=nev, sigma=-1e-6, which='LM')
    order = np.argsort(vals); vals, vecs = vals[order], vecs[:, order]
    sel = vals < tol * max(1.0, abs(G).max())
    if sel.all():
        raise RuntimeError("null space may be larger than nev; increase expect")
    Z, _ = np.linalg.qr(vecs[:, sel])
    nz = int(sel.sum())
    gap = (float(vals[nz - 1]) if nz else float('nan'), float(vals[nz]))
    return Z, gap


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--n', type=int, default=2); ap.add_argument('--p', type=int, default=2)
    ap.add_argument('--seed', type=int, default=3); ap.add_argument('--k', type=int, default=12)
    ap.add_argument('--dry-run', action='store_true'); ap.add_argument('--out', default=None)
    a = ap.parse_args()
    t0 = time.time()
    q = Quiver(a.n, a.p, a.seed)
    n, P, k = a.n, a.p, a.k
    zero = (0,) * (3 * n)
    letters = []
    for e, f in itertools.product(range(3), range(P)):
        m = q.mode(e, f, 0, 0)
        letters.append((('A', 'B', 'C')[e] + f"^{f}_11", m, tuple(int(x) for x in q.W[m])))
    blocks = {(kk, zero): q.block(kk, zero) for kk in (k - 3, k, k + 3)}
    for _, m, w in letters:
        for kk in (k - 2, k + 1, k + 4):
            blocks.setdefault((kk, w), q.block(kk, w))
    print(f"(n,p)=({n},{P}), seed {a.seed}, couplings C={q.C.tolist()}; enumeration {time.time()-t0:.0f}s, "
          f"{peak_gb():.2f} GB", flush=True)
    for key, b in sorted(blocks.items()):
        print(f"  block k={key[0]:2d} w={key[1]}: {len(b)} states", flush=True)
    big = max(len(b) for b in blocks.values())
    print(f"  largest block {big}; shift-invert on S^T S needs a sparse LU of that size (est < 2 GB)", flush=True)
    if a.dry_run:
        return
    del q.wv
    rec = dict(params=vars(a), couplings=q.C.tolist(), blocks={f"{kk}|{w}": len(b) for (kk, w), b in blocks.items()})

    # --- singlet BPS states at degree k, zero weight: ker of [Q_k; Q_{k-3}^T; all lowering operators]
    B9, B12, B15 = blocks[(k - 3, zero)], blocks[(k, zero)], blocks[(k + 3, zero)]
    rows = [q.Q(B12, B15), q.Q(B9, B12).T.tocsr()]
    for node in range(3):
        for l, kk_ in itertools.permutations(range(n), 2):
            w = list(zero); w[node * n + l] += 1; w[node * n + kk_] -= 1
            tgt = q.block(k, w) if False else None
            # targets: the weight block reached by E_{l k}; enumerate lazily from all masks of degree k
            sel = (q.kk == k)
            tm = q.masks[sel]
            # restrict to masks reachable: build with full degree-k list, then drop empty rows
            E = q.lowering(B12, tm, node, l, kk_)
            E = E[np.unique(E.nonzero()[0])]
            rows.append(E.tocsr())
    S = sp.vstack(rows).tocsr()
    Zs, gap = null_space(S, len(B12), expect=100)
    print(f"singlet BPS states at k={k}: {Zs.shape[1]} (expected 90); last zero / first nonzero eigenvalue of S^T S "
          f"{gap[0]:.1e} / {gap[1]:.3e}  [{time.time()-t0:.0f}s, {peak_gb():.2f} GB]",
          flush=True)
    rec['n_singlet_bps'] = int(Zs.shape[1])

    # --- obstruction matrices O = beta^T X z for each letter X_11
    obstructed = []
    Ostack = []
    for lab, m, w in letters:
        B10, B13, B16 = blocks[(k - 2, w)], blocks[(k + 1, w)], blocks[(k + 4, w)]
        S13 = sp.vstack([q.Q(B13, B16), q.Q(B10, B13).T.tocsr()]).tocsr()
        Zb, gapb = null_space(S13, len(B13))
        X = q.creation_chain(B12, B13, [m])
        # sanity: X z is Q-closed
        Xz = X @ Zs
        closed = float(np.abs(q.Q(B13, B16) @ Xz).max()) if B16.size else 0.0
        O = Zb.T @ Xz
        s = np.linalg.svd(O, compute_uv=False) if O.size else np.zeros(0)
        rank = int((s > 1e-8 * max(1.0, s.max() if s.size else 1.0)).sum())
        print(f"  letter {lab}: BPS states at k={k+1} in its block: {Zb.shape[1]}; |Q X z| = {closed:.1e}; "
              f"gap {gapb[0]:.1e}/{gapb[1]:.2e}; rank of <beta|X|z> = {rank} / {Zs.shape[1]}; top singular values {np.round(s[:4], 6).tolist()}",
              flush=True)
        obstructed.append(dict(letter=lab, n_bps_target=int(Zb.shape[1]), rank=rank, closed_residual=closed,
                               singular_values=[float(x) for x in s[:10]]))
        Ostack.append(O)
    Oall = np.vstack(Ostack)
    s_all = np.linalg.svd(Oall, compute_uv=False)
    r_all = int((s_all > 1e-8 * max(1.0, s_all.max())).sum())
    print(f"combined: {r_all} of {Zs.shape[1]} singlet classes are obstructed at leading order "
          f"(=> fortuitous); {Zs.shape[1] - r_all} pass the leading-order test  [{time.time()-t0:.0f}s, "
          f"{peak_gb():.2f} GB]", flush=True)
    rec.update(letters=obstructed, combined_rank=r_all, n_unobstructed=int(Zs.shape[1] - r_all),
               time_s=round(time.time() - t0, 1), peak_rss_gb=round(peak_gb(), 2))
    out = a.out or os.path.join(ROOT, 'results', 'data', f'quiver_fortuity_n{n}_p{P}_k{k}.json')
    json.dump(rec, open(out, 'w'), indent=1)
    print('saved', os.path.relpath(out, ROOT))


if __name__ == '__main__':
    main()
