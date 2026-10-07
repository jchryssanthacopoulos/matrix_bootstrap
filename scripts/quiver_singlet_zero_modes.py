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
from math import comb
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
U64 = np.uint64


def peak_gb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1e9 if sys.platform == 'darwin' else r / 1e6


def popcount(a):
    return np.bitwise_count(a).astype(np.int64)


class QuiverN2:
    def __init__(self, n, p, seed):
        assert n == 2, "Weyl reduction implemented for n = 2"
        self.n, self.p = n, p
        self.nm = 3 * p * n * n
        self.C = np.random.default_rng(seed).integers(1, 6, size=(p, p, p))   # as run_quiver_*.py
        # Weyl group: for each subset of nodes, swap index 0 <-> 1 at those nodes
        self.perms = []
        for sub in itertools.product([0, 1], repeat=3):
            sig = np.zeros(self.nm, dtype=np.int64)
            for e, f, x, y in itertools.product(range(3), range(p), range(n), range(n)):
                x2 = 1 - x if sub[e] else x
                y2 = 1 - y if sub[(e + 1) % 3] else y
                sig[self.mode(e, f, x, y)] = self.mode(e, f, x2, y2)
            self.perms.append(sig)
        # for the sign of U_g on a basis state: inversions = sum_{i occupied} popcount(mask & S_i),
        # S_i = {j > i : sigma(j) < sigma(i)}
        self.inv_masks = []
        for sig in self.perms:
            S = np.zeros(self.nm, dtype=U64)
            for i in range(self.nm):
                m = 0
                for j in range(i + 1, self.nm):
                    if sig[j] < sig[i]:
                        m |= 1 << j
                S[i] = U64(m)
            self.inv_masks.append(S)

    def mode(self, e, f, x, y):
        return ((e * self.p + f) * self.n + x) * self.n + y

    def weight_of_mode(self, e, x, y):
        w = [0] * (3 * self.n); w[e * self.n + x] += 1; w[((e + 1) % 3) * self.n + y] -= 1
        return w

    def block(self, k, w):
        """Sorted masks of degree k and U(2)^3 weight w (length 6)."""
        n, p = self.n, self.p
        cells = [(e, x, y) for e in range(3) for x in range(n) for y in range(n)]
        nc = len(cells)
        codes = np.arange((p + 1) ** nc, dtype=np.int64)
        occ = np.zeros((len(codes), nc), dtype=np.int8)
        rem = codes.copy()
        for c in range(nc):
            occ[:, c] = rem % (p + 1); rem //= (p + 1)
        del codes, rem
        occ = occ[occ.sum(1) == k]                                          # degree filter first (memory)
        W = np.array([self.weight_of_mode(*cl) for cl in cells], dtype=np.int8)   # (nc, 6)
        occ = occ[np.all(occ.astype(np.int16) @ W.astype(np.int16) == np.asarray(w, dtype=np.int16)[None, :], axis=1)]
        # expand flavour subsets cell by cell
        masks = np.zeros(len(occ), dtype=U64)
        rows = np.arange(len(occ))
        for c, (e, x, y) in enumerate(cells):
            o = occ[rows, c]
            subsets = {r: [sum(1 << self.mode(e, f, x, y) for f in fl) for fl in itertools.combinations(range(p), r)]
                       for r in range(p + 1)}
            cnt = np.array([comb(p, int(r)) for r in range(p + 1)])[o]
            rows = np.repeat(rows, cnt); masks = np.repeat(masks, cnt); o2 = np.repeat(o, cnt)
            # position within the repeat group
            starts = np.repeat(np.cumsum(cnt) - cnt, cnt)
            pos = np.arange(len(rows)) - starts
            add = np.zeros(len(rows), dtype=U64)
            for r in range(p + 1):
                mk = o2 == r
                if mk.any():
                    add[mk] = np.array(subsets[r], dtype=U64)[pos[mk]]
            masks |= add
        return np.sort(masks)

    def act(self, g, masks):
        """U_g on basis states: (new masks, signs)."""
        sig, S = self.perms[g], self.inv_masks[g]
        new = np.zeros_like(masks); inv = np.zeros(len(masks), dtype=np.int64)
        for i in range(self.nm):
            occ = (masks >> U64(i)) & U64(1)
            on = occ.astype(bool)
            new[on] |= U64(1) << U64(int(sig[i]))
            inv[on] += popcount(masks[on] & S[i])
        return new, 1.0 - 2.0 * (inv & 1)

    def canon(self, masks):
        """Orbit representative (minimum mask), sign of U_g e_s -> e_rep, orbit size, and 'dead' flag (stabiliser
        acting with a sign -1, so the orbit carries no invariant vector).  Two passes, no stored images (memory)."""
        best = masks.copy(); bsign = np.ones(len(masks), dtype=np.float32)
        for g in range(len(self.perms)):
            m2, s2 = self.act(g, masks)
            better = m2 < best
            best[better] = m2[better]; bsign[better] = s2[better]
        osize = np.zeros(len(masks), dtype=np.int8); dead = np.zeros(len(masks), dtype=bool)
        for g in range(len(self.perms)):
            m2, s2 = self.act(g, masks)
            hit = m2 == best
            osize += hit.astype(np.int8)
            dead |= hit & (s2 != bsign)
        orbit = (len(self.perms) // osize).astype(np.int8)
        return best, bsign, orbit, dead


def chain_apply(masks, chain):
    """c^dag_{chain[0]} ... c^dag_{chain[-1]} on masks (last acts first): (alive, new masks, signs)."""
    m = masks.copy(); s = np.ones(len(masks)); alive = np.ones(len(masks), bool)
    for mo in reversed(chain):
        bit = U64(1) << U64(mo)
        alive &= (m & bit) == 0
        s *= 1.0 - 2.0 * (popcount(m & (bit - U64(1))) & 1)
        m = m | bit
    return alive, m, s


def sym_basis(q, masks):
    best, bsign, orbit, dead = q.canon(masks)
    is_rep = (best == masks) & ~dead
    reps = masks[is_rep]                          # sorted (masks sorted)
    return dict(best=best, bsign=bsign, orbit=orbit, dead=dead, reps=reps, rep_orbit=orbit[is_rep])


def Q_sym(q, src, sb, tgt_sb):
    """Q_sym[o', o] = sqrt|Orb(o')| * sum_{s in orb o} b_o(s) Q_{rep(o'), s}, b_o(s) = sign(s->rep)/sqrt|Orb(o)|."""
    n, p = q.n, q.p
    alive_src = ~sb['dead']
    col = np.searchsorted(sb['reps'], sb['best'])
    bval = sb['bsign'] / np.sqrt(sb['orbit'].astype(np.float64))
    R, Cc, V = [], [], []
    treps = tgt_sb['reps']
    for a, b, c in itertools.product(range(p), repeat=3):
        for i, j, kk in itertools.product(range(n), repeat=3):
            ok, t, s = chain_apply(src, [q.mode(0, a, i, j), q.mode(1, b, j, kk), q.mode(2, c, kk, i)])
            ok &= alive_src
            idx = np.nonzero(ok)[0]
            if len(idx) == 0:
                continue
            tt = t[idx]
            pos = np.searchsorted(treps, tt)
            pos_c = np.minimum(pos, len(treps) - 1)
            hit = treps[pos_c] == tt
            idx, pos = idx[hit], pos[hit]
            R.append(pos.astype(np.int32)); Cc.append(col[idx].astype(np.int32))
            V.append((q.C[a, b, c] * s[idx] * bval[idx] * np.sqrt(tgt_sb['rep_orbit'][pos].astype(np.float64))).astype(np.float64))
    R, Cc, V = np.concatenate(R), np.concatenate(Cc), np.concatenate(V)
    return sp.csr_matrix((V, (R, Cc)), shape=(len(treps), len(sb['reps'])))


def raising_sym(q, src, sb, node, k):
    """E^{(node)}_{12} B : symmetric sector -> full block of weight (e_0 - e_1) at that node."""
    n, p = q.n, q.p
    w = [0] * 6; w[node * n + 0] += 1; w[node * n + 1] -= 1
    tgt = q.block(k, w)
    alive_src = ~sb['dead']
    col = np.searchsorted(sb['reps'], sb['best'])
    bval = sb['bsign'] / np.sqrt(sb['orbit'].astype(np.float64))
    R, Cc, V = [], [], []
    def hop(m_to, m_from, sgn):
        bf, bt = U64(1) << U64(m_from), U64(1) << U64(m_to)
        ok = ((src & bf) != 0) & alive_src
        idx = np.nonzero(ok)[0]
        m = src[idx]
        s = 1.0 - 2.0 * (popcount(m & (bf - U64(1))) & 1)
        m = m ^ bf
        ok2 = (m & bt) == 0
        idx, m, s = idx[ok2], m[ok2], s[ok2]
        s = s * (1.0 - 2.0 * (popcount(m & (bt - U64(1))) & 1)); m = m | bt
        r = np.searchsorted(tgt, m)
        assert np.all(tgt[np.minimum(r, len(tgt) - 1)] == m)
        R.append(r.astype(np.int32)); Cc.append(col[idx].astype(np.int32)); V.append(sgn * s * bval[idx])
    for e in range(3):
        for f in range(p):
            for z in range(n):
                if e == node:                       # fundamental index on node e: index 1 -> 0
                    hop(q.mode(e, f, 0, z), q.mode(e, f, 1, z), 1.0)
                if (e + 1) % 3 == node:             # antifundamental index: E_kl acts with - c^dag_(.,l) c_(.,k)
                    hop(q.mode(e, f, z, 1), q.mode(e, f, z, 0), -1.0)
    R, Cc, V = np.concatenate(R), np.concatenate(Cc), np.concatenate(V)
    M = sp.csr_matrix((V, (R, Cc)), shape=(len(tgt), len(sb['reps'])))
    return M[np.unique(R)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--n', type=int, default=2); ap.add_argument('--p', type=int, default=3)
    ap.add_argument('--k', type=int, required=True); ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('--mu', type=float, default=1.0); ap.add_argument('--nev', type=int, default=4)
    ap.add_argument('--dense-max', type=int, default=4000)
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
                               v0=np.random.default_rng(1).standard_normal(dim))
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
