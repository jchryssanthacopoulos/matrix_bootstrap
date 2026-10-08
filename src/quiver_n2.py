"""
The U(2)^3 fermionic quiver (p flavours per edge) on its Weyl-symmetric zero-weight sectors.

Model (research/notes/quiver_project.md): A in (2, 2bar, 1), B in (1, 2, 2bar), C in (2bar, 1, 2) of U(2)^3, p flavours
per edge, all fermionic creation operators; Q = sum_{abc} C_abc Tr(A^a B^b C^c), H = {Q, Q^dag}.  Mode index
((e p + f) n + x) n + y for edge e (0: A, 1: B, 2: C), flavour f, fundamental index x (on node e), antifundamental
index y (on node e+1).  Fock basis states are uint64 occupation masks; c_m, c^dag_m carry (-1)^(occupied modes below m).

Reduction used throughout:
  * zero U(2)^3 weight (contains every gauge singlet), enumerated through occupation patterns of the 4 cells per edge;
  * the Weyl group (S_2)^3 (index swap 0 <-> 1 at each node) acts by signed mode permutations U_g, a genuine gauge
    action, so every gauge-invariant operator O satisfies U_g O U_g^dag = O and every singlet is W-invariant.
    Basis of the W-invariant subspace: normalised orbit sums b_o; orbits whose stabiliser acts with a sign -1 carry no
    invariant vector and are dropped.  O_sym = B^T O B.
  * on W-invariant zero-weight vectors, sum_v ||E^{(v)}_{12} psi||^2 = <psi| sum_v J_v^2 |psi> (J_- J_+ = J^2 on
    J_z = 0), so the penalty P = sum_v (E^{(v)}_12 B)^T (E^{(v)}_12 B) is the sum of the three nodes' SU(2) Casimirs: it
    commutes with H, vanishes exactly on singlets and is >= 2 on every non-singlet.

Code moved here from scripts/quiver_singlet_zero_modes.py (2026-10-07); behaviour unchanged.
"""
import itertools
from math import comb

import numpy as np
import scipy.sparse as sp

U64 = np.uint64


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
        # sign of U_g on a basis state: inversions = sum_{i occupied} popcount(mask & S_i), S_i = {j > i : sig(j) < sig(i)}
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
        masks = np.zeros(len(occ), dtype=U64)
        rows = np.arange(len(occ))
        for c, (e, x, y) in enumerate(cells):
            o = occ[rows, c]
            subsets = {r: [sum(1 << self.mode(e, f, x, y) for f in fl) for fl in itertools.combinations(range(p), r)]
                       for r in range(p + 1)}
            cnt = np.array([comb(p, int(r)) for r in range(p + 1)])[o]
            rows = np.repeat(rows, cnt); masks = np.repeat(masks, cnt); o2 = np.repeat(o, cnt)
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
            on = ((masks >> U64(i)) & U64(1)).astype(bool)
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

    def Q_terms(self):
        """Q as a list of (coefficient, [(mode, dagger), ...]) in operator order."""
        n, p = self.n, self.p
        out = []
        for a, b, c in itertools.product(range(p), repeat=3):
            for i, j, kk in itertools.product(range(n), repeat=3):
                out.append((float(self.C[a, b, c]),
                            [(self.mode(0, a, i, j), True), (self.mode(1, b, j, kk), True), (self.mode(2, c, kk, i), True)]))
        return out


def word_apply(masks, ops):
    """o_1 o_2 ... o_L on masks (o_L acts first), o = (mode, dagger): c^dag_mode if dagger else c_mode.
    Returns (alive, new masks, signs)."""
    m = masks.copy(); s = np.ones(len(masks)); alive = np.ones(len(masks), bool)
    for mo, dag in reversed(ops):
        bit = U64(1) << U64(mo)
        occ = (m & bit) != 0
        alive &= ~occ if dag else occ
        s *= 1.0 - 2.0 * (popcount(m & (bit - U64(1))) & 1)
        m = m ^ bit
    return alive, m, s


def chain_apply(masks, chain):
    """c^dag_{chain[0]} ... c^dag_{chain[-1]} on masks (last acts first): (alive, new masks, signs)."""
    return word_apply(masks, [(mo, True) for mo in chain])


def sym_basis(q, masks):
    best, bsign, orbit, dead = q.canon(masks)
    is_rep = (best == masks) & ~dead
    reps = masks[is_rep]                          # sorted (masks sorted)
    return dict(best=best, bsign=bsign, orbit=orbit, dead=dead, reps=reps, rep_orbit=orbit[is_rep])


def op_sym(src, sb, tgt_sb, terms):
    """B^T O B for a gauge-invariant O = sum coef * word, from the full source block `src` (with its symmetric data
    `sb`) to the symmetric sector `tgt_sb`:  O_sym[o', o] = sqrt|Orb(o')| sum_{s in orb o} b_o(s) O_{rep(o'), s},
    b_o(s) = sign(s -> rep) / sqrt|Orb(o)|."""
    alive_src = ~sb['dead']
    col = np.searchsorted(sb['reps'], sb['best'])
    bval = sb['bsign'] / np.sqrt(sb['orbit'].astype(np.float64))
    R, Cc, V = [], [], []
    treps = tgt_sb['reps']
    tw = np.sqrt(tgt_sb['rep_orbit'].astype(np.float64))
    for coef, ops in terms:
        ok, t, s = word_apply(src, ops)
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
        V.append(coef * s[idx] * bval[idx] * tw[pos])
    if not R:
        return sp.csr_matrix((len(treps), len(sb['reps'])))
    R, Cc, V = np.concatenate(R), np.concatenate(Cc), np.concatenate(V)
    return sp.csr_matrix((V, (R, Cc)), shape=(len(treps), len(sb['reps'])))


def Q_sym(q, src, sb, tgt_sb):
    """Q restricted to the W-symmetric zero-weight sectors (degree k -> k+3)."""
    return op_sym(src, sb, tgt_sb, q.Q_terms())


def raising_sym(q, src, sb, node, k):
    """E^{(node)}_{12} B : symmetric sector -> full block of weight (e_0 - e_1) at that node (touched rows only)."""
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


def sector(q, k):
    """Zero-weight block of degree k and its W-symmetric data (or None outside 0..nm)."""
    if k < 0 or k > q.nm:
        return None
    m = q.block(k, [0] * 6)
    return dict(masks=m, sb=sym_basis(q, m))


def flavour_eps_map(q):
    """The flavour map F (p = 2): F A^a F^-1 = sum_a' eps_{a'a} A^{a'} on every edge, eps = [[0, 1], [-1, 0]], so that
    F Q(C) F^-1 = Q(eps.eps.eps C).  As a signed mode permutation: flavour 0 -> 1 with sign -1, flavour 1 -> 0 with +1.
    Returns (sigma, mode_sign, inversion masks) for perm_sym."""
    assert q.p == 2
    n, p = q.n, q.p
    sig = np.zeros(q.nm, dtype=np.int64); sgn = np.zeros(q.nm)
    for e in range(3):
        for f in range(p):
            for x in range(n):
                for y in range(n):
                    sig[q.mode(e, f, x, y)] = q.mode(e, 1 - f, x, y)
                    sgn[q.mode(e, f, x, y)] = -1.0 if f == 0 else 1.0
    S = np.zeros(q.nm, dtype=U64)
    for i in range(q.nm):
        m = 0
        for j in range(i + 1, q.nm):
            if sig[j] < sig[i]:
                m |= 1 << j
        S[i] = U64(m)
    return sig, sgn, S


def perm_apply(masks, sig, sgn, inv_masks, nm):
    """A signed mode permutation c^dag_m -> sgn[m] c^dag_{sig[m]} on basis states: (new masks, signs)."""
    new = np.zeros_like(masks); inv = np.zeros(len(masks), dtype=np.int64); s = np.ones(len(masks))
    for i in range(nm):
        on = ((masks >> U64(i)) & U64(1)).astype(bool)
        new[on] |= U64(1) << U64(int(sig[i]))
        inv[on] += popcount(masks[on] & inv_masks[i])
        if sgn[i] < 0:
            s[on] *= -1.0
    return new, s * (1.0 - 2.0 * (inv & 1))


def perm_sym(src, sb, tgt_sb, sig, sgn, inv_masks, nm):
    """B^T U B for a signed mode permutation U commuting with the Weyl action (same projection rule as op_sym)."""
    alive_src = ~sb['dead']
    col = np.searchsorted(sb['reps'], sb['best'])
    bval = sb['bsign'] / np.sqrt(sb['orbit'].astype(np.float64))
    t, s = perm_apply(src, sig, sgn, inv_masks, nm)
    treps = tgt_sb['reps']
    idx = np.nonzero(alive_src)[0]
    tt = t[idx]
    pos = np.searchsorted(treps, tt)
    pos_c = np.minimum(pos, len(treps) - 1)
    hit = treps[pos_c] == tt
    idx, pos = idx[hit], pos[hit]
    V = s[idx] * bval[idx] * np.sqrt(tgt_sb['rep_orbit'][pos].astype(np.float64))
    return sp.csr_matrix((V, (pos.astype(np.int32), col[idx].astype(np.int32))), shape=(len(treps), len(sb['reps'])))
