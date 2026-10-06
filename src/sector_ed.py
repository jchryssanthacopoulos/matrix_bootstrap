"""
Symmetry-resolved exact diagonalisation of Chen's three-matrix model in a fixed sector k = N_Psi and a fixed U(N)
weight block.  The default block is zero weight: every irrep occurring here has the zero weight (all weights lie in
the root lattice), so the lowest eigenvalue of the zero-weight block is the ground energy of the whole sector, and a
level of irrep lambda appears with multiplicity K_{lambda,0} (its zero-weight multiplicity) per copy.

H = E_c + H_1 + 9 Ydag^T Ydag  (D2.13 written as in D16.2), with two interchangeable builders:
  * 'python' : free_sectors.sector_operators on the tuple basis of cohomology.weight_basis (reference);
  * 'fast'   : the same operators from vectorised bit operations on int64 occupation masks.
Irrep labels come from <C_2> = (1/2) sum_{i != i'} |E_{i'i} v|^2 + (1/2)|mu|^2 with the one-body gauge generators
E_{kl} = sum_{a,j} c^dag_(a,k,j) c_(a,l,j) - sum_{a,i} c^dag_(a,i,l) c_(a,i,k) (D7.1 normalisation: adjoint c = N).

Conventions as in cohomology.py and free_sectors.py: mode m = a N^2 + i N + j, Psi^a_ij = c^dag_m, and c_m, c^dag_m
carry the sign (-1)^(number of occupied modes below m).
"""
import numpy as np
import scipy.sparse as sp

from cohomology import mode_index, kostant
import free_sectors as fs

ONE = np.int64(1)


# ----------------------------------------------------------------------------------------------- masks

def masks_from_tuples(basis):
    """int64 occupation masks of sorted mode tuples (same order as the input)."""
    arr = np.asarray(basis, dtype=np.int64)
    if arr.ndim == 1:
        arr = arr.reshape(len(basis), -1)
    return np.bitwise_or.reduce(np.left_shift(ONE, arr), axis=1) if arr.shape[1] else np.zeros(len(basis), np.int64)


def tuple_from_mask(m):
    m = int(m)
    out, b = [], 0
    while m:
        if m & 1:
            out.append(b)
        m >>= 1
        b += 1
    return tuple(out)


def _parity_below(masks, m):
    """(-1)^(number of set bits of each mask below bit m), as float64."""
    cnt = np.bitwise_count(masks & ((ONE << m) - ONE))
    return 1.0 - 2.0 * (cnt & 1)


def apply_hop(masks, alpha, beta):
    """c^dag_alpha c_beta on every mask.  Returns (source positions, new masks, signs) where it acts."""
    bb = ONE << beta
    src = np.nonzero(masks & bb)[0]
    m = masks[src]
    if alpha == beta:
        return src, m, np.ones(len(src))
    s1 = _parity_below(m, beta)
    m = m ^ bb
    aa = ONE << alpha
    ok = (m & aa) == 0
    src, m, s1 = src[ok], m[ok], s1[ok]
    return src, m | aa, s1 * _parity_below(m, alpha)


def apply_pair_annihilation(masks, m1, m2):
    """c_{m2} c_{m1} (m1 removed first) on every mask, m1 != m2.  Returns (source positions, new masks, signs)."""
    b1, b2 = ONE << m1, ONE << m2
    src = np.nonzero(((masks & b1) != 0) & ((masks & b2) != 0))[0]
    m = masks[src]
    s1 = _parity_below(m, m1)
    m = m ^ b1
    s2 = _parity_below(m, m2)
    return src, m ^ b2, s1 * s2


# ----------------------------------------------------------------------------------------------- fast builder

def build_fast(N, C, masks):
    """H_1 and Ydag of D16.2 on the block spanned by `masks` (sorted, closed under H).  Same operators as
    free_sectors.sector_operators, in the order of `masks`."""
    p = C.shape[0]
    dim = len(masks)
    assert np.all(np.diff(masks) > 0), "masks must be sorted and unique"
    M, L, _, _ = fs.flavour_tensors(C)
    h = {}
    for c in range(p):
        for d in range(p):
            for i in range(N):
                for j in range(N):
                    key = (mode_index(N, c, i, j), mode_index(N, d, i, j))
                    h[key] = h.get(key, 0.0) - 9 * N * M[c, d].real
                    key = (mode_index(N, c, i, i), mode_index(N, d, j, j))
                    h[key] = h.get(key, 0.0) + 9 * L[c, d].real
    rows, cols, vals = [], [], []
    for (al, be), coef in h.items():
        if coef == 0:
            continue
        src, new, sg = apply_hop(masks, al, be)
        if len(src) == 0:
            continue
        tgt = np.searchsorted(masks, new)
        tgt_c = np.minimum(tgt, dim - 1)
        assert np.all(masks[tgt_c] == new), "H_1 left the block"
        rows.append(tgt.astype(np.int32)); cols.append(src.astype(np.int32)); vals.append(coef * sg)
    H1 = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(dim, dim))
    del rows, cols, vals

    # Ydag: one row block per (a, i, kk), converted to CSR as soon as it is built so that the transient memory is
    # one block rather than the whole matrix; rows of a block are the distinct (k-2)-particle targets.
    Cc = np.conj(C).real
    blocks = []
    for a in range(p):
        for i in range(N):
            for kk in range(N):
                tm, ts, tv = [], [], []
                for b in range(p):
                    for c in range(p):
                        coef = Cc[a, b, c]
                        if coef == 0:
                            continue
                        for j in range(N):
                            m1, m2 = mode_index(N, b, i, j), mode_index(N, c, j, kk)
                            if m1 == m2:
                                continue
                            src, new, sg = apply_pair_annihilation(masks, m1, m2)
                            if len(src):
                                tm.append(new); ts.append(src.astype(np.int32)); tv.append(coef * sg)
                if not tm:
                    continue
                tm = np.concatenate(tm)
                uniq, inv = np.unique(tm, return_inverse=True)
                del tm
                blocks.append(sp.csr_matrix((np.concatenate(tv), (inv.astype(np.int32), np.concatenate(ts))),
                                            shape=(len(uniq), dim)))
                del uniq, inv, ts, tv
    if not blocks:                                       # k < 2: no pairs to annihilate
        return H1, sp.csr_matrix((0, dim))
    Yd = sp.vstack(blocks, format='csr')
    return H1, Yd


# ----------------------------------------------------------------------------------------------- labels

def casimir_expectations(N, p, masks, V, mu=None):
    """<v|C_2|v>/<v|v> for each column v of V, all of weight mu (default zero), basis order = `masks`."""
    V = V.reshape(len(masks), -1)
    out = np.zeros(V.shape[1])
    for l in range(N):
        for k in range(N):
            if l == k:
                continue
            tm, ts, tv = [], [], []
            for a in range(p):
                for j in range(N):                       # + c^dag_(a,l,j) c_(a,k,j)
                    src, new, sg = apply_hop(masks, mode_index(N, a, l, j), mode_index(N, a, k, j))
                    tm.append(new); ts.append(src); tv.append(sg)
                for i in range(N):                       # - c^dag_(a,i,k) c_(a,i,l)
                    src, new, sg = apply_hop(masks, mode_index(N, a, i, k), mode_index(N, a, i, l))
                    tm.append(new); ts.append(src); tv.append(-sg)
            tm = np.concatenate(tm)
            uniq, inv = np.unique(tm, return_inverse=True)
            W = sp.csr_matrix((np.concatenate(tv), (inv, np.concatenate(ts))), shape=(len(uniq), len(masks)))
            out += np.sum(np.abs(W @ V) ** 2, axis=0)
    out *= 0.5
    if mu is not None:
        out += 0.5 * float(np.dot(mu, mu)) * np.sum(np.abs(V) ** 2, axis=0)
    return out / np.sum(np.abs(V) ** 2, axis=0)


def sector_irreps(N, p, k):
    """Irreps of U(N) in Lambda^k(C^p (x) gl(N)): list of (lambda, multiplicity, c_lambda, dim, zero-weight mult)."""
    zero = (0,) * N
    poly = [dict() for _ in range(k + 1)]
    poly[0][zero] = 1
    for _ in range(p):
        for i in range(N):
            for j in range(N):
                r = [0] * N; r[i] += 1; r[j] -= 1
                for kk in range(k, 0, -1):
                    for w, m in poly[kk - 1].items():
                        w2 = tuple(x + y for x, y in zip(w, r))
                        poly[kk][w2] = poly[kk].get(w2, 0) + m
    irr = fs.decompose(N, poly[k])
    out = [(lam, m, fs.casimir2(lam) / 2, fs.weyl_dim(lam), kostant(N, lam, zero)) for lam, m in irr.items()]
    return sorted(out, key=lambda t: (t[2], t[0]))
