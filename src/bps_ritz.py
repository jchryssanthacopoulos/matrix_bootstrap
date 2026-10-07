"""
Phase 1 of docs/nearbps_bootstrap_plan.md: variational upper bounds on sector ground energies from BPS states
dressed by low-level words (Chen's three-matrix model, U(N), exact Fock-space vectors).

For BPS states B (H B = 0) in the zero-weight block of a window sector k_B and words O of charge q = k - k_B,

    E_0(k) <= min spec( H restricted to V ),      V = span{ O_ii B,  O_ij E_ji B (i != j) }  in sector k, weight 0.

Why this V.  An adjoint-valued word (L_1 ... L_n)_ij has weight e_i - e_j, so a zero-weight product O_mu B_nu needs
nu = e_j - e_i: a root or zero.  For every root alpha, E_alpha maps the zero-weight space of any finite-dimensional
representation onto its alpha-weight space (sl_2 strings through 0 and alpha), and the BPS space is G-invariant,
so BPS_alpha = E_alpha BPS_0.  Hence V is the zero-weight part of the gauge-covariant span {O_mu B_nu}; and since H
compressed to a G-invariant subspace commutes with G, its lowest Ritz value is attained at zero weight.  Any subset
of references or words gives a (weaker) valid bound.  This is the Rayleigh-Ritz analogue of the adjoint-gap bound of
Cho-Gabai-Lin-Yeh-Zheng 2025 with the BPS state as reference: <OB|H|OB> = |[Q,O}B|^2 + |[Qbar,O}B|^2.

Conventions as in sector_ed.py / free_sectors.py: mode m = a N^2 + i N + j, Psi^a_ij = c^dag_m, (Psibar^a)_ij =
c_(a,j,i); c_m and c^dag_m carry (-1)^(number of occupied modes below m).  Gauge generators
E_kl = sum_{a,j} c^dag_(a,k,j) c_(a,l,j) - sum_{a,i} c^dag_(a,i,l) c_(a,i,k)  (weight e_k - e_l).
All vectors are real (Chen's C is real).  Sparse vectors are pairs (masks, vals) that may repeat masks.
"""
import itertools

import numpy as np
import scipy.sparse as sp

from cohomology import weight_basis, mode_index
import sector_ed as se

ONE = np.int64(1)


# ----------------------------------------------------------------------------------------------- blocks

class Block:
    """Sorted occupation masks of the (sector k, weight mu) block of U(N) adjoint fermions with p flavours."""

    def __init__(self, N, p, k, mu=None):
        self.N, self.p, self.k = N, p, k
        self.mu = tuple(mu) if mu is not None else (0,) * N
        m = se.masks_from_tuples(weight_basis(N, p, k, self.mu))
        self.masks = np.sort(m)
        self.dim = len(self.masks)

    def positions(self, masks):
        pos = np.searchsorted(self.masks, masks)
        pos_c = np.minimum(pos, self.dim - 1)
        if not np.all(self.masks[pos_c] == masks):
            raise ValueError("vector leaves the block (wrong sector or weight)")
        return pos

    def dense(self, sv):
        """Sparse vector (masks, vals) -> dense vector on this block (duplicates summed)."""
        masks, vals = sv
        out = np.zeros(self.dim)
        if len(masks):
            np.add.at(out, self.positions(masks), vals)
        return out

    def sparse(self, v):
        nz = np.nonzero(v)[0]
        return self.masks[nz], v[nz].astype(float)


def compress(sv, tol=0.0):
    masks, vals = sv
    if len(masks) == 0:
        return sv
    u, inv = np.unique(masks, return_inverse=True)
    out = np.bincount(inv, weights=vals, minlength=len(u))
    keep = np.abs(out) > tol
    return u[keep], out[keep]


def sv_add(*svs):
    svs = [s for s in svs if len(s[0])]
    if not svs:
        return np.zeros(0, np.int64), np.zeros(0)
    return compress((np.concatenate([s[0] for s in svs]), np.concatenate([s[1] for s in svs])))


# ----------------------------------------------------------------------------------------------- fermion operators

def apply_mode(sv, m, create):
    """c^dag_m (create=True) or c_m on a sparse vector."""
    masks, vals = sv
    bit = ONE << m
    sel = ((masks & bit) == 0) if create else ((masks & bit) != 0)
    mm = masks[sel]
    return mm ^ bit, vals[sel] * se._parity_below(mm, m)


def letter_op(N, letter, x, y):
    """(mode, create) of the (x, y) matrix element of a letter ('P', a) = Psi^a or ('B', a) = Psibar^a."""
    kind, a = letter
    if kind == 'P':
        return mode_index(N, a, x, y), True
    return mode_index(N, a, y, x), False


def apply_word_component(N, word, i, j, sv):
    """(L_1 L_2 ... L_n)_ij sv for a word = [letter, ...] (L_n acts first)."""
    n = len(word)
    if n == 0:
        return sv if i == j else (np.zeros(0, np.int64), np.zeros(0))
    # S[x] = (L_m ... L_n)_{x j} sv, built from the right
    S = {}
    for x in range(N):
        m, cr = letter_op(N, word[-1], x, j)
        S[x] = apply_mode(sv, m, cr)
    for pos in range(n - 2, -1, -1):
        rows = [i] if pos == 0 else range(N)
        new = {}
        for y in rows:
            terms = []
            for x in range(N):
                if len(S[x][0]) == 0:
                    continue
                m, cr = letter_op(N, word[pos], y, x)
                terms.append(apply_mode(S[x], m, cr))
            new[y] = sv_add(*terms)
        S = new
    return compress(S[i]) if n > 1 else compress(S[i])


def apply_gauge(N, p, k_, l_, sv):
    """E_{k l} sv (weight e_k - e_l)."""
    terms = []
    for a in range(p):
        for j in range(N):        # + c^dag_(a,k,j) c_(a,l,j)
            masks, vals = sv
            src, new, sg = se.apply_hop(masks, mode_index(N, a, k_, j), mode_index(N, a, l_, j))
            terms.append((new, vals[src] * sg))
        for i in range(N):        # - c^dag_(a,i,l) c_(a,i,k)
            masks, vals = sv
            src, new, sg = se.apply_hop(masks, mode_index(N, a, i, l_), mode_index(N, a, i, k_))
            terms.append((new, -vals[src] * sg))
    return sv_add(*terms)


def words_of_charge(q, n, p=3):
    """All words of length n with (#Psi - #Psibar) = q, flavours unrestricted: list of tuples of letters."""
    if (n - q) % 2 or abs(q) > n:
        return []
    n_p = (n + q) // 2
    out = []
    for pos in itertools.combinations(range(n), n_p):
        kinds = ['P' if t in pos else 'B' for t in range(n)]
        for fl in itertools.product(range(p), repeat=n):
            out.append(tuple(zip(kinds, fl)))
    return out


def word_label(word):
    return ''.join(('' if k == 'P' else '~') + 'abc'[a] for k, a in word)


# ----------------------------------------------------------------------------------------------- Hamiltonian, Casimir

class SectorH:
    """H = E_c + H_1 + 9 Ydag^T Ydag on a zero-weight block (D2.13 as in sector_ed.build_fast)."""

    def __init__(self, N, C, block):
        import free_sectors as fs
        self.Ec = fs.vacuum_energy(N, C)
        self.H1, self.Yd = se.build_fast(N, C, block.masks)

    def __call__(self, X):
        return self.Ec * X + self.H1 @ X + 9 * (self.Yd.T @ (self.Yd @ X))

    def dense(self):
        return self.Ec * np.eye(self.H1.shape[0]) + self.H1.toarray() + 9 * (self.Yd.T @ self.Yd).toarray()


def casimir_W(N, p, block):
    """The maps W_lk = E_lk restricted to the zero-weight block (l != k); C_2 = (1/2) sum W^T W on that block."""
    Ws = []
    for l in range(N):
        for k in range(N):
            if l == k:
                continue
            tm, ts, tv = [], [], []
            for a in range(p):
                for j in range(N):
                    src, new, sg = se.apply_hop(block.masks, mode_index(N, a, l, j), mode_index(N, a, k, j))
                    tm.append(new); ts.append(src); tv.append(sg)
                for i in range(N):
                    src, new, sg = se.apply_hop(block.masks, mode_index(N, a, i, k), mode_index(N, a, i, l))
                    tm.append(new); ts.append(src); tv.append(-sg)
            tm = np.concatenate(tm)
            uniq, inv = np.unique(tm, return_inverse=True)
            Ws.append(sp.csr_matrix((np.concatenate(tv), (inv, np.concatenate(ts))), shape=(len(uniq), block.dim)))
    return Ws


def casimir_apply(Ws, X):
    return 0.5 * sum(W.T @ (W @ X) for W in Ws)


# ----------------------------------------------------------------------------------------------- BPS references

def bps_basis_dense(N, C, block, tol=1e-9):
    """Orthonormal zero modes of H on a small zero-weight block, resolved by the gauge Casimir.
    Returns (B (dim x n), c2 (n,)) with columns ordered by Casimir."""
    Hd = SectorH(N, C, block).dense()
    lam, U = np.linalg.eigh(Hd)
    Z = U[:, lam < tol]
    Ws = casimir_W(N, C.shape[0], block)
    Cz = Z.T @ casimir_apply(Ws, Z)
    c, V = np.linalg.eigh((Cz + Cz.T) / 2)
    return Z @ V, c, lam


# ----------------------------------------------------------------------------------------------- trial space and RR

def trial_vectors(N, p, refs, words, target, include_offdiag=True):
    """Columns O_ii B and O_ij E_ji B (i != j) on the target block, for each reference (sparse vector) and word.
    Returns (W dense dim x m, labels)."""
    cols, labels = [], []
    for r, B in enumerate(refs):
        lifted = {(i, j): (B if i == j else apply_gauge(N, p, j, i, B))
                  for i in range(N) for j in range(N) if include_offdiag or i == j}
        for word in words:
            for (i, j), Bij in lifted.items():
                if len(Bij[0]) == 0:
                    continue
                v = apply_word_component(N, word, i, j, Bij)
                if len(v[0]) == 0:
                    continue
                cols.append(target.dense(v)); labels.append((r, word_label(word), i, j))
    if not cols:
        return np.zeros((target.dim, 0)), labels
    return np.column_stack(cols), labels


def rayleigh_ritz(Hop, W, rcond=1e-10, n_ritz=6):
    """Lowest Ritz pairs of H on span(W).  Returns dict(values, vectors, residuals, rank).
    Memory: the Gram matrix is min(m, dim)^2 (never m x m when there are more trial vectors than states; the
    m x m form once requested 2.9 TB on an over-complete N=2 trial set and took the machine down, 2026-10-06)."""
    dim, m = W.shape
    if m > dim:
        G = W @ W.T                                  # dim x dim; its range is span(W)
        s, U = np.linalg.eigh((G + G.T) / 2)
        keep = s > rcond * max(s.max(), 1e-300)
        Qx = U[:, keep]
    else:
        G = W.T @ W
        s, U = np.linalg.eigh((G + G.T) / 2)
        keep = s > rcond * max(s.max(), 1e-300)
        X = W @ (U[:, keep] / np.sqrt(s[keep]))
        # one re-orthonormalisation pass (Gram-based orthonormalisation loses digits when W is ill-conditioned)
        Qx, _ = np.linalg.qr(X)
    del G
    HX = Hop(Qx)
    Hs = Qx.T @ HX
    e, c = np.linalg.eigh((Hs + Hs.T) / 2)
    n = min(n_ritz, len(e))
    V = Qx @ c[:, :n]
    HV = HX @ c[:, :n]
    res = np.linalg.norm(HV - V * e[:n], axis=0)
    return dict(values=e[:n], vectors=V, residuals=res, rank=int(keep.sum()), all_values=e)
