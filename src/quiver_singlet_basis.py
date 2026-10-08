"""
Singlet-adapted basis for the U(n)^3 fermionic quiver via Schur-Weyl duality (research/notes/quiver_project.md
section 11; docs/derivations.md D20).

Permutation picture.  A gauge singlet of degree k = 3m has m letters on every edge (edge e = v -> v+1 carries
A, B, C for e = 0, 1, 2).  At node v the node-v indices of the m out-letters (edge v -> v+1) are contracted with
those of the m in-letters (edge v-1 -> v) by a permutation sigma_v in S_m, so a singlet is a function x on S_m^3.
Letters carry flavours; with a fixed canonical order (flavour block 1 first, ...) the remaining freedom is the
relabelling group H_e = S_{r_1} x ... x S_{r_p} of each edge, and fermionic states satisfy |Y> = sgn(h) |h Y>.

Fourier picture.  x -> (X_v hat in End(V_{lambda_v}))_v with lambda_v |- m (Young orthogonal form, real orthogonal).
  * Finite n: only ell(lambda_v) <= n survives (Schur-Weyl: n^{c(pi)} = sum_lambda D_lambda(n) chi_lambda(pi)).
  * Relabelling: in-letters act on the row index of X_v hat, out-letters on the column index, both by rho(h).
    Edge e = v -> v+1 therefore lives on (column of node v) (x) (row of node v+1), and physical states are the
    sign-isotypic subspace E(lambda_v, lambda_{v+1}; r_e) of H_e acting diagonally there.
  * Singlets of degree 3m = (+)_{lambda's, r's} (x)_e E(lambda_v, lambda_{v+1}; r_e); dimension = n(3m).
  * Fock norm of a projected state: prod_e |H_e| * prod_v D_{lambda_v}(n) * (Euclidean norm of the Fourier data).
  * Q_{abc} = Tr(A^a B^b C^c) appends one letter per edge, closing a new trace: at each node sigma -> sigma (+) fixed
    point, i.e. X hat -> J X hat J^T with J : V_lambda -> V_mu (mu = lambda + box) the Gelfand-Tsetlin inclusion
    (append the new entry to the tableau).  The new letter (position m) is then relabelled to the end of its flavour
    block (sign (-1)^{#positions passed}).  With the canonical order A-block, B-block, C-block the operator ordering
    sign is +1.

Edge spaces.  The first flavour block is solved analytically: the GT basis restricts to S_{r_1} as V_nu (x) Skew,
and the sign-isotypic vector of V_nu (x) V_nu' is sum_U sgn(rowword U) U (x) U^t / sqrt(f_nu).  Further blocks are
handled numerically on the skew modules, by antisymmetrisers A_r = X_2 ... X_r / r!, X_k = 1 - sum_{j<k} (j k).
"""
import itertools
import math
from functools import lru_cache
import numpy as np
import scipy.sparse as sp


# ------------------------------------------------------------------------------------------------ partitions / tableaux
def partitions(m, maxrows):
    out = []

    def rec(rem, maxpart, cur):
        if rem == 0:
            out.append(tuple(cur))
            return
        if len(cur) == maxrows:
            return
        for k in range(min(rem, maxpart), 0, -1):
            rec(rem - k, k, cur + [k])
    rec(m, m, [])
    return out


def conjugate(shape):
    return tuple(sum(1 for r in shape if r > c) for c in range(shape[0])) if shape else ()


def contains(big, small):
    return len(small) <= len(big) and all(s <= b for s, b in zip(small, big))


def add_box(shape, row):
    s = list(shape) + [0]
    s[row] += 1
    if row > 0 and s[row] > s[row - 1]:
        return None
    return tuple(x for x in s if x)


def syt_words(shape):
    """Standard Young tableaux of a shape as Yamanouchi words (row index of entries 0..m-1)."""
    m = sum(shape)
    out = []
    counts = [0] * len(shape)

    def rec(word):
        if len(word) == m:
            out.append(tuple(word))
            return
        for r in range(len(shape)):
            if counts[r] < shape[r] and (r == 0 or counts[r] < counts[r - 1]):
                counts[r] += 1
                word.append(r)
                rec(word)
                word.pop()
                counts[r] -= 1
    rec([])
    return out


def word_cols(word):
    seen = {}
    cols = []
    for r in word:
        cols.append(seen.get(r, 0))
        seen[r] = seen.get(r, 0) + 1
    return cols


def perm_sign(seq):
    seq = list(seq)
    s = 1
    seen = [False] * len(seq)
    for i in range(len(seq)):
        if not seen[i]:
            j, L = i, 0
            while not seen[j]:
                seen[j] = True
                j = seq[j]
                L += 1
            if L % 2 == 0:
                s = -s
    return s


def udim(shape, n):
    """dim of the U(n) irrep S_shape(C^n) = prod (n + content) / hook."""
    if len(shape) > n:
        return 0
    conj = conjugate(shape)
    num, den = 1, 1
    for r, row in enumerate(shape):
        for c in range(row):
            num *= n + c - r
            den *= (row - c - 1) + (conj[c] - r - 1) + 1
    return num // den


class Shape:
    """GT basis of the S_m irrep `shape` with Young's orthogonal form for the adjacent transpositions."""

    def __init__(self, shape):
        self.shape = tuple(shape)
        self.m = sum(self.shape)
        self.words = syt_words(self.shape)
        self.index = {w: i for i, w in enumerate(self.words)}
        self.d = len(self.words)
        self._s = {}

    def s(self, i):
        """rho(s_i), s_i = (i, i+1) on entries 0..m-1, as a CSR matrix (columns = input tableaux)."""
        if i not in self._s:
            rows, cols, vals = [], [], []
            for t, w in enumerate(self.words):
                cw = word_cols(w)
                r1, r2, c1, c2 = w[i], w[i + 1], cw[i], cw[i + 1]
                if r1 == r2:
                    rows.append(t); cols.append(t); vals.append(1.0)
                elif c1 == c2:
                    rows.append(t); cols.append(t); vals.append(-1.0)
                else:
                    a = (c2 - r2) - (c1 - r1)
                    w2 = list(w); w2[i], w2[i + 1] = w2[i + 1], w2[i]
                    t2 = self.index[tuple(w2)]
                    rows += [t, t2]; cols += [t, t]; vals += [1.0 / a, math.sqrt(1.0 - 1.0 / a ** 2)]
            self._s[i] = sp.csr_matrix((vals, (rows, cols)), shape=(self.d, self.d))
        return self._s[i]

    def J(self, row, big):
        """GT inclusion V_shape -> V_big (big = shape + box in `row`): append entry m in that row."""
        rows = [big.index[w + (row,)] for w in self.words]
        return sp.csr_matrix((np.ones(self.d), (rows, np.arange(self.d))), shape=(big.d, self.d))


@lru_cache(maxsize=None)
def shape_obj(shape):
    return Shape(shape)


# ------------------------------------------------------------------------------------------------ edge spaces
def _apply_adjacent_pair(Sa, Sb, i, Y):
    """(rho_a(s_i) (x) rho_b(s_i)) applied to Y (d_a x d_b x K or d_a x d_b): rho_a Y rho_b^T."""
    if Y.ndim == 2:
        return Sa.s(i) @ Y @ Sb.s(i).T
    return np.einsum('ij,jkc->ikc', Sa.s(i).toarray(), np.einsum('jkc,lk->jlc', Y, Sb.s(i).toarray()))


def _skew_action(full_shape, nu_words0, skew_words, skew_index, i):
    """Matrix of rho(s_i) (i >= r1) on the skew index set, using the full tableau rules with a fixed prefix."""
    S = shape_obj(full_shape)
    u = nu_words0
    rows, cols, vals = [], [], []
    Si = S.s(i).tocsc()
    for j, sw in enumerate(skew_words):
        t = S.index[u + sw]
        lo, hi = Si.indptr[t], Si.indptr[t + 1]
        for tt, v in zip(Si.indices[lo:hi], Si.data[lo:hi]):
            w2 = S.words[tt]
            assert w2[:len(u)] == u
            rows.append(skew_index[w2[len(u):]]); cols.append(j); vals.append(v)
    n = len(skew_words)
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n))


def _antisym_apply(mats, positions, V):
    """Apply the product over blocks of the antisymmetriser A_r = X_2 X_3 ... X_r / r!, X_k = 1 - sum_{j<k} (j k),
    for the diagonal action; mats[i] = sparse matrix of s_i (dim x dim); positions = list of (start, length).
    The transposition (j k) = s_j s_{j+1} ... s_{k-2} s_{k-1} s_{k-2} ... s_j."""
    for start, length in positions:
        for k in range(start + length - 1, start, -1):          # X_r first (rightmost), X_2 last
            acc = V.copy()
            for j in range(k - 1, start - 1, -1):
                seq = list(range(j, k - 1)) + [k - 1] + list(range(k - 2, j - 1, -1))
                T = V
                for i in reversed(seq):
                    T = mats[i] @ T
                acc -= T
            V = acc / (k - start + 1)
    return V


@lru_cache(maxsize=None)
def edge_space(lam, mu, r):
    """Orthonormal basis (sparse, (d_lam*d_mu) x K, row index t_lam*d_mu + t_mu) of the sign-isotypic subspace of
    H_r = S_{r_1} x ... acting diagonally on V_lam (x) V_mu."""
    L, M = shape_obj(lam), shape_obj(mu)
    m = L.m
    assert M.m == m and sum(r) == m
    r1 = r[0]
    blocks = []
    o = r1
    for rf in r[1:]:
        if rf > 1:
            blocks.append((o, rf))
        o += rf
    cols = []
    for nu in partitions(r1, len(lam)) if r1 > 0 else [()]:
        if not contains(lam, nu) or not contains(mu, conjugate(nu)):
            continue
        nup = conjugate(nu)
        if r1 > 0:
            Nu = shape_obj(nu)
            U = Nu.words
            Ut = [tuple(word_cols(u)) for u in U]                       # transposed tableaux (rows = old columns)
            eps = []
            for u in U:
                rowread = [k for rr in range(len(nu)) for k in range(r1) if u[k] == rr]
                eps.append(perm_sign(rowread))
            fnu = len(U)
        else:
            U, Ut, eps, fnu = [()], [()], [1], 1
        SL = sorted({w[r1:] for w in L.words if w[:r1] == U[0]})
        SM = sorted({w[r1:] for w in M.words if w[:r1] == Ut[0]})
        if not SL or not SM:
            continue
        iL = {s: j for j, s in enumerate(SL)}
        iM = {s: j for j, s in enumerate(SM)}
        nL, nM = len(SL), len(SM)
        if blocks:
            mats = {}
            for (start, length) in blocks:
                for i in range(start, start + length - 1):
                    A = _skew_action(lam, U[0], SL, iL, i)
                    B = _skew_action(mu, Ut[0], SM, iM, i)
                    mats[i] = sp.kron(A, B, format='csr')
            dim = nL * nM
            rng = np.random.default_rng(12345)
            width = min(dim, 32)
            while True:                                   # adaptive width: grow until the range is not saturated
                R = rng.standard_normal((dim, width))
                Y = _antisym_apply(mats, blocks, R)
                Uq, sv, _ = np.linalg.svd(Y, full_matrices=False)
                # projected Gaussian block: non-zero singular values are O(sqrt(width)), the rest round-off
                rank = int((sv > 1e-6 * math.sqrt(width)).sum())
                if rank < width or width == dim:
                    break
                width = min(dim, 2 * width)
            Z = Uq[:, :rank]
        else:
            Z = np.eye(nL * nM)
        if Z.shape[1] == 0:
            continue
        # assemble columns in the full tableau basis
        rows_idx = []
        coef = []
        for a, (u, ut) in enumerate(zip(U, Ut)):
            for jl, sl in enumerate(SL):
                tL = L.index[u + sl]
                for jm, sm in enumerate(SM):
                    tM = M.index[ut + sm]
                    rows_idx.append(tL * M.d + tM)
                    coef.append((a, jl * nM + jm))
        rows_idx = np.array(rows_idx)
        avec = np.array([c[0] for c in coef])
        zrow = np.array([c[1] for c in coef])
        epsv = np.array(eps, float)[avec] / math.sqrt(fnu)
        for c in range(Z.shape[1]):
            cols.append((rows_idx, epsv * Z[zrow, c]))
    if not cols:
        return sp.csr_matrix((L.d * M.d, 0))
    data, ri, ci = [], [], []
    for c, (rr, vv) in enumerate(cols):
        keep = np.abs(vv) > 1e-14
        ri.append(rr[keep]); data.append(vv[keep]); ci.append(np.full(keep.sum(), c))
    return sp.csc_matrix((np.concatenate(data), (np.concatenate(ri), np.concatenate(ci))),
                         shape=(L.d * M.d, len(cols)))


def compositions(m, p, cap):
    """Compositions r of m into p parts with 0 <= r_f <= cap (cap = n^2: at most n^2 letters per flavour)."""
    out = []
    for r in itertools.product(range(min(m, cap) + 1), repeat=p):
        if sum(r) == m:
            out.append(r)
    return out


# ------------------------------------------------------------------------------------------------ the singlet complex
class SingletComplex:
    """Singlet sector of the U(n)^3 quiver with p flavours, degree by degree (k = 3m), in the Schur-Weyl basis."""

    def __init__(self, n, p, C):
        self.n, self.p = n, p
        self.C = np.asarray(C)
        self._edge = {}           # (lam, mu) -> list of (r, offset, K, basis)
        self._sectors = {}

    def edge(self, lam, mu):
        """Flavour-summed edge space F(lam, mu) = (+)_r E(lam, mu; r): list of (r, offset, K, basis), total dim."""
        key = (lam, mu)
        if key not in self._edge:
            m = sum(lam)
            blocks, off = [], 0
            for r in compositions(m, self.p, self.n * self.n):
                B = edge_space(lam, mu, r)
                if B.shape[1]:
                    blocks.append((r, off, B.shape[1], B))
                    off += B.shape[1]
            self._edge[key] = (blocks, off)
        return self._edge[key]

    def sectors(self, m):
        """List of (lam0, lam1, lam2, dims (R01, R12, R20), offset) for degree 3m, and the total dimension."""
        if m not in self._sectors:
            parts = partitions(m, self.n) if m > 0 else [()]
            out, off = [], 0
            for l0, l1, l2 in itertools.product(parts, repeat=3):
                R = (self.edge(l0, l1)[1], self.edge(l1, l2)[1], self.edge(l2, l0)[1])
                size = R[0] * R[1] * R[2]
                if size:
                    out.append((l0, l1, l2, R, off))
                    off += size
            self._sectors[m] = (out, off)
        return self._sectors[m]

    def dim(self, m):
        return self.sectors(m)[1]

    @lru_cache(maxsize=None)
    def edge_map(self, lam, lamp, mu, mup, f):
        """Flavour-summed edge map F(lam, lamp) -> F(mu, mup) adding one letter of flavour f (0-based), including the
        relabelling sign and the |H| normalisation factor sqrt(r_f + 1)."""
        src, Rs = self.edge(lam, lamp)
        dst, Rd = self.edge(mu, mup)
        L, Lp, Mu, Mup = (shape_obj(x) for x in (lam, lamp, mu, mup))
        m = L.m
        rowL = next(i for i in range(len(mu)) if (lam + (0,))[i] != (mu + (0,) * 2)[i]) if lam else 0
        rowLp = next(i for i in range(len(mup)) if (lamp + (0,))[i] != (mup + (0,) * 2)[i]) if lamp else 0
        JL = L.J(rowL, Mu) if lam else sp.csr_matrix(np.ones((1, 1)))
        JLp = Lp.J(rowLp, Mup) if lamp else sp.csr_matrix(np.ones((1, 1)))
        dst_by_r = {r: (off, K, B) for (r, off, K, B) in dst}
        out = np.zeros((Rd, Rs))
        for (r, off, K, B) in src:
            r2 = list(r); r2[f] += 1; r2 = tuple(r2)
            if r2 not in dst_by_r:
                continue
            off2, K2, B2 = dst_by_r[r2]
            t = sum(r[:f + 1])                      # target position (0-based) = end of block f in the new order
            Bd = B.toarray().reshape(L.d, Lp.d, K)
            Y = np.einsum('ia,abk->ibk', JL.toarray(), Bd)
            Y = np.einsum('ibk,jb->ijk', Y, JLp.toarray())          # (Mu.d, Mup.d, K)
            for i in range(m - 1, t - 1, -1):                       # rho(h) = rho(s_t) ... rho(s_{m-1}); s_{m-1} first
                Y = np.einsum('ij,jbk->ibk', Mu.s(i).toarray(), Y)
                Y = np.einsum('ibk,jb->ijk', Y, Mup.s(i).toarray())
            Y *= (-1) ** (m - t) * math.sqrt(r[f] + 1)
            out[off2:off2 + K2, off:off + K] = B2.T @ Y.reshape(Mu.d * Mup.d, K)
        return out

    def Q_blocks(self, m):
        """List of (source sector index, target sector index, node factor, (mu0, mu1, mu2), lam triple) for Q_m."""
        src, _ = self.sectors(m)
        dst, _ = self.sectors(m + 1)
        dindex = {(s[0], s[1], s[2]): i for i, s in enumerate(dst)}
        out = []
        for i, (l0, l1, l2, R, off) in enumerate(src):
            rows0 = range(len(l0) + 1) if l0 else [0]
            rows1 = range(len(l1) + 1) if l1 else [0]
            rows2 = range(len(l2) + 1) if l2 else [0]
            for a0, a1, a2 in itertools.product(rows0, rows1, rows2):
                mus = []
                ok = True
                for lam, a in ((l0, a0), (l1, a1), (l2, a2)):
                    mu = add_box(lam, a) if lam else (1,)
                    if mu is None or len(mu) > self.n:
                        ok = False
                        break
                    mus.append(mu)
                if not ok:
                    continue
                j = dindex.get(tuple(mus))
                if j is None:
                    continue
                fac = math.sqrt(np.prod([udim(mu, self.n) / max(udim(lam, self.n), 1) if lam else udim(mu, self.n)
                                         for lam, mu in zip((l0, l1, l2), mus)]))
                out.append((i, j, fac, tuple(mus)))
        return out

    def stacked_map(self, lam, lamp, mu, mup):
        """Edge maps for all flavours stacked: array (p, R(mu,mup), R(lam,lamp))."""
        key = ('stack', lam, lamp, mu, mup)
        if key not in self._sectors:
            self._sectors[key] = np.stack([self.edge_map(lam, lamp, mu, mup, f) for f in range(self.p)])
        return self._sectors[key]

    def apply_Q(self, m, x):
        """Q_m x (Fock-orthonormal singlet bases at degrees 3m and 3m+3)."""
        src, _ = self.sectors(m)
        dst, Nd = self.sectors(m + 1)
        y = np.zeros(Nd, dtype=np.result_type(x, self.C))
        for (i, j, fac, mus) in self._blocks(m):
            l0, l1, l2, R, off = src[i]
            m0, m1, m2, R2, off2 = dst[j]
            X = x[off:off + R[0] * R[1] * R[2]].reshape(R)
            Mc = self.stacked_map(l2, l0, m2, m0)                   # (p, R2'', R2)
            Mb = self.stacked_map(l1, l2, m1, m2)
            Ma = self.stacked_map(l0, l1, m0, m1)
            T = np.einsum('xyz,cwz->cxyw', X, Mc)                   # (c, R0, R1, R2')
            T = np.einsum('cxyw,bvy->bcxvw', T, Mb)                 # (b, c, R0, R1', R2')
            T = np.einsum('abc,bcxvw->axvw', self.C, T)            # (a, R0, R1', R2')
            acc = np.einsum('axvw,aux->uvw', T, Ma)                 # (R0', R1', R2')
            y[off2:off2 + acc.size] += fac * acc.ravel()
        return y

    def apply_Qdag(self, m, y):
        src, Ns = self.sectors(m)
        dst, _ = self.sectors(m + 1)
        x = np.zeros(Ns, dtype=np.result_type(y, self.C))
        Cc = np.conj(self.C)
        for (i, j, fac, mus) in self._blocks(m):
            l0, l1, l2, R, off = src[i]
            m0, m1, m2, R2, off2 = dst[j]
            Y = y[off2:off2 + R2[0] * R2[1] * R2[2]].reshape(R2)
            Mc = self.stacked_map(l2, l0, m2, m0)
            Mb = self.stacked_map(l1, l2, m1, m2)
            Ma = self.stacked_map(l0, l1, m0, m1)
            T = np.einsum('uvw,cwz->cuvz', Y, Mc)                   # (c, R0', R1', R2)
            T = np.einsum('cuvz,bvy->bcuyz', T, Mb)                 # (b, c, R0', R1, R2)
            T = np.einsum('abc,bcuyz->auyz', Cc, T)                # (a, R0', R1, R2)
            acc = np.einsum('auyz,aux->xyz', T, Ma)                 # (R0, R1, R2)
            x[off:off + acc.size] += fac * acc.ravel()
        return x

    def _blocks(self, m):
        key = ('blocks', m)
        if key not in self._sectors:
            self._sectors[key] = self.Q_blocks(m)
        return self._sectors[key]

    def Q_dense(self, m):
        Ns, Nd = self.dim(m), self.dim(m + 1)
        Q = np.zeros((Nd, Ns), dtype=np.result_type(self.C, float))
        for j in range(Ns):
            e = np.zeros(Ns); e[j] = 1.0
            Q[:, j] = self.apply_Q(m, e)
        return Q
