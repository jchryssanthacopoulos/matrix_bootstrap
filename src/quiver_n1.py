"""
The n = 1 member of the quiver family (research/notes/quiver_project.md section 10): gauge group U(1)^3, 3p complex
fermions a_f, b_f, c_f (f = 0..p-1) on the three edges, Q = sum_fgh C_fgh a_f^+ b_g^+ c_h^+, H = {Q, Q^+}.

Gauge singlets have N_a = N_b = N_c = m, so the singlet sector is a sum of blocks
    B_m = Lambda^m(C^p) (x) Lambda^m(C^p) (x) Lambda^m(C^p),   dim C(p,m)^3,   fermion number k = 3m,
and Q maps B_m -> B_{m+1}.  This is a three-species ("tripartite") N = 2 SYK model with q-hat = 3 restricted to its
charge-balanced sector; its large-p limit is the large-p limit of the whole family at n = 1 (docs/derivations.md D19).

Conventions.  |S_a, S_b, S_c> = (prod_{g in S_a} a_g^+)(prod b^+)(prod c^+)|0>, each product in increasing order.
a_f^+ b_g^+ c_h^+ acting on it gives (-1)^m s(f,S_a) s(g,S_b) s(h,S_c) |S_a+f, S_b+g, S_c+h>, s(f,S) = (-1)^{#{x in S:
x < f}}: c_h^+ passes 2m letters, b_g^+ passes m.  So Q_m = (-1)^m sum C_fgh E_f (x) E_g (x) E_h, with E_f the
single-species creation matrices Lambda^m -> Lambda^{m+1}.
"""
import itertools
import math
import numpy as np
import scipy.sparse as sp


def subsets(p, m):
    return list(itertools.combinations(range(p), m))


def creation_matrices(p, m):
    """E_f : Lambda^m(C^p) -> Lambda^{m+1}(C^p), |S> -> (-1)^{#{x in S: x < f}} |S u {f}> (zero if f in S)."""
    src, dst = subsets(p, m), subsets(p, m + 1)
    idx = {s: i for i, s in enumerate(dst)}
    out = []
    for f in range(p):
        rows, cols, vals = [], [], []
        for j, S in enumerate(src):
            if f in S:
                continue
            rows.append(idx[tuple(sorted(S + (f,)))])
            cols.append(j)
            vals.append((-1) ** sum(1 for x in S if x < f))
        out.append(sp.csr_matrix((np.array(vals, float), (rows, cols)), shape=(len(dst), len(src))))
    return out


class QuiverN1:
    """Q blocks of the n = 1 quiver; explicit sparse (small p) or matrix-free (large p) application."""

    def __init__(self, p, C):
        self.p = p
        self.C = np.asarray(C)
        assert self.C.shape == (p, p, p)
        self.E = {m: creation_matrices(p, m) for m in range(p)}
        self.D = [math.comb(p, m) for m in range(p + 1)]

    def dim(self, m):
        return self.D[m] ** 3

    # ---------------------------------------------------------------- explicit sparse matrices (p <= 8 or so)
    def Q_sparse(self, m):
        """Q_m : B_m -> B_{m+1} as a sparse matrix (row-major index S_a*D^2 + S_b*D + S_c)."""
        E = self.E[m]
        Q = None
        for f, g, h in itertools.product(range(self.p), repeat=3):
            c = self.C[f, g, h]
            if c == 0:
                continue
            term = sp.kron(E[f], sp.kron(E[g], E[h], format='csr'), format='csr') * c
            Q = term if Q is None else Q + term
        return ((-1) ** m) * Q.tocsr()

    # ---------------------------------------------------------------- matrix-free application (large p)
    def apply_Q(self, m, v):
        """Q_m v for v of length D_m^3 (returns length D_{m+1}^3).  The flavour contraction is one tensordot per f."""
        D, D1, p, E = self.D[m], self.D[m + 1], self.p, self.E[m]
        psi = np.asarray(v).reshape(D, D, D)
        dt = np.result_type(psi, self.C)
        out = np.zeros((D1, D1, D1), dtype=dt)
        flat = psi.reshape(D * D, D)
        Z = np.stack([np.asarray((E[h] @ flat.T).T).reshape(D, D, D1) for h in range(p)])   # (h, D, D, D1)
        for f in range(p):
            W = np.tensordot(self.C[f], Z, axes=([1], [0]))                                    # (g, D, D, D1)
            V = np.zeros((D, D1, D1), dtype=dt)
            for g in range(p):
                Wg = np.transpose(W[g], (1, 0, 2)).reshape(D, D * D1)                           # b-index first
                V += np.transpose((E[g] @ Wg).reshape(D1, D, D1), (1, 0, 2))                    # b_g^+ on axis 1
            out += (E[f] @ V.reshape(D, D1 * D1)).reshape(D1, D1, D1)                           # a_f^+ on axis 0
        return ((-1) ** m) * out.reshape(-1)

    def apply_Qdag(self, m, w):
        """Q_m^+ w for w of length D_{m+1}^3 (returns length D_m^3)."""
        D, D1, p, E = self.D[m], self.D[m + 1], self.p, self.E[m]
        phi = np.asarray(w).reshape(D1, D1, D1)
        dt = np.result_type(phi, self.C)
        Cc = np.conj(self.C)
        out = np.zeros((D, D, D), dtype=dt)
        flat = phi.reshape(D1 * D1, D1)
        Z = np.stack([np.asarray((E[h].T @ flat.T).T).reshape(D1, D1, D) for h in range(p)])  # (h, D1, D1, D)
        for f in range(p):
            W = np.tensordot(Cc[f], Z, axes=([1], [0]))                                       # (g, D1, D1, D)
            V = np.zeros((D1, D, D), dtype=dt)
            for g in range(p):
                Wg = np.transpose(W[g], (1, 0, 2)).reshape(D1, D1 * D)
                V += np.transpose((E[g].T @ Wg).reshape(D, D1, D), (1, 0, 2))
            out += (E[f].T @ V.reshape(D1, D * D)).reshape(D, D, D)
        return ((-1) ** m) * out.reshape(-1)


def dixon_index(p):
    """Singlet index of the n = 1 quiver: sum_m (-1)^m C(p,m)^3."""
    return sum((-1) ** m * math.comb(p, m) ** 3 for m in range(p + 1))
