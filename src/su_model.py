"""The traceless (su(N)) variant of the fermionic matrix models, and its exact Q-cohomology.

Motivation: for the u(N) model the BPS window has width N+1, and Chen 2025 (2.17)/(2.21) identifies the source as
the (1+q)^N from the N diagonal modes Psi_ii, each freely occupied.  su(N) has only N-1 Cartan directions, so the
window should shrink by one.  Since concentration needs the window no wider than the grading period q = 3
(note section 4j), that is exactly the margin that decides whether N = 3 concentrates.

Basis of su(N) (as a complex space): E_ij for i != j (weight e_i - e_j), and H_m = E_mm - E_{m+1,m+1} for
m = 1..N-1 (weight 0).  A mode is (flavour a, basis element alpha), so p(N^2-1) modes in all.  The supercharge is

    Q = sum_abc C_abc Tr[Psi^a Psi^b Psi^c] = sum g_{(a,alpha)(b,beta)(c,gamma)} c^dag c^dag c^dag,
    g = C_abc Tr[T_alpha T_beta T_gamma].

Q^2 = 0 for ANY g: the six-fermion product is totally antisymmetric while g (x) g is symmetric under exchanging
the two triples, which is an odd permutation of six objects.
"""
import itertools
from math import comb

import numpy as np

import cohomology as co


def su_basis(N):
    """(matrices, weights) for the su(N) basis: off-diagonal E_ij then the N-1 Cartan elements."""
    mats, wts = [], []
    for i in range(N):
        for j in range(N):
            if i != j:
                T = np.zeros((N, N)); T[i, j] = 1.0
                w = [0] * N; w[i] += 1; w[j] -= 1
                mats.append(T); wts.append(tuple(w))
    for m in range(N - 1):
        T = np.zeros((N, N)); T[m, m] = 1.0; T[m + 1, m + 1] = -1.0
        mats.append(T); wts.append(tuple([0] * N))
    return mats, wts


def modes(N, p):
    """Mode list [(a, alpha)] and their weights; index = a*(N^2-1) + alpha."""
    mats, wts = su_basis(N)
    nb = len(mats)
    out = [(a, al) for a in range(p) for al in range(nb)]
    w = [wts[al] for (a, al) in out]
    return out, w, mats


def cubic_terms(N, p, C, tol=1e-12):
    """Non-zero g_{m1 m2 m3} = C_abc Tr[T_al T_be T_ga] over ordered mode triples."""
    md, _, mats = modes(N, p)
    nb = len(mats)
    tr = np.einsum('aij,bjk,cki->abc', np.array(mats), np.array(mats), np.array(mats))
    terms = []
    for (i1, (a, al)) in enumerate(md):
        for (i2, (b, be)) in enumerate(md):
            for (i3, (c, ga)) in enumerate(md):
                g = C[a, b, c] * tr[al, be, ga]
                if abs(g) > tol:
                    terms.append((i1, i2, i3, g))
    return terms


def apply_Q(terms, state):
    """Q|state> as {new_state: coeff}, with Q given by its list of cubic terms."""
    out = {}
    for (m1, m2, m3, g) in terms:
        s, st = 1, state
        for m in (m3, m2, m1):               # rightmost acts first
            r = co._create(st, m)
            if r is None:
                s = 0; break
            sg, st = r; s *= sg
        if s:
            out[st] = out.get(st, 0) + s * g
    return {k: v for k, v in out.items() if abs(v) > 1e-9}


def offdiag_occupations(N, p, k_off, lam):
    """Off-diagonal occupations n_ij (i != j) in 0..p with total k_off and row-sums minus col-sums = lam.
    Pruned row by row, exactly as cohomology.occupation_matrices, but with the diagonal held at zero (in su(N)
    the diagonal is replaced by the N-1 zero-weight Cartan modes, handled separately)."""
    lam = tuple(lam)
    rows_choices = [c for c in itertools.product(range(p + 1), repeat=N)]
    out = []
    chosen, colsum = [], [0] * N

    def rec(i, used):
        if used > k_off or k_off - used > p * (N - 1) * (N - i):
            return
        for t in range(i):
            need = chosen[t][1] - lam[t] - colsum[t]
            if need < 0 or need > p * (N - 1 - (0 if i > t else 0)) * (N - i):
                return
        if i == N:
            if used == k_off and all(chosen[t][1] - lam[t] == colsum[t] for t in range(N)):
                out.append(tuple(r for r, _ in chosen))
            return
        for row in rows_choices:
            if row[i]:                       # diagonal entry must be zero
                continue
            r = sum(row)
            chosen.append((row, r))
            for j in range(N):
                colsum[j] += row[j]
            rec(i + 1, used + r)
            for j in range(N):
                colsum[j] -= row[j]
            chosen.pop()

    rec(0, 0)
    return out


def weight_basis(N, p, k, lam):
    """States of k modes with total weight lam.  Off-diagonal occupations n_ij (i != j) carry the weight; the
    remaining modes come from the pool of p(N-1) zero-weight Cartan modes -- that pool is the ONLY difference
    from u(N), where it has size pN."""
    nb = N * N - 1
    off = [(i, j) for i in range(N) for j in range(N) if i != j]
    pos = {ij: c for c, ij in enumerate(off)}
    zero_pool = [a * nb + (N * (N - 1) + m) for a in range(p) for m in range(N - 1)]
    out = []
    for k_off in range(min(k, p * len(off)) + 1):
        nz = k - k_off
        if nz < 0 or nz > len(zero_pool):
            continue
        for n in offdiag_occupations(N, p, k_off, lam):
            occ = [n[i][j] for (i, j) in off]
            choices = [[c for c in itertools.combinations(range(p), cnt)] for cnt in occ]
            for combo in itertools.product(*choices):
                base = []
                for (fl, (i, j)) in zip(combo, off):
                    for a in fl:
                        base.append(a * nb + (i * (N - 1) + (j if j < i else j - 1)))
                for zsub in itertools.combinations(zero_pool, nz):
                    out.append(tuple(sorted(base + list(zsub))))
    return out
