"""
Exact gauge-singlet multiplicities n(k) of the U(n)^3 fermionic quiver by dual Cauchy decomposition (moved from
scripts/quiver_singlet_index.py, 2026-10-07; see that script and research/notes/quiver_project.md section 4).

    Lambda(C^p (x) V_0 (x) Vbar_1) = (x)_f (+)_{lambda in n x n box} S_lambda(V_0) (x) S_lambda'(Vbar_1),
    n(t) = Tr[(D_t G T)^3],  G[P,Q] = int_{U(n)} s_P conj(s_Q) dU   (exact torus quadrature, rounded).
"""
import itertools, math
import numpy as np


def box_partitions(n):
    """Partitions with at most n parts, each <= n (the n x n box), as length-n tuples."""
    out = []
    def rec(prefix, maxpart):
        if len(prefix) == n:
            out.append(tuple(prefix)); return
        for x in range(maxpart, -1, -1):
            rec(prefix + [x], x)
    rec([], n)
    return out


def transpose(lam, n):
    return tuple(sum(1 for x in lam if x > j) for j in range(n))


def schur_on_grid(lam, X):
    """s_lam(x) for each grid point (rows of X, shape (pts, n)) via the bialternant."""
    n = X.shape[1]
    rho = np.arange(n - 1, -1, -1)
    num = np.linalg.det(X[:, None, :] ** (np.array(lam) + rho)[None, :, None])
    den = np.linalg.det(X[:, None, :] ** rho[None, :, None])
    return num / den, den


def singlet_series(n, p, M=None):
    parts = box_partitions(n)
    idx = {l: i for i, l in enumerate(parts)}
    tuples = list(itertools.product(range(len(parts)), repeat=p))
    tidx = {t: i for i, t in enumerate(tuples)}
    deg = np.array([sum(sum(parts[i]) for i in t) for t in tuples])
    tau = np.array([tidx[tuple(idx[transpose(parts[i], n)] for i in t)] for t in tuples])
    # torus grid: integrand exponents bounded by 2pn + (n-1) per variable
    M = M or (2 * p * n + n + 2)
    rng = np.random.default_rng(0)
    shifts = rng.uniform(0, 2 * np.pi / M, size=n)            # generic offsets avoid x_i = x_j
    grid = np.array(list(itertools.product(range(M), repeat=n)), dtype=float) * 2 * np.pi / M + shifts
    X = np.exp(1j * grid)
    S = []
    for lam in parts:
        s, den = schur_on_grid(lam, X)
        S.append(s)
    S = np.array(S)                                          # (parts, pts)
    w = np.abs(den) ** 2 / (math.factorial(n) * M ** n)   # Weyl measure, |a_rho|^2 / n!
    Sp = np.ones((len(tuples), X.shape[0]), dtype=complex)
    for f in range(p):
        Sp *= S[[t[f] for t in tuples]]
    Gc = (Sp * w) @ Sp.conj().T
    G = np.rint(Gc.real).astype(np.int64)
    err = float(np.abs(Gc - G).max())
    Tm = np.zeros((len(tuples), len(tuples)), dtype=np.int64)
    Tm[np.arange(len(tuples)), tau] = 1                      # (G T)[A, C] = G[A, tau C]
    GT = G @ Tm
    dmax = int(deg.max())
    Md = [GT * (deg == d)[:, None] for d in range(dmax + 1)]
    nk = np.zeros(3 * dmax + 1, dtype=object)
    # Tr(M_{d1} M_{d2} M_{d3}) for all degree triples (exact integers)
    pair = {}
    for d2 in range(dmax + 1):
        for d3 in range(dmax + 1):
            pair[(d2, d3)] = Md[d2] @ Md[d3]
    for d1 in range(dmax + 1):
        for d2 in range(dmax + 1):
            for d3 in range(dmax + 1):
                v = int(np.einsum('ij,ji->', Md[d1], pair[(d2, d3)]))
                if v:
                    nk[d1 + d2 + d3] += v
    return [int(x) for x in nk], err, len(tuples)


def refined_index(nk):
    I = [0, 0, 0]
    for k, m in enumerate(nk):
        I[k % 3] += (-1) ** ((k - k % 3) // 3) * m
    return I


