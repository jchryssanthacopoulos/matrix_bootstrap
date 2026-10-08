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




# ----------------------------------------------------------------------------------------------- exact index at t = -1

def partitions_in_box(rows, cols):
    """All partitions with at most `rows` parts, each <= cols, as length-`rows` tuples."""
    out = []
    def rec(prefix, maxpart):
        if len(prefix) == rows:
            out.append(tuple(prefix)); return
        for x in range(maxpart, -1, -1):
            rec(prefix + [x], x)
    rec([], cols)
    return out


def lr_product_coefficients(n, p, tuple_chunk=2000, verbose=False):
    """c[P, nu] = multiplicity of S_nu(V) in S_{lambda_1}(V) (x) ... (x) S_{lambda_p}(V), V = C^n, for every p-tuple
    P of partitions in the n x n box and nu in the n x (p n) box, by exact quadrature on the maximal torus of U(n):
    the integrand is a Laurent polynomial with exponents |e_i| <= p n + n - 1, so a uniform grid with
    M = p n + n points per angle (generically shifted) integrates it exactly.  Returns (c as int64, parts, nus, err)."""
    parts = box_partitions(n)
    nus = partitions_in_box(n, p * n)
    M = p * n + n + 1
    rng = np.random.default_rng(0)
    shifts = rng.uniform(0, 2 * np.pi / M, size=n)
    grid = np.array(list(itertools.product(range(M), repeat=n)), dtype=float) * 2 * np.pi / M + shifts
    X = np.exp(1j * grid)
    S = []
    for lam in parts:
        s, den = schur_on_grid(lam, X)
        S.append(s)
    S = np.array(S)
    w = np.abs(den) ** 2 / (math.factorial(n) * M ** n)
    Snu = np.array([schur_on_grid(nu, X)[0] for nu in nus])          # (n_nu, pts)
    B = (Snu.conj() * w).T                                             # (pts, n_nu)
    tuples = list(itertools.product(range(len(parts)), repeat=p))
    c = np.zeros((len(tuples), len(nus)), dtype=np.int64)
    err = 0.0
    T = np.array(tuples)
    for t0 in range(0, len(tuples), tuple_chunk):
        idx = T[t0:t0 + tuple_chunk]
        prod = np.ones((len(idx), X.shape[0]), dtype=complex)
        for f in range(p):
            prod *= S[idx[:, f]]
        cc = prod @ B
        rr = np.rint(cc.real)
        err = max(err, float(np.abs(cc - rr).max()))
        c[t0:t0 + len(idx)] = rr.astype(np.int64)
        if verbose:
            print(f"    LR coefficients: tuples {t0 + len(idx)}/{len(tuples)}, rounding error so far {err:.1e}", flush=True)
    return c, parts, nus, tuples, err


def singlet_index_exact(n, p, verbose=False):
    """Refined singlet index I_0 = sum_m (-1)^m n(3m) = n(t = -1) = Tr[R(-1)^3], with the edge transfer matrix
    R_{nu nu'}(t) = sum_P t^{|P|} c_nu(P) c_nu'(P^T) (dual Cauchy on each edge, Littlewood-Richardson at each node).
    Exact integer arithmetic for R(-1) and its cube."""
    c, parts, nus, tuples, err = lr_product_coefficients(n, p, verbose=verbose)
    idx = {l: i for i, l in enumerate(parts)}
    tidx = {t: i for i, t in enumerate(tuples)}
    tau = np.array([tidx[tuple(idx[transpose(parts[i], n)] for i in t)] for t in tuples])
    deg = np.array([sum(sum(parts[i]) for i in t) for t in tuples])
    sign = np.where(deg % 2 == 0, 1, -1).astype(np.int64)
    # R(-1) = c^T diag(sign) c[tau]   (entries checked to stay far below 2^63)
    R = (c * sign[:, None]).T @ c[tau]
    Ro = R.astype(object)
    R2 = Ro.dot(Ro)
    I = int(sum(R2[i, j] * Ro[j, i] for i in range(len(nus)) for j in range(len(nus))))
    return dict(index=I, n_tuples=len(tuples), n_nu=len(nus), lr_rounding_error=err, max_abs_R=int(np.abs(R).max()))


def dixon(N):
    """sum_m (-1)^m C(N, m)^3 (Dixon's identity: = (-1)^{N/2} (3N/2)! / ((N/2)!)^3 for even N, 0 for odd N)."""
    return sum((-1) ** m * math.comb(N, m) ** 3 for m in range(N + 1))


# ----------------------------------------------------------------------------------------------- flavour recursion

def lr_tensor(n, mu_cols, verbose=False):
    """L[mu, lam, nu] = c^nu_{mu lam} for mu in the n x mu_cols box, lam in the n x n box, nu in the n x (mu_cols + n)
    box, by exact torus quadrature on U(n) (exponents |e| <= mu_cols + n + n - 1; M = mu_cols + 2n points)."""
    mus = partitions_in_box(n, mu_cols)
    lams = box_partitions(n)
    nus = partitions_in_box(n, mu_cols + n)
    M = mu_cols + 2 * n
    rng = np.random.default_rng(1)
    shifts = rng.uniform(0, 2 * np.pi / M, size=n)
    grid = np.array(list(itertools.product(range(M), repeat=n)), dtype=float) * 2 * np.pi / M + shifts
    X = np.exp(1j * grid)
    _, den = schur_on_grid(tuple([0] * n), X)
    w = np.abs(den) ** 2 / (math.factorial(n) * M ** n)
    Sl = np.array([schur_on_grid(l, X)[0] for l in lams]) * w[None, :]
    Bn = np.array([schur_on_grid(v, X)[0] for v in nus]).conj().T             # (pts, n_nu)
    L = np.zeros((len(mus), len(lams), len(nus)), dtype=np.int64)
    err = 0.0
    for i, mu in enumerate(mus):
        smu = schur_on_grid(mu, X)[0]
        cc = (Sl * smu[None, :]) @ Bn
        rr = np.rint(cc.real)
        err = max(err, float(np.abs(cc - rr).max()))
        L[i] = rr.astype(np.int64)
        if verbose and (i % 50 == 0 or i == len(mus) - 1):
            print(f"    LR tensor n={n} mu_cols={mu_cols}: mu {i + 1}/{len(mus)}, rounding {err:.1e}", flush=True)
    return L, mus, lams, nus, err


def edge_transfer_minus1(n, p, verbose=False):
    """R^{(p)}(-1) over U(n) irreps nu in the n x pn box, by the flavour recursion
    R^{(q)} = sum_lam (-1)^{|lam|} L_lam^T R^{(q-1)} L_{lam^T},  R^{(1)}_{lam mu} = (-1)^{|lam|} delta_{mu, lam^T}."""
    lams = box_partitions(n)
    li = {l: i for i, l in enumerate(lams)}
    tr = np.array([li[transpose(l, n)] for l in lams])
    sg = np.array([(-1) ** sum(l) for l in lams], dtype=np.int64)
    R = np.zeros((len(lams), len(lams)), dtype=np.int64)
    R[np.arange(len(lams)), tr] = sg
    errs = []
    for q in range(2, p + 1):
        L, mus, _, nus, err = lr_tensor(n, (q - 1) * n, verbose=verbose)
        errs.append(err)
        Rn = np.zeros((len(nus), len(nus)), dtype=np.int64)
        for j in range(len(lams)):
            A = L[:, j, :]                                                   # (mu, nu)
            Bt = L[:, tr[j], :]
            Rn += sg[j] * (A.T @ (R @ Bt))
        R = Rn
        if verbose:
            print(f"    R^({q}): {R.shape[0]} irreps, max |entry| {int(np.abs(R).max())}", flush=True)
    return R, errs


def singlet_index_recursive(n, p, verbose=False):
    """I_0 = Tr[R^{(p)}(-1)^3] with exact integer arithmetic (int64 square if safe, Python-int final contraction)."""
    R, errs = edge_transfer_minus1(n, p, verbose=verbose)
    m = int(np.abs(R).max())
    if m * m * R.shape[0] < 2 ** 62:
        R2 = R @ R
        I = int(np.sum(R2.astype(object) * R.T.astype(object)))
    else:                                                                      # exact but slow fallback
        Ro = R.astype(object); R2 = Ro.dot(Ro)
        I = int(np.sum(R2 * Ro.T))
    return dict(index=I, n_nu=int(R.shape[0]), max_abs_R=m, lr_rounding_errors=errs)


def singlet_total_recursive(n, p, verbose=False):
    """Total number of U(n)^3 singlet states, n(t = +1) = Tr[R^{(p)}(+1)^3], with the same flavour recursion but all
    signs +1 (floating point: only the logarithm is used for large cases; exact below 2^53)."""
    lams = box_partitions(n)
    li = {l: i for i, l in enumerate(lams)}
    tr = np.array([li[transpose(l, n)] for l in lams])
    R = np.zeros((len(lams), len(lams)))
    R[np.arange(len(lams)), tr] = 1.0
    for q in range(2, p + 1):
        L, mus, _, nus, err = lr_tensor(n, (q - 1) * n, verbose=verbose)
        Lf = L.astype(float)
        Rn = np.zeros((len(nus), len(nus)))
        for j in range(len(lams)):
            Rn += Lf[:, j, :].T @ (R @ Lf[:, tr[j], :])
        R = Rn
    return float(np.trace(R @ R @ R))


def singlet_series_recursive(n, p, verbose=False):
    """
    Degree-resolved singlet counts n(k), k = 0..3pn^2, from the flavour transfer matrix with a fugacity t per fermion:
    R^{(q)}(t) = sum_lam t^{|lam|} L_lam^T R^{(q-1)}(t) L_{lam^T},  n(t) = Tr[R^{(p)}(t)^3].
    Evaluated at the K = 3pn^2 + 1 roots of unity and inverted by FFT (complex128), then rounded; returns
    (counts as Python ints, max rounding residual).  Cross-checks: sum = singlet_total_recursive, alternating sum over
    k = 3m with (-1)^m = singlet_index_recursive, palindromic.
    """
    lams = box_partitions(n)
    li = {l: i for i, l in enumerate(lams)}
    tr = np.array([li[transpose(l, n)] for l in lams])
    size = np.array([sum(l) for l in lams])
    K = 3 * p * n * n + 1
    ts = np.exp(2j * np.pi * np.arange(K) / K)
    Ls = [lr_tensor(n, (q - 1) * n, verbose=verbose) for q in range(2, p + 1)]
    vals = np.empty(K, dtype=complex)
    for it, t in enumerate(ts):
        w = t ** size
        R = np.zeros((len(lams), len(lams)), dtype=complex)
        R[np.arange(len(lams)), tr] = w
        for (L, mus, _, nus, err) in Ls:
            Lf = L.astype(float)
            Rn = np.zeros((len(nus), len(nus)), dtype=complex)
            for j in range(len(lams)):
                Rn += w[j] * (Lf[:, j, :].T @ (R @ Lf[:, tr[j], :]))
            R = Rn
        vals[it] = np.trace(R @ R @ R)
    coef = np.fft.fft(vals) / K                      # n(t) = sum_k c_k t^k  ->  c_k = (1/K) sum_j n(w^j) w^{-jk}
    counts = np.rint(coef.real)
    resid = float(max(np.abs(coef.real - counts).max(), np.abs(coef.imag).max()))
    return [int(c) for c in counts], resid
