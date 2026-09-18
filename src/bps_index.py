"""
Refined Witten index of the fermionic U(N) matrix models, per (N_Psi mod 3, U(N) irrep, Z_p flavor charge).

Fock space F = wedge^*(C^p (x) adj_{U(N)}); modes Psi^(m)_ij (Fourier flavor basis, Z_p charge m) with U(N) weight
u_i/u_j, N_Psi = 1.  Fock character
    Z(t, z; u) = prod_{m=0}^{p-1} prod_{i,j} (1 + t z^m u_i / u_j),
so the number of states with N_Psi = k, flavor charge w (mod p) and U(N) weight mu (an integer vector with sum 0)
is the coefficient of t^k z^w u^mu; irrep multiplicities follow by the Weyl alternating sum
    m(lambda) = sum_{sigma in S_N} sgn(sigma) W(lambda + delta - sigma(delta)),   delta = (N-1, ..., 1, 0).
Everything commutes with Q except N_Psi itself ([N_Psi, Q] = 3Q), so the index of the complex
(c = k mod 3, lambda, w) is the Euler characteristic
    I_{c,lambda,w} = sum_{k = c mod 3} (-1)^k n(k, lambda, w)
(sign convention (-1)^F with F = N_Psi; the overall sign of a complex is immaterial).  |I| is a lower bound on the
number of BPS states of the complex, with equality iff the cohomology sits in a single degree ("index saturation",
Chang-Chen-Sia-Yang 2024 sect. 5.1) -- which is what R-charge concentration asserts.

Exactness: weight multiplicities are built by exact integer arithmetic (the generating function is multiplied by
the p N^2 linear factors one at a time; exponents are stored modulo M > 4p(N-1), which is alias-free).  Feasible
for N <= 4 with int64 arrays (N=4: 49 x 3 x 37^3 entries); N=5 would need ~10 GB.  U(N) irreps are labelled by lambda = (lambda_1 >= ... >= lambda_N),
sum lambda_i = 0 (charge-0 irreps; for N=2, lambda = (j, -j) is spin j).
"""
import itertools
import numpy as np


def _weights_exact(N, p, M=None):
    """W[k, w, e_1, ..., e_{N-1}] (int64): number of Fock states with N_Psi = k, Z_p charge w and SU(N)-torus
    exponents e (exponent of u_i, i<N, after u_N = 1/(u_1...u_{N-1}); stored modulo M with M > 4p(N-1), so no
    aliasing).  Exact integer arithmetic: the generating function is built by multiplying the p N^2 linear factors
    (1 + t z^m u^v) one at a time, each being a shift of the array along (k, w, e)."""
    n_modes = p * N * N
    if M is None:
        M = 4 * p * (N - 1) + 1
    shape = (n_modes + 1, p) + (M,) * (N - 1)
    W = np.zeros(shape, dtype=np.int64); W[(0, 0) + (0,) * (N - 1)] = 1
    for m in range(p):
        for i in range(N):
            for j in range(N):
                mu = np.zeros(N, dtype=int); mu[i] += 1; mu[j] -= 1
                e = _exponents_of_weight(mu)
                shifted = np.roll(W, 1, axis=0)                 # t
                shifted[0] = 0
                shifted = np.roll(shifted, m, axis=1)           # z^m  (Z_p charge, cyclic)
                for ax, sh in enumerate(e):
                    if sh:
                        shifted = np.roll(shifted, sh, axis=2 + ax)   # u^e, cyclic = modulo M (no aliasing)
                W = W + shifted
    return W, M


_weights_by_dft = _weights_exact      # kept name for callers


def _exponents_of_weight(mu):
    """SU(N)-torus exponents of the U(N) weight mu (sum 0): u^mu = prod u_i^{mu_i}; with u_N = 1/(u_1...u_{N-1}),
    u^mu = prod_{i<N} u_i^{mu_i - mu_N}."""
    return tuple(int(m - mu[-1]) for m in mu[:-1])


def irrep_multiplicities(N, p=3, verbose=False):
    """Dict {(k, w, lambda): multiplicity} for all U(N) charge-0 irreps lambda appearing in the Fock space."""
    W, M = _weights_by_dft(N, p)
    n_modes = p * N * N
    delta = np.arange(N - 1, -1, -1)
    lam_max = p * (N - 1)
    # candidate highest weights: lambda_1 >= ... >= lambda_N, sum 0, lambda_1 <= lam_max
    lams = []
    for parts in itertools.product(range(-lam_max, lam_max + 1), repeat=N - 1):
        lam = list(parts) + [-sum(parts)]
        if all(lam[i] >= lam[i + 1] for i in range(N - 1)) and lam[0] <= lam_max and lam[-1] >= -lam_max:
            lams.append(tuple(lam))
    perms = list(itertools.permutations(range(N)))
    def sgn(perm):
        s = 1
        for a in range(N):
            for b in range(a + 1, N):
                if perm[a] > perm[b]:
                    s = -s
        return s
    def Wget(k, wq, mu):
        e = _exponents_of_weight(mu)
        idx = tuple(x % M for x in e)
        return W[(k, wq) + idx]
    out = {}
    for k in range(n_modes + 1):
        for wq in range(p):
            for lam in lams:
                m = 0
                for perm in perms:
                    mu = np.array(lam) + delta - delta[list(perm)]
                    m += sgn(perm) * Wget(k, wq, mu)
                if m != 0:
                    assert m > 0, (k, wq, lam, m)
                    out[(k, wq, lam)] = int(m)
    if verbose:
        tot = sum(m * irrep_dim(N, lam) for (k, wq, lam), m in out.items())
        print(f"  N={N} p={p}: {len(out)} (k,w,lambda) cells, total dimension {tot} (should be 2^{n_modes} = {2**n_modes})")
    return out


def irrep_dim(N, lam):
    """Weyl dimension formula for U(N)."""
    lam = np.array(lam, dtype=float); d = 1.0
    for i in range(N):
        for j in range(i + 1, N):
            d *= (lam[i] - lam[j] + j - i) / (j - i)
    return int(round(d))


def refined_index(mult, p=3):
    """I_{c, lambda, w} = sum_{k = c mod 3} (-1)^k n(k, lambda, w), from the multiplicities dict."""
    out = {}
    for (k, wq, lam), m in mult.items():
        key = (k % 3, wq, lam)
        out[key] = out.get(key, 0) + (-1) ** k * m
    return {key: v for key, v in out.items() if v != 0}


def class_totals(N, p=3):
    """Closed form for the total (all irreps, all flavors) index per class c = k mod 3, from
    sum_k (-1)^k C(n,k) z^k = (1 - z)^n with n = p N^2:  I_c = (1/3) sum_r w^{-rc} (1 - w^r)^n."""
    n = p * N * N; w = np.exp(2j * np.pi / 3)
    return [int(round(sum(w ** (-r * c) * (1 - w ** r) ** n for r in range(3)).real / 3)) for c in range(3)]
