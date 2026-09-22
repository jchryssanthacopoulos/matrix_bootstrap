"""
Exact Q-cohomology of the fermionic matrix models, computed weight space by weight space.

The BPS states are H_Q = ker Q / im Q in each degree k = N_Psi.  Working in the full sector is hopeless
(dim H_13 = 2 x 10^7 at N = 3), but Q is a gauge singlet and preserves the U(N) weight, so the complex splits
over weights:

    ... --Q--> W_lambda(k) --Q--> W_lambda(k+3) --Q--> ...,     W_lambda(k) = weight-lambda subspace of H_k,

and the weight spaces are small (9, 126, 36 for lambda = (6,0,-6) at N = 3, k = 10, 13, 16).  Two uses:

* For the **maximal** weight of the Fock space there is no higher weight, so the weight-space complex IS the
  multiplicity complex of that irrep (summed over the Z_p flavour charges): dim H^k(W_lambda) = sum_w m(k,lambda,w)
  restricted to cohomology, i.e. exactly the number of BPS multiplets of lambda at degree k.
* For lower weights, dim H^k(W_lambda) = sum_{mu >= lambda} K(mu, lambda) dim H^k_mu with K the Kostant
  multiplicity (dim of the weight-lambda space of the irrep mu); the triangular system is inverted from the top
  down, and the K's are read off from the same weight-space dimensions (dim W_lambda(k) = sum_mu m(k,mu) K(mu,lambda)).

Conventions match fermion_matrix_model / bps_index: modes are (a, i, j) with a the flavour and (i, j) the matrix
indices, weight e_i - e_j, and Q = sum_abc C_abc Tr[Psi^a Psi^b Psi^c] = sum_abc C_abc sum_ijk Psi^a_ij Psi^b_jk
Psi^c_ki.  States are sorted tuples of mode indices with the standard fermionic sign; ranks are computed exactly
by fraction-free (Bareiss-style) elimination over a large prime.
"""
import itertools
import numpy as np
from math import comb

PRIME = (1 << 31) - 1          # 2^31 - 1, Mersenne: exact ranks via elimination mod p (no integer overflow in int64)


def modes(N, p):
    """All modes (a, i, j) in a fixed order; index = a*N*N + i*N + j."""
    return [(a, i, j) for a in range(p) for i in range(N) for j in range(N)]


def mode_index(N, a, i, j):
    return a * N * N + i * N + j


def weight_of(N, p, idx):
    a, rem = divmod(idx, N * N); i, j = divmod(rem, N)
    w = [0] * N; w[i] += 1; w[j] -= 1
    return tuple(w)


def weight_basis(N, p, k, lam):
    """Sorted tuples of k mode indices whose weights sum to lam.  Enumerated through the occupation matrix
    n_ij (how many flavours of mode (.,i,j) are filled), which makes the search tiny."""
    lam = tuple(lam)
    out = []
    for flat in itertools.product(range(p + 1), repeat=N * N):
        n = np.array(flat).reshape(N, N)
        if n.sum() != k:
            continue
        if tuple(n.sum(axis=1) - n.sum(axis=0)) != lam:
            continue
        choices = []
        for i in range(N):
            for j in range(N):
                cnt = n[i, j]
                choices.append([tuple(sorted(c)) for c in itertools.combinations(range(p), cnt)])
        for combo in itertools.product(*choices):
            state = []
            pos = 0
            for i in range(N):
                for j in range(N):
                    for a in combo[pos]:
                        state.append(mode_index(N, a, i, j))
                    pos += 1
            out.append(tuple(sorted(state)))
    return out


def _create(state, m):
    """c^dag_m |state>: returns (sign, new_state) or None if already occupied."""
    if m in state:
        return None
    pos = 0
    while pos < len(state) and state[pos] < m:
        pos += 1
    return (-1) ** pos, state[:pos] + (m,) + state[pos:]


def apply_Q(N, p, C, state):
    """Q|state> as {new_state: coefficient}."""
    out = {}
    for a in range(p):
        for b in range(p):
            for c in range(p):
                coef = C[a, b, c]
                if coef == 0:
                    continue
                for i in range(N):
                    for j in range(N):
                        for kk in range(N):
                            ms = (mode_index(N, a, i, j), mode_index(N, b, j, kk), mode_index(N, c, kk, i))
                            s, st = 1, state
                            ok = True
                            for m in reversed(ms):          # rightmost operator acts first
                                r = _create(st, m)
                                if r is None:
                                    ok = False; break
                                sg, st = r; s *= sg
                            if ok:
                                out[st] = out.get(st, 0) + s * coef
    return {st: v for st, v in out.items() if v != 0}


def Q_matrix(N, p, C, k, lam, basis_k=None, basis_k3=None):
    """Matrix of Q : W_lambda(k) -> W_lambda(k+3) with integer entries (C must be integral up to a common factor)."""
    bk = basis_k if basis_k is not None else weight_basis(N, p, k, lam)
    b3 = basis_k3 if basis_k3 is not None else weight_basis(N, p, k + 3, lam)
    index = {s: r for r, s in enumerate(b3)}
    M = np.zeros((len(b3), len(bk)), dtype=np.int64)
    for col, st in enumerate(bk):
        for st2, v in apply_Q(N, p, C, st).items():
            M[index[st2], col] += int(round(v))
    return M, bk, b3


def rank_mod_p(M, prime=PRIME):
    """Exact rank of an integer matrix modulo a large prime (Gaussian elimination in int64)."""
    A = (np.asarray(M, dtype=np.int64) % prime).copy()
    rows, cols = A.shape
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i, c]:
                piv = i; break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), prime - 2, prime)
        A[r] = (A[r] * inv) % prime
        nz = np.nonzero(A[r + 1:, c])[0]
        if len(nz):
            A[r + 1:][nz] = (A[r + 1:][nz] - np.outer(A[r + 1:][nz, c], A[r]) ) % prime
        r += 1
        if r == rows:
            break
    return r


def complex_cohomology(N, p, C, lam, verbose=False):
    """Dimensions of H^k for the weight-lambda complex, all k.  Returns dict k -> (dim W, rank Q_k, dim H^k)."""
    n_modes = p * N * N
    bases = {}
    for k in range(n_modes + 1):
        b = weight_basis(N, p, k, lam)
        if b:
            bases[k] = b
    ranks = {}
    for k in sorted(bases):
        if k + 3 in bases:
            M, _, _ = Q_matrix(N, p, C, k, lam, bases[k], bases[k + 3])
            ranks[k] = rank_mod_p(M)
        else:
            ranks[k] = 0
        if verbose:
            print(f"    k={k}: dim W={len(bases[k])}, rank Q_k={ranks[k]}", flush=True)
    out = {}
    for k in sorted(bases):
        dim = len(bases[k])
        out[k] = (dim, ranks[k], dim - ranks[k] - ranks.get(k - 3, 0))
    return out


def flavour_permutation(N, p, state):
    """Z_p flavour rotation a -> a+1 on a state: returns (sign, new_state)."""
    new = sorted(mode_index(N, (idx // (N * N) + 1) % p, *divmod(idx % (N * N), N)) for idx in state)
    # fermionic sign of the permutation taking the image list (in the old order) to sorted order
    img = [mode_index(N, (idx // (N * N) + 1) % p, *divmod(idx % (N * N), N)) for idx in state]
    perm = sorted(range(len(img)), key=lambda t: img[t])
    sign, seen = 1, [False] * len(perm)
    for st in range(len(perm)):
        if seen[st]:
            continue
        ln, i = 0, st
        while not seen[i]:
            seen[i] = True; i = perm[i]; ln += 1
        if ln % 2 == 0:
            sign = -sign
    return sign, tuple(new)


def flavour_blocks(N, p, basis):
    """Projectors onto the Z_p flavour eigenspaces on a weight basis: returns list of (charge, P) with P a
    complex matrix whose columns span the eigenspace of the rotation with eigenvalue exp(2 pi i w / p)."""
    idx = {s: r for r, s in enumerate(basis)}
    U = np.zeros((len(basis), len(basis)), dtype=complex)
    for c, st in enumerate(basis):
        sg, st2 = flavour_permutation(N, p, st)
        U[idx[st2], c] = sg
    out = []
    for w in range(p):
        ev = np.exp(2j * np.pi * w / p)
        M = np.eye(len(basis), dtype=complex)
        for t in range(1, p):                       # projector (1/p) sum_t ev^{-t} U^t
            M = M + (ev ** (-t)) * np.linalg.matrix_power(U, t)
        M /= p
        rk = np.linalg.matrix_rank(M, tol=1e-8)
        if rk:
            u, sv, _ = np.linalg.svd(M)
            out.append((w, u[:, :rk]))
    return out


def complex_cohomology_flavour(N, p, C, lam, verbose=False):
    """Cohomology of the weight-lambda complex, refined by Z_p flavour charge."""
    n_modes = p * N * N
    bases, blocks = {}, {}
    for k in range(n_modes + 1):
        b = weight_basis(N, p, k, lam)
        if b:
            bases[k] = b; blocks[k] = dict(flavour_blocks(N, p, b))
    charges = sorted({w for k in blocks for w in blocks[k]})
    out = {}
    for w in charges:
        ranks, dims = {}, {}
        for k in sorted(bases):
            if w not in blocks[k]:
                dims[k] = 0; ranks[k] = 0; continue
            Pk = blocks[k][w]; dims[k] = Pk.shape[1]
            if k + 3 in bases and w in blocks[k + 3]:
                M, _, _ = Q_matrix(N, p, C, k, lam, bases[k], bases[k + 3])
                Mw = blocks[k + 3][w].conj().T @ M @ Pk
                ranks[k] = int(np.linalg.matrix_rank(Mw, tol=1e-7))
            else:
                ranks[k] = 0
        out[w] = {k: (dims[k], ranks[k], dims[k] - ranks[k] - ranks.get(k - 3, 0)) for k in sorted(bases)}
        if verbose:
            print(f"  flavour charge w={w}: H^k = { {k: v[2] for k, v in out[w].items() if v[2]} }", flush=True)
    return out


# ---------------------------------------------------------------------------------------------------------------
# Kostant weight multiplicities and the top-down inversion of the weight-space cohomology
# ---------------------------------------------------------------------------------------------------------------

def partition_function(v, roots, _memo=None):
    """Kostant partition function P(v): the number of ways to write v as a non-negative integer combination of the
    positive roots.  Recursion over the root list, pruned with the fact that every non-negative combination of
    positive roots e_i - e_j (i < j) has non-negative partial sums."""
    v = tuple(int(x) for x in v)
    if sum(v) != 0:
        return 0
    if _memo is None:
        _memo = {}

    def ok(u):
        s = 0
        for x in u:
            s += x
            if s < 0:
                return False
        return True

    def rec(u, i):
        if all(x == 0 for x in u):
            return 1
        if i == len(roots) or not ok(u):
            return 0
        key = (u, i)
        if key in _memo:
            return _memo[key]
        total, w = 0, u
        while True:
            total += rec(w, i + 1)
            w = tuple(a - b for a, b in zip(w, roots[i]))
            if not ok(w):
                break
        _memo[key] = total
        return total

    return rec(v, 0)


def positive_roots(N):
    out = []
    for i in range(N):
        for j in range(i + 1, N):
            r = [0] * N; r[i] = 1; r[j] = -1
            out.append(tuple(r))
    return out


def kostant(N, mu, lam):
    """Multiplicity of the weight lam in the U(N) irrep with highest weight mu (Kostant's formula)."""
    mu, lam = np.array(mu, dtype=int), np.array(lam, dtype=int)
    if mu.sum() != lam.sum():
        return 0
    rho = np.arange(N - 1, -1, -1)
    roots = positive_roots(N)
    total = 0
    for perm in itertools.permutations(range(N)):
        s = 1
        for a in range(N):
            for b in range(a + 1, N):
                if perm[a] > perm[b]:
                    s = -s
        v = (mu + rho)[list(perm)] - (lam + rho)
        if all(x == 0 for x in v):
            total += s
        elif sum(v) == 0:
            total += s * partition_function(tuple(int(x) for x in v), roots)
    return total


def peel(N, H_by_weight, order=None):
    """Invert dim H^k(W_lam) = sum_{mu >= lam} K(mu, lam) h^k_mu from the top down.

    H_by_weight maps a dominant weight lam -> {k: dim H^k(W_lam)}.  The set of weights must be downward-closed in
    the sense that every mu > lam appearing with non-zero cohomology is itself present (true when the weights are
    the highest ones of the Fock space, since higher weights have smaller weight spaces).  Returns {lam: {k: h^k}}.
    """
    if order is None:
        order = sorted(H_by_weight, key=lambda l: -sum(x * x for x in l))
    h = {}
    for lam in order:
        hh = dict(H_by_weight[lam])
        for mu in order:
            if mu == lam or mu not in h:
                continue
            kmu = kostant(N, mu, lam)
            if kmu:
                for k, v in h[mu].items():
                    hh[k] = hh.get(k, 0) - kmu * v
        h[lam] = {k: v for k, v in hh.items() if v}
    return h
