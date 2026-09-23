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


def occupation_matrices(N, p, k, lam):
    """All occupation matrices n_ij in 0..p with sum n = k and row-sums minus column-sums = lam.

    Enumerated row by row with pruning, which is what makes N >= 4 reachable: the naive product over all N^2
    cells is (p+1)^(N^2) = 4.3e9 at N=4, p=3, while the pruned search visits only what survives the constraints
    (the maximal weight at N=4 has just 4096 basis states in total)."""
    lam = tuple(lam)
    rows_choices = [c for c in itertools.product(range(p + 1), repeat=N)]
    out = []
    chosen, colsum = [], [0] * N

    def rec(i, used):
        if used > k or k - used > p * N * (N - i):
            return
        # a fixed row i' < i needs its column to finish at r_{i'} - lam_{i'}; later rows can add 0..p*(N-i)
        for t in range(i):
            need = chosen[t][1] - lam[t] - colsum[t]
            if need < 0 or need > p * (N - i):
                return
        # a row not yet chosen will have r_t = colsum_t(final) + lam_t, which must be a legal row sum
        for t in range(i, N):
            if colsum[t] + lam[t] > p * N:
                return
        if i == N:
            if used == k and all(chosen[t][1] - lam[t] == colsum[t] for t in range(N)):
                out.append(tuple(r for r, _ in chosen))
            return
        for row in rows_choices:
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
    """Sorted tuples of k mode indices whose weights sum to lam.  The search runs over the occupation matrix
    n_ij (how many flavours of mode (.,i,j) are filled), then over which flavours fill each cell."""
    out = []
    for n in occupation_matrices(N, p, k, lam):
        choices = [[tuple(sorted(c)) for c in itertools.combinations(range(p), n[i][j])]
                   for i in range(N) for j in range(N)]
        for combo in itertools.product(*choices):
            state, pos = [], 0
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


def exact_rank(M, cross_check=False):
    """Exact rank of an integer matrix: blocked (BLAS) elimination once the matrix is big enough for it to pay,
    plain elimination otherwise.  With cross_check, the rank is recomputed over a second prime and the two must
    agree (a rank over F_p can only under-estimate the rank over Q)."""
    n = min(M.shape)
    if n == 0:
        return 0
    if n < 1200:
        r = rank_mod_p(M)
        if cross_check:
            assert r == rank_mod_p(M, 2147483629), "rank disagreed between primes"
        return r
    r = rank_mod_p_blocked(M, SMALL_PRIMES[0], block=128)
    if cross_check:
        assert r == rank_mod_p_blocked(M, SMALL_PRIMES[1], block=128), "rank disagreed between primes"
    return r


def complex_cohomology(N, p, C, lam, verbose=False, cross_check=False):
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
            ranks[k] = exact_rank(M, cross_check=cross_check)
            del M
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


# ---------------------------------------------------------------------------------------------------------------
# Blocked exact rank: the same elimination, but with the trailing update done by BLAS
# ---------------------------------------------------------------------------------------------------------------

# Primes small enough that block * (p-1)^2 < 2^53, so a float64 matmul of inner dimension <= block is EXACT
# integer arithmetic and can be reduced mod p afterwards.  Two of them, to cross-check ranks.
SMALL_PRIMES = (1048573, 1048571)


def rank_mod_p_blocked(M, prime=SMALL_PRIMES[0], block=256, copy=True):
    """Exact rank of an integer matrix over F_prime, by right-looking blocked Gaussian elimination (LU with
    partial pivoting, skipping rank-deficient columns).

    Identical mathematics to `rank_mod_p`, but the O(n^3) work is moved into one `L21 @ U12` matmul per panel,
    which numpy hands to BLAS.  That is the difference between hours and minutes at n ~ 10^4.  Exactness of the
    float64 arithmetic requires block * (prime - 1)^2 < 2^53 (asserted): each product of two reduced residues is
    < (prime-1)^2 and at most `block` of them are summed before the reduction mod prime.

    A rank over F_prime is a lower bound on the rank over Q, with equality unless the prime divides the relevant
    minors; callers should confirm with a second prime (see SMALL_PRIMES).

    With copy=False and a float64 input whose entries are already reduced mod `prime`, the elimination runs in
    place and DESTROYS the input; this avoids two full-size temporaries, which matters at n ~ 10^4 (1.3 GB each).
    """
    assert block * (prime - 1) ** 2 < 2 ** 53, "float64 matmul would not be exact"
    A = np.asarray(M)
    if A.dtype == np.float64 and not copy:
        pass                                   # already reduced mod prime; eliminate in place, no extra copies
    elif A.dtype == np.float64:
        A = A.copy()
    else:
        A = (A.astype(np.int64) % prime).astype(np.float64)
    m, n = A.shape
    r, j = 0, 0
    while r < m and j < n:
        b = min(block, n - j)
        r0, pcols = r, []
        # --- panel: eliminate inside columns [j, j+b) only, leaving the unit multipliers in place of the zeros ---
        for c in range(j, j + b):
            if r >= m:
                break
            nz = np.nonzero(A[r:, c])[0]
            if len(nz) == 0:
                continue                      # rank-deficient column: no pivot, nothing to do
            piv = r + int(nz[0])
            if piv != r:
                A[[r, piv]] = A[[piv, r]]     # full-row swap, so multipliers already stored travel with the row
            inv = pow(int(A[r, c]), prime - 2, prime)
            col = A[r + 1:, c]                # becomes the column of unit multipliers, stored in place
            np.multiply(col, inv, out=col)
            np.mod(col, prime, out=col)
            if c + 1 < j + b:                 # eliminate, inside the panel only; no fancy indexing (it is slow)
                blk = A[r + 1:, c + 1:j + b]
                blk -= np.outer(col, A[r, c + 1:j + b])
                np.mod(blk, prime, out=blk)
            pcols.append(c)
            r += 1
        k = r - r0
        # --- trailing update: U12 = L11^{-1} A12, then A22 -= L21 @ U12 ---
        if k and j + b < n:
            U = A[r0:r0 + k, j + b:].copy()
            L11 = A[r0:r0 + k, pcols]                 # unit lower triangular in its strict lower part
            for t in range(k - 1):
                col = L11[t + 1:, t]
                nzt = np.nonzero(col)[0]
                if len(nzt):
                    U[t + 1:][nzt] = (U[t + 1:][nzt] - np.outer(col[nzt], U[t])) % prime
            A[r0:r0 + k, j + b:] = U
            if r < m:                         # the one big matmul; in place, to avoid n^2 temporaries
                A22 = A[r:, j + b:]
                A22 -= A[r:, pcols] @ U
                np.mod(A22, prime, out=A22)
        j += b
    return r


# ---------------------------------------------------------------------------------------------------------------
# Z_p flavour blocking: splits each weight space into p blocks of ~n/p, i.e. p^2-fold less rank work
# ---------------------------------------------------------------------------------------------------------------

def _is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def blocking_primes(p=3, count=2, below=1 << 20):
    """Primes P just under `below` with P = 1 (mod p), so that F_P contains a primitive p-th root of unity and
    the flavour eigenbasis can be written down exactly over F_P."""
    out, n = [], below - 1
    while len(out) < count and n > 2:
        if n % p == 1 and _is_prime(n):
            out.append(n)
        n -= 1
    return out


def root_of_unity(prime, p=3):
    """A primitive p-th root of unity in F_prime (requires prime = 1 mod p)."""
    assert prime % p == 1
    for g in range(2, prime):
        w = pow(g, (prime - 1) // p, prime)
        if w != 1 and pow(w, p, prime) == 1:
            return w
    raise RuntimeError("no primitive root found")


def flavour_orbits(N, p, basis):
    """Orbits of the Z_p flavour rotation on a weight basis, carrying the fermionic signs.

    Returns a list of orbits; each is a list of (index_into_basis, d_t) of length 1 or p, where d_t is the sign
    defined by sigma^t |s_0> = d_t |s_t> (so d_0 = 1).  A fixed point must have d = +1, since sigma^p = 1."""
    index = {s: i for i, s in enumerate(basis)}
    seen = [False] * len(basis)
    orbits = []
    for i, s in enumerate(basis):
        if seen[i]:
            continue
        orb, st, d = [], s, 1
        for _ in range(p):
            orb.append((index[st], d))
            seen[index[st]] = True
            sg, st = flavour_permutation(N, p, st)
            d = d * sg
            if st == s:
                break
        assert d == 1, "sigma^p must act as the identity with sign +1"
        orbits.append(orb)
    return orbits


def Q_matrix_flavour(N, p, C, k, lam, basis_k, basis_k3, prime, orbits_k=None, orbits_k3=None):
    """The p flavour blocks of Q : W_lambda(k) -> W_lambda(k+3), as matrices over F_prime.

    Q commutes with the flavour rotation, so it is block diagonal in the eigenbasis.  For a source orbit with
    representative s_0, writing Q|s_0> = sum_u a_u |u>, the entry into a target orbit with signs d_r is

        (Q_w)_{O', O} = sum_r a_{u_r} d_r omega^{w r},

    derived from |u_r> = (d_r / p) sum_w omega^{w r} v'_w and Q v_w = sum_t omega^{-w t} sigma^t (Q|s_0>).
    Orbits of size 1 exist only in the w = 0 block.  Returns a list of p matrices (int64, entries in F_prime).
    """
    ok = orbits_k if orbits_k is not None else flavour_orbits(N, p, basis_k)
    o3 = orbits_k3 if orbits_k3 is not None else flavour_orbits(N, p, basis_k3)
    omega = root_of_unity(prime, p)
    idx3 = {}                                   # basis index -> (orbit number, position r, sign d_r)
    for oi, orb in enumerate(o3):
        for r, (bi_, d) in enumerate(orb):
            idx3[bi_] = (oi, r, d)
    # column / row numbering inside each block
    cols = [{} for _ in range(p)]
    for oi, orb in enumerate(ok):
        for w in (range(p) if len(orb) == p else [0]):
            cols[w][oi] = len(cols[w])
    rows = [{} for _ in range(p)]
    for oi, orb in enumerate(o3):
        for w in (range(p) if len(orb) == p else [0]):
            rows[w][oi] = len(rows[w])
    mats = [np.zeros((len(rows[w]), len(cols[w])), dtype=np.float64) for w in range(p)]
    for oi, orb in enumerate(ok):
        s0 = basis_k[orb[0][0]]
        a = apply_Q(N, p, C, s0)
        contrib = {}                            # (target orbit, w) -> accumulated entry
        for st, v in a.items():
            oj, r, d = idx3[index_of(basis_k3, st)]
            for w in (range(p) if len(o3[oj]) == p else [0]):
                key = (oj, w)
                contrib[key] = (contrib.get(key, 0) + int(v) * d * pow(omega, w * r, prime)) % prime
        for (oj, w), val in contrib.items():
            if val and oi in cols[w] and oj in rows[w]:
                mats[w][rows[w][oj], cols[w][oi]] = val
    return mats, cols, rows


_index_cache = {}


def index_of(basis, state):
    """Index of a state in a basis list (cached per basis object)."""
    key = id(basis)
    d = _index_cache.get(key)
    if d is None or len(d) != len(basis):
        d = {s: i for i, s in enumerate(basis)}
        _index_cache.clear()
        _index_cache[key] = d
    return d[state]


def complex_cohomology_blocked(N, p, C, lam, prime=None, verbose=False):
    """Cohomology of the weight-lambda complex, computed one flavour block at a time.

    Same answer as `complex_cohomology` (summed over the blocks) at ~p^2 less cost and ~p^2 less memory, and it
    additionally resolves the Z_p flavour charge.  Returns (per_charge, totals) where per_charge[w][k] = dim H^k
    in that block and totals[k] = (dim W, rank Q_k, dim H^k)."""
    prime = prime or blocking_primes(p, 1)[0]
    bases, orbits = {}, {}
    for k in range(p * N * N + 1):
        b = weight_basis(N, p, k, lam)
        if b:
            bases[k] = b
            orbits[k] = flavour_orbits(N, p, b)
    ranks = {w: {} for w in range(p)}
    dims = {w: {} for w in range(p)}
    for k in sorted(bases):
        nblk = [sum(1 for o in orbits[k] if len(o) == p or w == 0) for w in range(p)]
        for w in range(p):
            dims[w][k] = nblk[w]
        if k + 3 in bases:
            mats, _, _ = Q_matrix_flavour(N, p, C, k, lam, bases[k], bases[k + 3], prime,
                                          orbits[k], orbits[k + 3])
            for w in range(p):                # rank in place, releasing each block as it is consumed
                ranks[w][k] = rank_mod_p_blocked(mats[w], prime, block=128, copy=False) if mats[w].size else 0
                mats[w] = None
            del mats
        else:
            for w in range(p):
                ranks[w][k] = 0
        if verbose:
            print(f"    k={k}: dim W={len(bases[k])} -> blocks {nblk}, ranks {[ranks[w][k] for w in range(p)]}",
                  flush=True)
    per_charge = {w: {k: dims[w][k] - ranks[w][k] - ranks[w].get(k - 3, 0) for k in sorted(bases)}
                  for w in range(p)}
    totals = {k: (len(bases[k]),
                  sum(ranks[w][k] for w in range(p)),
                  sum(per_charge[w][k] for w in range(p))) for k in sorted(bases)}
    return per_charge, totals


# ---------------------------------------------------------------------------------------------------------------
# Screening: supercharges of arbitrary odd degree q,  Q = sum C_{a1..aq} Tr[Psi^{a1} ... Psi^{aq}]
# ---------------------------------------------------------------------------------------------------------------

def apply_Q_degree(N, p, C, state, q):
    """Q|state> for a degree-q supercharge.  C is a rank-q flavour tensor (C = 1 for p = 1).

    Only ODD q is admissible: Tr[Psi^q] picks up (-1)^{q-1} under a cyclic shift, so it vanishes identically for
    even q; and Q^2 = 0 is automatic exactly when q is odd, since exchanging the two q-tuples in C (x) C is a
    permutation of 2q fermions with q^2 transpositions, odd iff q is odd, against a totally antisymmetric product.
    """
    assert q % 2 == 1, "only odd q gives a non-zero nilpotent Tr[Psi^q]"
    out = {}
    for flav in itertools.product(range(p), repeat=q):
        coef = C[flav] if p > 1 else 1.0
        if coef == 0:
            continue
        for loop in itertools.product(range(N), repeat=q):
            ms = tuple(mode_index(N, flav[t], loop[t], loop[(t + 1) % q]) for t in range(q))
            s, st, ok = 1, state, True
            for m in reversed(ms):
                r = _create(st, m)
                if r is None:
                    ok = False; break
                sg, st = r; s *= sg
            if ok:
                out[st] = out.get(st, 0) + s * coef
    return {st: v for st, v in out.items() if v != 0}


def maximal_weight_window(N, p, C, q, prime=PRIME):
    """Cohomology of the maximal-weight complex for a degree-q supercharge.  Returns (window, h, dims).

    The grading is by k mod q (since [N_Psi, Q] = q Q), so concentration on this complex requires the window to
    contain at most one degree per residue class.
    """
    lam = tuple(p * (N + 1 - 2 * i) for i in range(1, N + 1))
    n_modes = p * N * N
    bases = {}
    for k in range(n_modes + 1):
        b = weight_basis(N, p, k, lam)
        if b:
            bases[k] = b
    ranks = {}
    for k in sorted(bases):
        if k + q in bases:
            idx = {s: r for r, s in enumerate(bases[k + q])}
            M = np.zeros((len(bases[k + q]), len(bases[k])), dtype=np.int64)
            for col, st in enumerate(bases[k]):
                for st2, v in apply_Q_degree(N, p, C, st, q).items():
                    M[idx[st2], col] += int(round(v))
            r1 = rank_mod_p(M, prime)
            assert r1 == rank_mod_p(M, 2147483629), "rank disagreed between primes"
            ranks[k] = r1
        else:
            ranks[k] = 0
    h = {k: len(bases[k]) - ranks[k] - ranks.get(k - q, 0) for k in sorted(bases)}
    h = {k: v for k, v in h.items() if v}
    return sorted(h), h, {k: len(b) for k, b in bases.items()}
