"""Build delta = beta ^ - explicitly on the cohomology of alpha, and measure its rank.

The mapping-cone recursion for a generic q-form omega_{n+1} = alpha + e_{n+1} ^ beta reads

    dim H^m_{n+1} = dim coker( delta : H^{m-q}_n -> H^{m-1}_n ) + dim ker( delta : H^m_n -> H^{m+q-1}_n ),

with delta = beta ^ - the induced map on H(alpha).  Everything in the recursion is forced except rank(delta),
so this script computes that rank directly rather than inferring it from measured Betti numbers.

    rank(delta : H^j -> H^{j+q-1}) = rank[ B_j N_j | A_{j-1} ] - rank A_{j-1},

where N_j is a basis of ker(alpha ^ -) on Lambda^j, B_j is beta ^ - on Lambda^j, and A_{j-1} is the matrix of
alpha ^ - into Lambda^{j+q-1} (whose image is what we quotient by).  All linear algebra is exact over F_P.

    python scripts/run_delta_corank.py --q 3 --nmax 13
"""

import argparse
import itertools
from math import comb
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cohomology import PRIME  # noqa: E402

ALT_PRIME = 2147483629


def _rref(M, p):
    M = np.array(M % p, dtype=np.int64, copy=True)
    rows, cols = M.shape
    piv, r = [], 0
    for c in range(cols):
        if r >= rows:
            break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0:
            continue
        i = r + int(nz[0])
        if i != r:
            M[[r, i]] = M[[i, r]]
        M[r] = (M[r] * pow(int(M[r, c]), p - 2, p)) % p
        col = M[:, c].copy()
        col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr):
            M[nzr] = (M[nzr] - np.outer(col[nzr], M[r])) % p
        piv.append(c)
        r += 1
    return M, piv


def _rank(M, p):
    if M.size == 0 or min(M.shape) == 0:
        return 0
    return len(_rref(M, p)[1])


def _nullspace(M, p):
    """Basis of ker M as columns."""
    cols = M.shape[1]
    if M.size == 0 or M.shape[0] == 0:
        return np.eye(cols, dtype=np.int64)
    R, piv = _rref(M, p)
    pivset = set(piv)
    free = [c for c in range(cols) if c not in pivset]
    N = np.zeros((cols, len(free)), dtype=np.int64)
    for k, f in enumerate(free):
        N[f, k] = 1
        for i, pc in enumerate(piv):
            N[pc, k] = (-int(R[i, f])) % p
    return N


def wedge_matrix(n, deg, form, src_k):
    """Matrix of (form ^ -) : Lambda^{src_k} C^n -> Lambda^{src_k+deg} C^n."""
    src = list(itertools.combinations(range(n), src_k))
    tgt = list(itertools.combinations(range(n), src_k + deg))
    idx = {t: i for i, t in enumerate(tgt)}
    M = np.zeros((len(tgt), len(src)), dtype=np.int64)
    for col, s in enumerate(src):
        ss = set(s)
        for t, c in form.items():
            if ss & set(t):
                continue
            merged = t + s
            order = sorted(range(len(merged)), key=lambda i: merged[i])
            sign, lst = 1, list(order)
            for a in range(len(lst)):
                for b in range(a + 1, len(lst)):
                    if lst[a] > lst[b]:
                        sign = -sign
            M[idx[tuple(sorted(merged))], col] += sign * c
    return M


def generic_form(n, deg, rng):
    return {t: int(rng.integers(1, 9)) for t in itertools.combinations(range(n), deg)}


def analyse(n, q, seed, p):
    rng = np.random.default_rng(seed)
    alpha = generic_form(n, q, rng)
    beta = generic_form(n, q - 1, rng)
    A = {k: wedge_matrix(n, q, alpha, k) for k in range(n + 1) if k + q <= n}
    B = {k: wedge_matrix(n, q - 1, beta, k) for k in range(n + 1) if k + q - 1 <= n}
    rkA = {k: _rank(A[k], p) for k in A}
    dim = lambda k: comb(n, k) if 0 <= k <= n else 0
    h = {k: dim(k) - rkA.get(k, 0) - rkA.get(k - q, 0) for k in range(n + 1)}
    h = {k: v for k, v in h.items() if v > 0}
    out = []
    for j in sorted(h):
        tgt = j + q - 1
        if tgt not in h:
            continue
        Nj = _nullspace(A[j], p) if j in A else np.eye(dim(j), dtype=np.int64)
        img = (B[j] @ Nj) % p if j in B else np.zeros((dim(tgt), Nj.shape[1]), dtype=np.int64)
        prev = A.get(tgt - q)
        if prev is None or prev.size == 0:
            r = _rank(img, p)
            base = 0
        else:
            base = _rank(prev, p)
            r = _rank(np.concatenate([img, prev % p], axis=1), p) - base
        out.append((j, tgt, h[j], h[tgt], r))
    return h, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, default=3)
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--nmax", type=int, default=13)
    ap.add_argument("--seed", type=int, default=5)
    a = ap.parse_args()
    print("delta = beta ^ - on H(alpha), built explicitly.  q = %d\n" % a.q)
    print("%3s %4s %-22s | %s" % ("n", "w_s", "H(alpha)", "delta maps: H^j -> H^{j+q-1}  rank/min(dim)  corank"))
    for n in range(a.nmin, a.nmax + 1):
        if n < a.q:
            continue
        h, maps = analyse(n, a.q, a.seed, PRIME)
        h2, maps2 = analyse(n, a.q, a.seed, ALT_PRIME)
        assert [m[4] for m in maps] == [m[4] for m in maps2], "rank disagreed between primes at n=%d" % n
        if not maps:
            cells = "none (no two occupied degrees are q-1 apart)"
        else:
            cells = "   ".join("H^%d->H^%d %d/%d corank %d" % (j, t, r, min(dj, dt), min(dj, dt) - r)
                               for j, t, dj, dt, r in maps)
        print("%3d %4d %-22s | %s" % (n, len(h), str(sorted(h.items())), cells))
        sys.stdout.flush()


if __name__ == "__main__":
    main()
