"""Test  Z_lambda(t) = m_lambda * t^{k_lambda} * (1+t)^{pN}  irrep by irrep, for general odd q.

For p < q there is no self-loop, so the Cartan differential vanishes and z_slot = (1+t)^p.  The claim under
test is that the cohomology at every weight is a FREE Lambda^bullet(Z)-module of rank m_lambda generated in a
single degree k_lambda.  At q=3 we found k_lambda constant, = p*C(N,2); this script checks whether the
factorisation survives at larger q with k_lambda allowed to vary.

    python scripts/run_irrep_factorisation.py --p 1 --N 3 --q 5 --seed 7
"""

import argparse
import itertools
import os
import sys
from math import comb, prod

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cohomology import apply_Q_degree, exact_rank, kostant, weight_basis  # noqa: E402


def cohomology_q(N, p, C, lam, q):
    bases = {}
    for k in range(p * N * N + 1):
        b = weight_basis(N, p, k, lam)
        if b:
            bases[k] = b
    ranks = {}
    for k in sorted(bases):
        if k + q not in bases:
            ranks[k] = 0
            continue
        idx = {s: r for r, s in enumerate(bases[k + q])}
        M = np.zeros((len(bases[k + q]), len(bases[k])), dtype=np.int64)
        for col, st in enumerate(bases[k]):
            for st2, v in apply_Q_degree(N, p, C, st, q).items():
                M[idx[st2], col] += int(round(v))
        ranks[k] = exact_rank(M, cross_check=True)
        del M
    return {k: len(bases[k]) - ranks[k] - ranks.get(k - q, 0) for k in sorted(bases)}


def peel_general(N, H, order):
    out = {}
    for lam in order:
        hh = dict(H[lam])
        for mu in order:
            if mu == lam or mu not in out:
                continue
            kmu = kostant(N, mu, lam)
            if kmu:
                for k, v in out[mu].items():
                    hh[k] = hh.get(k, 0) - kmu * v
        out[lam] = {k: v for k, v in hh.items() if v}
    return out


def factorise(h, pN):
    """Is h(t) = m * t^k * (1+t)^pN ?  Return (m, k) or None."""
    if not h:
        return (0, None)
    k = min(h)
    m = h[k]
    target = {k + j: m * comb(pN, j) for j in range(pN + 1)}
    return (m, k) if h == target else None


def dimU(N, l):
    return round(prod((l[i] - l[j] + j - i) / (j - i) for i in range(N) for j in range(i + 1, N)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--p", type=int, required=True)
    ap.add_argument("--q", type=int, default=3)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    N, p, q = a.N, a.p, a.q
    assert q % 2 == 1 and q <= 2 * N - 1, "need odd q <= 2N-1"

    rng = np.random.default_rng(a.seed)
    C = rng.integers(-3, 4, size=(p,) * q).astype(np.int64) if p > 1 else np.ones((1,) * q)

    weights = []
    for lam in itertools.product(range(-p * (N - 1), p * (N - 1) + 1), repeat=N):
        if sum(lam) or any(lam[i] < lam[i + 1] for i in range(N - 1)):
            continue
        if any(weight_basis(N, p, k, lam) for k in range(p * N * N + 1)):
            weights.append(lam)
    weights.sort(key=lambda l: -sum(x * x for x in l))

    H = {lam: {k: v for k, v in cohomology_q(N, p, C, lam, q).items() if v} for lam in weights}
    irreps = peel_general(N, H, weights)

    print("(p,N,q) = (%d,%d,%d)  seed %d   pN = %d,  p*C(N,2) = %d" % (p, N, q, a.seed, p * N, p * N * (N - 1) // 2))
    ok, bad, ks, tot = True, [], {}, 0
    for lam in weights:
        hh = {k: v for k, v in irreps[lam].items() if v}
        if not hh:
            continue
        f = factorise(hh, p * N)
        if f is None:
            ok = False
            bad.append((lam, hh))
            print("  %-16s NOT of the form m t^k (1+t)^%d :  %s" % (str(lam), p * N, hh))
        else:
            m, k = f
            ks[lam] = (m, k)
            tot += m * dimU(N, lam)
            print("  %-16s m=%-3d k=%-3d W=%d  %s" % (str(lam), m, k, p * N + 1, [hh[x] for x in sorted(hh)]))
    print("  ---")
    print("  factorises in every irrep : %s" % ok)
    if ok:
        kk = sorted({k for _, k in ks.values()})
        print("  generator degrees k_lambda: %s  (%s)" % (kk, "constant" if len(kk) == 1 else "varies"))
        print("  widths                    : all %d = pN+1" % (p * N + 1))
        print("  sum m*dim                 : %d   q^{p C(N,2)} = %d   %s"
              % (tot, q ** (p * N * (N - 1) // 2), "MATCH" if tot == q ** (p * N * (N - 1) // 2) else "differs"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
