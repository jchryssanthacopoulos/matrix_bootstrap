"""Test the general prediction   supp Z_lambda = k_0 + supp( z_slot^N )   irrep by irrep.

z_slot is the cohomology of a generic q-form on Lambda^bullet C^p (the rank-one / SYK problem), computed here
rather than assumed.  The prediction is that at EVERY weight the BPS window is the same set of degrees,
namely k_0 = p*C(N,2) shifted by the support of z_slot^N.  We also record whether the stronger factorisation
Z_lambda = m_lambda t^{k_0} z_slot(t)^N holds, which it need not once p >= q.
"""

import argparse, itertools, os, sys
from math import prod
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cohomology import complex_cohomology, kostant, weight_basis, exact_rank  # noqa: E402
from multitrace import apply_omega  # noqa: E402


def z_slot(p, q, seed=3):
    rng = np.random.default_rng(seed)
    om = {t: int(rng.integers(1, 9)) for t in itertools.combinations(range(p), q)}
    if not om:
        return {j: 1 for j in range(p + 1)}          # Lambda^q C^p = 0: whole exterior algebra survives
    st = {k: list(itertools.combinations(range(p), k)) for k in range(p + 1)}
    rk = {}
    for k in range(p + 1):
        if k + q > p:
            rk[k] = 0; continue
        idx = {s: r for r, s in enumerate(st[k + q])}
        M = np.zeros((len(st[k + q]), len(st[k])), dtype=np.int64)
        for col, s in enumerate(st[k]):
            for s2, v in apply_omega(om, s).items():
                M[idx[s2], col] += int(v)
        rk[k] = exact_rank(M)
    h = {k: len(st[k]) - rk[k] - rk.get(k - q, 0) for k in range(p + 1)}
    return {k: v for k, v in h.items() if v}


def polypow(poly, n):
    out = {0: 1}
    for _ in range(n):
        new = {}
        for a, x in out.items():
            for b, y in poly.items():
                new[a + b] = new.get(a + b, 0) + x * y
        out = new
    return out


def dimU(N, l):
    return round(prod((l[i] - l[j] + j - i) / (j - i) for i in range(N) for j in range(i + 1, N)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--p", type=int, required=True)
    ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args()
    N, p, q = a.N, a.p, 3
    k0 = p * N * (N - 1) // 2
    zs = z_slot(p, q)
    zN = polypow(zs, N)
    pred_supp = tuple(sorted(k0 + d for d in zN))
    print("(p,N,q)=(%d,%d,%d)   z_slot = %s   w_s = %d" % (p, N, q, sorted(zs.items()), len(zs)))
    print("  k_0 = %d,  predicted window = %d..%d  (W = %d)" % (k0, pred_supp[0], pred_supp[-1], len(pred_supp)))
    print("  predicted profile shape (m=1): %s" % [zN[d] for d in sorted(zN)])

    rng = np.random.default_rng(a.seed)
    C = rng.integers(-3, 4, size=(p,) * 3).astype(np.int64)
    weights = []
    for lam in itertools.product(range(-p * (N - 1), p * (N - 1) + 1), repeat=N):
        if sum(lam) or any(lam[i] < lam[i + 1] for i in range(N - 1)):
            continue
        if any(weight_basis(N, p, k, lam) for k in range(p * N * N + 1)):
            weights.append(lam)
    weights.sort(key=lambda l: -sum(x * x for x in l))
    H = {}
    for lam in weights:
        out = complex_cohomology(N, p, C, lam, cross_check=True)
        H[lam] = {k: v[2] for k, v in out.items() if v[2]}
    irr = {}
    for lam in weights:
        hh = dict(H[lam])
        for mu in weights:
            if mu == lam or mu not in irr:
                continue
            km = kostant(N, mu, lam)
            if km:
                for k, v in irr[mu].items():
                    hh[k] = hh.get(k, 0) - km * v
        irr[lam] = {k: v for k, v in hh.items() if v}

    supp_ok = fact_ok = True
    print("\n  %-16s %-14s %-6s %s" % ("irrep", "window", "factors?", "profile"))
    for lam in weights:
        hh = {k: v for k, v in irr[lam].items() if v}
        if not hh:
            continue
        sup = tuple(sorted(hh))
        m = hh[min(hh)] // zN[min(zN)] if hh[min(hh)] % zN[min(zN)] == 0 else None
        f = m is not None and hh == {k0 + d: m * zN[d] for d in zN}
        if sup != pred_supp:
            supp_ok = False
        if not f:
            fact_ok = False
        print("  %-16s %-14s %-6s %s" % (str(lam), "%d..%d" % (sup[0], sup[-1]),
                                         ("m=%d" % m) if f else "NO", [hh[k] for k in sorted(hh)]))
    print("\n  support == predicted in every irrep : %s" % supp_ok)
    print("  full factorisation m t^k0 z^N       : %s" % fact_ok)


if __name__ == "__main__":
    main()
