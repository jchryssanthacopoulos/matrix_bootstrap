"""BPS cohomology of the three-node quiver  Q = Tr(A B C)  with bifundamentals of U(n)^3.

A_{ij} in (n, nbar, 1), B_{jk} in (1, n, nbar), C_{ki} in (nbar, 1, n); 3n^2 complex fermions; Tr(ABC) is the
unique cubic invariant.  Q = omega ^ - with omega = sum_{ijk} a_ij b_jk c_ki.  Every mode has non-zero weight
(e^{(A)}_i - e^{(B)}_j with nodes A != B), so V_0 = 0: none of our zero-weight obstructions applies.

Cohomology is computed per (degree, weight) block with exact ranks over two primes; U(n)^3 irrep multiplicities
per degree follow from the Weyl-denominator formula  m_lambda = sum_{w in W} eps(w) d(lambda + rho - w rho).

    python scripts/run_quiver_cohomology.py --n 2
"""
import argparse, itertools, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cohomology import rank_mod_p, PRIME  # noqa: E402
from multitrace import apply_omega, _add_term  # noqa: E402


def setup(n):
    m = lambda e, x, y: e * n * n + x * n + y
    om = {}
    for i, j, k in itertools.product(range(n), repeat=3):
        _add_term(om, (m(0, i, j), m(1, j, k), m(2, k, i)), 1)
    wt = []
    for e in range(3):
        for x in range(n):
            for y in range(n):
                w = [0] * (3 * n)
                w[e * n + x] += 1                    # node e, fundamental index x
                w[((e + 1) % 3) * n + y] -= 1        # node e+1, antifundamental index y
                wt.append(tuple(w))
    return om, wt


def states_by_block(n, wt):
    nm = 3 * n * n
    blocks = {}
    for mask in range(1 << nm):
        s = tuple(i for i in range(nm) if mask >> i & 1)
        w = [0] * (3 * n)
        for i in s:
            for a, v in enumerate(wt[i]):
                w[a] += v
        blocks.setdefault((len(s), tuple(w)), []).append(s)
    return blocks


def cohomology(om, blocks):
    rk = {}
    for (k, w), b in blocks.items():
        tgt = blocks.get((k + 3, w))
        if not tgt:
            rk[(k, w)] = 0; continue
        idx = {s: r for r, s in enumerate(tgt)}
        ent = {}
        for c, s in enumerate(b):
            for s2, v in apply_omega(om, s).items():
                ent[(idx[s2], c)] = ent.get((idx[s2], c), 0) + v
        rs = []
        for P in (PRIME, 2147483629):
            M = np.zeros((len(tgt), len(b)), dtype=np.int64)
            for (r_, c_), v in ent.items():
                M[r_, c_] = v % P
            rs.append(rank_mod_p(M, P))
        assert rs[0] == rs[1]
        rk[(k, w)] = rs[0]
    return {key: len(b) - rk[key] - rk.get((key[0] - 3, key[1]), 0) for key, b in blocks.items()}


def irreps(n, H):
    rho = tuple(n - 1 - i for i in range(n))
    perms = list(itertools.permutations(range(n)))
    def sgn(p):
        s = 1
        for a in range(n):
            for b in range(a + 1, n):
                if p[a] > p[b]: s = -s
        return s
    out = {}
    keys = {(k, w) for (k, w), v in H.items() if v}
    cand = {(k, w) for (k, w) in keys
            if all(all(w[e*n+a] >= w[e*n+a+1] for a in range(n - 1)) for e in range(3))}
    for (k, lam) in cand:
        tot = 0
        for p3 in itertools.product(perms, repeat=3):
            eps = 1; mu = list(lam)
            for e, p in enumerate(p3):
                eps *= sgn(p)
                for a in range(n):
                    mu[e*n+a] += rho[a] - rho[p[a]]
            tot += eps * H.get((k, tuple(mu)), 0)
        if tot:
            out[(k, lam)] = tot
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--n", type=int, default=2); a = ap.parse_args()
    n = a.n
    om, wt = setup(n)
    blocks = states_by_block(n, wt)
    H = cohomology(om, blocks)
    irr = irreps(n, H)
    byl = {}
    for (k, lam), m in irr.items():
        byl.setdefault(lam, {})[k] = m
    tot = sum(v for v in H.values())
    print("n=%d: %d modes, total BPS states %d, BPS degrees %s" % (n, 3*n*n, tot, sorted({k for (k, w), v in H.items() if v})))
    print("%-28s %-12s %s" % ("irrep (nodes A|B|C)", "window", "multiplicities by degree"))
    bad = 0
    for lam in sorted(byl, key=lambda l: -sum(x*x for x in l)):
        hk = byl[lam]; ks = sorted(hk)
        classes = {}
        for k in ks: classes.setdefault(k % 3, []).append(k)
        conc = all(len(v) == 1 for v in classes.values())
        if not conc: bad += 1
        lab = "|".join(",".join(str(x) for x in lam[e*n:(e+1)*n]) for e in range(3))
        print("%-28s %-12s %s   %s" % (lab, "%d..%d" % (ks[0], ks[-1]), [hk[k] for k in ks], "" if conc else "<- two degrees in one Z3 class"))
    print("irreps with BPS states: %d, of which NOT concentrated: %d" % (len(byl), bad))


if __name__ == "__main__":
    main()
