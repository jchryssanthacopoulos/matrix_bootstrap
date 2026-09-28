"""Per-irrep refined Witten index of the U(n)^3 quiver Q = Tr(ABC), from the Fock character.

The Fock character is prod_modes (1 + t x^{w}).  For each dominant lambda and degree k the Fock multiplicity is
m_lambda(k) = sum_{w in W} eps(w) d(k, lambda + rho - w rho).  Since Q raises k by 3, the Euler characteristic of
the strand c = k mod 3 in irrep lambda is I_c(lambda) = sum_{k = c mod 3} (-1)^{(k-c)/3} m_lambda(k), and
|I_c(lambda)| is an exact lower bound on the BPS multiplicity of that complex.

    python scripts/quiver_index.py --n 3
"""
import argparse, itertools
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--n", type=int, default=3); ap.add_argument("--p", type=int, default=1); a = ap.parse_args()
    n, P = a.n, a.p
    wts = []
    for e in range(3):
        for x in range(n):
            for y in range(n):
                w = [0] * (3 * n); w[e*n + x] += 1; w[((e+1) % 3)*n + y] -= 1
                wts.extend([tuple(w)] * P)                      # p flavours per edge
    Z = {(0, tuple([0] * (3*n))): 1}
    for w in wts:
        new = defaultdict(int)
        for (k, v), c in Z.items():
            new[(k, v)] += c
            new[(k+1, tuple(a_ + b_ for a_, b_ in zip(v, w)))] += c
        Z = new
    rho = tuple(n - 1 - i for i in range(n))
    perms = list(itertools.permutations(range(n)))
    def sgn(p):
        s = 1
        for i in range(n):
            for j in range(i+1, n):
                if p[i] > p[j]: s = -s
        return s
    dom = {v for (k, v) in Z if all(all(v[e*n+i] >= v[e*n+i+1] for i in range(n-1)) for e in range(3))}
    ks = range(3*P*n*n + 1)
    idx = {}; fock = {}
    for lam in dom:
        m = {}
        for k in ks:
            tot = 0
            for p3 in itertools.product(perms, repeat=3):
                eps = 1; mu = list(lam)
                for e, p in enumerate(p3):
                    eps *= sgn(p)
                    for i in range(n): mu[e*n+i] += rho[i] - rho[p[i]]
                tot += eps * Z.get((k, tuple(mu)), 0)
            if tot: m[k] = tot
        if not m: continue
        fock[lam] = m
        I = {c: sum((-1) ** ((k - c)//3) * v for k, v in m.items() if k % 3 == c) for c in range(3)}
        idx[lam] = I
    nz = {lam: I for lam, I in idx.items() if any(I.values())}
    mx = max(max(abs(x) for x in I.values()) for I in nz.values())
    maxfock = max(sum(m.values()) for m in fock.values())
    sing = tuple([0]*(3*n))
    print("n=%d p=%d (%d modes): %d irreps occur; %d have non-zero index" % (n, P, 3*P*n*n, len(fock), len(nz)))
    print("  max |I_c(lambda)| over all irreps and strands = %d" % mx)
    print("  max total Fock multiplicity of any irrep     = %d" % maxfock)
    print("  sum over irreps of sum_c |I_c| (index lower bound on BPS multiplets) = %d" % sum(sum(abs(x) for x in I.values()) for I in nz.values()))
    print("  singlet: Fock multiplicity %s, index %s" % (sum(fock.get(sing, {}).values()), idx.get(sing)))
    hist = defaultdict(int)
    for I in nz.values():
        for x in I.values():
            if x: hist[abs(x)] += 1
    print("  histogram of |I_c| over non-zero strands:", dict(sorted(hist.items())))


if __name__ == "__main__":
    main()
