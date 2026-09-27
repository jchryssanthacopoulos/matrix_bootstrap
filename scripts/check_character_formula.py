"""Check the conjectured closed form against a stored irrep-resolved spectrum.

    Z_BPS(t,x) = t^{p C(N,2)} (1+t)^{pN} prod_{i<j} ( x_i/x_j + 1 + x_j/x_i )^p        (q = 3)

Checks, irrep by irrep: (i) Z_lambda = m_lambda t^{k0} (1+t)^{pN}; (ii) m_lambda equals the multiplicity of
chi_lambda in the character product, computed by multiplying by the Weyl denominator.
"""

import argparse
import itertools
import json
import os
import sys
from math import comb, prod


def character_mults(N, p):
    W = {tuple([0] * N): 1}
    for i in range(N):
        for j in range(i + 1, N):
            for _ in range(p):
                new = {}
                for w, c in W.items():
                    for m in (-1, 0, 1):
                        v = list(w); v[i] += m; v[j] -= m
                        t = tuple(v); new[t] = new.get(t, 0) + c
                W = new
    delta = tuple(range(N - 1, -1, -1))
    P = {}
    for sg in itertools.permutations(range(N)):
        s, l = 1, list(sg)
        for a in range(N):
            for b in range(a + 1, N):
                if l[a] > l[b]:
                    s = -s
        shift = tuple(delta[sg[i]] for i in range(N))
        for w, c in W.items():
            k = tuple(w[i] + shift[i] for i in range(N))
            P[k] = P.get(k, 0) + s * c
    out = {}
    for k, c in P.items():
        if c and all(k[i] > k[i + 1] for i in range(N - 1)):
            out[tuple(k[i] - delta[i] for i in range(N))] = c
    return out


def dimU(N, l):
    return round(prod((l[i] - l[j] + j - i) / (j - i) for i in range(N) for j in range(i + 1, N)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", required=True)
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--p", type=int, required=True)
    a = ap.parse_args()
    N, p = a.N, a.p
    k0 = p * N * (N - 1) // 2
    d = json.load(open(a.file))
    meas = {}
    shape_ok = True
    for k, v in d.items():
        lam = eval(k)
        h = {int(x): y for x, y in v.items() if y}
        if not h:
            continue
        m = h[min(h)]
        want = {k0 + j: m * comb(p * N, j) for j in range(p * N + 1)}
        if h != want:
            shape_ok = False
            print("  shape FAILS at %s: %s" % (str(lam), h))
        meas[lam] = m
    pred = character_mults(N, p)
    mult_ok = pred == meas
    print("(p,N)=(%d,%d)  irreps with BPS: %d" % (p, N, len(meas)))
    print("  shape  Z_lam = m t^%d (1+t)^%d : %s" % (k0, p * N, shape_ok))
    print("  m_lam  == character multiplicity : %s" % mult_ok)
    if not mult_ok:
        for lam in sorted(set(pred) | set(meas), key=lambda l: -sum(x * x for x in l)):
            if pred.get(lam, 0) != meas.get(lam, 0):
                print("    %-16s pred %d  meas %d" % (str(lam), pred.get(lam, 0), meas.get(lam, 0)))
    print("  total BPS states = %d,  2^{pN} 3^{k0} = %d  %s"
          % (sum(m * dimU(N, l) for l, m in meas.items()) * 2 ** (p * N),
             2 ** (p * N) * 3 ** k0,
             "MATCH" if sum(m * dimU(N, l) for l, m in meas.items()) == 3 ** k0 else "differs"))
    return 0 if (shape_ok and mult_ok) else 1


if __name__ == "__main__":
    sys.exit(main())
