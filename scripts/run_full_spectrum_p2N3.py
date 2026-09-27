"""Complete BPS spectrum of the (p,N)=(2,3) model, resolved into U(3) irreps.

The point is to test, on a spectrum that is complete rather than truncated to high Casimir, whether the BPS
window is the same in every irrep: same width, and indeed the same absolute degrees.  The rank law predicts
W = N(w_s-1)+1 = 7 at the maximal weight, since at p=2 a self-loop needs three distinct flavours at one site
and none exists, so z_slot = (1+t)^2 and w_s = 3.
"""

import argparse
import itertools
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cohomology import complex_cohomology, peel, weight_basis  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=3)
    ap.add_argument("--p", type=int, default=2)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--out", default="results/data/cohomology_N3_p2_full.json")
    a = ap.parse_args()
    N, p = a.N, a.p

    rng = np.random.default_rng(a.seed)
    C = rng.integers(-3, 4, size=(p,) * 3).astype(np.int64)

    weights = []
    for lam in itertools.product(range(-p * (N - 1), p * (N - 1) + 1), repeat=N):
        if sum(lam) or any(lam[i] < lam[i + 1] for i in range(N - 1)):
            continue                                  # dominant weights only
        if any(weight_basis(N, p, k, lam) for k in range(p * N * N + 1)):
            weights.append(lam)
    weights.sort(key=lambda l: -sum(x * x for x in l))
    print("dominant weights: %d" % len(weights), flush=True)

    H = {}
    for lam in weights:
        out = complex_cohomology(N, p, C, lam, cross_check=True)
        H[lam] = {k: v[2] for k, v in out.items() if v[2]}
        print("  lam=%-14s dimW=%-7d window=%s" % (str(lam), sum(v[0] for v in out.values()),
                                                   sorted(H[lam])), flush=True)

    irreps = peel(N, H, order=weights)
    print("\n%-14s %8s %8s %6s %-18s %s" % ("irrep", "C2", "total", "W", "window", "profile"))
    c2 = lambda l: sum(x * x for x in l)
    rows = []
    for lam in weights:
        hh = {k: v for k, v in irreps[lam].items() if v}
        if not hh:
            print("%-14s %8d %8s %6s %-18s %s" % (str(lam), c2(lam), "0", "-", "-", "(no BPS states)"))
            continue
        w = sorted(hh)
        rows.append((lam, w))
        print("%-14s %8d %8d %6d %-18s %s" % (str(lam), c2(lam), sum(hh.values()), len(w),
                                              "%d..%d" % (w[0], w[-1]), [hh[k] for k in w]))
    widths = {len(w) for _, w in rows}
    spans = {(w[0], w[-1]) for _, w in rows}
    print("\ndistinct widths across irreps with BPS states: %s" % sorted(widths))
    print("distinct absolute windows:                      %s" % sorted(spans))
    json.dump({str(k): {str(a): b for a, b in v.items()} for k, v in irreps.items()},
              open(a.out, "w"), indent=1)
    print("written to %s" % a.out)


if __name__ == "__main__":
    main()
