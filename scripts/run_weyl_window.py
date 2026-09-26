"""Is the Weyl group, rather than the trace structure, what keeps the BPS window open?

At the maximal weight the whole complex is the Cartan sector: p*N modes (flavour x site) carrying a 3-form.
Gauge invariance leaves that form invariant under the Weyl group, which for u(N)/su(N) is S_N permuting sites.
This script compares, at the same number of modes n = p*N:

  generic     an unconstrained 3-form            (= N=2 SYK; the concentrating benchmark)
  weyl        a generic S_N-invariant 3-form     (everything gauge invariance can possibly allow)
  single      the single-trace Cartan form       Tr[H H H] = sum_i H_ii H_ii H_ii
  multi       single + double + triple trace     the delocalised multi-trace model

If 'weyl' grows with N while 'generic' stays bounded, no choice of gauge-invariant couplings -- multi-trace
included -- can concentrate with adjoint matter, and the obstruction is the Weyl group itself.

  python scripts/run_weyl_window.py --p 1 2 3 --nmax 15 --seed 11
"""

import argparse
import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cohomology import exact_rank                      # noqa: E402
from multitrace import apply_omega, _add_term          # noqa: E402


def _sort_sign(t):
    """Parity of the permutation sorting t, with the sorted tuple."""
    s, lst = 1, list(t)
    for a in range(len(lst)):
        for b in range(a + 1, len(lst)):
            if lst[a] > lst[b]:
                s = -s
    return s, tuple(sorted(t))


def weyl_form(p, N, seed):
    """A generic S_N-invariant 3-form on modes (a, i) -> a*N + i.

    A 3-subset touches at most three sites, so its S_N-orbit is labelled by the pattern obtained after
    relabelling those sites canonically -- only S_3 matters, not S_N.  An orbit whose stabiliser reverses
    orientation is forced to have coefficient zero.
    """
    rng = np.random.default_rng(seed)
    n = p * N
    coeff, dead = {}, set()
    canon_of = {}
    for t in itertools.combinations(range(n), 3):
        sites = sorted({m % N for m in t})
        seen = {}
        for perm in itertools.permutations(range(len(sites))):
            relab = {s: perm[k] for k, s in enumerate(sites)}
            img = tuple((m // N) * N + relab[m % N] for m in t)
            sg, srt = _sort_sign(img)
            seen.setdefault(srt, set()).add(sg)
        best = min(seen)
        if len(seen[best]) > 1:                  # orientation-reversing stabiliser: the orbit must vanish
            dead.add(best)
        canon_of[t] = (best, next(iter(seen[best])))
        if best not in coeff:
            coeff[best] = int(rng.integers(1, 9))
    om = {}
    for t, (best, sg) in canon_of.items():
        if best in dead:
            continue
        om[t] = sg * coeff[best]
    return om


def trace_cartan_form(p, N, seed, multi, _tries=20):
    """The Cartan restriction of the gauge-invariant traces, on the same p*N modes.

    Only the totally antisymmetric part of each coupling survives the wedge, and for small p a random integer
    tensor can have that part vanish by accident, so redraw until the form is non-zero.
    """
    rng = np.random.default_rng(seed)
    for _ in range(_tries):
        C1 = rng.integers(-4, 5, size=(p, p, p))
        C2 = rng.integers(-4, 5, size=(p, p, p))
        C3 = rng.integers(-4, 5, size=(p, p, p))
        om = _build_trace_form(p, N, C1, C2, C3, multi)
        if om:
            return om
    return om


def _build_trace_form(p, N, C1, C2, C3, multi):
    m = lambda a, i: a * N + i
    om = {}
    for a, b, c in itertools.product(range(p), repeat=3):
        for i in range(N):                                       # Tr[H H H] : one site
            _add_term(om, (m(a, i), m(b, i), m(c, i)), int(C1[a, b, c]))
        if not multi:
            continue
        for i, j in itertools.product(range(N), repeat=2):       # Tr[H] Tr[H H] : two sites
            _add_term(om, (m(a, i), m(b, j), m(c, j)), int(C2[a, b, c]))
        for i, j, k in itertools.product(range(N), repeat=3):    # Tr[H] Tr[H] Tr[H] : three sites
            _add_term(om, (m(a, i), m(b, j), m(c, k)), int(C3[a, b, c]))
    return om


def generic_form(n, seed):
    rng = np.random.default_rng(seed)
    return {t: int(rng.integers(1, 9)) for t in itertools.combinations(range(n), 3)}


def window(n, om):
    if not om:
        return None, None
    st = {k: list(itertools.combinations(range(n), k)) for k in range(n + 1)}
    rk = {}
    for k in range(n + 1):
        if k + 3 > n:
            rk[k] = 0
            continue
        idx = {s: r for r, s in enumerate(st[k + 3])}
        M = np.zeros((len(st[k + 3]), len(st[k])), dtype=np.int64)
        for col, s in enumerate(st[k]):
            for s2, v in apply_omega(om, s).items():
                M[idx[s2], col] += int(v)
        rk[k] = exact_rank(M)
    h = {k: len(st[k]) - rk[k] - rk.get(k - 3, 0) for k in range(n + 1)}
    h = {k: v for k, v in h.items() if v}
    return len(h), h


def check_invariance(p, N, om, seed=0, trials=3):
    """Verify the constructed form really is S_N-invariant."""
    rng = np.random.default_rng(seed)
    for _ in range(trials):
        perm = rng.permutation(N)
        img = {}
        for t, c in om.items():
            _add_term(img, tuple((m // N) * N + perm[m % N] for m in t), c)
        if img != om:
            return False
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p", type=int, nargs="+", default=[1, 2, 3])
    ap.add_argument("--nmax", type=int, default=15, help="cap on n = p*N")
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args()

    print("n = p*N modes in the Cartan sector; q = 3, so concentration needs W <= 3\n")
    hdr = "%2s %3s %4s | %-7s %-7s %-7s %-7s" % ("p", "N", "n", "generic", "weyl", "single", "multi")
    print(hdr); print("-" * len(hdr))
    for p in a.p:
        for N in range(2, a.nmax + 1):
            n = p * N
            if n < a.nmin or n > a.nmax:
                continue
            forms = {
                "generic": generic_form(n, a.seed),
                "weyl": weyl_form(p, N, a.seed),
                "single": trace_cartan_form(p, N, a.seed, multi=False),
                "multi": trace_cartan_form(p, N, a.seed, multi=True),
            }
            assert check_invariance(p, N, forms["weyl"], a.seed), "weyl_form is not S_N-invariant"
            assert check_invariance(p, N, forms["multi"], a.seed), "multi-trace form is not S_N-invariant"
            row, detail = [], {}
            for key in ("generic", "weyl", "single", "multi"):
                W, h = window(n, forms[key])
                row.append("W=%s" % (W if W is not None else "-"))
                detail[key] = h
            print("%2d %3d %4d | %-7s %-7s %-7s %-7s" % (p, N, n, *row))
            sys.stdout.flush()
        print()


if __name__ == "__main__":
    main()
