"""BPS cohomology of the Biggs-Lin-Maldacena SU(2) model, Q ~ sum C^j_{m1 m2 m3} psi_m1 psi_m2 psi_m3.

N = 2j+1 complex fermions psi_m (m = -j..j) in the spin-j irrep, j odd so that the 3j symbol (j j j) is
antisymmetric.  Q = omega ^ - with omega a 3-form, so H = ker Q / im Q is computed per block of fixed
fermion number k and J3 = M, and SU(2) multiplicities follow from h_l(k) = h(k, M=l) - h(k, M=l+1).

Exactness.  The 3j symbols are square roots of rationals.  By the Racah formula,
    (j j j; m1 m2 m3) = sqrt(Delta) * prod_i sqrt((j+m_i)!(j-m_i)!) * S(m),   S rational,
so the diagonal rescaling psi_m -> psi_m / sqrt((j+m)!(j-m)!) together with dropping the common sqrt(Delta)
makes every coefficient rational.  A diagonal GL(V) change of basis preserves J3 and all cohomology
dimensions, so exact integer ranks over F_P apply.

BLM conventions: N_psi = k - (2j+1)/2, R = N_psi/3, so R = -1/6, +1/6 <-> k = j, j+1;
R = -1/2, +1/2 <-> k = j-1, j+2; R = -5/6, +5/6 <-> k = j-2, j+3.

    python scripts/run_su2_3j_model.py --jmax 9
"""
import argparse, itertools, os, sys
from fractions import Fraction
from math import factorial, lcm
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cohomology import rank_mod_p, PRIME   # noqa: E402
from multitrace import apply_omega, _add_term   # noqa: E402


def racah_S(j, m1, m2, m3):
    """Rational part S(m) of the 3j symbol (j j j; m1 m2 m3), sign included."""
    if m1 + m2 + m3 != 0:
        return Fraction(0)
    s = Fraction(0)
    for k in range(0, 3 * j + 1):
        args = [k, j - j + k + m1, j - j + k - m2, j + j - j - k, j - k - m1, j - k + m2]
        if min(args) < 0:
            continue
        den = 1
        for a in args:
            den *= factorial(a)
        s += Fraction((-1) ** k, den)
    return (1 if (j - j - m3) % 2 == 0 else -1) * s


def omega(j):
    N = 2 * j + 1
    raw = {}
    for m1, m2, m3 in itertools.product(range(-j, j + 1), repeat=3):
        c = racah_S(j, m1, m2, m3)
        if c:
            raw[(m1 + j, m2 + j, m3 + j)] = c
    L = 1
    for c in raw.values():
        L = lcm(L, c.denominator)
    om = {}
    for ms, c in raw.items():
        _add_term(om, ms, int(c * L))
    from math import gcd
    g = 0
    for v in om.values():
        g = gcd(g, abs(v))
    return {k: v // g for k, v in om.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jmin", type=int, default=3)
    ap.add_argument("--jmax", type=int, default=9)
    a = ap.parse_args()
    for j in range(a.jmin, a.jmax + 1, 2):
        N = 2 * j + 1
        om = omega(j)
        mval = [i - j for i in range(N)]
        blocks = {}
        for k in range(N + 1):
            for s in itertools.combinations(range(N), k):
                M = sum(mval[i] for i in s)
                if M >= 0:
                    blocks.setdefault((k, M), []).append(s)
        rk = {}
        for (k, M), b in blocks.items():
            tgt = blocks.get((k + 3, M))
            if not tgt:
                rk[(k, M)] = 0
                continue
            idx = {s: r for r, s in enumerate(tgt)}
            entries = {}
            for col, s in enumerate(b):
                for s2, v in apply_omega(om, s).items():
                    entries[(idx[s2], col)] = entries.get((idx[s2], col), 0) + v   # exact Python ints
            ranks = []
            for P in (PRIME, 2147483629):
                Mt = np.zeros((len(tgt), len(b)), dtype=np.int64)
                for (r_, c_), v in entries.items():
                    Mt[r_, c_] = v % P
                ranks.append(rank_mod_p(Mt, P))
            assert ranks[0] == ranks[1], "prime disagreement"
            rk[(k, M)] = ranks[0]
        h = {key: len(b) - rk[key] - rk.get((key[0] - 3, key[1]), 0) for key, b in blocks.items()}
        # SU(2) multiplets: h_l(k) = h(k,l) - h(k,l+1)
        mult = {}
        for (k, M), v in h.items():
            d = v - h.get((k, M + 1), 0)
            if d:
                mult.setdefault(k, {})[M] = d
        print("j=%d  (N=%d fermions, maximal block %d)" % (j, N, max(len(b) for b in blocks.values())))
        for k in sorted(mult):
            R = Fraction(2 * k - N, 6)
            states = sum(d * (2 * l + 1) for l, d in mult[k].items())
            ls = sorted(mult[k].items())
            short = ", ".join("l=%d:%d" % (l, d) for l, d in ls) if len(ls) <= 6 else "%d spins, l=%d..%d" % (len(ls), ls[0][0], ls[-1][0])
            print("   k=%-3d R=%-6s BPS states=%-7d multiplets=%-5d  %s" % (k, str(R), states, sum(mult[k].values()), short))
        print("   3^j = %d" % 3 ** j)
        sys.stdout.flush()


if __name__ == "__main__":
    main()
