"""Harness for src/cohomology.py: weight-space Q-cohomology and the Kostant inversion.
Run:  cd src && PYTHONPATH=. ../.venv/bin/python ../tests/test_cohomology.py  (or pytest)."""
import sys, os, itertools
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from fermion_matrix_model import chen_C
import cohomology as co
import bps_index as bi


def _C(p=3):
    return (3 * chen_C(p)).real.round().astype(np.int64)


def test_Q_squared_is_zero():
    """Q^2 = 0 on every basis state of a few weight spaces (the premise of the whole computation)."""
    for N, p, lam, ks in [(2, 3, (2, -2), [3, 4, 5]), (3, 3, (6, 0, -6), [9, 10, 11, 12])]:
        C = _C(p)
        for k in ks:
            for st in co.weight_basis(N, p, k, lam):
                out = {}
                for st1, v1 in co.apply_Q(N, p, C, st).items():
                    for st2, v2 in co.apply_Q(N, p, C, st1).items():
                        out[st2] = out.get(st2, 0) + v1 * v2
                assert all(v == 0 for v in out.values()), (N, lam, k, st)
    print("  Q^2 = 0 on all tested weight spaces: OK")


def test_weight_basis_counts():
    """Summing the weight spaces over all weights at fixed k must give binom(pN^2, k)."""
    from math import comb
    N, p = 2, 3
    for k in range(0, p * N * N + 1):
        tot = 0
        for a in range(-2 * p, 2 * p + 1):
            tot += len(co.weight_basis(N, p, k, (a, -a)))
        assert tot == comb(p * N * N, k), (k, tot, comb(p * N * N, k))
    print("  weight bases partition each N_Psi sector: OK")


def test_kostant_against_weyl_dimension():
    """sum_lambda K(mu, lambda) = dim mu for every irrep used in the N=3 peel."""
    for mu in [(1, 0, -1), (2, 0, -2), (3, 0, -3), (5, 0, -5), (6, 0, -6),
               (6, -1, -5), (5, 1, -6), (6, -2, -4), (4, 2, -6), (6, -3, -3), (3, 3, -6)]:
        tot = sum(co.kostant(3, mu, (a, b, -a - b)) for a in range(-9, 10) for b in range(-9, 10))
        assert tot == bi.irrep_dim(3, mu), (mu, tot, bi.irrep_dim(3, mu))
    assert co.kostant(3, (1, 0, -1), (0, 0, 0)) == 2          # the adjoint has a 2-dimensional zero-weight space
    assert co.kostant(3, (1, 0, -1), (1, 0, -1)) == 1         # highest weight is simple
    assert co.kostant(3, (1, 0, -1), (3, -1, -2)) == 0        # outside the weight diagram
    print("  Kostant multiplicities vs Weyl dimensions: OK")


def test_rank_is_prime_independent():
    """Ranks mod two different primes and in floating point must agree."""
    N, p, lam = 3, 3, (6, 0, -6)
    C = _C(p)
    bases = {k: b for k in range(p * N * N + 1) if (b := co.weight_basis(N, p, k, lam))}
    for k in sorted(bases):
        if k + 3 not in bases:
            continue
        M, _, _ = co.Q_matrix(N, p, C, k, lam, bases[k], bases[k + 3])
        r1 = co.rank_mod_p(M, (1 << 31) - 1)
        r2 = co.rank_mod_p(M, 2147483629)
        r3 = int(np.linalg.matrix_rank(M.astype(float)))
        assert r1 == r2 == r3, (k, r1, r2, r3)
    print("  ranks independent of the prime: OK")


def test_N2_cohomology_matches_index_and_ED():
    """N=2: cohomology sits at k=5,6,7 only, and the peeled counts reproduce the known multiplet numbers."""
    N, p = 2, 3
    C = _C(p)
    lams = [(3, -3), (2, -2), (1, -1), (0, 0)]
    Hw = {}
    for lam in lams:
        res = co.complex_cohomology(N, p, C, lam)
        Hw[lam] = {k: v[2] for k, v in res.items() if v[2]}
        assert set(Hw[lam]) <= {5, 6, 7}, (lam, Hw[lam])
    h = co.peel(N, Hw, order=lams)
    expected = {(3, -3): 9, (2, -2): 18, (1, -1): 27, (0, 0): 9}     # multiplets at k=5, from ED
    for lam, n in expected.items():
        assert h[lam].get(5, 0) == n, (lam, h[lam], n)
        assert h[lam].get(6, 0) == 2 * n and h[lam].get(7, 0) == n, (lam, h[lam])   # 1:2:1 within each irrep
    per_k = {k: sum(h[lam].get(k, 0) * bi.irrep_dim(N, lam) for lam in lams) for k in (5, 6, 7)}
    assert per_k == {5: 243, 6: 486, 7: 243}, per_k
    assert sum(per_k.values()) == 972
    print("  N=2 peeled cohomology matches ED; BPS states 243:486:243 = 1:2:1 over k=5,6,7: OK")


def test_euler_characteristic_equals_refined_index():
    """For every peeled irrep and class c, sum_{k=c mod 3} (-1)^k h^k = I_{c,lambda}."""
    N, p = 3, 3
    C = _C(p)
    lams = [(6, 0, -6), (6, -1, -5), (5, 1, -6)]
    Hw = {lam: {k: v[2] for k, v in co.complex_cohomology(N, p, C, lam).items() if v[2]} for lam in lams}
    h = co.peel(N, Hw, order=lams)
    ridx = bi.refined_index(bi.irrep_multiplicities(N, p))
    for lam in lams:
        assert all(v > 0 for v in h[lam].values()), (lam, h[lam])          # non-negativity
        assert set(h[lam]) <= {12, 13, 14, 15}, (lam, h[lam])
        assert h[lam].get(12) == h[lam].get(15) and h[lam].get(13) == h[lam].get(14)   # particle-hole
        for c in range(3):
            chi = sum((-1) ** k * v for k, v in h[lam].items() if k % 3 == c)
            I = sum(v for (cc, w, l), v in ridx.items() if cc == c and l == lam)
            assert chi == I, (lam, c, chi, I)
    print("  Euler characteristics equal the refined index: OK")


def test_blocked_rank_matches_reference():
    """The BLAS-blocked elimination must give the same rank as the plain one, for every block size and prime."""
    rng = np.random.default_rng(0)
    for (m, n, rk) in [(80, 60, 25), (200, 200, 137), (300, 150, 90), (150, 300, 120),
                       (513, 400, 301), (1000, 700, 512), (77, 77, 0)]:
        X = rng.integers(-5, 6, size=(m, rk)); Y = rng.integers(-5, 6, size=(rk, n))
        A = (X @ Y).astype(np.int64) if rk else np.zeros((m, n), dtype=np.int64)
        ref = co.rank_mod_p(A)
        assert ref == int(np.linalg.matrix_rank(A.astype(float))) == (rk if rk else 0), (m, n, rk, ref)
        for prime in co.SMALL_PRIMES:
            for block in (16, 64, 256):
                assert co.rank_mod_p_blocked(A, prime, block=block) == ref, (m, n, prime, block)
    print("  blocked rank == reference rank for all sizes/blocks/primes: OK")


def test_flavour_blocking_matches_dense():
    """Z_3 blocking must reproduce the dense cohomology, and the block ranks must sum to the dense rank."""
    C = _C(3)
    for N, lam, expected in [(2, (3, -3), {5: 9, 6: 18, 7: 9}),
                             (3, (6, 0, -6), {12: 27, 13: 81, 14: 81, 15: 27})]:
        pc, tot = co.complex_cohomology_blocked(N, 3, C, lam)
        assert {k: v[2] for k, v in tot.items() if v[2]} == expected, (N, lam, tot)
        dense = co.complex_cohomology(N, 3, C, lam)
        for k in dense:
            assert tot[k][1] == dense[k][1], (N, lam, k, tot[k][1], dense[k][1])   # ranks agree
            assert tot[k][0] == dense[k][0]
        # the three charges are equal for these weights, and sum to the total
        for k, v in expected.items():
            assert sum(pc[w][k] for w in range(3)) == v, (N, lam, k)
    # the primes used for blocking must admit a cube root of unity
    for prime in co.blocking_primes(3, 2):
        w = co.root_of_unity(prime, 3)
        assert w != 1 and pow(w, 3, prime) == 1 and prime % 3 == 1
    print("  Z_3 flavour blocking reproduces the dense cohomology and ranks: OK")


def test_N4_maximal_weight():
    """N=4 maximal weight: the enumerator must match the D6 generating function, and the cohomology must be
    81*binom(4,j) over k=22..26 with Euler characteristics equal to the refined index."""
    N, p, lam = 4, 3, (9, 3, -3, -9)
    C = _C(p)
    W, M = bi._weights_exact(N, p)
    e = tuple(int(x) % M for x in bi._exponents_of_weight(np.array(lam)))
    bases = {}
    for k in range(p * N * N + 1):
        b = co.weight_basis(N, p, k, lam)
        assert len(b) == int(W[(k, slice(None)) + e].sum()), (k, len(b))
        if b:
            bases[k] = b
    ranks = {}
    for k in sorted(bases):
        if k + 3 in bases:
            Mx, _, _ = co.Q_matrix(N, p, C, k, lam, bases[k], bases[k + 3])
            ranks[k] = co.rank_mod_p(Mx)
        else:
            ranks[k] = 0
    h = {k: len(bases[k]) - ranks[k] - ranks.get(k - 3, 0) for k in sorted(bases)}
    h = {k: v for k, v in h.items() if v}
    assert h == {22: 81, 23: 324, 24: 486, 25: 324, 26: 81}, h        # 81 * binom(4, j)
    ridx = bi.refined_index(bi.irrep_multiplicities(N, p))
    for c in range(3):
        chi = sum((-1) ** k * v for k, v in h.items() if k % 3 == c)
        I = sum(v for (cc, w, l), v in ridx.items() if cc == c and l == lam)
        assert chi == I, (c, chi, I)
    # saturation fails in c=1,2 even though their index is non-zero
    assert sorted(k for k in h if k % 3 == 0) == [24]
    assert sorted(k for k in h if k % 3 == 1) == [22, 25]
    assert sorted(k for k in h if k % 3 == 2) == [23, 26]
    print("  N=4 maximal weight: h = 81*binom(4,j) over k=22..26, Euler = index: OK")


if __name__ == '__main__':
    for fn in [test_Q_squared_is_zero, test_weight_basis_counts, test_kostant_against_weyl_dimension,
               test_rank_is_prime_independent, test_N2_cohomology_matches_index_and_ED,
               test_euler_characteristic_equals_refined_index, test_blocked_rank_matches_reference,
               test_flavour_blocking_matches_dense, test_N4_maximal_weight]:
        print(f"{fn.__name__} ...")
        fn()
    print("\nall cohomology tests passed")
