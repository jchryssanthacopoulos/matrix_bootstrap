"""Multi-trace supercharges for the (p, q, N) fermionic matrix model.

The single-trace model has Q = C_abc Tr[Psi^a Psi^b Psi^c].  Its restriction to the Cartan (which is the whole
story at the maximal weight, and controls the window bound W >= 2 alpha_0 - p r + 1 at every weight) is

    Tr[H_1 H_2 H_3] = sum_i (H_1)_ii (H_2)_ii (H_3)_ii ,

i.e. it LOCALISES to a single site.  That is what makes the Cartan interaction hypergraph sparse, gives
alpha_0 = min(p, q-1) N, and forces W >= N + 1.

Multi-trace terms are equally gauge invariant and do NOT localise:

    Tr[Psi^a] Tr[Psi^b Psi^c] |_Cartan = sum_{i,j} Psi^a_ii Psi^b_jj Psi^c_jj      (two sites)
    Tr[Psi^a] Tr[Psi^b] Tr[Psi^c] |_Cartan = sum_{i,j,k} Psi^a_ii Psi^b_jj Psi^c_kk (three sites)

so they add Cartan triangles that straddle sites and drive alpha_0 down.  This module builds Q as an explicit
3-form omega (a dict from a sorted mode triple to a coefficient), which makes any mixture of trace structures
uniform and, since Q = omega ^ - with omega of ODD degree, gives Q^2 = 0 for free.

Conventions (modes, weights, bases) are taken unchanged from cohomology.py.
"""

import itertools
import numpy as np

from cohomology import (PRIME, mode_index, weight_basis, rank_mod_p, _create)


# --------------------------------------------------------------------------------------------------- the form

def _add_term(omega, ms, coeff):
    """Accumulate coeff * e_{ms[0]} ^ e_{ms[1]} ^ e_{ms[2]} into omega, in sorted order with the sign."""
    if len(set(ms)) < len(ms):
        return                                   # repeated mode: the wedge vanishes
    order = sorted(range(len(ms)), key=lambda t: ms[t])
    sign, seen = 1, list(order)                  # parity of the sorting permutation
    for a in range(len(seen)):
        for b in range(a + 1, len(seen)):
            if seen[a] > seen[b]:
                sign = -sign
    key = tuple(sorted(ms))
    v = omega.get(key, 0) + sign * coeff
    if v:
        omega[key] = v
    else:
        omega.pop(key, None)


def omega_form(N, p, C1=None, C2=None, C3=None):
    """The 3-form of  C1 Tr[PPP]  +  C2 Tr[P] Tr[PP]  +  C3 Tr[P] Tr[P] Tr[P].

    Each C is a generic (p, p, p) tensor; the accumulation projects it onto whichever symmetry type the
    corresponding trace structure actually supports, so no symmetrisation is needed on input.
    """
    om = {}
    if C1 is not None:                                                   # single trace: one index loop i->j->k->i
        for a, b, c in itertools.product(range(p), repeat=3):
            if C1[a, b, c] == 0:
                continue
            for i, j, k in itertools.product(range(N), repeat=3):
                _add_term(om, (mode_index(N, a, i, j), mode_index(N, b, j, k), mode_index(N, c, k, i)), C1[a, b, c])
    if C2 is not None:                                                   # double trace: a self-loop and a 2-loop
        for a, b, c in itertools.product(range(p), repeat=3):
            if C2[a, b, c] == 0:
                continue
            for i, j, k in itertools.product(range(N), repeat=3):
                _add_term(om, (mode_index(N, a, i, i), mode_index(N, b, j, k), mode_index(N, c, k, j)), C2[a, b, c])
    if C3 is not None:                                                   # triple trace: three self-loops
        for a, b, c in itertools.product(range(p), repeat=3):
            if C3[a, b, c] == 0:
                continue
            for i, j, k in itertools.product(range(N), repeat=3):
                _add_term(om, (mode_index(N, a, i, i), mode_index(N, b, j, j), mode_index(N, c, k, k)), C3[a, b, c])
    return om


def random_couplings(p, seed=0, lo=-4, hi=5):
    rng = np.random.default_rng(seed)
    return [rng.integers(lo, hi, size=(p, p, p)).astype(np.int64) for _ in range(3)]


# ---------------------------------------------------------------------------------------------- the complexes

def apply_omega(omega, state):
    """(omega ^ -)|state> as {new_state: coefficient}."""
    out = {}
    for ms, coeff in omega.items():
        s, st, ok = 1, state, True
        for m in reversed(ms):                    # rightmost operator acts first, as in cohomology.apply_Q
            r = _create(st, m)
            if r is None:
                ok = False
                break
            sg, st = r
            s *= sg
        if ok:
            out[st] = out.get(st, 0) + s * coeff
    return {st: v for st, v in out.items() if v != 0}


def cohomology_at_weight(N, p, omega, lam, prime=PRIME, cross_check=True):
    """Window, Betti numbers and dimensions of the weight-lam complex under Q = omega ^ -."""
    bases = {}
    for k in range(p * N * N + 1):
        b = weight_basis(N, p, k, lam)
        if b:
            bases[k] = b
    ranks = {}
    for k in sorted(bases):
        if k + 3 not in bases:
            ranks[k] = 0
            continue
        idx = {s: r for r, s in enumerate(bases[k + 3])}
        M = np.zeros((len(bases[k + 3]), len(bases[k])), dtype=np.int64)
        for col, st in enumerate(bases[k]):
            for st2, v in apply_omega(omega, st).items():
                M[idx[st2], col] += int(v)
        r1 = rank_mod_p(M, prime)
        if cross_check:
            assert r1 == rank_mod_p(M, 2147483629), "rank disagreed between primes"
        ranks[k] = r1
    h = {k: len(bases[k]) - ranks[k] - ranks.get(k - 3, 0) for k in sorted(bases)}
    h = {k: v for k, v in h.items() if v}
    return sorted(h), h, {k: len(b) for k, b in bases.items()}


def maximal_weight(N, p):
    return tuple(p * (N + 1 - 2 * i) for i in range(1, N + 1))


def check_nilpotent(N, p, omega, lam, kmax=None):
    """Q^2 = 0 on the weight-lam complex.  Automatic for an odd form, but cheap insurance on the bookkeeping."""
    for k in range(p * N * N + 1):
        if kmax is not None and k > kmax:
            break
        for st in weight_basis(N, p, k, lam):
            acc = {}
            for st1, v1 in apply_omega(omega, st).items():
                for st2, v2 in apply_omega(omega, st1).items():
                    acc[st2] = acc.get(st2, 0) + v1 * v2
            if any(v != 0 for v in acc.values()):
                return False
    return True


# -------------------------------------------------------------------------------- the Cartan triangle hypergraph

def cartan_triangles(N, p, omega):
    """Triangles all of whose modes are diagonal, indexed by (flavour, site) in 0..p*N-1."""
    diag = {mode_index(N, a, i, i): a * N + i for a in range(p) for i in range(N)}
    out = []
    for ms in omega:
        if all(m in diag for m in ms):
            out.append(frozenset(diag[m] for m in ms))
    return sorted(set(out), key=sorted)


def alpha_triangle_free(n, triangles):
    """Largest triangle-free subset, by exhaustive search over subsets (n <= ~22)."""
    tri = [tuple(sorted(t)) for t in triangles]
    best, best_set = 0, ()
    for size in range(n, 0, -1):
        if size <= best:
            break
        for S in itertools.combinations(range(n), size):
            s = set(S)
            if all(not s.issuperset(t) for t in tri):
                return size, S
    return best, best_set
