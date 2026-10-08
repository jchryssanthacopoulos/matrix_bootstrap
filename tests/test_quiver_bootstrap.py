#!/usr/bin/env python3
"""
Validity checks of the quiver singlet bootstrap (src/quiver_bootstrap.py): every equality row, every reality relation
and every cone must hold on functionals of exact states (multiplet-averaged), and the SDP bounds must respect the
exact values.  Run as a script: python3 tests/test_quiver_bootstrap.py

  1. (n,p) = (1,2), m = 1 and (1,3), m = 1: the averaged BPS functional satisfies all rows including the BPS rows;
     the BPS margin is ~0 (feasible) and the energy bound is <= 0.
  2. (n,p) = (2,1), m = 1, 3: the averaged singlet ground-state functional satisfies all rows including the Casimir,
     Gauss, sandwiched-Casimir (cone irrep) and finite-n rows; the energy lower bound is <= the exact singlet E_0.
"""
import itertools, os, sys
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
import quiver_trace as qt                                                         # noqa: E402
from quiver_bootstrap import QuiverSDP                                           # noqa: E402


def check(name, ok, detail=''):
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}", flush=True)
    if not ok:
        sys.exit(1)


# ------------------------------------------------------------------------------------------------ exact states
def explicit_Q(X, C):
    n, p = X.n, X.p
    Q = sp.csr_matrix((X.dim, X.dim))
    for a, b, c in itertools.product(range(p), repeat=3):
        if C[a, b, c] == 0:
            continue
        for i, j, k in itertools.product(range(n), repeat=3):
            Q = Q + C[a, b, c] * (X.letter_op(qt.P(0, a), i, j) @ X.letter_op(qt.P(1, b), j, k) @ X.letter_op(qt.P(2, c), k, i))
    return Q.tocsr()


def explicit_casimir(X, v):
    n, p = X.n, X.p
    e_in = (v - 1) % 3
    G = {}
    for i, j in itertools.product(range(n), repeat=2):
        g = sp.csr_matrix((X.dim, X.dim))
        for f in range(p):
            for k in range(n):
                g = g + X.letter_op(qt.P(v, f), i, k) @ X.letter_op(qt.B(v, f), k, j) \
                      + X.letter_op(qt.B(e_in, f), i, k) @ X.letter_op(qt.P(e_in, f), k, j)
        G[i, j] = g
    S = sum(G[i, j] @ G[j, i] for i, j in itertools.product(range(n), repeat=2))
    T = sum(G[i, i] for i in range(n))
    return (n * S - T @ T).tocsr()


def edge_numbers(X):
    per = X.p * X.n * X.n
    return np.array([[sum(1 for mo in s if mo // per == e) for e in range(3)] for s in X.states])


def sector_states(X, C, m, bps=False, tol=1e-8):
    """Orthonormal real basis of the singlet ground space (or the singlet BPS space) with N_e = m."""
    Ne = edge_numbers(X)
    sel = np.where((Ne == m).all(axis=1))[0]
    Q = explicit_Q(X, C)
    H = (Q @ Q.T + Q.T @ Q).tocsr()[sel][:, sel].toarray()
    K = sum(explicit_casimir(X, v) for v in range(3)).tocsr()[sel][:, sel].toarray()
    w, U = np.linalg.eigh(K)
    S = U[:, w < tol * max(1.0, abs(w).max())]                         # singlet subspace
    Hs = S.T @ H @ S
    e, V = np.linalg.eigh(Hs)
    pick = (e < tol * max(1.0, e.max())) if bps else (e < e[0] + 1e-8 * max(1.0, abs(e[0])))
    vecs = np.zeros((X.dim, int(pick.sum())))
    vecs[sel] = S @ V[:, pick]
    return vecs, e


def trace_apply(X, word, vec):
    """Tr[word] vec with an O(L n^3) index sweep (letters applied right to left)."""
    n = X.n
    if len(word) == 0:
        return n * vec
    cur = {(i0, j): X.letter_op(word[-1], j, i0) @ vec for i0 in range(n) for j in range(n)}
    for t in range(len(word) - 2, -1, -1):
        cur = {(i0, j): sum(X.letter_op(word[t], j, k) @ cur[(i0, k)] for k in range(n)) for i0 in range(n) for j in range(n)}
    return sum(cur[(i0, i0)] for i0 in range(n))


def functional(X, vecs, monos):
    phi = {}
    for mo in monos:
        val = 0.0
        for c in range(vecs.shape[1]):
            v = vecs[:, c]
            u = v
            for w in reversed(mo):
                u = trace_apply(X, w, u)
            val += float(v @ u)
        phi[mo] = val / vecs.shape[1]
    return phi


def validate(name, S, X, vecs, exact_E):
    phi = functional(X, vecs, S.monos)
    r = S.check_functional(phi)
    ok = r['worst_row'] < 1e-10 and r['worst_reality'] < 1e-10 and r['min_cone_eig'] > -1e-10 and abs(r['energy'] - exact_E) < 1e-8 * max(1, exact_E)
    check(name, ok, f"rows {r['worst_row']:.1e}, reality {r['worst_reality']:.1e}, min cone eig {r['min_cone_eig']:.1e}, "
                    f"E {r['energy']:.6f} (exact {exact_E:.6f}); {len(S.monos)} monomials, {len(S.rows)} rows")


def main():
    rng = np.random.default_rng(5)
    # ---- 1. n = 1: BPS functionals
    for p, m in ((2, 1), (3, 1)):
        C = rng.integers(1, 6, size=(p, p, p)).astype(float)
        X = qt.ExplicitQuiver(1, p)
        vecs, e = sector_states(X, C, m, bps=True)
        check(f'(1,{p}) m={m}: BPS singlets exist', vecs.shape[1] > 0, f'{vecs.shape[1]} states')
        S = QuiverSDP(C, 1, m, L_adj=2, L_sing=3, bps=True, gauss=(2, 2), cone_casimir=True, workers=4)
        validate(f'(1,{p}) m={m}: all rows incl. BPS on the BPS functional', S, X, vecs, 0.0)
        r = S.margin(eps=1e-9)
        check(f'(1,{p}) m={m}: BPS cell not excluded', r['margin'] > -1e-6, f"t* = {r['margin']:.2e} [{r['status']}]")
    # ---- 2. (2,1): singlet ground states, Casimir/Gauss/finite-n rows
    C = np.ones((1, 1, 1))
    X = qt.ExplicitQuiver(2, 1)
    for m in (1, 3):                    # (2,1) has no singlets at m = 2
        vecs, e = sector_states(X, C, m)
        S = QuiverSDP(C, 2, m, L_adj=2, L_sing=4, gauss=(2, 2), finite_n_len=6, cone_casimir=True, workers=4)
        validate(f'(2,1) m={m}: all rows on the singlet ground state', S, X, vecs, e[0])
        r = S.energy(eps=1e-9)
        check(f'(2,1) m={m}: energy bound <= exact singlet E_0', r['obj'] <= e[0] + 1e-6 * max(1, e[0]),
              f"bound {r['obj']:.6f}, exact {e[0]:.6f} [{r['status']}]")
    print('all quiver_bootstrap checks passed')


if __name__ == '__main__':
    main()
