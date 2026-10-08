#!/usr/bin/env python3
"""
Checks of the quiver trace-word algebra (src/quiver_trace.py) against explicit Fock-space operators.
Run as a script: python3 tests/test_quiver_trace.py
"""
import itertools, os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
import trace_algebra as ta                                                       # noqa: E402
import quiver_trace as qt                                                        # noqa: E402


def check(name, ok, detail=''):
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}")
    if not ok:
        sys.exit(1)


def explicit_Q(X, C):
    n, p = X.n, X.p
    Q = None
    for a, b, c in itertools.product(range(p), repeat=3):
        if C[a, b, c] == 0:
            continue
        for i, j, k in itertools.product(range(n), repeat=3):
            t = X.letter_op(qt.P(0, a), i, j) @ X.letter_op(qt.P(1, b), j, k) @ X.letter_op(qt.P(2, c), k, i)
            Q = C[a, b, c] * t if Q is None else Q + C[a, b, c] * t
    return Q.tocsr()


def explicit_casimir(X, v):
    n, p = X.n, X.p
    G = {}
    e_in = (v - 1) % 3
    for i, j in itertools.product(range(n), repeat=2):
        g = None
        for f in range(p):
            for k in range(n):
                t1 = X.letter_op(qt.P(v, f), i, k) @ X.letter_op(qt.B(v, f), k, j)
                t2 = X.letter_op(qt.B(e_in, f), i, k) @ X.letter_op(qt.P(e_in, f), k, j)
                g = t1 + t2 if g is None else g + t1 + t2
        G[i, j] = g
    S = sum(G[i, j] @ G[j, i] for i, j in itertools.product(range(n), repeat=2))
    T = sum(G[i, i] for i in range(n))
    return n * S - T @ T


def compare(name, X, expr, M, vecs):
    worst = 0.0
    for v in vecs:
        a = X.expr_apply(expr, v)
        b = M @ v
        worst = max(worst, np.abs(a - b).max() / max(1.0, np.abs(b).max()))
    check(name, worst < 1e-10, f'max rel dev {worst:.1e}')


def main():
    rng = np.random.default_rng(0)
    # (n, p) = (2, 1): full Fock space
    X = qt.ExplicitQuiver(2, 1)
    C = np.ones((1, 1, 1))
    vecs = [rng.standard_normal(X.dim) for _ in range(3)]
    Qe = explicit_Q(X, C)
    compare('Q at (2,1)', X, qt.supercharge(C), Qe, vecs)
    H = (Qe @ Qe.T + Qe.T @ Qe).tocsr()
    compare('H = {Q,Q^+} at (2,1)', X, qt.hamiltonian(C), H, vecs)
    import scipy.sparse as sp
    compare('N_tot at (2,1)', X, qt.number_operator(1), sp.diags(X.degree().astype(float)), vecs)
    for v in range(3):
        compare(f'node-{v} Casimir at (2,1)', X, qt.casimir(1, v), explicit_casimir(X, v), vecs)
    Q0 = Qe @ np.eye(X.dim)[:, X.index[()]]
    cas = [np.abs(X.expr_apply(qt.casimir(1, v), Q0)).max() for v in range(3)]
    check('Casimirs vanish on the singlet Q|0>', max(cas) < 1e-12, f'{cas}')
    blocks = [(qt.P(0, 0), qt.B(0, 0))] * 3
    rel = qt.finite_n_relation(blocks, 2)
    worst = max(np.abs(X.expr_apply(rel, v)).max() for v in vecs) if rel else 0.0
    check('finite-n relation (3 blocks, node 0) vanishes at n=2', worst < 1e-10, f'{len(rel)} monomials, max {worst:.1e}')
    # (n, p) = (1, 2): flavour structure, random couplings
    X1 = qt.ExplicitQuiver(1, 2)
    C1 = rng.standard_normal((2, 2, 2))
    v1 = [rng.standard_normal(X1.dim) for _ in range(3)]
    Q1 = explicit_Q(X1, C1)
    compare('H at (1,2), random couplings', X1, qt.hamiltonian(C1), (Q1 @ Q1.T + Q1.T @ Q1).tocsr(), v1)
    # (n, p) = (2, 2) restricted to degree <= 6: H on degree-3 vectors
    X2 = qt.ExplicitQuiver(2, 2, kmax=6)
    C2 = rng.standard_normal((2, 2, 2))
    deg = X2.degree()
    v2 = []
    for _ in range(2):
        v = np.zeros(X2.dim); v[deg == 3] = rng.standard_normal((deg == 3).sum()); v2.append(v)
    Q2 = explicit_Q(X2, C2)
    compare('H at (2,2) on degree-3 vectors, random couplings', X2, qt.hamiltonian(C2), (Q2 @ Q2.T + Q2.T @ Q2).tocsr(), v2)
    print('all quiver_trace checks passed')


if __name__ == '__main__':
    main()
