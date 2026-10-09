#!/usr/bin/env python3
"""
Checks for the large-p Schwinger-Dyson analysis (src/n2syk_sd.py; docs/derivations.md D22).
Run as a script: python3 tests/test_n2syk.py

  1. the three Wick pairings of <Q Q^+ Q Q^+>_infinity for the tripartite (n = 1) supercharge
     Q = sum C_I T_I, T_I = a_f^+ b_g^+ c_h^+, equal (p^2+p)^3/64, (p^2+p)^3/64 and p^3/8 (explicit Fock operators,
     p = 2, 3), so that E_C Var(H) = (3 sigma^4/16) p^3 (p^2 + p + 1);
  2. Q^2 = 0 and <H>_infinity = sum_I C_I^2 / 4 for one coupling draw;
  3. SD at high temperature: E -> J/24 and dE/dbeta -> -J^2/32 per Majorana (qh = 3), i.e. J/12 and J^2/16 per
     complex fermion, the p -> infinity limits of item 1;
  4. SD at low temperature approaches the conformal solution (FGMS (2.29)-(2.30)) at tau = beta/2.
"""
import itertools, os, sys
import numpy as np
import scipy.sparse as sp

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
from n2syk_sd import solve, energy, conformal_G                                  # noqa: E402


def check(name, ok, detail=''):
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}")
    if not ok:
        sys.exit(1)


def fock_creators(nmodes):
    """Jordan-Wigner creation operators c_i^+ on 2^nmodes states (sparse)."""
    I, Z = sp.identity(2, format='csr'), sp.csr_matrix(np.diag([1.0, -1.0]))
    up = sp.csr_matrix(np.array([[0.0, 0.0], [1.0, 0.0]]))       # |0> -> |1>
    ops = []
    for i in range(nmodes):
        M = None
        for j in range(nmodes):
            f = Z if j < i else (up if j == i else I)
            M = f if M is None else sp.kron(M, f, format='csr')
        ops.append(M)
    return ops


def tripartite_T(p):
    cr = fock_creators(3 * p)
    a, b, c = cr[:p], cr[p:2 * p], cr[2 * p:]
    return {(f, g, h): (a[f] @ b[g] @ c[h]).tocsr() for f, g, h in itertools.product(range(p), repeat=3)}


def test_pairings():
    for p in (2, 3):
        T = tripartite_T(p)
        D = 2 ** (3 * p)
        tr = lambda M: M.diagonal().sum() / D
        P = {I: (t @ t.T).tocsr() for I, t in T.items()}            # T_I T_I^+ (real operators: ^+ = ^T)
        Pb = {I: (t.T @ t).tocsr() for I, t in T.items()}           # T_I^+ T_I
        s1 = sum(tr(P[I] @ P[K]) for I in T for K in T)
        s2 = sum(tr(T[I] @ Pb[K] @ T[I].T) for I in T for K in T)
        s3 = sum(tr(T[I] @ T[K].T @ T[I] @ T[K].T) for I in T for K in T)
        e12, e3 = (p * p + p) ** 3 / 64, p ** 3 / 8
        check(f'Wick pairings p={p}', abs(s1 - e12) < 1e-12 and abs(s2 - e12) < 1e-12 and abs(s3 - e3) < 1e-12,
              f'({s1:.6f}, {s2:.6f}, {s3:.6f}) vs ({e12:.6f}, {e12:.6f}, {e3:.6f})')


def test_draw():
    p = 2
    T = tripartite_T(p)
    C = np.random.default_rng(7).standard_normal((p, p, p)) / p
    Q = sum(C[I] * t for I, t in T.items())
    H = (Q @ Q.T + Q.T @ Q).toarray()
    D = H.shape[0]
    check('Q^2 = 0', abs(Q @ Q).max() < 1e-14 if (Q @ Q).nnz else True)
    check('<H>_inf = sum C^2 / 4', abs(np.trace(H) / D - (C ** 2).sum() / 4) < 1e-14)


def test_sd_high_T():
    E = {}
    for bJ in (0.005, 0.01):
        E[bJ] = energy(solve(bJ, qh=3, M=2 ** 10, tol=1e-14))
    slope = (E[0.01] - E[0.005]) / 0.005
    e0 = E[0.005] - 0.005 * slope
    check('SD: E(beta -> 0) = J/24', abs(e0 - 1 / 24) < 1e-6, f'{e0:.8f} vs {1/24:.8f}')
    check('SD: dE/dbeta(0) = -J^2/32', abs(slope + 1 / 32) < 2e-4, f'{slope:.6f} vs {-1/32:.6f}')


def test_sd_conformal():
    bJ = 100.0
    r = solve(bJ, qh=3, M=2 ** 14, tol=1e-12)
    k = r['M'] // 2
    tau = r['tau'][k - 1:k + 1].mean()                              # tau = beta/2 (between the two midpoints)
    ratio = r['G'][k - 1:k + 1].mean() / conformal_G(tau, r['beta'])
    check('SD -> conformal at beta J = 100, tau = beta/2', abs(ratio - 1) < 0.01, f'G/G_conf = {ratio:.4f}')


if __name__ == '__main__':
    test_pairings()
    test_draw()
    test_sd_high_T()
    test_sd_conformal()
    print('all checks passed')
