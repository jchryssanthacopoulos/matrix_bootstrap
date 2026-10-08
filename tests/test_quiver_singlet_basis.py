#!/usr/bin/env python3
"""
Checks for the Schur-Weyl singlet basis of the quiver (src/quiver_singlet_basis.py; docs/derivations.md D20).
Run as a script: python3 tests/test_quiver_singlet_basis.py

  1. edge spaces: sign condition rho(s_i) (x) rho(s_i) B = -B for every generator of H_r, orthonormal columns;
  2. singlet dimensions per degree equal the exact counts n(k) (transfer matrix) at (2,2), (2,3) and (3,2);
  3. Q^2 = 0;
  4. the (2,2) pair spectra equal the stored Fock-space spectra level by level (seed-3 couplings);
  5. 90 BPS singlets at (2,2), k = 12.
"""
import json, os, sys
import numpy as np
import scipy.sparse as sp

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
from quiver_singlet_basis import SingletComplex, edge_space, shape_obj, partitions, compositions   # noqa: E402
from quiver_index import singlet_series_recursive                                               # noqa: E402


def check(name, ok, detail=''):
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}")
    if not ok:
        sys.exit(1)


def test_edge_spaces():
    worst_sign, worst_orth, count = 0.0, 0.0, 0
    for m in range(2, 6):
        for lam in partitions(m, 3):
            for mu in partitions(m, 3):
                for r in compositions(m, 2, 9):
                    B = edge_space(lam, mu, r).toarray()
                    if B.shape[1] == 0:
                        continue
                    L, M = shape_obj(lam), shape_obj(mu)
                    worst_orth = max(worst_orth, np.abs(B.T @ B - np.eye(B.shape[1])).max())
                    o = 0
                    for rf in r:
                        for i in range(o, o + rf - 1):
                            G = sp.kron(L.s(i), M.s(i)).toarray()
                            worst_sign = max(worst_sign, np.abs(G @ B + B).max())
                        o += rf
                    count += 1
    check('edge spaces: sign condition and orthonormality', worst_sign < 1e-10 and worst_orth < 1e-10,
          f'{count} spaces, sign {worst_sign:.1e}, orth {worst_orth:.1e}')


def test_dims():
    for n, p, mmax in ((2, 2, 8), (2, 3, 12), (3, 2, 7)):
        ref, _ = singlet_series_recursive(n, p)
        sc = SingletComplex(n, p, np.ones((p, p, p)))
        dims = [sc.dim(m) for m in range(mmax + 1)]
        ok = all(d == ref[3 * m] for m, d in enumerate(dims))
        check(f'singlet dimensions ({n},{p}), m <= {mmax}', ok, f'{dims}')


def test_Q_and_spectra():
    C = np.random.default_rng(3).integers(1, 6, size=(2, 2, 2)).astype(float)
    sc = SingletComplex(2, 2, C)
    ref = json.load(open(os.path.join(ROOT, 'results', 'data', 'quiver_singlet_pairs_n2_p2.json')))['pairs']
    Qs = {m: sc.Q_dense(m) for m in range(8)}
    nil = max(np.abs(Qs[m + 1] @ Qs[m]).max() for m in range(7))
    check('Q^2 = 0 at (2,2)', nil < 1e-10, f'max {nil:.1e}')
    worst = 0.0
    for m in range(8):
        ev = np.linalg.eigvalsh(Qs[m].T @ Qs[m])
        lev = np.sort(ev[ev > 1e-8 * max(1.0, ev.max())])
        rl = np.sort(np.array(ref[f'{3 * m},{3 * m + 3}']['levels']))
        assert len(rl) == len(lev)
        worst = max(worst, np.abs(rl - lev).max() / max(1.0, rl.max()))
    check('(2,2) pair spectra equal the stored Fock-space spectra', worst < 1e-12, f'max relative deviation {worst:.1e}')
    H = Qs[4].T @ Qs[4] + Qs[3] @ Qs[3].T
    w = np.linalg.eigvalsh(H)
    nb = int((w < 1e-8 * w.max()).sum())
    check('90 BPS singlets at (2,2), k = 12', nb == 90, f'{nb}')


if __name__ == '__main__':
    test_edge_spaces()
    test_dims()
    test_Q_and_spectra()
    print('all quiver_singlet_basis checks passed')
