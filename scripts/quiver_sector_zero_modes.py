#!/usr/bin/env python3
"""
Zero modes of H = {Q, Q^dag} in the full Fock sector N_A = N_B = N_C = m of the (n, p) quiver, ALL gauge irreps
included (no singlet projection).  Q is built by bit-vectorised fermion algebra on the three sectors m-1, m, m+1 and
the lowest eigenvalues of H on sector m are found by Lanczos.

Purpose: diagnose the bootstrap frontier (research/notes/quiver_project.md section 12).  The singlet bootstrap
excludes a degree only if its relaxation can rule out every functional obeying the rows; where non-singlet BPS
states exist in the same N_e sector, the singlet rows must do all the work.  At (2,2) the level-2 frontier sits
exactly at the onset of non-singlet BPS states (k = 9).

    python scripts/quiver_sector_zero_modes.py --n 2 --p 2 --m 3 --seed 3 --nev 12

Memory: the three sectors have binom(p n^2, m)^3 states each; Q_hi has about dim(m) * p^3 n^3 * (fill) nonzeros.
(2,2), m = 3: 0.4 GB.
"""
import argparse, itertools, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla


def sector_states(n, p, m):
    per = p * n * n
    subsets = [sum(1 << i for i in c) for c in itertools.combinations(range(per), m)]
    out = [s0 | (s1 << per) | (s2 << 2 * per) for s0 in subsets for s1 in subsets for s2 in subsets]
    return np.array(sorted(out), dtype=np.int64)


def parity_below(s, j):
    """(-1)^(number of occupied modes below j), vectorised over states s."""
    x = s & ((np.int64(1) << j) - 1)
    c = np.zeros_like(x)
    while np.any(x):
        c += x & 1
        x >>= 1
    return 1 - 2 * (c & 1)


def build_Q(n, p, C, src, dst):
    """Q = sum C_abc sum_ijk A^a_ij B^b_jk C^c_ki from sector src to sector dst (mode order: edge, flavour, row, col)."""
    per = p * n * n
    mode = lambda e, f, x, y: e * per + (f * n + x) * n + y
    rows, cols, vals = [], [], []
    for a, b, c in itertools.product(range(p), repeat=3):
        if C[a, b, c] == 0:
            continue
        for i, j, k in itertools.product(range(n), repeat=3):
            s = src.copy()
            sign = np.ones(len(s), dtype=np.int64)
            ok = np.ones(len(s), dtype=bool)
            for mo in (mode(2, c, k, i), mode(1, b, j, k), mode(0, a, i, j)):      # C acts first
                bit = np.int64(1) << mo
                ok &= (s & bit) == 0
                sign *= parity_below(s, mo)
                s = s | bit
            idx = np.searchsorted(dst, s[ok])
            assert np.all(dst[idx] == s[ok])
            rows.append(idx); cols.append(np.where(ok)[0]); vals.append(C[a, b, c] * sign[ok])
    return sp.csr_matrix((np.concatenate(vals).astype(float), (np.concatenate(rows), np.concatenate(cols))),
                         shape=(len(dst), len(src)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--p', type=int, required=True)
    ap.add_argument('--m', type=int, required=True, help='fermions per edge (degree k = 3m)')
    ap.add_argument('--seed', type=int, default=3, help='integer couplings 1..5 from default_rng(seed)')
    ap.add_argument('--nev', type=int, default=12)
    a = ap.parse_args()
    n, p, m = a.n, a.p, a.m
    C = np.random.default_rng(a.seed).integers(1, 6, size=(p, p, p)).astype(float)
    t0 = time.time()
    Vm1, V0, V1 = sector_states(n, p, m - 1), sector_states(n, p, m), sector_states(n, p, m + 1)
    Qlo, Qhi = build_Q(n, p, C, Vm1, V0), build_Q(n, p, C, V0, V1)
    print(f'(n,p,m)=({n},{p},{m}): dims {len(Vm1)}, {len(V0)}, {len(V1)}; Q nnz {Qlo.nnz}, {Qhi.nnz}; '
          f'max |Q^2| {abs(Qhi @ Qlo).max():.1e}  [{time.time() - t0:.0f}s]', flush=True)
    H = sla.LinearOperator((len(V0), len(V0)), matvec=lambda v: Qhi.T @ (Qhi @ v) + Qlo @ (Qlo.T @ v), dtype=float)
    w = sla.eigsh(H, k=a.nev, which='SA', tol=1e-9, ncv=max(60, 3 * a.nev), return_eigenvectors=False)
    print('lowest eigenvalues of H on the full N_e sector:', np.sort(w).round(6), f'[{time.time() - t0:.0f}s]')


if __name__ == '__main__':
    main()
