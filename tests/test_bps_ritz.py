"""
Checks of src/bps_ritz.py at N=2 (12 modes, full Fock space 4096) against independent full-space sparse matrices.

    .venv/bin/python tests/test_bps_ritz.py      (or pytest)
"""
import os
import sys

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
from fermion_matrix_model import chen_C                                       # noqa: E402
from cohomology import mode_index                                             # noqa: E402
import free_sectors as fs                                                     # noqa: E402
import bps_ritz as br                                                         # noqa: E402

N, P = 2, 3
C = chen_C()
NM = P * N * N


def full_c(m):
    """c_m on the full 2^NM space, basis index = occupation mask, sign (-1)^(occupied below m)."""
    dim = 1 << NM
    masks = np.arange(dim, dtype=np.int64)
    sel = (masks >> m) & 1 == 1
    src = masks[sel]
    sg = 1.0 - 2.0 * (np.bitwise_count(src & ((1 << m) - 1)) & 1)
    return sp.csr_matrix((sg, (src ^ (1 << m), src)), shape=(dim, dim))


CM = [full_c(m) for m in range(NM)]


def full_letter(letter, x, y):
    kind, a = letter
    return CM[mode_index(N, a, x, y)].T.tocsr() if kind == 'P' else CM[mode_index(N, a, y, x)]


def full_word(word, i, j):
    """(L_1 ... L_n)_ij as a full-space matrix by explicit index sums."""
    dim = 1 << NM
    acc = sp.csr_matrix((dim, dim))
    n = len(word)
    for mid in np.ndindex(*([N] * (n - 1))):
        idx = (i,) + tuple(mid) + (j,)
        M = sp.identity(dim, format='csr')
        for t in range(n):
            M = M @ full_letter(word[t], idx[t], idx[t + 1])
        acc = acc + M
    return acc


def embed(block, v):
    out = np.zeros(1 << NM)
    out[block.masks] = v
    return out


def test_word_action_and_gauge():
    rng = np.random.default_rng(0)
    blk = br.Block(N, P, 5)
    v = rng.standard_normal(blk.dim)
    sv = blk.sparse(v)
    for word in [(('B', 0),), (('P', 1), ('B', 2)), (('B', 0), ('P', 2), ('B', 1)), (('B', 2), ('B', 2), ('P', 0))]:
        for i in range(N):
            for j in range(N):
                ref = full_word(word, i, j) @ embed(blk, v)
                got = np.zeros(1 << NM)
                m, x = br.apply_word_component(N, word, i, j, sv)
                np.add.at(got, m, x)
                assert np.abs(got - ref).max() < 1e-12, (word, i, j)
    # gauge generator E_kl = sum c^dag_(a,k,j) c_(a,l,j) - sum c^dag_(a,i,l) c_(a,i,k)
    for k_, l_ in [(0, 1), (1, 0)]:
        E = sp.csr_matrix(((1 << NM), (1 << NM)))
        for a in range(P):
            for j in range(N):
                E = E + CM[mode_index(N, a, k_, j)].T @ CM[mode_index(N, a, l_, j)]
            for i in range(N):
                E = E - CM[mode_index(N, a, i, l_)].T @ CM[mode_index(N, a, i, k_)]
        ref = E @ embed(blk, v)
        got = np.zeros(1 << NM)
        m, x = br.apply_gauge(N, P, k_, l_, sv)
        np.add.at(got, m, x)
        assert np.abs(got - ref).max() < 1e-12


def test_bps_and_ritz():
    blk5 = br.Block(N, P, 5)
    B, c2, lam = br.bps_basis_dense(N, C, blk5)
    # cohomology_N2_full.json: h^5 = 9, 27, 18, 9 in irreps C2 = 0, 2, 6, 12 (zero-weight multiplicity 1 each)
    assert B.shape[1] == 63
    vals, counts = np.unique(np.round(c2, 6), return_counts=True)
    assert list(vals) == [0, 2, 6, 12] and list(counts) == [9, 27, 18, 9]
    # gauge rotations of BPS states are BPS
    blk_r = br.Block(N, P, 5, (1, -1))
    Hr = br.SectorH(N, C, blk_r)
    for n in (10, 40):
        w = blk_r.dense(br.apply_gauge(N, P, 0, 1, blk5.sparse(B[:, n])))
        assert np.linalg.norm(Hr(w)) < 1e-9 * max(1, np.linalg.norm(w))
    # trial vectors in sector 4 from all BPS references and the three charge -1 letters, and the RR bound
    blk4 = br.Block(N, P, 4)
    H4 = br.SectorH(N, C, blk4)
    refs = [blk5.sparse(B[:, n]) for n in range(B.shape[1])]
    W, labels = br.trial_vectors(N, P, refs, br.words_of_charge(-1, 1), blk4)
    rr = br.rayleigh_ritz(H4, W)
    E0 = np.linalg.eigvalsh(H4.dense())[0]
    assert rr['values'][0] >= E0 - 1e-10
    # <v|H|v> against |Q v|^2 + |Q^dag v|^2 on one trial vector
    x = W[:, 3] / np.linalg.norm(W[:, 3])
    psi = {br.se.tuple_from_mask(m): x[n] for n, m in enumerate(blk4.masks) if x[n] != 0}
    assert abs(x @ H4(x) - fs.energy_from_Q(N, P, C, psi)) < 1e-9


if __name__ == '__main__':
    for fn in [test_word_action_and_gauge, test_bps_and_ritz]:
        print(f"{fn.__name__} ...", flush=True)
        fn()
    print("all bps_ritz tests passed")
