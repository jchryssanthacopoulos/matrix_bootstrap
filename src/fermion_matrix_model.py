"""
General purely fermionic matrix model with a cubic supercharge

    Q = C_abc Tr[Psi^a Psi^b Psi^c],   H = {Q, Qbar},   a,b,c = 0..p-1,

built explicitly (sparse Jordan-Wigner) at finite N.  Conventions follow Chen 2025:
Psi^a_ij are creation operators, Psibar^a_ij = (Psi^a_ji)^dagger, {Psi^a_ij, Psibar^b_kl} = d^ab d_il d_jk.

Mode index: m = a*N*N + i*N + j.  NOTE: in this JW representation a creation operator *clears*
a bit, so basis index 0 is the FILLED state and index 2^n-1 is the Fock vacuum.  Use
`charge_sectors` (built from the N_Psi operator) rather than popcounts.

Provides: chen_C(p), cyclic(C), build_model(N, p, C) -> dict with letters, Q, Qbar, H, Npsi,
gauge generators Ja and Casimir C2, and helpers for open-index words.
"""
import itertools
import numpy as np
import scipy.sparse as sp
from msyk_model import jw_creation_sparse, su_generators


def cyclic(C):
    return (C + np.transpose(C, (1, 2, 0)) + np.transpose(C, (2, 0, 1))) / 3


def chen_C(p=3):
    """Cyclically symmetrised coefficients of Q = sum_{a<=b<=c} Tr[Psi^a Psi^b Psi^c]."""
    Ct = np.zeros((p, p, p))
    for a, b, c in itertools.product(range(p), repeat=3):
        if a <= b <= c:
            Ct[a, b, c] = 1
    return cyclic(Ct)


def build_model(N, p, C, fourier=False):
    """fourier=True: letters are the Z_p Fourier combinations Psi^(m) = p^{-1/2} sum_c w^{mc} Psi^c (m = 0..p-1),
    which carry definite Z_p flavor charge (+m for Psi^(m), -m for its conjugate); requires C cyclically symmetric."""
    n = p * N * N
    dim = 2 ** n
    cd = jw_creation_sparse(n)
    an = [x.getH().tocsr() for x in cd]

    def Psi(a, i, j):            # creation
        return cd[a * N * N + i * N + j]

    def Psibar(a, i, j):         # annihilation, Psibar_ij = (Psi_ji)^dagger
        return an[a * N * N + j * N + i]

    Z = lambda: sp.csr_matrix((dim, dim), dtype=complex)
    Q = Z()
    for a, b, c in itertools.product(range(p), repeat=3):
        if abs(C[a, b, c]) < 1e-14:
            continue
        for i, j, k in itertools.product(range(N), repeat=3):
            Q = Q + C[a, b, c] * (Psi(a, i, j) @ Psi(b, j, k) @ Psi(c, k, i))
    Q = Q.tocsr()
    Qbar = Q.getH().tocsr()
    H = (Q @ Qbar + Qbar @ Q).tocsr()
    Npsi = sum((cd[m] @ an[m] for m in range(n)), Z()).tocsr()

    T = su_generators(N)
    Id = np.eye(N)
    Ja = []
    for Ta in T:
        Ga = np.kron(Ta, Id) - np.kron(Id, Ta.T)
        J = Z()
        for a in range(p):
            base = a * N * N
            for u in range(N * N):
                for v in range(N * N):
                    if abs(Ga[u, v]) > 1e-14:
                        J = J + Ga[u, v] * (cd[base + u] @ an[base + v])
        Ja.append(J.tocsr())
    C2 = sum((J @ J for J in Ja), Z()).tocsr()

    charge = np.rint(Npsi.diagonal().real).astype(int)
    sectors = {k: np.where(charge == k)[0] for k in range(n + 1)}
    return dict(N=N, p=p, n=n, dim=dim, C=C, Psi=Psi, Psibar=Psibar, Q=Q, Qbar=Qbar, H=H,
                Npsi=Npsi, Ja=Ja, C2=C2, charge=charge, sectors=sectors, fourier=fourier)


def letter_op(model, let, i, j):
    """Sparse operator for matrix element (i,j) of letter let=(a, 'P'|'B').  In the Fourier basis a is the Z_p charge."""
    a, t = let
    p = model['p']
    if not model.get('fourier') or p == 1:
        return model['Psi'](a, i, j) if t == 'P' else model['Psibar'](a, i, j)
    w = np.exp(2j * np.pi / p)
    if t == 'P':
        return sum((w ** (a * c) * model['Psi'](c, i, j) for c in range(p)), sp.csr_matrix((model['dim'], model['dim']), dtype=complex)) / np.sqrt(p)
    return sum((w ** (-a * c) * model['Psibar'](c, i, j) for c in range(p)), sp.csr_matrix((model['dim'], model['dim']), dtype=complex)) / np.sqrt(p)


def word_zcharge(model, word):
    """Z_p flavor charge of a word in the Fourier basis (0 if not in the Fourier basis)."""
    if not model.get('fourier') or model['p'] == 1:
        return 0
    return sum((a if t == 'P' else -a) for a, t in word) % model['p']


# ---------------------------------------------------------------- open-index words
# A "letter" is (flavor a, kind) with kind 'P' (Psi) or 'B' (Psibar); its charge is +1 / -1.
# An open word w is an operator-valued N x N matrix: (w)_ij = sum_m (l1)_im (l2)_mj ... .

def letters(model):
    return [(a, t) for a in range(model['p']) for t in ('P', 'B')]


def letter_matrix(model, let):
    N = model['N']
    return [[letter_op(model, let, i, j) for j in range(N)] for i in range(N)]


def word_matrix(model, word):
    """Operator-valued N x N matrix of the open word (tuple of letters); word () is the identity."""
    N, dim = model['N'], model['dim']
    if len(word) == 0:
        I = sp.identity(dim, dtype=complex, format='csr')
        Z = sp.csr_matrix((dim, dim), dtype=complex)
        return [[I if i == j else Z for j in range(N)] for i in range(N)]
    mats = [letter_matrix(model, l) for l in word]
    cur = mats[0]
    for M in mats[1:]:
        new = [[None] * N for _ in range(N)]
        for i in range(N):
            for j in range(N):
                acc = sp.csr_matrix((dim, dim), dtype=complex)
                for m in range(N):
                    acc = acc + cur[i][m] @ M[m][j]
                new[i][j] = acc.tocsr()
        cur = new
    return cur


def word_charge(word):
    return sum(1 if t == 'P' else -1 for _, t in word)


def all_words(model, Lmax, normal_ordered=False):
    """All open words up to length Lmax.  normal_ordered=True keeps only words with every Psi left of
    every Psibar.  WARNING: this is NOT a lossless reduction -- reordering fermion operators preserves the
    index-contraction pattern of the original word, and the resulting operators are not matrix products
    of the reordered letters (e.g. (Psi Psibar Psi)_ij = N Psi_ij - Psi_im Psi_nj Psibar_mn).  Verified
    2026-09-16: it weakens the N=4, k=5 bound from 11.4 to 0.  Kept only as an option for experiments."""
    out = [()]
    lets = letters(model)
    for L in range(1, Lmax + 1):
        for w in itertools.product(lets, repeat=L):
            kinds = [t for _, t in w]
            if normal_ordered and any(kinds[i] == 'B' and kinds[j] == 'P' for i in range(L) for j in range(i + 1, L)):
                continue
            out.append(w)
    return out


def word_label(word):
    if not word:
        return '1'
    return ''.join(('Ψ' if t == 'P' else 'Ψ̄') + str(a + 1) for a, t in word)
