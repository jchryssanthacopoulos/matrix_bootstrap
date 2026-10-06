"""
Free sectors of Chen's three-matrix model: the objects of derivation D16 and their numerical checks.

D16 (docs/derivations.md).  With H = H_4 + H_1 + E_c (D2.13), H_4 = 9 sum_{a,i,k} (X_a)_ik (X_a)_ik^dag >= 0 and H_1
the one-body operator whose lowest level, -38N, is the flavour-symmetric gauge-traceless band, one has
E_0(k;N) >= 16N^3 - (15+38k)N for k <= N^2-1, with equality iff some psi in Lambda^k sl(N) (built from that band) is
annihilated by the pair-contraction map delta, equivalently has gauge Casimir N k.  Products
prod_m phi^dag(a_m)|0> of pairwise commuting, linearly independent traceless a_m are such states.

Conventions follow cohomology.py: mode index m = a N^2 + i N + j, Psi^a_ij creates mode (a,i,j), a basis state is
the sorted tuple of occupied modes, and the Jordan-Wigner sign of c^dag_m / c_m is (-1)^(number of occupied modes
below m).  Psibar^a_ij annihilates mode (a,j,i).

Provides
  * Fock tools on superpositions {sorted tuple: amplitude}: creation by a linear combination of modes,
    Q and Q^dag (energy <H> = |Q psi|^2 + |Q^dag psi|^2, independent of the D2.13 rewriting of H);
  * free_state(N, mats): prod_m phi^dag(a_m)|0> with phi = Psi^s = (Psi^1+Psi^2+Psi^3)/sqrt3;
  * sector_operators(N, C, basis): sparse H_1 and the pair-annihilation matrix Ydag with H_4 = 9 Ydag^T Ydag (D2.13);
  * wedge_sl_characters / decompose / casimir2 / weyl_dim: irreducible content of Lambda^k sl(N);
  * delta_kernel_dims(N): exact dim ker(delta) on Lambda^k sl(N), weight-blocked, over F_P.
"""
import itertools
from array import array

import numpy as np
import scipy.sparse as sp

from cohomology import _create, apply_Q, mode_index, rank_mod_p, PRIME


# ----------------------------------------------------------------------------------------------- flavour data

def flavour_tensors(C):
    """M, L, |C|^2, <C, C^rev> of D2.4 (C cyclic).  For Chen's C: M = 5/9 1 + 11/9 J, L = 4/9 1 + 11/9 J, 16/3, 5."""
    M = np.einsum('abc,abd->cd', C, C.conj())
    L = np.einsum('abc,bad->cd', C, C.conj())
    n2 = float(np.sum(C * C.conj()).real)
    rev = float(np.einsum('abc,cba->', C, C.conj()).real)
    return M, L, n2, rev


def vacuum_energy(N, C):
    """H_0 = 3 N^3 |C|^2 - 3 N <C, C^rev>  (D2.10); 16 N^3 - 15 N for Chen's C."""
    _, _, n2, rev = flavour_tensors(C)
    return 3 * N ** 3 * n2 - 3 * N * rev


def free_energy(N, k):
    """The D16 value 16 N^3 - (15 + 38 k) N."""
    return 16 * N ** 3 - (15 + 38 * k) * N


# ----------------------------------------------------------------------------------------------- Fock tools

def _annihilate(state, m):
    """c_m |state>: (sign, new_state) or None."""
    try:
        pos = state.index(m)
    except ValueError:
        return None
    return (-1) ** pos, state[:pos] + state[pos + 1:]


def apply_creation(coeffs, psi):
    """(sum_m coeffs[m] c^dag_m) psi for a superposition psi = {state: amp}."""
    out = {}
    for st, amp in psi.items():
        for m, c in coeffs.items():
            r = _create(st, m)
            if r is not None:
                s, new = r
                out[new] = out.get(new, 0) + s * c * amp
    return {st: v for st, v in out.items() if v != 0}


def apply_Q_super(N, p, C, psi):
    """Q psi, using cohomology.apply_Q on each basis state."""
    out = {}
    for st, amp in psi.items():
        for new, v in apply_Q(N, p, C, st).items():
            out[new] = out.get(new, 0) + v * amp
    return out


def apply_Qdag_super(N, p, C, psi):
    """Q^dag psi with Q^dag = sum conj(C_abc) (Psi^c_ki)^dag (Psi^b_jk)^dag (Psi^a_ij)^dag: only ordered triples of
    distinct occupied modes (a,i,j), (b,j,k), (c,k,i) contribute, so we iterate over those."""
    NN = N * N
    out = {}
    for st, amp in psi.items():
        dec = [(m // NN, (m % NN) // N, m % N) for m in st]
        n = len(st)
        for x in range(n):
            a, i, j = dec[x]
            for y in range(n):
                if y == x or dec[y][1] != j:
                    continue
                b, _, k = dec[y]
                for z in range(n):
                    if z == x or z == y or dec[z][1] != k or dec[z][2] != i:
                        continue
                    c = dec[z][0]
                    coef = np.conj(C[a, b, c])
                    if coef == 0:
                        continue
                    s1, t1 = _annihilate(st, st[x])
                    s2, t2 = _annihilate(t1, st[y])
                    s3, t3 = _annihilate(t2, st[z])
                    out[t3] = out.get(t3, 0) + coef * s1 * s2 * s3 * amp
    return out


def norm2(psi):
    return float(sum(abs(v) ** 2 for v in psi.values()))


def energy_from_Q(N, p, C, psi):
    """<psi|H|psi>/<psi|psi> with H = {Q, Q^dag}, i.e. (|Q psi|^2 + |Q^dag psi|^2)/|psi|^2."""
    return (norm2(apply_Q_super(N, p, C, psi)) + norm2(apply_Qdag_super(N, p, C, psi))) / norm2(psi)


def free_state(N, mats, p=3):
    """prod_m phi^dag(a_m) |0>, phi^dag(a) = sum_ij a_ij Psi^s_ij, Psi^s = p^{-1/2} sum_c Psi^c (the -38N band).
    The leftmost matrix in `mats` is the leftmost operator."""
    psi = {(): 1.0}
    for a in reversed(mats):
        coeffs = {}
        for i in range(N):
            for j in range(N):
                if a[i, j] != 0:
                    for c in range(p):
                        coeffs[mode_index(N, c, i, j)] = a[i, j] / np.sqrt(p)
        psi = apply_creation(coeffs, psi)
    return psi


def block_matrices(N):
    """The lfloor N^2/4 rfloor matrix units E_ij, i < r <= j, r = lfloor N/2 rfloor: pairwise products vanish, so
    they span an abelian subalgebra of sl(N) of the maximal dimension."""
    r = N // 2
    out = []
    for i in range(r):
        for j in range(r, N):
            E = np.zeros((N, N)); E[i, j] = 1.0
            out.append(E)
    return out


def cartan_matrices(N):
    """h_m = E_mm - E_{m+1,m+1}, m < N-1."""
    out = []
    for m in range(N - 1):
        h = np.zeros((N, N)); h[m, m] = 1.0; h[m + 1, m + 1] = -1.0
        out.append(h)
    return out


# ----------------------------------------------------------------------------------------------- D2.13 on a block

def sector_operators(N, C, basis):
    """Sparse H_1 (one-body part of D2.13) and Ydag on a list of k-particle basis states closed under H.

    Ydag has one row per (a, i, k, (k-2)-particle state) reached by (X_a)_ik^dag, (X_a)_ik = sum_{b,c,j} C_abc
    Psi^b_ij Psi^c_jk, so that H_4 = 9 Ydag^T Ydag on the block.  Returns (H1, Ydag)."""
    p = C.shape[0]
    NN = N * N
    M, L, _, _ = flavour_tensors(C)
    index = {s: n for n, s in enumerate(basis)}
    r1, c1, v1 = array('l'), array('l'), array('d')
    ry, cy, vy = array('l'), array('l'), array('d')
    yrow = {}
    Cc = np.conj(C)
    for col, st in enumerate(basis):
        dec = [(m // NN, (m % NN) // N, m % N) for m in st]
        for pos, (d, i, j) in enumerate(dec):
            s1 = (-1) ** pos
            rest = st[:pos] + st[pos + 1:]
            for c in range(p):
                coef = -9 * N * M[c, d].real
                if coef != 0:
                    r = _create(rest, mode_index(N, c, i, j))
                    if r is not None:
                        s2, new = r
                        r1.append(index[new]); c1.append(col); v1.append(coef * s1 * s2)
                if i == j:
                    coef = 9 * L[c, d].real
                    if coef != 0:
                        for ii in range(N):
                            r = _create(rest, mode_index(N, c, ii, ii))
                            if r is not None:
                                s2, new = r
                                r1.append(index[new]); c1.append(col); v1.append(coef * s1 * s2)
        n = len(st)
        for x in range(n):
            b, i, j = dec[x]
            sx, t1 = _annihilate(st, st[x])
            for y in range(n):
                if y == x or dec[y][1] != j:
                    continue
                c, _, kk = dec[y]
                sy, t2 = _annihilate(t1, st[y])
                for a in range(p):
                    coef = Cc[a, b, c].real
                    if coef == 0:
                        continue
                    key = (a, i, kk, t2)
                    row = yrow.setdefault(key, len(yrow))
                    ry.append(row); cy.append(col); vy.append(coef * sx * sy)
    dim = len(basis)
    H1 = sp.csr_matrix((np.frombuffer(v1, dtype=float), (np.frombuffer(r1, dtype=int), np.frombuffer(c1, dtype=int))),
                       shape=(dim, dim))
    Yd = sp.csr_matrix((np.frombuffer(vy, dtype=float), (np.frombuffer(ry, dtype=int), np.frombuffer(cy, dtype=int))),
                       shape=(len(yrow), dim))
    return H1, Yd


# ----------------------------------------------------------------------------------------------- characters

def wedge_sl_characters(N, kmax):
    """Weight multiplicities of Lambda^k sl(N), k = 0..kmax, as dicts {U(N) weight tuple: mult}."""
    zero = (0,) * N
    poly = [dict() for _ in range(kmax + 1)]
    poly[0][zero] = 1

    def times(shift):                      # multiply by (1 + t e^shift); k descending keeps poly[k-1] old
        for k in range(kmax, 0, -1):
            src, dst = poly[k - 1], poly[k]
            for w, m in src.items():
                w2 = tuple(x + y for x, y in zip(w, shift))
                dst[w2] = dst.get(w2, 0) + m

    for i in range(N):
        for j in range(N):
            if i != j:
                r = [0] * N; r[i] += 1; r[j] -= 1
                times(tuple(r))
    for _ in range(N - 1):
        times(zero)
    return poly


def decompose(N, char):
    """Irreducible multiplicities {lambda: m} of a Weyl-invariant character {weight: mult}.  Multiplying by the Weyl
    denominator prod_{alpha>0} (1 - e^{-alpha}) leaves m_lambda as the coefficient of every dominant e^lambda."""
    cur = dict(char)
    for i in range(N):
        for j in range(i + 1, N):
            new = dict(cur)
            for w, m in cur.items():
                w2 = list(w); w2[i] -= 1; w2[j] += 1
                w2 = tuple(w2)
                new[w2] = new.get(w2, 0) - m
            cur = {w: m for w, m in new.items() if m != 0}
    return {w: m for w, m in cur.items() if all(w[t] >= w[t + 1] for t in range(N - 1))}


def casimir2(lam):
    """2 c_lambda = sum_i lambda_i (lambda_i + N + 1 - 2i) (D7.1 normalisation: adjoint has c = N)."""
    N = len(lam)
    return sum(l * (l + N + 1 - 2 * (i + 1)) for i, l in enumerate(lam))


def weyl_dim(lam):
    N = len(lam)
    num = den = 1
    for i in range(N):
        for j in range(i + 1, N):
            num *= lam[i] - lam[j] + j - i
            den *= j - i
    return num // den


# ----------------------------------------------------------------------------------------------- the map delta

def sl_basis(N):
    """Integer basis of sl(N): matrix units E_ij (i != j) and h_m = E_mm - E_{m+1,m+1}.  Each entry is
    (dict {(i,j): coeff}, weight tuple)."""
    out = []
    for i in range(N):
        for j in range(N):
            if i != j:
                w = [0] * N; w[i] += 1; w[j] -= 1
                out.append(({(i, j): 1}, tuple(w)))
    for m in range(N - 1):
        out.append(({(m, m): 1, (m + 1, m + 1): -1}, (0,) * N))
    return out


def _bracket(A, B):
    out = {}
    for (i, j), x in A.items():
        for (j2, k), y in B.items():
            if j == j2:
                out[(i, k)] = out.get((i, k), 0) + x * y
    for (i, j), y in B.items():
        for (j2, k), x in A.items():
            if j == j2:
                out[(i, k)] = out.get((i, k), 0) - y * x
    return {key: v for key, v in out.items() if v != 0}


def delta_kernel_dims(N, prime=PRIME):
    """dim ker(delta) on Lambda^k sl(N) for every k, where delta(a_1 ^ ... ^ a_k) = sum_{l<m} (-1)^{l+m}
    [a_l, a_m] (x) (a_1 ^ ..^ a_l-hat ..^ a_m-hat .. ^ a_k), computed exactly over F_prime one weight block at a
    time (delta preserves weights).  This is the map whose kernel is the D16 ground space."""
    basis = sl_basis(N)
    n = len(basis)
    brackets = {(s, t): _bracket(basis[s][0], basis[t][0]) for s in range(n) for t in range(s + 1, n)}
    blocks = {}
    for k in range(n + 1):
        for S in itertools.combinations(range(n), k):
            w = tuple(sum(basis[s][1][t] for s in S) for t in range(N))
            blocks.setdefault((k, w), []).append(S)
    ker = {k: 0 for k in range(n + 1)}
    for (k, w), cols in blocks.items():
        if k < 2:
            ker[k] += len(cols)
            continue
        rows = {}
        entries = []
        for c, S in enumerate(cols):
            for l in range(k):
                for m in range(l + 1, k):
                    br = brackets[(S[l], S[m])]
                    if not br:
                        continue
                    rest = S[:l] + S[l + 1:m] + S[m + 1:]
                    sgn = (-1) ** (l + m)
                    for ij, x in br.items():
                        r = rows.setdefault((ij, rest), len(rows))
                        entries.append((r, c, sgn * x))
        if not rows:
            ker[k] += len(cols)
            continue
        Mx = np.zeros((len(rows), len(cols)), dtype=np.int64)
        for r, c, x in entries:
            Mx[r, c] += x
        ker[k] += len(cols) - rank_mod_p(Mx % prime, prime)
    return ker
