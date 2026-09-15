"""
mqm_fermion_ed.py  --  matrix-free, R-charge-resolved ground-state energy for the
fermionic matrix models, so it scales past dense/full-sparse ED.

Models (select with p):
   p = 1 : Chen's SINGLE-matrix model      Q = Tr[Psi^3]
   p = 3 : the 3-matrix "matrix SYK" model Q = sum_{i<=j<=k} Tr[Psi_i Psi_j Psi_k]

Method: work in the fixed-R-charge occupation basis. Q is a sum of 3-creation monomials.
Within charge sector k,  H_k = A A^dag + B^dag B  with  A = Q_{k-3->k},  B = Q_{k->k+3},
applied matrix-free (vectorized bit operations); smallest eigenvalue via Lanczos.
No 2^{pN^2} matrix is ever formed.

Self-contained (numpy + scipy). See run_*.py for the user-facing drivers.
"""
import numpy as np
from math import comb
from itertools import combinations
import scipy.sparse.linalg as sla


# ---------- popcount on int64 arrays ----------
def popcount(x):
    x = x.astype(np.int64)
    x = x - ((x >> 1) & 0x5555555555555555)
    x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333)
    x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f
    return ((x * 0x0101010101010101) >> 56) & 0x7f


# ---------- supercharge as a list of 3-creation monomials (mode triples) ----------
def Q_monomials(N, p):
    """Return list of (a,b,c) mode indices; Q = sum_monomials  c^dag_a c^dag_b c^dag_c."""
    def mode(f, i, j): return f * N * N + i * N + j
    mons = []
    flavor_triples = [(f, g, h) for f in range(p) for g in range(f, p) for h in range(g, p)]
    for (f, g, h) in flavor_triples:
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    mons.append((mode(f, i, j), mode(g, j, k), mode(h, k, i)))
    return mons


# ---------- basis of a fixed-charge sector ----------
def charge_states(M, k):
    """Sorted int64 array of bitmasks with popcount k (the charge-k basis)."""
    if k < 0 or k > M:
        return np.zeros(0, dtype=np.int64)
    states = np.fromiter((sum(1 << m for m in combo) for combo in combinations(range(M), k)),
                         dtype=np.int64, count=comb(M, k))
    states.sort()
    return states


# ---------- vectorized single creation c^dag_p on an array of states ----------
def _create(x, p):
    bit = np.int64(1) << p
    occupied = (x & bit) != 0
    lower = x & (bit - 1)
    sign = 1 - 2 * (popcount(lower) & 1)          # +/-1
    xnew = x | bit
    return xnew, sign.astype(np.float64), ~occupied


def _apply_Q(v, src_states, tgt_states, mons):
    """(Q v) restricted to charge j -> j+3.  v defined on src_states, output on tgt_states."""
    out = np.zeros(len(tgt_states), dtype=v.dtype)
    for (a, b, c) in mons:
        x = src_states.copy(); sgn = np.ones(len(x)); ok = np.ones(len(x), bool)
        for p in (c, b, a):                        # right-to-left operator order
            x, s, val = _create(x, p)
            sgn = sgn * s; ok = ok & val
        idx = np.where(ok)[0]
        if idx.size == 0:
            continue
        pos = np.searchsorted(tgt_states, x[idx])
        np.add.at(out, pos, sgn[idx] * v[idx])
    return out


def _apply_Qdag(v, src_states, tgt_states, mons):
    """(Q^dag v) restricted to charge j -> j-3 (adjoint of _apply_Q from j-3 to j)."""
    out = np.zeros(len(tgt_states), dtype=v.dtype)
    # Q^dag = sum c_c c_b c_a ; annihilation = reverse of creation with matching signs.
    for (a, b, c) in mons:
        x = src_states.copy(); sgn = np.ones(len(x)); ok = np.ones(len(x), bool)
        for p in (a, b, c):                        # adjoint order
            bit = np.int64(1) << p
            occ = (x & bit) != 0
            lower = x & (bit - 1)
            s = 1 - 2 * (popcount(lower) & 1)
            x = x & ~bit
            sgn = sgn * s; ok = ok & occ
        idx = np.where(ok)[0]
        if idx.size == 0:
            continue
        pos = np.searchsorted(tgt_states, x[idx])
        np.add.at(out, pos, sgn[idx] * v[idx].astype(np.float64))
    return out


def E0_charge(N, p, k, which='SA', tol=1e-6):
    """Smallest eigenvalue of H in the charge-k sector (matrix-free Lanczos)."""
    M = p * N * N
    mons = Q_monomials(N, p)
    Sk = charge_states(M, k)
    Skm = charge_states(M, k - 3)
    Skp = charge_states(M, k + 3)
    dk = len(Sk)
    if dk == 0:
        return None, 0
    if dk == 1:
        # single state: H = <s|H|s>
        v = np.array([1.0])
        Av = _apply_Qdag(v, Sk, Skm, mons) if len(Skm) else np.zeros(0)
        Bv = _apply_Q(v, Sk, Skp, mons) if len(Skp) else np.zeros(0)
        return float(Av @ Av + Bv @ Bv), 1

    def matvec(v):
        v = np.asarray(v, dtype=np.float64).ravel()
        # A A^dag v : A = Q_{k-3->k};  A^dag v maps k->k-3, then A maps k-3->k
        if len(Skm):
            AdV = _apply_Qdag(v, Sk, Skm, mons)          # k -> k-3
            AAdV = _apply_Q(AdV, Skm, Sk, mons)          # k-3 -> k
        else:
            AAdV = np.zeros(dk)
        # B^dag B v : B = Q_{k->k+3}
        if len(Skp):
            BV = _apply_Q(v, Sk, Skp, mons)              # k -> k+3
            BdBV = _apply_Qdag(BV, Skp, Sk, mons)        # k+3 -> k
        else:
            BdBV = np.zeros(dk)
        return AAdV + BdBV

    # tiny sectors: build the dense block via matvec on basis vectors (robust; ARPACK
    # is unreliable for very small dimensions). Larger sectors use matrix-free Lanczos.
    if dk <= 80:
        cols = np.empty((dk, dk))
        e = np.zeros(dk)
        for i in range(dk):
            e[:] = 0.0; e[i] = 1.0
            cols[:, i] = matvec(e)
        cols = (cols + cols.T) / 2
        e0 = float(np.linalg.eigvalsh(cols)[0])
        return max(e0, 0.0), dk
    # large sectors: matrix-free. E_0=0 sectors are highly degenerate, so 'SA' stalls;
    # get lambda_max, then largest eigenvalue of (lambda_max*I - H).
    H = sla.LinearOperator((dk, dk), matvec=matvec, dtype=np.float64)
    lam_max = float(sla.eigsh(H, k=1, which='LM', tol=tol,
                              return_eigenvectors=False, maxiter=5000)[0])
    c = lam_max * 1.01 + 1.0
    Hs = sla.LinearOperator((dk, dk), matvec=lambda v: c * np.asarray(v).ravel() - matvec(v),
                            dtype=np.float64)
    w = sla.eigsh(Hs, k=min(6, dk - 1), which='LM', tol=tol,
                  return_eigenvectors=False, maxiter=5000)
    e0 = c - float(np.max(w))
    return max(e0, 0.0), dk


def profile(N, p, charges=None, verbose=True):
    M = p * N * N
    if charges is None:
        charges = range(M + 1)
    out = {}
    for k in charges:
        e0, dk = E0_charge(N, p, k)
        out[k] = dict(dim=int(comb(M, k)), E0=(None if e0 is None else round(e0, 6)))
        if verbose and e0 is not None:
            print(f"   N_Psi={k:3d}  sector dim={comb(M,k):>12}  E_0={e0:.5f}"
                  f"{'   <-- BPS' if e0 < 1e-4 else ''}")
    return out


if __name__ == "__main__":
    # self-test against known results
    print("SELF-TEST  single matrix (p=1):")
    print(" N=2 (expect window [1,3], E0=0 there; 18 outside):")
    profile(2, 1)
    print(" N=3 (expect window [3,6]; 72,45,18 outside):")
    profile(3, 1)
    print("\nSELF-TEST  3-matrix SYK (p=3), N=2 (expect BPS only at N_Psi=5,6,7):")
    profile(2, 3)
