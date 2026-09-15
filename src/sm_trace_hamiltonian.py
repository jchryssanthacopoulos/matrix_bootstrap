"""
Derive the single-/double-trace form of  H = {Tr[Psi^3], Tr[Psibar^3]}
for the single-matrix model, valid at all N, by matching to explicit H at N=2,3.

Letters:  P = Psi (creation matrix),  B = Psibar,  with Psibar_ab = Psi^dag_ba.
Gauge-invariant, R-neutral (2 P, 2 B) candidate operators:
  single-trace:  Tr[PPBB], Tr[PBPB]
  double-trace:  Tr[PB]Tr[PB], Tr[PB]Tr[BP], Tr[BP]Tr[BP], Tr[PP]Tr[BB]
  quadratic:     Tr[PB] (=N_Psi), Tr[BP], and identity
We solve H = sum_i c_i O_i on the operator vector space and read off c_i(N).
"""
import numpy as np
from mm_model import build_model


def letter_matrices(M):
    N = M['N']
    # operator-valued NxN matrices P_ab = Psi_ab, B_ab = Psibar_ab = Psi^dag_ba
    P = [[M['C'](a, b) for b in range(N)] for a in range(N)]           # creation
    B = [[M['A'](b, a) for b in range(N)] for a in range(N)]           # Psibar_ab = A(b,a)
    return P, B


def opmatmul(X, Y, N, dim):
    Z = [[np.zeros((dim, dim), complex) for _ in range(N)] for _ in range(N)]
    for a in range(N):
        for c in range(N):
            acc = np.zeros((dim, dim), complex)
            for b in range(N):
                acc = acc + X[a][b] @ Y[b][c]
            Z[a][c] = acc
    return Z


def optrace(X, N):
    return sum(X[a][a] for a in range(N))


def word_op(word, P, B, N, dim):
    mats = {'P': P, 'B': B}
    cur = mats[word[0]]
    for ch in word[1:]:
        cur = opmatmul(cur, mats[ch], N, dim)
    return optrace(cur, N)


def build_ops(M):
    N, dim = M['N'], M['dim']
    P, B = letter_matrices(M)
    I = np.eye(dim, dtype=complex)
    TrPB = word_op('PB', P, B, N, dim)
    TrBP = word_op('BP', P, B, N, dim)
    TrPP = word_op('PP', P, B, N, dim)
    TrBB = word_op('BB', P, B, N, dim)
    ops = {}
    # all single-trace words with 2 P and 2 B (cyclicity carries a sign for fermions,
    # so keep every ordering; the fit sorts out the relations)
    from itertools import permutations
    for w in sorted(set(''.join(p) for p in permutations('PPBB'))):
        ops['Tr[' + w + ']'] = word_op(w, P, B, N, dim)
    # double-traces
    ops['Tr[PB]Tr[PB]'] = TrPB @ TrPB
    ops['Tr[PB]Tr[BP]'] = TrPB @ TrBP
    ops['Tr[BP]Tr[BP]'] = TrBP @ TrBP
    ops['Tr[PP]Tr[BB]'] = TrPP @ TrBB
    # quadratic + identity
    ops['Tr[PB]'] = TrPB
    ops['I'] = I
    return ops


def fit_H(M):
    """Confirm H lies in the operator span (min-norm) AND give unique coefficients
    in a greedily-chosen independent sub-basis."""
    import scipy.linalg as sla
    H = M['H']
    ops = build_ops(M)
    names = list(ops.keys())
    A = np.array([ops[nm].ravel() for nm in names]).T       # dim^2 x nops
    b = H.ravel()
    # (a) completeness check
    coef_all, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    resid = np.linalg.norm(A @ coef_all - b)
    # (b) unique coefficients in an independent sub-basis
    _, R, piv = sla.qr(A, mode='economic', pivoting=True)
    r = int(np.sum(np.abs(np.diag(R)) > 1e-9))
    keep = sorted(piv[:r])
    Ai = A[:, keep]
    ci, _, _, _ = np.linalg.lstsq(Ai, b, rcond=None)
    canon = {names[keep[t]]: ci[t] for t in range(r)}
    return names, resid, canon


if __name__ == "__main__":
    print("Confirm H is a quartic single/double-trace operator, and read canonical coefficients:\n")
    canons = {}
    for N in (2, 3):
        M = build_model(N)
        names, resid, canon = fit_H(M)
        canons[N] = canon
        print(f"  N={N}:  ||H - (trace operators)|| = {resid:.2e}   (0 => H is in the span)")
        for nm, c in canon.items():
            print(f"        c[{nm:16s}] = {c.real:+.6f}")
        print()
    print("Canonical coefficients vs N (independent basis is the same set at N=2,3):")
    common = [nm for nm in canons[2] if nm in canons[3]]
    for nm in common:
        print(f"   {nm:16s}:  {canons[2][nm].real:+.5f} (N=2)   {canons[3][nm].real:+.5f} (N=3)")
