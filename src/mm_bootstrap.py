"""
Truncated eigenstate bootstrap for Q = Tr[Psi^3], formulation 2.3(a).

Linear functional  phi(O) = Tr[rho O],  rho Hermitian (NOT assumed PSD).
Truncation ingredients:
  - positivity of the moment matrix  M_IJ = phi(B_I^dag B_J) >= 0  on a chosen
    operator set B (linearly-independent fermion words up to level L),
  - equations of motion  phi([H, C]) = 0  for words C up to level Lc,
  - normalization phi(1)=1,
Objective: min phi(H) (rigorous lower bound on E_0), max phi(C2), max/min phi(H).

The only physics input is *which* positivity / EOM constraints are imposed; rho is a
general Hermitian matrix (the most general linear functional), so a value returned by
'minH' is a genuine truncated *bound*, not exact diagonalization.
"""
import numpy as np
import scipy.linalg as sla
import cvxpy as cp
from mm_model import build_model


# ---------------- operator-word enumeration ----------------
def enumerate_words(M, L):
    n, dim = M['n'], M['dim']
    letters = []
    for m in range(n):
        letters.append(M['cdag'][m]); letters.append(M['c'][m])
    I = np.eye(dim, dtype=complex)
    seen, basis = {}, []
    def key(mat):
        flat = mat.ravel(); idx = np.argmax(np.abs(flat))
        if abs(flat[idx]) < 1e-9: return None
        return np.round(flat / flat[idx], 6).tobytes()
    basis.append(I); seen[key(I)] = True; frontier = [I]
    for _ in range(L):
        newf = []
        for w in frontier:
            for let in letters:
                nw = let @ w
                if np.linalg.norm(nw) < 1e-9: continue
                k = key(nw)
                if k is None or k in seen: continue
                seen[k] = True; basis.append(nw); newf.append(nw)
        frontier = newf
        if not newf: break
    return basis


def independent_subset(mats, tol=1e-8):
    """Greedy linearly-independent subset (column-pivoted QR on flattened operators)."""
    A = np.array([m.ravel() for m in mats]).T          # dim^2 x K
    _, R, piv = sla.qr(A, mode='economic', pivoting=True)
    r = int(np.sum(np.abs(np.diag(R)) > tol))
    keep = sorted(piv[:r])
    return [mats[i] for i in keep]


# ---------------- bootstrap SDP ----------------
def bootstrap(M, L, Lc=None, objective='minH', charge=None, casimir=None,
              extra_ops=None, solver='SCS'):
    if Lc is None: Lc = L
    dim = M['dim']
    H, C2, Npsi = M['H'], M['C2'], M['Npsi']

    words = enumerate_words(M, L)
    if extra_ops: words = words + list(extra_ops)       # e.g. supercharge Q, Qbar
    B = independent_subset(words)                        # moment-matrix operators
    m = len(B)
    Bdag = [b.conj().T for b in B]

    rho = cp.Variable((dim, dim), hermitian=True)
    y = cp.vec(rho)                                     # column-major vec

    def coord(O):                                       # phi(O) = coord(O) . vec(rho)
        return (O.T).ravel(order='F')

    # vectorized moment matrix:  Mflat[I*m+J] = phi(B_I^dag B_J)
    Amat = np.empty((m * m, dim * dim), complex)
    for I in range(m):
        for J in range(m):
            Amat[I * m + J, :] = coord(Bdag[I] @ B[J])
    Mexpr = cp.reshape(Amat @ y, (m, m), order='C')

    cons = [cp.real(cp.trace(rho)) == 1.0, Mexpr >> 0]

    # equations of motion
    nEOM = 0
    for Cw in enumerate_words(M, Lc):
        comm = H @ Cw - Cw @ H
        if np.linalg.norm(comm) < 1e-9: continue
        cons.append(coord(comm) @ y == 0); nEOM += 1

    if charge is not None:
        cons.append(coord(Npsi) @ y == charge)
        cons.append(coord(Npsi @ Npsi) @ y == charge ** 2)
    if casimir is not None:
        cons.append(coord(C2) @ y == casimir)
        cons.append(coord(C2 @ C2) @ y == casimir ** 2)

    obj = {'minH': cp.Minimize(cp.real(coord(H) @ y)),
           'maxH': cp.Maximize(cp.real(coord(H) @ y)),
           'maxC2': cp.Maximize(cp.real(coord(C2) @ y)),
           'feasible': cp.Minimize(0)}[objective]

    prob = cp.Problem(obj, cons)
    prob.solve(solver=solver)
    return dict(status=prob.status, value=prob.value, m=m, nEOM=nEOM,
                rho=rho.value)


if __name__ == "__main__":
    N = 2
    M = build_model(N)
    Q, Qbar = M['Q'], M['Qbar']
    const = 3 * N * (N * N - 1); c2max = N * (N * N - 1) / 3.0
    print(f"################  N={N}  (Fock dim {M['dim']})  ################")
    print(f"exact:  E_0=0,  E_top={const},  C2_max={c2max}\n")

    print("[Exp 1] lower bound on E_0 = min phi(H)  vs moment level L  (words only)")
    for L in (1, 2):
        r = bootstrap(M, L=L, Lc=max(L, 3), objective='minH')
        print(f"   L={L}          : m={r['m']:3d}  #EOM={r['nEOM']:3d}  ->  E_lower = {r['value']:+.5f}  [{r['status']}]")

    print("\n[Exp 1b] add the supercharge Q, Qbar to the operator set (SUSY sum-of-squares)")
    for L in (1, 2):
        r = bootstrap(M, L=L, Lc=max(L, 3), objective='minH', extra_ops=[Q, Qbar])
        print(f"   L={L} + {{Q,Qbar}}: m={r['m']:3d}  #EOM={r['nEOM']:3d}  ->  E_lower = {r['value']:+.5e}  [{r['status']}]")

    print("\n[Exp 2] max phi(H)  vs level  (exact top =", const, ")")
    for L in (1, 2):
        r = bootstrap(M, L=L, Lc=max(L, 3), objective='maxH', extra_ops=[Q, Qbar])
        print(f"   L={L} + {{Q,Qbar}}: m={r['m']:3d}  ->  max<H> = {r['value']:+.5f}  [{r['status']}]")

    print("\n[Exp 3] max phi(C2)  (exact C2_max =", c2max, ";  <C2>=(const-<H>)/9)")
    for L in (1, 2):
        r = bootstrap(M, L=L, Lc=max(L, 3), objective='maxC2', extra_ops=[Q, Qbar])
        print(f"   L={L} + {{Q,Qbar}}: m={r['m']:3d}  ->  max<C2> = {r['value']:+.5f}  [{r['status']}]")

    print("\n[Exp 4] feasibility of the exact |lambda> functional (are the constraints correct?)")
    lam = M['lam']; rl = np.outer(lam, lam.conj())
    B = independent_subset(enumerate_words(M, 2) + [Q, Qbar]); Bd = [b.conj().T for b in B]
    Mm = np.array([[np.trace(rl @ (Bd[i] @ B[j])) for j in range(len(B))] for i in range(len(B))])
    eig = np.linalg.eigvalsh((Mm + Mm.conj().T) / 2).real
    eom = max(abs(np.trace(rl @ (M['H'] @ Cw - Cw @ M['H']))) for Cw in enumerate_words(M, 3))
    print(f"   min eig M[rho_lambda] = {eig.min():.2e} (>=0 expected)   max|phi([H,C])| = {eom:.2e} (0 expected)")

    print("\n[Exp 5] state-(a) target |lambda>: min phi(H) with sharp q*=N(N-1)/2 and C2=C2max")
    r = bootstrap(M, L=2, Lc=3, objective='minH', charge=N * (N - 1) // 2, casimir=c2max, extra_ops=[Q, Qbar])
    print(f"   -> phi(H) = {r['value']:+.3e}  [{r['status']}]  (0 expected; selects a BPS state)")
    print("\nDone.")
