"""
Covariant single-trace bootstrap of E_0(N_Psi) for the single-matrix model.

Idea: use ONLY gauge-invariant single-trace operators Tr[word] as the bootstrap
operators (their algebra is what abstracts to large N).  The functional
phi(O)=Tr[rho O] is a general Hermitian functional (NO master-field/factorization
assumption).  We reduce to the operator span actually touched (SVD) so the number
of SDP variables is small and independent of the 2^{N^2} Hilbert space -> scalable.

Test: does gauge-invariant data alone bound the (non-singlet!) E_0(N_Psi) tightly?
Compared to the exact baseline E_0(N_Psi) from sm_charge_ED.
"""
import numpy as np
import scipy.linalg as sla
import cvxpy as cp
from mm_model import build_model
from sm_trace_hamiltonian import letter_matrices, word_op


# ---- gauge-invariant single-trace operators Tr[word], words in {P=Psi, B=Psibar} ----
def single_trace_ops(M, Lmax):
    N, dim = M['N'], M['dim']
    P, B = letter_matrices(M)
    ops, meta = [], []
    from itertools import product
    seen = {}
    for L in range(1, Lmax + 1):
        for w in product('PB', repeat=L):
            word = ''.join(w)
            mat = word_op(word, P, B, N, dim)
            if np.linalg.norm(mat) < 1e-9:
                continue
            # dedup by normalized matrix
            fl = mat.ravel(); idx = np.argmax(np.abs(fl))
            key = np.round(fl / fl[idx], 6).tobytes()
            if key in seen:
                continue
            seen[key] = True
            charge = word.count('P') - word.count('B')
            ops.append(mat); meta.append((word, charge))
    return ops, meta


def reduced_basis(op_list, tol=1e-9):
    """Orthonormal operator basis {E_a} (rows Vh) spanning the given operators."""
    A = np.array([o.ravel() for o in op_list])          # nops x dim^2
    # SVD to get row space
    U, s, Vh = np.linalg.svd(A, full_matrices=False)
    r = int(np.sum(s > tol))
    return Vh[:r]                                        # r x dim^2  (orthonormal rows)


_SETUP_CACHE = {}

def _setup(M, Lmom, Leom):
    """k-independent precompute: reduced basis + all coordinate vectors."""
    key = (id(M), Lmom, Leom)
    if key in _SETUP_CACHE:
        return _SETUP_CACHE[key]
    dim = M['dim']; H, Npsi = M['H'], M['Npsi']
    I = np.eye(dim, dtype=complex)
    st_ops, meta = single_trace_ops(M, Lmom)
    Bops = [I] + st_ops
    m = len(Bops)
    P = [[Bops[a].conj().T @ Bops[b] for b in range(m)] for a in range(m)]
    eom_src, _ = single_trace_ops(M, Leom)
    eom_ops = [c for c in (H @ O - O @ H for O in eom_src) if np.linalg.norm(c) > 1e-9]
    # ground-state positivity operators  Ngs_ab = <O_a^dag [H,O_b]>
    HP = [[Bops[a].conj().T @ (H @ Bops[b] - Bops[b] @ H) for b in range(m)] for a in range(m)]
    allops = ([I, H, Npsi, Npsi @ Npsi]
              + [P[a][b] for a in range(m) for b in range(m)]
              + [HP[a][b] for a in range(m) for b in range(m)] + eom_ops)
    Vh = reduced_basis(allops); r = Vh.shape[0]
    def coord(O): return Vh.conj() @ O.ravel()
    R = np.array([coord(Vh[a].reshape(dim, dim).conj().T) for a in range(r)])
    Cmat = np.array([[coord(P[a][b]) for b in range(m)] for a in range(m)]).reshape(m * m, r)
    Gmat = np.array([[coord(HP[a][b]) for b in range(m)] for a in range(m)]).reshape(m * m, r)
    setup = dict(dim=dim, m=m, r=r, Vh=Vh, R=R, Cmat=Cmat, Gmat=Gmat,
                 cI=coord(I), cH=coord(H), cN=coord(Npsi), cN2=coord(Npsi @ Npsi),
                 cEOM=[coord(c) for c in eom_ops])
    _SETUP_CACHE[key] = setup
    return setup


def covariant_bootstrap(M, k, Lmom=3, Leom=4, gs=True, solver='SCS'):
    s = _setup(M, Lmom, Leom)
    m, r = s['m'], s['r']
    y = cp.Variable(r, complex=True)
    cons = [s['R'] @ y == cp.conj(y), s['cI'] @ y == 1.0,
            s['cN'] @ y == k, s['cN2'] @ y == k * k]
    for c in s['cEOM']:
        cons.append(c @ y == 0)
    cons.append(cp.reshape(s['Cmat'] @ y, (m, m), order='C') >> 0)
    if gs:   # ground-state positivity  <O^dag [H,O]> >= 0
        Gm = cp.reshape(s['Gmat'] @ y, (m, m), order='C')
        cons.append((Gm + cp.conj(Gm.T)) / 2 >> 0)
    prob = cp.Problem(cp.Minimize(cp.real(s['cH'] @ y)), cons)
    prob.solve(solver=solver)
    return dict(status=prob.status, E0=prob.value, m=m, r=r)


if __name__ == "__main__":
    for N in (2, 3):
        M = build_model(N)
        const = 3 * N * (N * N - 1)
        # exact baseline
        pc = np.array([bin(i).count("1") for i in range(M['dim'])])
        Hd = M['H']
        print(f"\n===== N={N}: covariant single-trace bootstrap of E_0(N_Psi) =====")
        print(f"   {'N_Psi':>5} {'E0_exact':>9} {'E0_boot(gauge-inv)':>18} {'m':>4} {'r':>4} status")
        for k in range(N * N + 1):
            ix = np.where(pc == k)[0]
            e0ex = float(np.linalg.eigvalsh(Hd[np.ix_(ix, ix)])[0].real)
            try:
                out = covariant_bootstrap(M, k, Lmom=3, Leom=4)
                eb = out['E0']; st = out['status']; mm = out['m']; rr = out['r']
            except Exception as e:
                eb, st, mm, rr = float('nan'), str(e)[:20], 0, 0
            print(f"   {k:>5} {e0ex:>9.3f} {eb:>18.3f} {mm:>4} {rr:>4} {st}")
