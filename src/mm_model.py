"""
Single-matrix fermionic model  Q = Tr[Psi^3]  (Chen, arXiv:2511.00790).
Matrix-index fermions Psi_ij built by Jordan-Wigner on N^2 modes.

Exposes: build_model(N) -> dict with Q, Qbar, H, Npsi, C2, |lambda>, ED spectrum,
and an oracle ev(rho, O) for expectation values.
"""
import numpy as np
import itertools

# ---------- su(N) generators (generalized Gell-Mann), Tr(T^a T^b)=1/2 delta^ab ----------
def su_generators(N):
    gens = []
    for i in range(N):
        for j in range(i + 1, N):
            S = np.zeros((N, N), complex); S[i, j] = 1; S[j, i] = 1
            gens.append(S / 2.0)
            A = np.zeros((N, N), complex); A[i, j] = -1j; A[j, i] = 1j
            gens.append(A / 2.0)
    for k in range(1, N):
        d = np.zeros(N, complex)
        for m in range(k):
            d[m] = 1
        d[k] = -k
        D = np.diag(d) * np.sqrt(1.0 / (k * (k + 1)))
        gens.append(D / np.sqrt(2))
    return gens  # length N^2-1

# ---------- Jordan-Wigner fermions ----------
def jw_creation(n):
    I = np.eye(2); Z = np.array([[1, 0], [0, -1]], float); sp = np.array([[0, 1], [0, 0]], float)
    cdag = []
    for a in range(n):
        op = np.array([[1.0]])
        for k in range(n):
            if k < a:   op = np.kron(op, Z)
            elif k == a: op = np.kron(op, sp)
            else:        op = np.kron(op, I)
        cdag.append(op.astype(complex))
    return cdag

def build_model(N):
    n = N * N                        # number of fermionic modes (matrix-index, full U(N))
    dim = 2 ** n
    cdag = jw_creation(n)            # c^dag_m,  m = i*N + j
    c = [x.conj().T for x in cdag]

    def C(i, j): return cdag[i * N + j]   # Psi_ij  (creation)
    def A(i, j): return c[i * N + j]       # (Psi_ij)^dagger = c_{ij} (annihilation)
    # barPsi_ij = Psi^dag_ji = c_{ji} = A(j,i)

    # supercharge Q = Tr[Psi^3] = sum_ijk Psi_ij Psi_jk Psi_ki
    Q = np.zeros((dim, dim), complex)
    for i in range(N):
        for j in range(N):
            for k in range(N):
                Q += C(i, j) @ C(j, k) @ C(k, i)
    Qbar = Q.conj().T
    H = Q @ Qbar + Qbar @ Q

    # fermion number / R-charge  N_Psi = Tr[Psi barPsi] = sum_m c^dag_m c_m
    Npsi = sum(cdag[m] @ c[m] for m in range(n))

    # SU(N) generators J^a = sum_{mn} (G^a)_{mn} c^dag_m c_n,  G^a = T^a (x) I - I (x) (T^a)^T
    T = su_generators(N)
    Id = np.eye(N)
    Ja = []
    for Ta in T:
        Ga = np.kron(Ta, Id) - np.kron(Id, Ta.T)   # acts on vec(M) row-major
        J = np.zeros((dim, dim), complex)
        for m in range(n):
            for nn in range(n):
                if abs(Ga[m, nn]) > 1e-14:
                    J += Ga[m, nn] * (cdag[m] @ c[nn])
        Ja.append(J)
    C2 = sum(J @ J for J in Ja)

    # Fock vacuum: state annihilated by all c_m
    vac = np.zeros(dim, complex)
    # find it: it's the common 0-eigenvector of number ops with all empty -> index where all spins "down"
    Ntot = Npsi
    # empty state has Npsi eigenvalue 0
    w = np.round(np.diag(Ntot).real)
    idx0 = int(np.where(w == 0)[0][0])
    vac[idx0] = 1.0

    # reference state |lambda> = prod_{i>j} Psi_ij |0>
    lam = vac.copy()
    for i in range(N):
        for j in range(N):
            if i > j:
                lam = C(i, j) @ lam
    nrm = np.linalg.norm(lam)
    lam = lam / nrm if nrm > 0 else lam

    return dict(N=N, n=n, dim=dim, cdag=cdag, c=c, C=C, A=A,
                Q=Q, Qbar=Qbar, H=H, Npsi=Npsi, Ja=Ja, C2=C2, vac=vac, lam=lam)


if __name__ == "__main__":
    for N in (2, 3):
        M = build_model(N)
        H, Q, Qbar, C2, Npsi, lam = M['H'], M['Q'], M['Qbar'], M['C2'], M['Npsi'], M['lam']
        dim = M['dim']
        const = 3 * N * (N * N - 1)
        idn = np.eye(dim)
        print(f"\n===== N={N}  (modes={M['n']}, Fock dim={dim}) =====")
        print(" Q^2 = 0 :", np.allclose(Q @ Q, 0))
        print(" H = {Q,Qbar} Hermitian :", np.allclose(H, H.conj().T))
        print(" operator identity H = 3N(N^2-1) - 9 C2 :", np.allclose(H, const * idn - 9 * C2))
        ev = np.linalg.eigvalsh(H).real
        vals, cts = np.unique(np.round(ev, 4), return_counts=True)
        print(" spectrum (E:deg):", dict(zip(vals, cts)))
        # reference state checks
        def q(O): return float(np.real(lam.conj() @ (O @ lam)))
        print(" |lambda> : <H> =", round(q(H), 8), " ||Q|lam>|| =", round(np.linalg.norm(Q @ lam), 8),
              " ||Qbar|lam>|| =", round(np.linalg.norm(Qbar @ lam), 8))
        print(" |lambda> : <N_Psi> =", round(q(Npsi), 6), " expected N(N-1)/2 =", N * (N - 1) // 2)
        print(" |lambda> : <C2> =", round(q(C2), 6), " expected max N(N^2-1)/3 =", N * (N * N - 1) / 3)
