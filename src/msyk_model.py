"""
Matrix SYK model (Chen arXiv:2511.00790, sec 4):  p fermion matrices Psi_f,
cubic supercharge   Q = sum_{1<=i<=j<=k<=p} Tr[Psi_i Psi_j Psi_k],  H = {Q, Qbar}.

Purely fermionic; modes m=(f,i,j), f=0..p-1, i,j=0..N-1  ->  p*N^2 modes.
Sparse Jordan-Wigner. Default p=3, N=2  ->  12 complex fermions, dim 4096.
"""
import numpy as np
import scipy.sparse as sp
import itertools


def su_generators(N):
    gens = []
    for i in range(N):
        for j in range(i + 1, N):
            S = np.zeros((N, N), complex); S[i, j] = 1; S[j, i] = 1; gens.append(S / 2)
            A = np.zeros((N, N), complex); A[i, j] = -1j; A[j, i] = 1j; gens.append(A / 2)
    for k in range(1, N):
        d = np.zeros(N, complex); d[:k] = 1; d[k] = -k
        gens.append(np.diag(d) * np.sqrt(1.0 / (k * (k + 1))) / np.sqrt(2))
    return gens


def jw_creation_sparse(n):
    Z = sp.csr_matrix(np.array([[1, 0], [0, -1]], complex))
    sp_ = sp.csr_matrix(np.array([[0, 1], [0, 0]], complex))
    I2 = sp.identity(2, format='csr', dtype=complex)
    cdag = []
    for a in range(n):
        op = sp.identity(1, format='csr', dtype=complex)
        for k in range(n):
            op = sp.kron(op, Z if k < a else (sp_ if k == a else I2), format='csr')
        cdag.append(op)
    return cdag


def build_msyk(N=2, p=3):
    n = p * N * N
    dim = 2 ** n
    cdag = jw_creation_sparse(n)
    c = [x.getH().tocsr() for x in cdag]

    def C(f, i, j): return cdag[f * N * N + i * N + j]     # Psi^f_ij (creation)

    # supercharge Q = sum_{f<=g<=h} sum_{ijk} Psi^f_ij Psi^g_jk Psi^h_ki
    Q = sp.csr_matrix((dim, dim), dtype=complex)
    for f in range(p):
        for g in range(f, p):
            for h in range(g, p):
                for i in range(N):
                    for j in range(N):
                        for k in range(N):
                            Q = Q + C(f, i, j) @ C(g, j, k) @ C(h, k, i)
    Q = Q.tocsr()
    Qbar = Q.getH().tocsr()
    H = (Q @ Qbar + Qbar @ Q).tocsr()

    Npsi = sp.csr_matrix((dim, dim), dtype=complex)
    for m in range(n):
        Npsi = Npsi + cdag[m] @ c[m]
    Npsi = Npsi.tocsr()

    # gauge SU(N): J^a = sum_f sum_{mn} (G^a)_{mn} c^dag_{f,m} c_{f,n},  G^a=T^a⊗I - I⊗T^aT
    T = su_generators(N); Id = np.eye(N)
    Ja = []
    for Ta in T:
        Ga = np.kron(Ta, Id) - np.kron(Id, Ta.T)
        J = sp.csr_matrix((dim, dim), dtype=complex)
        for f in range(p):
            base = f * N * N
            for a in range(N * N):
                for b in range(N * N):
                    if abs(Ga[a, b]) > 1e-14:
                        J = J + Ga[a, b] * (cdag[base + a] @ c[base + b])
        Ja.append(J.tocsr())
    C2gauge = sum((J @ J for J in Ja), sp.csr_matrix((dim, dim), dtype=complex)).tocsr()

    # flavor number operators  Nf[f] = sum_ij n_{f,i,j}
    Nf = []
    for f in range(p):
        Nfop = sp.csr_matrix((dim, dim), dtype=complex)
        for i in range(N):
            for j in range(N):
                m = f * N * N + i * N + j
                Nfop = Nfop + cdag[m] @ c[m]
        Nf.append(Nfop.tocsr())

    return dict(N=N, p=p, n=n, dim=dim, cdag=cdag, c=c, C=C,
                Q=Q, Qbar=Qbar, H=H, Npsi=Npsi, Ja=Ja, C2gauge=C2gauge, Nf=Nf)


if __name__ == "__main__":
    M = build_msyk(N=2, p=3)
    dim, Q, Qbar, H, Npsi = M['dim'], M['Q'], M['Qbar'], M['H'], M['Npsi']
    print(f"p={M['p']}, N={M['N']}, modes={M['n']}, dim={dim}")
    print(" Q^2 = 0 :", abs(Q @ Q).max() < 1e-9)
    print(" H Hermitian :", abs((H - H.getH())).max() < 1e-9)
    # [H, Npsi] = 0
    print(" [H,N_Psi]=0 :", abs((H @ Npsi - Npsi @ H)).max() < 1e-9)
    # [H, J^a]=0 (gauge)
    comm = max(abs((H @ J - J @ H)).max() for J in M['Ja'])
    print(" [H,J^a]=0 (gauge SU(N)) :", comm < 1e-9)
    # [Q, Npsi] = 3 Q  (R-charge of Q)
    print(" [N_Psi,Q]=3Q :", abs((Npsi @ Q - Q @ Npsi) - 3 * Q).max() < 1e-9)
    # flavor S3 check: does permuting flavors leave H invariant? test one transposition via relabel
    print(" H nnz:", H.nnz, " ||H||_max:", round(abs(H).max(), 3))
