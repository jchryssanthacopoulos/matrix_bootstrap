"""
Symmetry reduction for the sector bootstrap (docs/derivations.md D3, D4).

For a fermion-number sector k we restrict the bootstrap functional to be invariant under
    G = SU(2)_gauge (only implemented for N=2)  x  Z_p (cyclic flavor permutation, p>1).
An invariant rho_k is block diagonal over the isotypic components of the k-block,
    rho_k = (+)_{R} sigma_R (x) 1_{dim R},   R = (j, z),
so any operator X enters only through its reduced matrices
    Xtilde^R_{ba} = sum_mu <R,b,mu| X |R,a,mu>,
and phi(X) = sum_R Tr[sigma_R Xtilde^R].  The reduced operator space has dimension sum_R m_R^2 << d_k^2.

Bases |R,a,mu> are built from an orthonormal basis of the highest-weight space of each (j,z)
(kernel of J_+ inside the J_z = j, C_2 = j(j+1), R = omega^z eigenspace) by lowering with J_-;
the normalisation is copy-independent, so the intertwiners between copies are the identity.
"""
import numpy as np
import scipy.sparse as sp
from msyk_model import su_generators


def cyclic_flavor_unitary(model):
    """Unitary U on Fock space implementing Psi^a -> Psi^{a+1 mod p} (verified: U c^dag_m U^dag = c^dag_{perm(m)})."""
    n, N, p, dim = model['n'], model['N'], model['p'], model['dim']
    perm = [((m // (N * N) + 1) % p) * N * N + m % (N * N) for m in range(n)]
    rows, cols, vals = [], [], []
    for s in range(dim):
        occ = [k for k in range(n) if (s >> (n - 1 - k)) & 1]
        new = [perm[k] for k in occ]
        sgn = 1; arr = new[:]
        for i in range(len(arr)):
            for j in range(len(arr) - 1 - i):
                if arr[j] > arr[j + 1]:
                    arr[j], arr[j + 1] = arr[j + 1], arr[j]; sgn = -sgn
        rows.append(sum(1 << (n - 1 - k) for k in arr)); cols.append(s); vals.append(sgn)
    U = sp.csr_matrix((vals, (rows, cols)), shape=(dim, dim), dtype=complex)
    # verify the intertwining property on a few modes
    cd = [model['Psi'](a, i, j) for a in range(p) for i in range(N) for j in range(N)]
    for m in (0, n // 2, n - 1):
        assert abs(U @ cd[m] @ U.getH() - cd[perm[m]]).max() < 1e-12, "flavor unitary convention error"
    return U


class IsotypicReducer:
    """Isotypic decomposition of the k-block under SU(2)_gauge (N=2) x Z_p and the map X -> reduced vector."""

    def __init__(self, model, k, use_gauge=True, use_flavor=True, tol=1e-8):
        N, p = model['N'], model['p']
        ix = model['sectors'][k]; d = len(ix)
        self.d = d
        self.use_gauge = use_gauge and N == 2
        self.use_flavor = use_flavor and p > 1
        R = lambda A: A[ix][:, ix].tocsr()
        # commuting operators to diagonalise jointly: J_z, C2 (gauge), U (flavor)
        ops = []
        if self.use_gauge:
            Ja = [R(J) for J in model['Ja']]                 # su(2) generators in the T^a = sigma^a/2 normalisation
            Jz = Ja[2]; Jp = (Ja[0] + 1j * Ja[1]).tocsr(); Jm = Jp.getH().tocsr()
            C2 = R(model['C2'])
            ops += [Jz.toarray(), C2.toarray()]
        if self.use_flavor:
            Uf = R(cyclic_flavor_unitary(model))
            Hf = (Uf + Uf.getH()).toarray() / 2; Af = ((Uf - Uf.getH()) / 2j).toarray()   # Hermitian parts of U
            ops += [Hf, Af]
        if not ops:
            self.blocks = [dict(label='all', mult=d, dim=1, V=[np.eye(d)])]
            self.D = d * d
            return
        # joint eigen-decomposition by a generic real combination (all ops commute)
        rng = np.random.default_rng(0)
        A = sum(rng.normal() * (O + O.conj().T) / 2 for O in ops)
        w, v = np.linalg.eigh(A)
        # group eigenvectors by their (rounded) eigenvalues under each operator
        labels = []
        for i in range(d):
            vec = v[:, i]; lab = []
            for O in ops:
                lab.append(np.round(np.real(vec.conj() @ (O @ vec)), 6))
            labels.append(tuple(lab))
        blocks = []
        if self.use_gauge:
            # highest-weight vectors: Jz = j with C2 = j(j+1)
            groups = {}
            for i, lab in enumerate(labels):
                jz, c2 = lab[0], lab[1]; j = (-1 + np.sqrt(1 + 4 * c2)) / 2
                if abs(jz - j) < 1e-5:
                    key = (round(2 * j) / 2,) + tuple(lab[2:])
                    groups.setdefault(key, []).append(i)
            Jm_d = Jm.toarray()
            for key, idx in sorted(groups.items()):
                j = key[0]; m = len(idx)
                # orthonormal highest-weight basis inside this (j, flavor) space
                Hw = v[:, idx]; Hw, _ = np.linalg.qr(Hw)
                V = [Hw]
                cur = Hw
                for step in range(int(round(2 * j))):
                    cur = Jm_d @ cur
                    nrm = np.linalg.norm(cur[:, 0])           # copy-independent
                    cur = cur / nrm
                    V.append(cur)
                zlab = key[1:] if self.use_flavor else ()
                blocks.append(dict(label=(j,) + zlab, mult=m, dim=len(V), V=V))
        else:
            groups = {}
            for i, lab in enumerate(labels):
                groups.setdefault(lab, []).append(i)
            for lab, idx in sorted(groups.items()):
                Vb, _ = np.linalg.qr(v[:, idx])
                blocks.append(dict(label=lab, mult=len(idx), dim=1, V=[Vb]))
        self.blocks = blocks
        self.D = int(sum(b['mult'] ** 2 for b in blocks))
        # completeness check
        tot = sum(b['mult'] * b['dim'] for b in blocks)
        assert tot == d, f"isotypic decomposition incomplete: {tot} != {d}"

    def reduce(self, X):
        """Reduced vector (concatenated Xtilde^R, row-major) of a sparse d x d operator X."""
        out = np.empty(self.D, dtype=complex); pos = 0
        for b in self.blocks:
            m = b['mult']; Xt = np.zeros((m, m), dtype=complex)
            for Vmu in b['V']:
                Xt += Vmu.conj().T @ (X @ Vmu)
            out[pos:pos + m * m] = Xt.ravel(); pos += m * m
        return out

    def adjoint(self, v):
        """Reduced vector of X^dagger from that of X (blockwise conjugate transpose)."""
        out = np.empty_like(v); pos = 0
        for b in self.blocks:
            m = b['mult']
            out[pos:pos + m * m] = v[pos:pos + m * m].reshape(m, m).conj().T.ravel(); pos += m * m
        return out

    def lift(self, v):
        """Dense d x d invariant operator with reduced vector v (inverse of reduce on invariant operators):
        X = sum_R (1/dim R) sum_mu V_mu vtilde^R V_mu^dag."""
        X = np.zeros((self.d, self.d), dtype=complex); pos = 0
        for b in self.blocks:
            m = b['mult']; vt = v[pos:pos + m * m].reshape(m, m); pos += m * m
            for Vmu in b['V']:
                X += (Vmu @ vt @ Vmu.conj().T) / b['dim']
        return X

    def summary(self):
        return [(b['label'], b['mult'], b['dim']) for b in self.blocks]
