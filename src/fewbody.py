"""
Few-body exact diagonalisation in a fixed N_Psi = k sector without the full Fock space.

Modes: (a, i, j) for flavor a and matrix indices i, j;  Psi^a_ij = c^dag_(a,i,j),  Psibar^a_ij = c_(a,j,i), so that
{Psi^a_ij, Psibar^b_kl} = delta^ab delta_il delta_jk, N_Psi = Tr[Psi Psibar] counts particles and the vacuum is the
state annihilated by all Psibar.  A k-particle basis state is a sorted k-tuple of modes; operators given as
trace_algebra Exprs (sums of products of trace words) are applied by explicit index sums.  Dimension binom(p N^2, k):
N=3: k<=4 (17 550); N=4: k<=3 (17 296); N=10: k=2 (44 850).  Used to cross-check the trace bootstrap (plan V4/V5).
"""
import itertools
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla


def _apply_mode(state, mode, create):
    """Apply c^dag_mode (create=True) or c_mode to a basis state (sorted tuple).  Returns (sign, new_state) or None."""
    if create:
        if mode in state:
            return None
        pos = 0
        while pos < len(state) and state[pos] < mode:
            pos += 1
        return (-1) ** pos, state[:pos] + (mode,) + state[pos:]
    if mode not in state:
        return None
    pos = state.index(mode)
    return (-1) ** pos, state[:pos] + state[pos + 1:]


def _apply_word_entry(state, word, i, j, N):
    """Apply the (i,j) entry of the open word (matrix product, operator order = matrix order) to a basis state.
    Returns dict new_state -> amplitude.  Rightmost letter acts first."""
    L = len(word)
    out = {}
    # index chains i = i_0, i_1, ..., i_L = j
    for inner in itertools.product(range(N), repeat=L - 1):
        chain = (i,) + inner + (j,)
        st, amp = state, 1
        ok = True
        for t in range(L - 1, -1, -1):          # rightmost first
            a, typ = word[t]
            r, c = chain[t], chain[t + 1]
            mode = (a, r, c) if typ == 'P' else (a, c, r)
            res = _apply_mode(st, mode, typ == 'P')
            if res is None:
                ok = False; break
            s, st = res; amp *= s
        if ok:
            out[st] = out.get(st, 0) + amp
    return out


def apply_trace_word(state_dict, word, N):
    """Apply Tr[word] to a superposition {state: amp}."""
    out = {}
    for st, amp in state_dict.items():
        for i in range(N):
            for st2, a2 in _apply_word_entry(st, word, i, i, N).items():
                out[st2] = out.get(st2, 0) + amp * a2
    return out


def apply_expr(expr, state_dict, N):
    """Apply a trace_algebra Expr (coefficients evaluated at N) to a superposition."""
    out = {}
    for mono, coeff in expr.items():
        c = coeff.at(N)
        if abs(c) < 1e-15:
            continue
        cur = dict(state_dict)
        for word in reversed(mono):
            cur = apply_trace_word(cur, word, N) if word else {s: N * a for s, a in cur.items()}
            if not cur:
                break
        for st, a in cur.items():
            out[st] = out.get(st, 0) + c * a
    return out


def sector_basis(N, p, k):
    modes = [(a, i, j) for a in range(p) for i in range(N) for j in range(N)]
    return [tuple(s) for s in itertools.combinations(modes, k)]


def sector_matrix(expr, N, p, k, basis=None, verbose=False):
    """Sparse matrix of the operator `expr` on the k-particle sector."""
    basis = basis or sector_basis(N, p, k)
    index = {s: n for n, s in enumerate(basis)}
    rows, cols, vals = [], [], []
    for n, s in enumerate(basis):
        for st, a in apply_expr(expr, {s: 1.0}, N).items():
            if abs(a) > 1e-14:
                rows.append(index[st]); cols.append(n); vals.append(a)
        if verbose and n % 2000 == 0:
            print(f"    fewbody: column {n}/{len(basis)}", flush=True)
    return sp.csr_matrix((vals, (rows, cols)), shape=(len(basis), len(basis)))


def sector_ground_energy(H_expr, N, p, k, n_eig=1, verbose=False):
    """Lowest eigenvalue(s) of H_expr on the k-particle sector (dense for dim <= 2000, else Lanczos)."""
    basis = sector_basis(N, p, k)
    Hk = sector_matrix(H_expr, N, p, k, basis, verbose)
    herm = abs(Hk - Hk.getH()).max()
    if herm > 1e-8:
        raise RuntimeError(f"sector Hamiltonian not Hermitian: {herm}")
    if len(basis) <= 2000:
        lam = np.linalg.eigvalsh(Hk.toarray())
        return lam[:n_eig], len(basis)
    lam = sla.eigsh(Hk, k=max(n_eig, 2), which='SA', return_eigenvectors=False)
    return np.sort(lam)[:n_eig], len(basis)
