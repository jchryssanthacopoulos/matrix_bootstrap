"""
R-charge-sector bootstrap with adjoint-valued (open-index) operators, at finite N.

For a sector of fixed fermion number k the functional is phi(O) = Tr[rho_k P_k O P_k] with rho_k a
general Hermitian matrix on the k-block (NOT assumed positive).  Constraints (see docs/derivations.md D3):
  * adjoint channel   phi(Tr[w^dag w']) >= 0  (as a matrix over open words w of equal charge)
  * singlet channel   phi(Tr[w]^dag Tr[w']) >= 0  (over traced words, incl. Q and Qbar)
  * EOM               phi([H, X]) = 0 for every charge-neutral operator X appearing above (+ short traces)
  * optional ground-state positivity  phi(Tr[w^dag [H, w']]) >= 0 for charge-preserving words only
  * optional gauge Ward identities    phi([J^x, X]) = 0
  * normalisation phi(1) = 1.   Objective: minimise phi(H)  ->  rigorous lower bound on E_0(k).
The words are explicit finite-N operators, so this is a validation tool, not the large-N engine.
"""
import itertools, time
import numpy as np
import scipy.sparse as sp
import cvxpy as cp
from fermion_matrix_model import (build_model, chen_C, all_words, word_charge, word_label)


def _restrict(A, rows, cols):
    return A[rows][:, cols].tocsr()


class MemoryBudgetExceeded(RuntimeError):
    pass


def _letter_block(model, let, k):
    """Rectangular blocks (N x N list) of letter `let` mapping sector k -> sector k+q(let)."""
    a, t = let
    sec = model['sectors']; N = model['N']
    q = 1 if t == 'P' else -1
    if not (0 <= k + q <= model['n']):
        return None
    src, tgt = sec[k], sec[k + q]
    f = model['Psi'] if t == 'P' else model['Psibar']
    return [[_restrict(f(a, i, j), tgt, src) for j in range(N)] for i in range(N)]


def word_blocks(model, word, k, cache):
    """B^w_ij = P_{k+q(w)} w_ij P_k as rectangular sparse blocks, built letter by letter from the RIGHT
    (the word acts on sector k first with its last letter) so that nothing is ever formed on the full space.
    Returns (blocks, q) or (None, q) if the charge leaves the Hilbert space or the block vanishes."""
    N = model['N']
    q = word_charge(word)
    if not (0 <= k + q <= model['n']):
        return None, q
    d = len(model['sectors'][k])
    if len(word) == 0:
        I = sp.identity(d, dtype=complex, format='csr'); Z = sp.csr_matrix((d, d), dtype=complex)
        return [[I if i == j else Z for j in range(N)] for i in range(N)], 0
    key = (word, k)
    if key in cache:
        return cache[key], q
    prefix, last = word[:-1], word[-1]
    Lb = _letter_block(model, last, k)                 # sector k -> k + q(last)
    if Lb is None:
        return None, q
    kk = k + (1 if last[1] == 'P' else -1)
    Pb, _ = word_blocks(model, prefix, kk, cache)       # sector kk -> k + q
    if Pb is None:
        return None, q
    out = [[None] * N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            acc = None
            for m in range(N):
                term = Pb[i][m] @ Lb[m][j]
                acc = term if acc is None else acc + term
            acc = acc.tocsr(); acc.eliminate_zeros()
            out[i][j] = acc
    if max(b.nnz for row in out for b in row) == 0:
        cache[key] = None
        return None, q
    cache[key] = out
    return out, q


def sector_operators(model, k, L_adj=2, L_sing=2, L_eom=3, add_Q=True, gs=True, verbose=False):
    """Build all block-restricted operators entering the SDP for sector k (never touching the full space)."""
    sec = model['sectors']; ix = sec[k]; d = len(ix); n = model['n']; N = model['N']
    Hk = _restrict(model['H'], ix, ix)
    cache = {}
    words = []
    for w in all_words(model, max(L_adj, L_sing)):
        B, q = word_blocks(model, w, k, cache)
        if B is None:
            continue
        words.append(dict(word=w, q=q, L=len(w), B=B))
    Zd = lambda: sp.csr_matrix((d, d), dtype=complex)
    adj_groups = {}
    for wd in words:
        if wd['L'] <= L_adj:
            adj_groups.setdefault(wd['q'], []).append(wd)
    def gram_adj(wa, wb):
        acc = Zd()
        for i in range(N):
            for j in range(N):
                acc = acc + wa['B'][i][j].getH() @ wb['B'][i][j]
        return acc.tocsr()
    adj_blocks = []
    for q, grp in sorted(adj_groups.items()):
        m = len(grp)
        X = [[None] * m for _ in range(m)]
        for a in range(m):
            for b in range(a, m):
                X[a][b] = gram_adj(grp[a], grp[b])
                if b != a:
                    X[b][a] = X[a][b].getH().tocsr()
        adj_blocks.append(dict(kind='adj', q=q, labels=[word_label(g['word']) for g in grp], X=X))
    sing = []
    for wd in words:
        if wd['L'] <= L_sing:
            S = sum((wd['B'][i][i] for i in range(N)), sp.csr_matrix((len(sec[k + wd['q']]), d), dtype=complex)).tocsr()
            S.eliminate_zeros()
            if S.nnz:
                sing.append(dict(label='Tr' + word_label(wd['word']), q=wd['q'], S=S))
    if add_Q and L_sing < 3:
        for lab, Op, q in (('Q', model['Q'], 3), ('Qbar', model['Qbar'], -3)):
            if 0 <= k + q <= n:
                S = _restrict(Op, sec[k + q], ix); S.eliminate_zeros()
                if S.nnz:
                    sing.append(dict(label=lab, q=q, S=S))
    sing_groups = {}
    for s_ in sing:
        sing_groups.setdefault(s_['q'], []).append(s_)
    sing_blocks = []
    for q, grp in sorted(sing_groups.items()):
        m = len(grp)
        X = [[None] * m for _ in range(m)]
        for a in range(m):
            for b in range(a, m):
                X[a][b] = (grp[a]['S'].getH() @ grp[b]['S']).tocsr()
                if b != a:
                    X[b][a] = X[a][b].getH().tocsr()
        sing_blocks.append(dict(kind='sing', q=q, labels=[g['label'] for g in grp], X=X))
    gs_blocks = []
    if gs:
        for src, name in ((adj_groups.get(0, []), 'adj'), (sing_groups.get(0, []), 'sing')):
            m = len(src)
            if m == 0:
                continue
            Y = [[None] * m for _ in range(m)]
            for a in range(m):
                for b in range(m):
                    if name == 'adj':
                        acc = Zd()
                        for i in range(N):
                            for j in range(N):
                                Bb = src[b]['B'][i][j]
                                acc = acc + src[a]['B'][i][j].getH() @ (Hk @ Bb - Bb @ Hk)
                        Y[a][b] = acc.tocsr()
                    else:
                        Sa, Sb = src[a]['S'], src[b]['S']
                        Y[a][b] = (Sa.getH() @ (Hk @ Sb - Sb @ Hk)).tocsr()
            gs_blocks.append(dict(kind='gs-' + name, q=0, labels=[], X=Y))
    eom_ops = []
    for blk in adj_blocks + sing_blocks:
        for row in blk['X']:
            eom_ops += row
    for wd in words:
        if wd['q'] == 0 and wd['L'] <= L_eom:
            S = sum((wd['B'][i][i] for i in range(N)), Zd()).tocsr(); S.eliminate_zeros()
            if S.nnz:
                eom_ops.append(S)
    if verbose:
        print(f"  sector k={k}: d={d}, open words={len(words)}, adjoint blocks="
              f"{[len(b['labels']) for b in adj_blocks]}, singlet blocks={[len(b['labels']) for b in sing_blocks]}, "
              f"EOM ops={len(eom_ops)}", flush=True)
    return dict(k=k, d=d, ix=ix, Hk=Hk, adj_blocks=adj_blocks, sing_blocks=sing_blocks, gs_blocks=gs_blocks,
                eom_ops=eom_ops)


def _flat(X, d):
    """flattened (row-major) sparse row vector of X."""
    return sp.csr_matrix(X.reshape((1, d * d), order='C'))


def _hermitian_parts(oplist, d):
    """Distinct, unit-normalised Hermitian/anti-Hermitian parts of oplist as flattened sparse rows."""
    rows, seen = [], set()
    for X in oplist:
        for Y in ((X + X.getH()) * 0.5, (X - X.getH()) * (-0.5j)):
            Y = Y.tocsr(); Y.eliminate_zeros()
            if Y.nnz == 0:
                continue
            nrm = np.sqrt(abs(Y.multiply(Y.conj())).sum().real)
            if nrm < 1e-12:
                continue
            Yn = Y / nrm
            # cheap duplicate detection (up to sign) via a hash of the sparsity pattern and leading values
            Yn.sort_indices()
            v = Yn.data[:8]
            sgn = 1 if (v[0].real if abs(v[0].real) > 1e-12 else v[0].imag) > 0 else -1
            key = (Yn.nnz, tuple(Yn.indices[:8]), tuple(np.round(sgn * v, 8).view(float)))
            if key in seen:
                continue
            seen.add(key)
            rows.append(_flat(Yn, d))
    return sp.vstack(rows).tocsr()


def estimate_memory_bytes(n_h, nnzF, d):
    """Rough peak memory of the reduction.  Two routes: Gram of the n_h touched operators (n_h x n_h) or
    Gram of their 2d^2 real components (2d^2 x 2d^2); the smaller is used.  Dense symmetric matrix plus
    eigh workspace ~4x its size, plus sparse F copies."""
    m = min(n_h, 2 * d * d)
    return 8 * m * m * 4 + 3 * 24 * nnzF


def _reduced_basis(oplist, d, tol=1e-10, budget_bytes=6e9, verbose=False):
    """Orthonormal Hermitian basis {E_alpha} of span{Hermitian/anti-Hermitian parts of oplist}.
    Returns (F, W) with F the (n_h x d^2) sparse matrix of flattened Hermitian operators and W (n_h x r)
    such that E_alpha = sum_a W[a,alpha] O_a."""
    F = _hermitian_parts(oplist, d)
    n_h = F.shape[0]
    est = estimate_memory_bytes(n_h, F.nnz, d)
    route = 'ops' if n_h <= 2 * d * d else 'components'
    if verbose:
        print(f"    reduction: n_h={n_h} touched Hermitian ops, d^2={d*d}, nnz(F)={F.nnz}, route={route}, "
              f"est. peak memory {est/1e9:.2f} GB", flush=True)
    if est > budget_bytes:
        raise MemoryBudgetExceeded(f"estimated {est/1e9:.1f} GB > budget {budget_bytes/1e9:.1f} GB (n_h={n_h}, d={d})")
    Fr, Fi = F.real.tocsr(), F.imag.tocsr()
    if route == 'ops':
        G = np.empty((n_h, n_h))
        step = max(1, int(2e8 // max(n_h, 1)))
        for s0 in range(0, n_h, step):
            s1 = min(n_h, s0 + step)
            G[s0:s1] = (Fr[s0:s1] @ Fr.T).toarray() + (Fi[s0:s1] @ Fi.T).toarray()
        G = (G + G.T) / 2
        lam, U = np.linalg.eigh(G)
        keep = lam > tol * lam.max()
        W = U[:, keep] / np.sqrt(lam[keep])                 # E_alpha = sum_a W[a,alpha] O_a
        return F, W
    # components route: R = [Fr | Fi] is n_h x 2d^2 real; orthonormal basis of its row space from R^T R
    R = sp.hstack([Fr, Fi]).tocsr()
    M = (R.T @ R).toarray()
    M = (M + M.T) / 2
    lam, V = np.linalg.eigh(M)
    keep = lam > tol * lam.max()
    Vk = V[:, keep]                                         # orthonormal directions in component space
    # express E_alpha = sum_a W[a,alpha] O_a : need W with R^T W = Vk ; W = R Vk / lam (least squares, exact on range)
    W = (R @ Vk) / lam[keep]
    return F, W


def _independent_rows(C, tol=1e-9):
    """Orthonormal basis of the row space of a complex coefficient matrix C (rows = linear functionals of y),
    as a real matrix acting on real y.  Removes the (typically overwhelming) linear dependence among EOM rows."""
    R = np.vstack([C.real, C.imag])
    R = R[np.abs(R).sum(1) > 0]
    if R.shape[0] == 0:
        return np.zeros((0, C.shape[1]))
    U, sv, Vt = np.linalg.svd(R, full_matrices=False)
    return Vt[sv > tol * sv.max()]


class SectorSDP:
    """Builds the sector SDP once (operators, reduction, constraint coefficients); then solve() with
    different objectives / extra constraints.  See docs/derivations.md D3."""

    def __init__(self, model, k, L_adj=2, L_sing=2, L_eom=3, add_Q=True, gs=True, ward=False,
                 use_adj=True, use_sing=True, observables=None, eigenstate=False, verbose=False, budget_gb=6.0):
        t0 = time.time()
        self.model, self.k = model, k
        ops = sector_operators(model, k, L_adj, L_sing, L_eom, add_Q, gs, verbose)
        d, Hk, ix = ops['d'], ops['Hk'], ops['ix']
        self.d, self.Hk, self.ix = d, Hk, ix
        I = sp.identity(d, dtype=complex, format='csr')
        blocks = []
        if use_adj: blocks += ops['adj_blocks']
        if use_sing: blocks += ops['sing_blocks']
        blocks += ops['gs_blocks']
        eom = [Hk @ X - X @ Hk for X in ops['eom_ops']]
        if ward:
            Jk = [_restrict(J, ix, ix) for J in model['Ja']]
            eom += [J @ X - X @ J for X in ops['eom_ops'] for J in Jk]
        eom = [E.tocsr() for E in eom if E.nnz and abs(E).max() > 1e-12]
        # observables: dict name -> block-restricted operator
        obs = dict(observables or {})
        touched = [I, Hk] + [X for blk in blocks for row in blk['X'] for X in row] + eom + list(obs.values())
        # eigenstate constraints phi(O H) = E phi(O): need the operators O H_k for the Gram-entry operators O
        eig_pairs = []
        if eigenstate:
            for X in ops['eom_ops']:
                XH = (X @ Hk).tocsr()
                if XH.nnz:
                    eig_pairs.append((X, XH))
            touched += [XH for _, XH in eig_pairs]
        F, W = _reduced_basis(touched, d, budget_bytes=budget_gb * 1e9, verbose=verbose)
        self.r = W.shape[1]
        Fc = F.conj().tocsr()
        def coords(Xs, chunk=256):
            """rows c(X) with phi(X) = c(X) . y, computed in chunks to bound memory (n_h x chunk dense)."""
            out = np.empty((len(Xs), W.shape[1]), dtype=complex)
            for s0 in range(0, len(Xs), chunk):
                Xf = sp.vstack([_flat(X, d) for X in Xs[s0:s0 + chunk]]).tocsr()
                T = (Fc @ Xf.T).toarray()                 # n_h x chunk
                out[s0:s0 + chunk] = (W.T @ T).T
            return out
        self.c1 = coords([I])[0]
        self.cH = coords([Hk])[0]
        self.cobs = {name: coords([X])[0] for name, X in obs.items()}
        self.block_coeffs = []
        for blk in blocks:
            m = len(blk['X'])
            self.block_coeffs.append((m, coords([blk['X'][a][b] for a in range(m) for b in range(m)])))
        self.eom_rows = _independent_rows(coords(eom)) if eom else np.zeros((0, self.r))
        self.n_eom_raw = len(eom)
        if eigenstate:
            self.cO = coords([X for X, _ in eig_pairs]); self.cOH = coords([XH for _, XH in eig_pairs])
        self.eigenstate = eigenstate
        self.exact = float(np.linalg.eigvalsh(Hk.toarray())[0]) if d <= 3000 else float('nan')
        self.build_time = time.time() - t0
        if verbose:
            print(f"    SDP: r={self.r} variables, {self.n_eom_raw} EOM rows -> {self.eom_rows.shape[0]} independent, "
                  f"blocks {[m for m, _ in self.block_coeffs]}, build {self.build_time:.0f}s", flush=True)

    def _problem(self, objective_coeffs, sense='min', fixE=None):
        y = cp.Variable(self.r)
        cons = [cp.real(self.c1 @ y) == 1]
        for m, Cm in self.block_coeffs:
            M = cp.reshape(Cm @ y, (m, m), order='C')
            cons.append((M + cp.conj(M.T)) / 2 >> 0)
        if self.eom_rows.shape[0]:
            cons.append(self.eom_rows @ y == 0)
        if fixE is not None:
            cons.append(cp.real(self.cH @ y) == fixE)
            if self.eigenstate:
                A = _independent_rows(self.cOH - fixE * self.cO)
                if A.shape[0]:
                    cons.append(A @ y == 0)
        obj_expr = cp.real(objective_coeffs @ y) if objective_coeffs is not None else 0
        obj = cp.Minimize(obj_expr) if sense == 'min' else cp.Maximize(obj_expr)
        return cp.Problem(obj, cons), y

    def solve(self, objective='H', sense='min', fixE=None, solver='CLARABEL'):
        """objective: 'H' or an observable name (or None for a feasibility test at fixed E)."""
        c = None if objective is None else (self.cH if objective == 'H' else self.cobs[objective])
        prob, y = self._problem(c, sense, fixE)
        try:
            prob.solve(solver=solver)
        except cp.error.SolverError:
            # CLARABEL occasionally fails numerically on these (not strictly feasible) SDPs; retry with SCS,
            # bounded in iterations so a hard instance cannot stall a scan.
            try:
                prob.solve(solver='SCS', eps=1e-6, max_iters=20000)
            except Exception as e2:
                return dict(value=float('nan'), status='error: ' + str(e2)[:60])
        return dict(value=prob.value, status=prob.status)


def solve_sector(model, k, L_adj=2, L_sing=2, L_eom=3, add_Q=True, gs=True, ward=False,
                 use_adj=True, use_sing=True, solver='CLARABEL', verbose=False, budget_gb=6.0, dry_run=False):
    """Rigorous lower bound on E_0(k).  dry_run=True returns the memory estimate without forming the Gram matrix."""
    t0 = time.time()
    if dry_run:
        ops = sector_operators(model, k, L_adj, L_sing, L_eom, add_Q, gs, verbose)
        d, Hk, ix = ops['d'], ops['Hk'], ops['ix']
        I = sp.identity(d, dtype=complex, format='csr')
        blocks = ([] if not use_adj else ops['adj_blocks']) + ([] if not use_sing else ops['sing_blocks']) + ops['gs_blocks']
        eom = [E for E in (Hk @ X - X @ Hk for X in ops['eom_ops']) if E.nnz and abs(E).max() > 1e-12]
        Fd = _hermitian_parts([I, Hk] + [X for blk in blocks for row in blk['X'] for X in row] + eom, d)
        est = estimate_memory_bytes(Fd.shape[0], Fd.nnz, d)
        return dict(k=k, d=d, n_h=Fd.shape[0], nnzF=Fd.nnz, est_gb=est / 1e9, nEOM=len(eom),
                    block_sizes=[len(b['X']) for b in blocks], time=time.time() - t0)
    S = SectorSDP(model, k, L_adj, L_sing, L_eom, add_Q, gs, ward, use_adj, use_sing, verbose=verbose, budget_gb=budget_gb)
    res = S.solve('H', 'min', solver=solver)
    return dict(k=k, d=S.d, r=S.r, E_bound=res['value'], E_exact=S.exact, status=res['status'],
                nEOM=S.eom_rows.shape[0], nEOM_raw=S.n_eom_raw, time=time.time() - t0)


if __name__ == '__main__':
    import sys
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    model = build_model(N, 1, np.ones((1, 1, 1)))
    print(f"single-matrix model, N={N}: sector-resolved lower bounds on E_0(k)")
    for k in range(model['n'] + 1):
        r = solve_sector(model, k, L_adj=2, L_sing=2, L_eom=3, gs=True, verbose=True)
        print(f"   k={k:2d} d={r['d']:4d} r={r['r']:4d}: bound {r['E_bound']:9.4f}  exact {r['E_exact']:9.4f}  [{r['status']}, {r['nEOM']} EOM, {r['time']:.1f}s]")
