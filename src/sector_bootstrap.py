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
Symmetry-reduced mode (symmetry=True, docs/derivations.md D4): rho_k is restricted to SU(2)_gauge x Z_p invariant
operators (isotypic decomposition from symmetry_reduction.IsotypicReducer, N=2 only for the gauge part); every
operator is reduced to its D = sum_R mult_R^2 invariant components as soon as it is built and the touched/EOM
spans are accumulated as 2D x 2D Gram matrices, so memory no longer scales with d^2.  With an invariant
functional the Cho et al. singlet/adjoint split of the Gram matrix (adjoint_projected=True) is available and is
what makes the three-matrix k=2 sector exact at level 3.  Use Fourier flavor letters (build_model(..., fourier=True))
so that the cones are graded by the Z_p charge.
The words are explicit finite-N operators, so this is a validation tool, not the large-N engine.
"""
import itertools, time
import numpy as np
import scipy.sparse as sp
import cvxpy as cp
from fermion_matrix_model import (build_model, chen_C, all_words, word_charge, word_label, letter_op, word_zcharge)
from symmetry_reduction import IsotypicReducer


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
    return [[_restrict(letter_op(model, let, i, j), tgt, src) for j in range(N)] for i in range(N)]


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


def sector_words(model, k, L_adj, L_sing, L_eom, add_Q=True):
    """Open words (block form) for sector k and their grouping into Gram cones; cheap (no products)."""
    sec = model['sectors']; ix = sec[k]; d = len(ix); n = model['n']; N = model['N']
    cache = {}
    words = []
    for w in all_words(model, max(L_adj, L_sing, L_eom)):
        B, q = word_blocks(model, w, k, cache)
        if B is None:
            continue
        words.append(dict(word=w, q=q, z=word_zcharge(model, w), L=len(w), B=B,
                          S=sum((B[i][i] for i in range(N)), sp.csr_matrix((len(sec[k + q]), d), dtype=complex)).tocsr()))
    adj_groups = {}
    for wd in words:
        if wd['L'] <= L_adj:
            adj_groups.setdefault((wd['q'], wd['z']), []).append(wd)
    sing = []
    for wd in words:
        if wd['L'] <= L_sing:
            S = wd['S']; S.eliminate_zeros()
            if S.nnz:
                sing.append(dict(label='Tr' + word_label(wd['word']), q=wd['q'], z=wd['z'], S=S))
    if add_Q and L_sing < 3:
        for lab, Op, q in (('Q', model['Q'], 3), ('Qbar', model['Qbar'], -3)):
            if 0 <= k + q <= n:
                S = _restrict(Op, sec[k + q], ix); S.eliminate_zeros()
                if S.nnz:
                    sing.append(dict(label=lab, q=q, z=0, S=S))     # Q is Z_p neutral (C cyclic)
    sing_groups = {}
    for s_ in sing:
        sing_groups.setdefault((s_['q'], s_['z']), []).append(s_)
    eom_words = [wd for wd in words if wd['q'] == 0 and wd['z'] == 0 and wd['L'] <= L_eom and wd['S'].nnz]
    return dict(d=d, ix=ix, words=words, adj_groups=adj_groups, sing_groups=sing_groups, eom_words=eom_words)


def plan_sector(model, k, L_adj=2, L_sing=2, L_eom=3, add_Q=True, gs=True, symmetry=False):
    """Cone sizes, number of Gram-entry operators and a conservative memory estimate, without forming any product."""
    sw = sector_words(model, k, L_adj, L_sing, L_eom, add_Q)
    d = sw['d']
    adj = [len(g) for _, g in sorted(sw['adj_groups'].items())]
    sng = [len(g) for _, g in sorted(sw['sing_groups'].items())]
    gsz = ([len(sw['adj_groups'].get((0, 0), [])), len(sw['sing_groups'].get((0, 0), []))] if gs else [])
    n_entries = sum(m * (m + 1) // 2 for m in adj + sng) + sum(m * m for m in gsz)
    n_eom = n_entries + len(sw['eom_words'])
    n_ops = 2 * n_entries + n_eom + 2                      # both triangles stored, plus commutators, I, H
    sum_m2 = sum(m * m for m in adj + sng + gsz)
    if symmetry:
        Rd = IsotypicReducer(model, k, use_gauge=True, use_flavor=bool(model.get('fourier')))
        D = Rd.D
        # streaming build: upper-triangle Gram entries (D complex each) + two 2D x 2D Gram accumulators
        # + the SDP coefficient arrays (sum m^2 x r complex, r <= D) with a factor 3 for cvxpy's canonicalisation
        est = 16 * n_entries * D + 2 * 8 * (2 * D) ** 2 + 3 * 16 * sum_m2 * D
    else:
        D = None
        est = 16 * n_ops * d * d                            # dense worst case for the stored operators alone
    return dict(k=k, d=d, D=D, cones=adj + sng + gsz, sum_m2=sum_m2, n_entries=n_entries, n_eom_raw=n_eom, n_ops=n_ops,
                est_gb=est / 1e9)


def sector_operators(model, k, L_adj=2, L_sing=2, L_eom=3, add_Q=True, gs=True, verbose=False, adjoint_projected=False,
                     ward=False, eigenstate=False, transform=None, sinks=None):
    """Build all block-restricted operators entering the SDP for sector k (never touching the full space).
    transform: optional map applied to every operator as soon as it is formed (e.g. IsotypicReducer.reduce);
    the d x d matrix is then discarded, so memory scales with the reduced dimension D, not d^2.
    sinks=(touched, eom): optional _SpanAccumulator pair; then every (transformed) operator is streamed into
    `touched`, EOM/Ward rows into `eom`, EOM rows are NOT returned and only upper triangles of the Gram cones
    are kept (the lower triangle is the adjoint, phi(X^dag) = phi(X)^*)."""
    sw = sector_words(model, k, L_adj, L_sing, L_eom, add_Q)
    d, ix, words = sw['d'], sw['ix'], sw['words']
    N = model['N']
    Hk = _restrict(model['H'], ix, ix)
    T = (lambda X: X.tocsr()) if transform is None else transform
    Zd = lambda: sp.csr_matrix((d, d), dtype=complex)
    adj_groups, sing_groups = sw['adj_groups'], sw['sing_groups']
    def gram_adj(wa, wb):
        acc = Zd()
        for i in range(N):
            for j in range(N):
                acc = acc + wa['B'][i][j].getH() @ wb['B'][i][j]
        if adjoint_projected:      # Cho et al. (3.5): adjoint tensor structure = Tr[w^dag w'] - (1/N) Tr w^dag Tr w'
            acc = acc - (wa['S'].getH() @ wb['S']) / N
        return acc.tocsr()
    def comm(X):
        return (Hk @ X - X @ Hk).tocsr()
    Jk = [_restrict(J, ix, ix) for J in model['Ja']] if ward else []
    eom, eig_pairs = [], []
    streaming = sinks is not None
    if streaming:
        s_touched, s_eom = sinks
        keep_eom = lambda v: (s_touched.add(v), s_eom.add(v))
        keep_op = lambda v: (s_touched.add(v), v)[1]
    else:
        keep_eom = eom.append
        keep_op = lambda v: v
    def register(X):
        """EOM (and optional Ward / eigenstate) rows for a charge-neutral operator X; everything transformed."""
        E = comm(X)
        if E.nnz and abs(E).max() > 1e-12:
            keep_eom(T(E))
        for J in Jk:
            E = (J @ X - X @ J).tocsr()
            if E.nnz and abs(E).max() > 1e-12:
                keep_eom(T(E))
        if eigenstate:
            XH = (X @ Hk).tocsr()
            if XH.nnz:
                eig_pairs.append((keep_op(T(X)), keep_op(T(XH))))
    def gram_block(grp, prod, kind, key, labels):
        m = len(grp)
        X = [[None] * m for _ in range(m)]
        for a in range(m):
            for b in range(a, m):
                Xab = prod(grp[a], grp[b])
                register(Xab)
                X[a][b] = keep_op(T(Xab))
                if b != a:
                    Xba = Xab.getH().tocsr()
                    register(Xba)            # [H, X^dag] is not implied by [H, X] for a real functional? it is:
                    if not streaming:        # phi([H,X^dag]) = -phi([H,X])^*, so the row is redundant; kept for the
                        X[b][a] = T(Xba)     # unreduced path only to leave its numerics unchanged
        return dict(kind=kind, q=key, labels=labels, X=X)
    adj_blocks = [gram_block(grp, gram_adj, 'adj', key, [word_label(g['word']) for g in grp])
                  for key, grp in sorted(adj_groups.items())]
    sing_blocks = [gram_block(grp, lambda a, b: (a['S'].getH() @ b['S']).tocsr(), 'sing', key, [g['label'] for g in grp])
                   for key, grp in sorted(sing_groups.items())]
    gs_blocks = []
    if gs:
        for src, name in ((adj_groups.get((0, 0), []), 'adj'), (sing_groups.get((0, 0), []), 'sing')):
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
                        Y[a][b] = keep_op(T(acc.tocsr()))
                    else:
                        Sa, Sb = src[a]['S'], src[b]['S']
                        Y[a][b] = keep_op(T((Sa.getH() @ (Hk @ Sb - Sb @ Hk)).tocsr()))
            gs_blocks.append(dict(kind='gs-' + name, q=0, labels=[], X=Y))
    for wd in sw['eom_words']:
        register(wd['S'].tocsr())
    if verbose:
        print(f"  sector k={k}: d={d}, open words={len(words)}, adjoint blocks="
              f"{[len(b['labels']) for b in adj_blocks]}, singlet blocks={[len(b['labels']) for b in sing_blocks]}, "
              f"EOM rows={len(eom)}", flush=True)
    I = sp.identity(d, dtype=complex, format='csr')
    return dict(k=k, d=d, ix=ix, Hk=Hk, I_t=keep_op(T(I)), Hk_t=keep_op(T(Hk)), adj_blocks=adj_blocks,
                sing_blocks=sing_blocks, gs_blocks=gs_blocks, eom=eom, eig_pairs=eig_pairs)


def _flat(X, d):
    """flattened (row-major) sparse row vector of X."""
    return sp.csr_matrix(X.reshape((1, d * d), order='C'))


def _hermitian_parts(oplist, d, flat=None):
    """Distinct, unit-normalised Hermitian/anti-Hermitian parts of oplist, as rows of a matrix F.
    flat=None: sparse flattened d^2 rows.  flat=callable: dense rows flat(Y) (e.g. symmetry-reduced vectors)."""
    rows, seen = [], set()
    for X in oplist:
        cut = 1e-9 * np.sqrt(abs(X.multiply(X.conj())).sum().real)
        for Y in ((X + X.getH()) * 0.5, (X - X.getH()) * (-0.5j)):
            Y = Y.tocsr(); Y.eliminate_zeros()
            if Y.nnz == 0:
                continue
            nrm = np.sqrt(abs(Y.multiply(Y.conj())).sum().real)
            if nrm <= max(cut, 1e-14):              # round-off part, see _SpanAccumulator.add
                continue
            Yn = Y / nrm
            Yn.sort_indices()
            v = Yn.data[:8]
            sgn = 1 if (v[0].real if abs(v[0].real) > 1e-12 else v[0].imag) > 0 else -1
            key = (Yn.nnz, tuple(Yn.indices[:8]), tuple(np.round(sgn * v, 8).view(float)))
            if key in seen:
                continue
            seen.add(key)
            if flat is None:
                rows.append(_flat(Yn, d))
            else:
                r = flat(Yn)
                nr = np.linalg.norm(r)
                if nr > 1e-12:                      # operators with vanishing invariant part drop out
                    rows.append(r / nr)
    if flat is None:
        return sp.vstack(rows).tocsr()
    return np.array(rows) if rows else np.zeros((0, 1), dtype=complex)


class _SpanAccumulator:
    """Streaming Gram matrix M = sum_rows e e^T over the real embedding emb(v) = [Re v, Im v] in R^{2D} of the
    Hermitian and anti-Hermitian parts of reduced operator vectors.  An orthonormal basis of the accumulated
    span is read off from the eigenvectors of M.  Nothing but the 2D x 2D matrix is stored."""
    def __init__(self, reducer, chunk=2048):
        self.R, self.D, self.chunk = reducer, reducer.D, chunk
        self.M = np.zeros((2 * self.D, 2 * self.D))
        self.buf = []; self.count = 0
    @staticmethod
    def emb(v):
        return np.concatenate([v.real, v.imag])
    def parts(self, v):
        va = self.R.adjoint(v)
        return (v + va) * 0.5, (v - va) * (-0.5j)
    def add(self, v, rel_tol=1e-9):
        # a part below rel_tol * |v| is round-off (e.g. the Hermitian part of an anti-Hermitian commutator);
        # normalising it would inject a spurious unit row, so it is dropped
        cut = rel_tol * np.linalg.norm(v)
        for Y in self.parts(v):
            nrm = np.linalg.norm(Y)
            if nrm > max(cut, 1e-14):
                self.buf.append(self.emb(Y / nrm)); self.count += 1
        if len(self.buf) >= self.chunk:
            self.flush()
    def flush(self):
        if self.buf:
            B = np.array(self.buf); self.M += B.T @ B; self.buf = []
    def basis(self, tol=1e-10):
        """Orthonormal basis (columns, 2D x r) of the span of everything added."""
        self.flush()
        if self.count == 0:
            return np.zeros((2 * self.D, 0))
        lam, V = np.linalg.eigh((self.M + self.M.T) / 2)
        return V[:, lam > tol * lam.max()]


def _hermitian_parts_vec(vecs, adjoint):
    """Same as _hermitian_parts but for operators given as symmetry-reduced vectors (dense rows)."""
    rows, seen = [], set()
    for v in vecs:
        va = adjoint(v); cut = 1e-9 * np.linalg.norm(v)
        for Y in ((v + va) * 0.5, (v - va) * (-0.5j)):
            nrm = np.linalg.norm(Y)
            if nrm <= max(cut, 1e-14):              # round-off part, see _SpanAccumulator.add
                continue
            Yn = Y / nrm
            nz = np.flatnonzero(np.abs(Yn) > 1e-9)
            v0 = Yn[nz[0]]
            sgn = 1 if (v0.real if abs(v0.real) > 1e-9 else v0.imag) > 0 else -1
            key = np.round(sgn * Yn, 7).tobytes()
            if key in seen:
                continue
            seen.add(key)
            rows.append(Yn)
    return np.array(rows) if rows else np.zeros((0, len(vecs[0])), dtype=complex)


def estimate_memory_bytes(n_h, nnzF, d, D=None):
    """Rough peak memory of the reduction step.  Unreduced: Gram over min(n_h, 2d^2) with sparse F.
    Reduced (D given): dense F (n_h x D complex) plus Gram over min(n_h, 2D)."""
    if D is None:
        m = min(n_h, 2 * d * d)
        return 8 * m * m * 4 + 3 * 24 * nnzF
    m = min(n_h, 2 * D)
    return 8 * m * m * 4 + 3 * 16 * n_h * D


def _reduced_basis(F, d, tol=1e-10, budget_bytes=6e9, verbose=False, D=None):
    """Orthonormal Hermitian basis {E_alpha} of the span of the rows of F (touched Hermitian operators,
    flattened sparse or symmetry-reduced dense).  Returns W (n_h x r) with E_alpha = sum_a W[a,alpha] O_a."""
    n_h = F.shape[0]
    dense = D is not None
    nnz = int(F.nnz) if not dense else n_h * F.shape[1]
    est = estimate_memory_bytes(n_h, nnz, d, D if dense else None)
    ncomp = 2 * (D if dense else d * d)
    route = 'ops' if n_h <= ncomp else 'components'
    if verbose:
        print(f"    reduction: n_h={n_h} touched Hermitian ops, components={ncomp//2}{' (symmetry-reduced)' if dense else ''}, "
              f"route={route}, est. peak memory {est/1e9:.2f} GB", flush=True)
    if est > budget_bytes:
        raise MemoryBudgetExceeded(f"estimated {est/1e9:.1f} GB > budget {budget_bytes/1e9:.1f} GB (n_h={n_h})")
    if dense:
        Fr, Fi = F.real, F.imag
        if route == 'ops':
            G = Fr @ Fr.T + Fi @ Fi.T
            G = (G + G.T) / 2
            lam, U = np.linalg.eigh(G); keep = lam > tol * lam.max()
            return U[:, keep] / np.sqrt(lam[keep])
        R = np.hstack([Fr, Fi]); M = R.T @ R; M = (M + M.T) / 2
        lam, V = np.linalg.eigh(M); keep = lam > tol * lam.max()
        return (R @ V[:, keep]) / lam[keep]
    Fr, Fi = F.real.tocsr(), F.imag.tocsr()
    if route == 'ops':
        G = np.empty((n_h, n_h))
        step = max(1, int(2e8 // max(n_h, 1)))
        for s0 in range(0, n_h, step):
            s1 = min(n_h, s0 + step)
            G[s0:s1] = (Fr[s0:s1] @ Fr.T).toarray() + (Fi[s0:s1] @ Fi.T).toarray()
        G = (G + G.T) / 2
        lam, U = np.linalg.eigh(G); keep = lam > tol * lam.max()
        return U[:, keep] / np.sqrt(lam[keep])
    R = sp.hstack([Fr, Fi]).tocsr()
    M = (R.T @ R).toarray(); M = (M + M.T) / 2
    lam, V = np.linalg.eigh(M); keep = lam > tol * lam.max()
    return (R @ V[:, keep]) / lam[keep]


def _independent_rows(C, tol=1e-9):
    """Orthonormal basis of the row space of a complex coefficient matrix C (rows = linear functionals of y),
    as a real matrix acting on real y.  Removes the (typically overwhelming) linear dependence among EOM rows."""
    R = np.vstack([C.real, C.imag])
    R = R[np.abs(R).sum(1) > 0]
    if R.shape[0] == 0:
        return np.zeros((0, C.shape[1]))
    U, sv, Vt = np.linalg.svd(R, full_matrices=False)
    return Vt[sv > tol * sv.max()]


def solver_options(solver, eps=1e-7, max_iters=50000):
    """Solver settings that work on these SDPs.  CLARABEL: cvxpy's real embedding of complex Hermitian cones makes the
    KKT system degenerate; with chordal decomposition and equilibration switched off CLARABEL converges (k=2 level 3:
    'optimal' in 22 s), with either enabled it stops with InsufficientProgress.  SCS: always bounded."""
    if solver == 'SCS':
        return dict(eps=eps, max_iters=max_iters)
    if solver == 'CLARABEL':
        return dict(chordal_decomposition_enable=False, equilibrate_enable=False)
    return {}


class SectorSDP:
    """Builds the sector SDP once (operators, reduction, constraint coefficients); then solve() with
    different objectives / extra constraints.  See docs/derivations.md D3."""

    def __init__(self, model, k, L_adj=2, L_sing=2, L_eom=3, add_Q=True, gs=True, ward=False,
                 use_adj=True, use_sing=True, observables=None, eigenstate=False, verbose=False, budget_gb=6.0,
                 symmetry=False, adjoint_projected=False):
        """symmetry=True: restrict to G-invariant functionals (SU(2) gauge for N=2, Z_p flavor if the model
        is in the Fourier basis) via the isotypic reduction of symmetry_reduction.IsotypicReducer.
        adjoint_projected=True: use the sharper adjoint tensor structure Tr[w^dag w'] - Tr w^dag Tr w'/N
        (valid only together with symmetry=True or ward=True)."""
        t0 = time.time()
        self.model, self.k = model, k
        if adjoint_projected and not (symmetry or ward):
            raise ValueError("adjoint_projected requires a gauge-invariant functional (symmetry=True or ward=True)")
        plan = plan_sector(model, k, L_adj, L_sing, L_eom, add_Q, gs, symmetry)
        if verbose:
            print(f"    plan: cones {plan['cones']}, {plan['n_entries']} Gram entries, ~{plan['n_eom_raw']} EOM rows, "
                  f"conservative build estimate {plan['est_gb']:.2f} GB", flush=True)
        if plan['est_gb'] * 1e9 > budget_gb * 1e9:
            raise MemoryBudgetExceeded(f"build estimate {plan['est_gb']:.1f} GB > budget {budget_gb:.1f} GB "
                                       f"(n_ops={plan['n_ops']}, d={plan['d']}, D={plan['D']})")
        self.reducer = None
        transform = None
        if symmetry:
            self.reducer = IsotypicReducer(model, k, use_gauge=True, use_flavor=bool(model.get('fourier')))
            transform = self.reducer.reduce
            if verbose:
                print(f"    symmetry reduction: {len(self.reducer.blocks)} isotypic blocks, operator components "
                      f"d^2={plan['d']**2} -> D={self.reducer.D}", flush=True)
        sinks = None
        if symmetry:
            sinks = (_SpanAccumulator(self.reducer), _SpanAccumulator(self.reducer))
        ops = sector_operators(model, k, L_adj, L_sing, L_eom, add_Q, gs, verbose, adjoint_projected=adjoint_projected,
                               ward=ward, eigenstate=eigenstate, transform=transform, sinks=sinks)
        d, Hk, ix = ops['d'], ops['Hk'], ops['ix']
        self.d, self.Hk, self.ix = d, Hk, ix
        I_t, Hk_t = ops['I_t'], ops['Hk_t']
        blocks = []
        if use_adj: blocks += ops['adj_blocks']
        if use_sing: blocks += ops['sing_blocks']
        blocks += ops['gs_blocks']
        eom, eig_pairs = ops['eom'], ops['eig_pairs']
        T = (lambda X: X.tocsr()) if transform is None else transform
        obs = {name: T(X) for name, X in (observables or {}).items()}
        if transform is None:
            touched = [I_t, Hk_t] + [X for blk in blocks for row in blk['X'] for X in row] + eom + list(obs.values()) \
                      + [XH for _, XH in eig_pairs]
            F = _hermitian_parts(touched, d)
            W = _reduced_basis(F, d, budget_bytes=budget_gb * 1e9, verbose=verbose)
            Fc = F.conj().tocsr()
            def coords(Xs, chunk=256):
                """rows c(X) with phi(X) = c(X) . y, computed in chunks to bound memory (n_h x chunk dense)."""
                out = np.empty((len(Xs), W.shape[1]), dtype=complex)
                for s0 in range(0, len(Xs), chunk):
                    Xf = sp.vstack([_flat(X, d) for X in Xs[s0:s0 + chunk]]).tocsr()
                    T_ = (Fc @ Xf.T).toarray()                 # n_h x chunk
                    out[s0:s0 + chunk] = (W.T @ T_).T
                return out
            self.r = W.shape[1]
            def block_coeffs(blk):
                m = len(blk['X'])
                return coords([blk['X'][a][b] for a in range(m) for b in range(m)])
            eom_rows = _independent_rows(coords(eom)) if eom else np.zeros((0, self.r))
            self.n_eom_raw = len(eom)
        else:
            s_touched, s_eom = sinks
            for v in obs.values():
                s_touched.add(v)
            E = s_touched.basis()                              # 2D x r orthonormal, columns span the touched space
            self.r = E.shape[1]; self._E = E
            Dr = self.reducer.D
            emb, parts = _SpanAccumulator.emb, s_touched.parts
            if verbose:
                print(f"    streaming reduction: {s_touched.count} touched rows, {s_eom.count} EOM rows -> r={self.r} "
                      f"(max 2D={2*Dr})", flush=True)
            def coords(Xs, chunk=1024):
                out = np.empty((len(Xs), self.r), dtype=complex)
                for s0 in range(0, len(Xs), chunk):
                    Hs, As = zip(*(parts(v) for v in Xs[s0:s0 + chunk]))
                    out[s0:s0 + chunk] = (np.array([emb(h) for h in Hs]) @ E) + 1j * (np.array([emb(a) for a in As]) @ E)
                return out
            def block_coeffs(blk):
                m = len(blk['X'])
                C = np.empty((m, m, self.r), dtype=complex)
                for a in range(m):
                    C[a, a:] = coords(blk['X'][a][a:])
                    C[a:, a] = C[a, a:].conj()             # phi(X^dag) = phi(X)^*  ->  c(X^dag) = conj c(X)
                return C.reshape(m * m, self.r)
            Fe = s_eom.basis(); self._Fe = Fe                 # orthonormal EOM directions, already in span(E)
            eom_rows = (Fe.T @ E) if Fe.shape[1] else np.zeros((0, self.r))
            self.n_eom_raw = s_eom.count
        I, Hk = I_t, Hk_t
        self.c1 = coords([I])[0]
        self.cH = coords([Hk])[0]
        self.cobs = {name: coords([X])[0] for name, X in obs.items()}
        self.block_coeffs = []
        for blk in blocks:
            m = len(blk['X']); Cm = block_coeffs(blk)
            # diagonal congruence M -> S^{-1/2} M S^{-1/2} (rescale each word so that its diagonal coefficient row
            # has unit norm): PSD-ness is unchanged, but the solver sees O(1) entries instead of O(1)..O(10^3)
            sc = np.linalg.norm(Cm.reshape(m, m, -1)[np.arange(m), np.arange(m)], axis=1)
            sc = np.where(sc > 1e-12, 1 / np.sqrt(sc), 1.0)
            Cm = (Cm.reshape(m, m, -1) * sc[:, None, None] * sc[None, :, None]).reshape(m * m, -1)
            self.block_coeffs.append((m, Cm))
        self.eom_rows = eom_rows
        if eigenstate:
            self.cO = coords([X for X, _ in eig_pairs]); self.cOH = coords([XH for _, XH in eig_pairs])
        self.eigenstate = eigenstate
        self.exact = float(np.linalg.eigvalsh(self.Hk.toarray())[0]) if d <= 3000 else float('nan')
        self.build_time = time.time() - t0
        if verbose:
            print(f"    SDP: r={self.r} variables, {self.n_eom_raw} EOM rows -> {self.eom_rows.shape[0]} independent, "
                  f"blocks {[m for m, _ in self.block_coeffs]}, build {self.build_time:.0f}s", flush=True)

    def exact_coordinates(self, which='ground'):
        """Coordinates y of the exact ground-state functional phi(X) = Tr[P X]/Tr P (P = projector on the exact
        ground eigenspace of H_k, which is symmetry invariant) -- diagnostic for feasibility of the SDP."""
        lam, U = np.linalg.eigh(self.Hk.toarray())
        sel = lam < lam[0] + 1e-7
        P = U[:, sel] @ U[:, sel].conj().T / sel.sum()
        if self.reducer is None:
            raise NotImplementedError
        E = self._E
        y = np.empty(self.r)
        for al in range(self.r):
            v = E[:self.reducer.D, al] + 1j * E[self.reducer.D:, al]
            y[al] = np.trace(P @ self.reducer.lift(v)).real
        return y

    def check_feasibility(self, y):
        """Residuals of all constraints at coordinates y: (phi(1)-1, max |EOM row . y|, min eigenvalue per cone)."""
        out = dict(norm=float((self.c1 @ y).real - 1), eom=float(np.abs(self.eom_rows @ y).max()) if self.eom_rows.shape[0] else 0.0)
        mins = []
        for m, Cm in self.block_coeffs:
            M = (Cm @ y).reshape(m, m); M = (M + M.conj().T) / 2
            mins.append(float(np.linalg.eigvalsh(M)[0]))
        out['cone_min_eig'] = mins; out['H'] = float((self.cH @ y).real)
        return out

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

    def solve(self, objective='H', sense='min', fixE=None, solver='CLARABEL', eps=1e-7, max_iters=50000):
        """objective: 'H' or an observable name (or None for a feasibility test at fixed E).
        eps/max_iters apply to SCS (always bounded, so a hard instance cannot stall a scan)."""
        c = None if objective is None else (self.cH if objective == 'H' else self.cobs[objective])
        prob, y = self._problem(c, sense, fixE)
        try:
            prob.solve(solver=solver, **solver_options(solver, eps, max_iters))
        except cp.error.SolverError:
            # CLARABEL occasionally fails numerically on these (not strictly feasible) SDPs; retry with SCS,
            # bounded in iterations so a hard instance cannot stall a scan.
            try:
                prob.solve(solver='SCS', eps=1e-6, max_iters=20000)
            except Exception as e2:
                return dict(value=float('nan'), status='error: ' + str(e2)[:60])
        return dict(value=prob.value, status=prob.status)


def solve_sector(model, k, L_adj=2, L_sing=2, L_eom=3, add_Q=True, gs=True, ward=False,
                 use_adj=True, use_sing=True, solver='CLARABEL', verbose=False, budget_gb=6.0, dry_run=False,
                 symmetry=False, adjoint_projected=False):
    """Rigorous lower bound on E_0(k).  dry_run=True returns the memory estimate without forming the Gram matrix."""
    t0 = time.time()
    if dry_run:
        plan = plan_sector(model, k, L_adj, L_sing, L_eom, add_Q, gs, symmetry)
        plan['time'] = time.time() - t0
        return plan
    S = SectorSDP(model, k, L_adj, L_sing, L_eom, add_Q, gs, ward, use_adj, use_sing, verbose=verbose, budget_gb=budget_gb,
                  symmetry=symmetry, adjoint_projected=adjoint_projected)
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
