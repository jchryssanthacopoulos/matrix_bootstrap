"""
Exact finite-N trace bootstrap of the R-charge sector ground energy E_0(k; N) (docs/trace_bootstrap_plan.md, M2).

Variables: phi(m) for every canonical neutral monomial m (product of trace words) touched by the constraints;
real unknowns z = (Re phi(m), Im phi(m)).  Constraints (all exact at the given integer N, see plan 2.2):
  * normalisation phi(1) = 1;  reality phi(m^dag) = conj phi(m);
  * EOM  phi([H, X]) = 0 and sector rows phi((N_Psi - k) X) = phi(X (N_Psi - k)) = 0 for X in the Gram entries and
    the neutral trace words up to L_eom;
  * optional finite-N trace relations (antisymmetriser over N+1 indices) up to total length finiteN_len;
  * cones: singlet  phi(T_a^dag T_b) >= 0 over trace words of equal (q,z);  adjoint-projected
    phi(Tr[w_a^dag w_b]) - phi(Tr w_a^dag Tr w_b)/N >= 0 over open words of equal (q,z)  (plain adjoint if
    adjoint_projected=False).
The conic program is assembled in standard form  A z + s = b, s in {0} x PSD x ...  and handed directly to
Clarabel or SCS (no cvxpy).  Objective: minimise Re phi(H)  ->  rigorous lower bound on E_0(k; N).
"""
import itertools, time
import numpy as np
import scipy.sparse as sp
import trace_algebra as ta
from trace_algebra import NPoly, Expr


class MemoryBudgetExceeded(RuntimeError):
    pass


def _all_words(letters, Lmax, Lmin=1):
    for L in range(Lmin, Lmax + 1):
        for w in itertools.product(letters, repeat=L):
            yield tuple(w)


def _expr_key(e):
    """Key identifying an Expr up to a complex scalar factor (for de-duplication of constraint rows)."""
    items = sorted(e.items(), key=lambda kv: ta.mono_key(kv[0]))
    m0, c0 = items[0]
    k0 = max(c0.items(), key=lambda kv: abs(kv[1]))          # normalise by the largest coefficient of the leading monomial
    norm = k0[1]
    return tuple((m, tuple(sorted((k, np.round(v / norm, 9)) for k, v in c.items()))) for m, c in items)


def _dedupe_exprs(exprs):
    seen, out = set(), []
    for e in exprs:
        key = _expr_key(e)
        if key not in seen:
            seen.add(key); out.append(e)
    return out


def _independent_rows(A_eq, b_eq, max_dense=4e7, tol=1e-10):
    """Drop linearly dependent equality rows (dense pivoted QR on the transpose) when affordable."""
    m, n = A_eq.shape
    if m * n > max_dense:
        return A_eq, b_eq, False
    from scipy.linalg import qr
    M = A_eq.toarray()
    _, R, piv = qr(M.T, mode='economic', pivoting=True)
    d = np.abs(np.diag(R)); r = int((d > tol * d[0]).sum())
    keep = np.sort(piv[:r])
    return A_eq[keep], b_eq[keep], True


# ---- parallel fixed-N build: module-level state for forked workers -------------------------------------------------
_W = {}


def _w_init(C, p, k, N, L_eom_gram, adjoint_projected):
    _W['H'] = ta.hamiltonian(C); _W['Nk'] = ta.number_operator(p) + ta.scalar(-k)
    _W['N'] = N; _W['L_eom_gram'] = L_eom_gram; _W['adj'] = adjoint_projected


def _w_adj_entry(pair):
    """Adjoint(-projected) Gram entry for open words (wa, wb) and, if short enough, its EOM and sector rows."""
    wa, wb = pair
    wad = ta.word_dagger(wa)
    A = ta.trace(wad + wb)
    if _W['adj']:
        A = A - ta.canonical((wad, wb)) * NPoly.N(-1)
    N = _W['N']
    rows = []
    if _W['L_eom_gram'] is None or len(wa) + len(wb) <= _W['L_eom_gram']:
        rows = [ta.commutator(_W['H'], A).at(N), ta.mul(_W['Nk'], A).at(N)]
    return A.at(N), rows


def _w_sing_entry(pair):
    Ea, Eb = pair                      # Exprs of the two traced operators
    S = ta.mul(ta.dagger(Ea), Eb)
    N = _W['N']
    rows = []
    if _W['L_eom_gram'] is None or S.total_length() <= _W['L_eom_gram']:
        rows = [ta.commutator(_W['H'], S).at(N), ta.mul(_W['Nk'], S).at(N)]
    return S.at(N), rows


def _w_word_rows(e):
    N = _W['N']
    return [ta.commutator(_W['H'], e).at(N), ta.mul(_W['Nk'], e).at(N)]


def _w_dagger(m):
    return m, ta.canonical(tuple(ta.word_dagger(w) for w in reversed(m))).at(_W['N'])


class TraceSDP:
    def __init__(self, C, p, k, L_adj=2, L_sing=2, L_eom=3, finiteN_len=0, N_for_relations=None,
                 adjoint_projected=True, fourier=True, gs=False, verbose=False):
        """C: cyclic couplings in the FLAVOR basis (p x p x p).  The algebra runs in the Z_p Fourier basis when
        fourier=True (cones graded by the Z_p charge).  All Exprs carry N-polynomial coefficients; the numerical
        value of N enters only in assemble(N) -- except the finite-N relations, which are generated for the integer
        N_for_relations (default: none)."""
        t0 = time.time()
        self.p, self.k = p, k
        self.fourier = fourier and p > 1
        self.C = ta.fourier_C(np.asarray(C, dtype=complex)) if self.fourier else np.asarray(C, dtype=complex)
        self.letters = [(a, t) for a in range(p) for t in ('P', 'B')]
        self.H = ta.hamiltonian(self.C)
        self.Npsi = ta.number_operator(p)
        self.Nk = self.Npsi + ta.scalar(-k)
        zch = (lambda w: ta.word_z(w, p)) if self.fourier else (lambda w: 0)
        # ---- cones -------------------------------------------------------------------------------------------
        self.cones = []                 # dict(name, labels, entries: {(a,b): Expr} for a<=b)
        groups = {}
        for w in _all_words(self.letters, L_adj):
            groups.setdefault((ta.word_q(w), zch(w)), []).append(w)
        for key, ws in sorted(groups.items()):
            ent = {}
            for a in range(len(ws)):
                wa = ta.word_dagger(ws[a])
                for b in range(a, len(ws)):
                    A = ta.trace(wa + ws[b])
                    if adjoint_projected:
                        A = A - ta.canonical((wa, ws[b])) * NPoly.N(-1)
                    ent[(a, b)] = A
            self.cones.append(dict(name=f"adj q={key[0]} z={key[1]}", labels=ws, entries=ent))
        # singlet cone over trace words: keep words with distinct non-zero canonical traces
        sing = {}
        seen = {}
        for w in _all_words(self.letters, L_sing):
            e = ta.trace(w)
            if not e:
                continue
            key = tuple(sorted((m, tuple(sorted(c.items()))) for m, c in e.items()))
            if key in seen:
                continue
            seen[key] = w
            sing.setdefault((ta.word_q(w), zch(w)), []).append(w)
        if L_sing < 3:                  # Q, Qbar explicitly (charge +-3, neutral under Z_p)
            for lab, e in (('Q', ta.supercharge(self.C)), ('Qbar', ta.dagger(ta.supercharge(self.C)))):
                sing.setdefault(('Q' if lab == 'Q' else 'Qbar', 0), []).append(e)
        for key, ws in sorted(sing.items(), key=lambda kv: str(kv[0])):
            ent = {}
            exprs = [w if isinstance(w, Expr) else None for w in ws]
            for a in range(len(ws)):
                for b in range(a, len(ws)):
                    if exprs[a] is None and exprs[b] is None:
                        ent[(a, b)] = ta.canonical((ta.word_dagger(ws[a]), ws[b]))
                    else:
                        Ea = exprs[a] if exprs[a] is not None else ta.trace(ws[a])
                        Eb = exprs[b] if exprs[b] is not None else ta.trace(ws[b])
                        ent[(a, b)] = ta.mul(ta.dagger(Ea), Eb)
            self.cones.append(dict(name=f"sing q={key[0]} z={key[1]}", labels=ws, entries=ent))
        if gs:
            raise NotImplementedError("ground-state positivity not implemented in the trace engine yet")
        # ---- equality rows -------------------------------------------------------------------------------------
        Xs = [e for cone in self.cones for e in cone['entries'].values()]
        for w in _all_words(self.letters, L_eom):
            if ta.word_q(w) == 0 and zch(w) == 0:
                e = ta.trace(w)
                if e:
                    Xs.append(e)
        self.rows = []                  # list of Expr that must vanish
        for X in Xs:
            self.rows.append(ta.commutator(self.H, X))
            self.rows.append(ta.mul(self.Nk, X))
            self.rows.append(ta.mul(X, self.Nk))
        self.rows.append(self.Nk)       # X = 1
        self.rows = _dedupe_exprs([r for r in self.rows if r])
        self.n_eom_rows = len(self.rows)
        # ---- variables -----------------------------------------------------------------------------------------
        monos = set([()])
        for cone in self.cones:
            for e in cone['entries'].values():
                monos.update(e.keys())
        for r in self.rows:
            monos.update(r.keys())
        monos.update(self.H.keys())
        for m in monos:
            assert ta.mono_q(m) == 0 and (not self.fourier or ta.mono_z(m, p) == 0), m
        # finite-N relations restricted to touched monomials
        self.n_rel = 0
        if finiteN_len and N_for_relations:
            Nr = N_for_relations
            for rel in self._finiteN_rows(Nr, finiteN_len, zch):
                if rel and all(m in monos for m in rel):
                    self.rows.append(rel); self.n_rel += 1
        # close the variable set under the dagger (needed for the reality rows), then build those rows
        daggers = {}
        todo = list(monos)
        while todo:
            m = todo.pop()
            md = ta.canonical(tuple(ta.word_dagger(w) for w in reversed(m)))
            daggers[m] = md
            for m2 in md:
                if m2 not in monos:
                    monos.add(m2); todo.append(m2)
        self.monos = sorted(monos, key=ta.mono_key)
        self.index = {m: i for i, m in enumerate(self.monos)}
        self.reality = []               # conj phi(m) - phi(m^dag) = 0
        for m in self.monos:
            md = daggers[m]
            if len(md) == 1 and m in md and abs(md[m].at(1) - 1) < 1e-12 and abs(md[m].at(2) - 1) < 1e-12:
                self.reality.append((m, None))          # Hermitian monomial: Im phi(m) = 0
            else:
                self.reality.append((m, md))
        self.build_time = time.time() - t0
        if verbose:
            print(f"  TraceSDP k={k} level ({L_adj},{L_sing},{L_eom}) finiteN_len={finiteN_len}: {len(self.monos)} monomials, "
                  f"cones {[len(c['labels']) for c in self.cones]}, {self.n_eom_rows} EOM/sector rows, {self.n_rel} finite-N rows, "
                  f"build {self.build_time:.1f}s", flush=True)

    @classmethod
    def fixed_N(cls, C, p, k, N, L_adj=2, L_sing=2, L_eom=3, L_eom_gram=None, adjoint_projected=True, fourier=True,
                workers=8, verbose=False):
        """Build directly at numerical N with all coefficients evaluated (no N-polynomials stored) and the algebra
        run in `workers` forked processes.  L_eom_gram: EOM/sector rows only for Gram entries of total word length
        <= L_eom_gram (None = all); at level 4 the length-8 entries' commutators create ~10^6 length-10 monomials that
        carry almost no information (Hilbert-space experience: 9 new EOM directions out of 3x10^5 rows)."""
        import multiprocessing as mp
        Pool = mp.get_context('fork').Pool        # fork: no re-import of the caller's script (macOS defaults to spawn)
        self = cls.__new__(cls)
        t0 = time.time()
        self.p, self.k, self.N_fixed = p, k, N
        self.fourier = fourier and p > 1
        self.C = ta.fourier_C(np.asarray(C, dtype=complex)) if self.fourier else np.asarray(C, dtype=complex)
        self.letters = [(a, t) for a in range(p) for t in ('P', 'B')]
        self.H = ta.hamiltonian(self.C); self.Npsi = ta.number_operator(p); self.Nk = self.Npsi + ta.scalar(-k)
        zch = (lambda w: ta.word_z(w, p)) if self.fourier else (lambda w: 0)
        groups = {}
        for w in _all_words(self.letters, L_adj):
            groups.setdefault((ta.word_q(w), zch(w)), []).append(w)
        sing, seen = {}, {}
        for w in _all_words(self.letters, L_sing):
            e = ta.trace(w)
            if not e:
                continue
            key = tuple(sorted((m, tuple(sorted(c.items()))) for m, c in e.items()))
            if key in seen:
                continue
            seen[key] = w; sing.setdefault((ta.word_q(w), zch(w)), []).append(e)
        if L_sing < 3:
            Q = ta.supercharge(self.C)
            sing.setdefault(('Q', 0), []).append(Q); sing.setdefault(('Qbar', 0), []).append(ta.dagger(Q))
        word_ops = [ta.trace(w) for w in _all_words(self.letters, L_eom) if ta.word_q(w) == 0 and zch(w) == 0]
        word_ops = [e for e in word_ops if e]
        self.cones, self.rows = [], []
        with Pool(workers, initializer=_w_init, initargs=(self.C, p, k, N, L_eom_gram, adjoint_projected)) as pool:
            for key, ws in sorted(groups.items()):
                pairs = [(ws[a], ws[b]) for a in range(len(ws)) for b in range(a, len(ws))]
                ent = {}
                for (a, b), (A, rows) in zip(((a, b) for a in range(len(ws)) for b in range(a, len(ws))),
                                             pool.imap(_w_adj_entry, pairs, chunksize=32)):
                    ent[(a, b)] = A; self.rows += rows
                self.cones.append(dict(name=f"adj q={key[0]} z={key[1]}", labels=ws, entries=ent))
                if verbose:
                    print(f"    cone {key}: {len(ws)} words, {time.time()-t0:.0f}s", flush=True)
            for key, es in sorted(sing.items(), key=lambda kv: str(kv[0])):
                pairs = [(es[a], es[b]) for a in range(len(es)) for b in range(a, len(es))]
                ent = {}
                for (a, b), (S, rows) in zip(((a, b) for a in range(len(es)) for b in range(a, len(es))),
                                             pool.imap(_w_sing_entry, pairs, chunksize=8)):
                    ent[(a, b)] = S; self.rows += rows
                self.cones.append(dict(name=f"sing q={key[0]} z={key[1]}", labels=es, entries=ent))
            for rows in pool.imap(_w_word_rows, word_ops, chunksize=8):
                self.rows += rows
            self.rows.append(self.Nk.at(N))
            self.rows = [r for r in self.rows if r]
            self.n_eom_rows = len(self.rows); self.n_rel = 0
            monos = set([()])
            for cone in self.cones:
                for e in cone['entries'].values():
                    monos.update(e.keys())
            for r in self.rows:
                monos.update(r.keys())
            monos.update(self.H.keys())
            # dagger closure + reality rows
            daggers = {}
            todo = list(monos)
            while todo:
                batch, todo = todo, []
                for m, md in pool.imap_unordered(_w_dagger, batch, chunksize=64):
                    daggers[m] = md
                    for m2 in md:
                        if m2 not in monos:
                            monos.add(m2); todo.append(m2)
        self.monos = sorted(monos, key=ta.mono_key)
        self.index = {m: i for i, m in enumerate(self.monos)}
        self.reality = []
        for m in self.monos:
            md = daggers[m]
            if len(md) == 1 and m in md and abs(md[m] - 1) < 1e-12:
                self.reality.append((m, None))
            else:
                self.reality.append((m, md))
        self._evaluated = True
        self.build_time = time.time() - t0
        if verbose:
            print(f"  TraceSDP(fixed N={N}) k={k} level ({L_adj},{L_sing},{L_eom}) L_eom_gram={L_eom_gram}: {len(self.monos)} monomials, "
                  f"cones {[len(c['labels']) for c in self.cones]}, {self.n_eom_rows} EOM/sector rows, build {self.build_time:.0f}s", flush=True)
        return self

    def _finiteN_rows(self, N, Lmax, zch):
        """Antisymmetriser relations with N+1 blocks of total length <= Lmax, neutral in (q,z)."""
        out = []
        for total in range(N + 1, Lmax + 1):
            for comp in itertools.product(range(1, total - N + 1), repeat=N + 1):
                if sum(comp) != total:
                    continue
                for lets in itertools.product(self.letters, repeat=total):
                    if ta.word_q(lets) != 0 or zch(lets) != 0:
                        continue
                    blocks, pos = [], 0
                    for c in comp:
                        blocks.append(lets[pos:pos + c]); pos += c
                    e = ta.finiteN_relations(blocks, N)
                    if e:
                        out.append(e)
        return out

    # ------------------------------------------------------------------------------------------------------------
    def _ev(self, obj, N):
        """Coefficients at N: Expr -> dict, or an already-evaluated dict."""
        if getattr(self, '_evaluated', False):
            return obj
        return obj.at(N)

    def assemble(self, N, solver='clarabel', prune=True):
        """Standard-form data at numerical N.  Real unknowns z[2i] = Re phi(m_i), z[2i+1] = Im phi(m_i)."""
        if getattr(self, '_evaluated', False):
            assert N == self.N_fixed, "this instance was built at a fixed N"
        n = 2 * len(self.monos)
        rows_i, rows_j, rows_v, b = [], [], [], []
        r = 0
        def add_complex_row(coeffs):        # sum_m c_m phi(m) = 0  ->  two real rows
            nonlocal r
            # a row whose coefficients are all round-off (cancellations leave ~1e-16) must not be normalised into a
            # spurious O(1) constraint; drop such rows and the round-off entries of genuine rows (genuine
            # coefficients are >= 1/9 in magnitude)
            coeffs = {m: c for m, c in coeffs.items() if abs(c) > 1e-9}     # absolute floor: genuine |c| >~ 1e-3
            if not coeffs:
                return
            for part in (0, 1):
                for m, c in coeffs.items():
                    i = self.index[m]
                    if part == 0:   # real part: Re c Re x - Im c Im x
                        rows_i.extend([r, r]); rows_j.extend([2 * i, 2 * i + 1]); rows_v.extend([c.real, -c.imag])
                    else:           # imag part: Im c Re x + Re c Im x
                        rows_i.extend([r, r]); rows_j.extend([2 * i, 2 * i + 1]); rows_v.extend([c.imag, c.real])
                b.append(0.0); r += 1
        # normalisation
        i0 = self.index[()]
        rows_i += [r]; rows_j += [2 * i0]; rows_v += [1.0]; b.append(1.0); r += 1
        rows_i += [r]; rows_j += [2 * i0 + 1]; rows_v += [1.0]; b.append(0.0); r += 1
        # reality
        for m, md in self.reality:
            i = self.index[m]
            if md is None:
                rows_i += [r]; rows_j += [2 * i + 1]; rows_v += [1.0]; b.append(0.0); r += 1
            else:
                # conj x_m - sum_j d_j x_j = 0 :  real: u_m - sum(Re d u_j - Im d v_j);  imag: -v_m - sum(Im d u_j + Re d v_j)
                cj = {mj: d for mj, d in self._ev(md, N).items() if abs(d) > 1e-10}
                rows_i += [r]; rows_j += [2 * i]; rows_v += [1.0]
                for mj, d in cj.items():
                    j = self.index[mj]; rows_i += [r, r]; rows_j += [2 * j, 2 * j + 1]; rows_v += [-d.real, d.imag]
                b.append(0.0); r += 1
                rows_i += [r]; rows_j += [2 * i + 1]; rows_v += [-1.0]
                for mj, d in cj.items():
                    j = self.index[mj]; rows_i += [r, r]; rows_j += [2 * j, 2 * j + 1]; rows_v += [-d.imag, -d.real]
                b.append(0.0); r += 1
        # EOM / sector / finite-N
        for e in self.rows:
            add_complex_row(self._ev(e, N))
        n_eq = r
        A_eq = sp.csc_matrix((rows_v, (rows_i, rows_j)), shape=(n_eq, n)); b_eq = np.array(b)
        # normalise every equality row to unit 2-norm (rows mix coefficients from N^0 to N^4; without this the
        # rank-revealing pruning and the solver's equilibration lose the small rows at large N)
        rn = np.sqrt(np.asarray(A_eq.multiply(A_eq).sum(axis=1)).ravel()); rn[rn == 0] = 1.0
        A_eq = sp.diags(1.0 / rn) @ A_eq; b_eq = b_eq / rn
        A_eq = A_eq.tocsc()
        if prune:
            A_eq, b_eq, pruned = _independent_rows(A_eq, b_eq)
            n_eq = A_eq.shape[0]
        else:
            pruned = False
        rows_i, rows_j, rows_v, b = [], [], [], list(b_eq); r = n_eq
        # PSD cones: real embedding [[X, -Y], [Y, X]] of the Hermitian M = X + iY, M_ab = sum_m c_m phi(m)
        psd_dims = []
        for cone in self.cones:
            d = len(cone['labels']); D = 2 * d
            psd_dims.append(D)
            # entries of the embedding as linear functions of z: collect (row_index_in_svec, col z, value)
            # s = b - A z must equal svec(R(z)) -> A = -svec-map, b = 0
            entries = {}       # (I,J) with I<=J -> dict z_col -> value
            def acc(I, J, col, val):
                if I > J:
                    I, J = J, I
                entries.setdefault((I, J), {}); entries[(I, J)][col] = entries[(I, J)].get(col, 0.0) + val
            for (a, bb), e in cone['entries'].items():
                for m, c in self._ev(e, N).items():
                    if abs(c) < 1e-12:
                        continue
                    i = self.index[m]
                    # M_ab = c x = (c_r + i c_i)(u + i v) = (c_r u - c_i v) + i (c_i u + c_r v)
                    Xr = [(2 * i, c.real), (2 * i + 1, -c.imag)]
                    Yi = [(2 * i, c.imag), (2 * i + 1, c.real)]
                    for col, val in Xr:      # X_ab = X_ba
                        acc(a, bb, col, val); acc(a + d, bb + d, col, val)
                    if a != bb:
                        for col, val in Yi:  # Y_ab = -Y_ba: R[a, b+d] = -Y_ab, R[b, a+d] = -Y_ba = Y_ab
                            acc(a, bb + d, col, -val); acc(bb, a + d, col, val)
                    # a == bb: Y_aa = 0 by hermiticity (enforced by reality rows; not added)
            order = []
            if solver == 'clarabel':
                for J in range(D):
                    for I in range(J + 1):
                        order.append((I, J))
            else:
                for J in range(D):
                    for I in range(J, D):
                        order.append((J, I))
            for (I, J) in order:
                for col, val in entries.get((I, J), {}).items():
                    rows_i.append(r); rows_j.append(col); rows_v.append(-val * (1.0 if I == J else np.sqrt(2)))
                b.append(0.0); r += 1
        A = sp.vstack([A_eq, sp.csc_matrix((rows_v, (rows_i, rows_j)), shape=(r, n))[n_eq:]]).tocsc()
        c = np.zeros(n)
        for m, cf in self.H.at(N).items():
            i = self.index[m]; c[2 * i] += cf.real; c[2 * i + 1] += -cf.imag     # Re(cf x)
        # 't Hooft rescaling of the unknowns: phi(m) = N^{w(m)} y_m with w = sum_traces (1 + L_i/2) (the natural
        # upper bound on the size of a trace word), objective divided by N^3.  Undone in solve().
        w = np.array([sum(1 + len(word) / 2 for word in m) for m in self.monos])
        colscale = np.repeat(np.power(float(N), w), 2)
        A = (A @ sp.diags(colscale)).tocsc()
        c = c * colscale / float(N) ** 3
        # reality at the variable level: z = T y with y the independent real unknowns (Hermitian monomials keep
        # only Re; a monomial whose dagger is a single monomial times a phase shares its unknowns with the partner)
        T = self._reality_substitution(N)
        A = (A @ T).tocsr(); c = T.T @ c; b = np.array(b)
        # equality rows made empty by the substitution (the reality rows it absorbed) would break the solvers'
        # row scaling: drop them
        rown = np.sqrt(np.asarray(A.multiply(A).sum(axis=1)).ravel())
        nz = rown > 1e-8                          # rows were unit-normalised before the substitution
        nz[n_eq:] = True
        assert np.all(np.abs(b[:n_eq][~nz[:n_eq]]) < 1e-8)
        A = A[nz].tocsc(); b = b[nz]; n_eq = int(nz[:n_eq].sum())
        A.data[np.abs(A.data) < 1e-13] = 0.0; A.eliminate_zeros()
        return dict(A=A, b=b, c=c, n_eq=n_eq, psd_dims=psd_dims, n=T.shape[1], n_full=n, pruned=pruned,
                    colscale=colscale, objscale=float(N) ** 3, T=T)

    def _reality_substitution(self, N):
        """Sparse T (2*n_monos x n_red) with z = T y.  Uses the reality data: Hermitian monomials -> Im = 0;
        single-term daggers  conj x_m = c x_r  ->  x_m = conj(c) conj(x_r)  eliminate (u_m, v_m) in favour of the
        representative r (the smaller monomial); everything else stays independent (and is constrained by the
        reality rows already in A)."""
        n_m = len(self.monos)
        rep = {}                                  # monomial index -> (r, conj(c)) meaning x_m = conj(c) * conj(x_r)
        herm = set()
        for m, md in self.reality:
            i = self.index[m]
            if md is None:
                herm.add(i); continue
            mdv = self._ev(md, N)
            if len(mdv) == 1:
                (mr, cval), = mdv.items()
                if abs(abs(cval) - 1) < 1e-12 and mr != m:
                    j = self.index[mr]
                    if j < i:                     # eliminate the larger index in favour of the smaller
                        rep[i] = (j, np.conj(cval))
        keep_cols = []                            # list of (z_col) that remain independent, in order
        col_of = {}
        for i in range(n_m):
            if i in rep:
                continue
            col_of[(i, 0)] = len(keep_cols); keep_cols.append(2 * i)
            if i not in herm:
                col_of[(i, 1)] = len(keep_cols); keep_cols.append(2 * i + 1)
        rows, cols, vals = [], [], []
        for i in range(n_m):
            if i in rep:
                j, cc = rep[i]
                assert j not in rep               # dagger is an involution and rep only points to smaller indices
                # x_i = cc * conj(x_j):  u_i = Re(cc) u_j + Im(cc) v_j ;  v_i = Im(cc) u_j - Re(cc) v_j
                ju = col_of[(j, 0)]; jv = col_of.get((j, 1))
                rows += [2 * i]; cols += [ju]; vals += [cc.real]
                if jv is not None:
                    rows += [2 * i]; cols += [jv]; vals += [cc.imag]
                rows += [2 * i + 1]; cols += [ju]; vals += [cc.imag]
                if jv is not None:
                    rows += [2 * i + 1]; cols += [jv]; vals += [-cc.real]
            else:
                rows += [2 * i]; cols += [col_of[(i, 0)]]; vals += [1.0]
                if i not in herm:
                    rows += [2 * i + 1]; cols += [col_of[(i, 1)]]; vals += [1.0]
        return sp.csc_matrix((vals, (rows, cols)), shape=(2 * n_m, len(keep_cols)))

    def solve(self, N, solver='clarabel', eps=1e-8, max_iters=100000, verbose=False, budget_gb=10.0, prune=None, indirect=False,
              adaptive_scale=True):
        """Clarabel: pruned equality rows, chordal decomposition off, equilibration on (settings found by sweep,
        2026-09-18).  SCS: unpruned rows (first-order methods tolerate redundancy; the pruned system converged worse)."""
        if prune is None:
            prune = (solver == 'clarabel')
        data = self.assemble(N, solver, prune=prune)
        A, b, c = data['A'], data['b'], data['c']
        est = (A.nnz * 12 + A.shape[0] * 16) * 3 / 1e9
        if est > budget_gb:
            raise MemoryBudgetExceeded(f"estimated {est:.1f} GB > {budget_gb} GB")
        t0 = time.time()
        if solver == 'clarabel':
            import clarabel
            P = sp.csc_matrix((data['n'], data['n']))
            cones = [clarabel.ZeroConeT(data['n_eq'])] + [clarabel.PSDTriangleConeT(D) for D in data['psd_dims']]
            st = clarabel.DefaultSettings(); st.verbose = verbose
            st.chordal_decomposition_enable = False
            st.tol_gap_abs = eps; st.tol_gap_rel = eps; st.tol_feas = eps; st.max_iter = 500
            st.static_regularization_constant = 1e-7 if data['pruned'] else 1e-5   # rank-deficient equality block
            sol = clarabel.DefaultSolver(P, c, A, b, cones, st).solve()
            status, x, val = str(sol.status), np.array(sol.x), float(sol.obj_val)
        else:
            import scs
            # rho_x = 1e-3 (default 1e-6): with the redundant, rank-deficient equality block SCS stalls at the default
            # (found 2026-09-18 on N=2, k=1: 20 000 iterations without convergence vs 'solved' in 5 925)
            # linear solver: Apple Accelerate LDL when available (level-4 KKT: 120 s setup, 0.7 s/iteration, 7 GB,
            # where QDLDL had not finished its factorisation after 30 min); the CG 'indirect' solver is hopeless here
            import platform
            lin = {"linear_solver": "cpu_indirect"} if indirect else ({"linear_solver": "accelerate"} if platform.system() == 'Darwin' else {})
            # adaptive_scale=False: each rescaling re-factorises the KKT system (~2 min at level 4), which turned a
            # 0.7 s/iteration solve into a >4 h one; without it the scale is fixed after setup
            sol = scs.SCS(dict(A=A, b=b, c=c), dict(z=data['n_eq'], s=data['psd_dims']), eps_abs=eps, eps_rel=eps,
                          max_iters=max_iters, verbose=verbose, rho_x=1e-3, adaptive_scale=adaptive_scale, **lin).solve()
            status, x, val = sol['info']['status'], np.array(sol['x']), float(sol['info']['pobj'])
        val = val * data['objscale']; x = (data['T'] @ x) * data['colscale']
        self.last = dict(status=status, value=val, x=x, solve_time=time.time() - t0, n_rows=A.shape[0], nnz=A.nnz,
                         n=data['n'], n_eq=data['n_eq'], psd_dims=data['psd_dims'], pruned=data['pruned'])
        return self.last

    # ------------------------------------------------------------------------------------------------------------
    def exact_functional(self, model):
        """phi_GS(m) for every monomial from the exact sector ground state (multiplet-averaged) of `model`
        (must be the same N, p, C, basis).  Diagnostic for small N."""
        ix = model['sectors'][self.k]
        Hk = model['H'][ix][:, ix].toarray()
        lam, U = np.linalg.eigh(Hk); sel = lam < lam[0] + 1e-7
        vecs = U[:, sel]
        full = np.zeros((model['dim'], vecs.shape[1]), dtype=complex); full[ix] = vecs
        from fermion_matrix_model import word_matrix
        Nn = model['N']; cache = {}
        def T(word):
            if word not in cache:
                W = word_matrix(model, word)
                cache[word] = sum(W[i][i] for i in range(Nn)).tocsr() if word else Nn * sp.identity(model['dim'], format='csr')
            return cache[word]
        phi = {}
        for m in self.monos:
            val = 0j
            for j in range(full.shape[1]):
                v = full[:, j]
                for word in reversed(m):
                    v = T(word) @ v
                val += np.vdot(full[:, j], v)
            phi[m] = val / full.shape[1]
        return phi, float(lam[0])

    def check_exact(self, model):
        """Residuals of all constraints on the exact ground-state functional (rows, reality, cones)."""
        N = model['N']
        phi, E0 = self.exact_functional(model)
        row_res = max((abs(sum(c * phi[m] for m, c in self._ev(e, N).items())) for e in self.rows), default=0.0)
        real_res = 0.0
        for m, md in self.reality:
            if md is None:
                real_res = max(real_res, abs(phi[m].imag))
            else:
                real_res = max(real_res, abs(np.conj(phi[m]) - sum(d * phi[mj] for mj, d in self._ev(md, N).items())))
        cone_min = []
        for cone in self.cones:
            d = len(cone['labels']); M = np.zeros((d, d), dtype=complex)
            for (a, b), e in cone['entries'].items():
                M[a, b] = sum(c * phi[m] for m, c in self._ev(e, N).items()); M[b, a] = np.conj(M[a, b])
            cone_min.append(float(np.linalg.eigvalsh(M)[0]) if d else 0.0)
        EH = sum(c * phi[m] for m, c in self.H.at(N).items())
        return dict(E0=E0, phi_H=EH, max_row_residual=row_res, max_reality_residual=real_res, cone_min_eig=cone_min)
