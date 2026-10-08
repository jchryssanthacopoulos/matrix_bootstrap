"""
Singlet-sector trace bootstrap for the U(n)^3 fermionic quiver (docs/derivations.md D21; research/notes/
quiver_project.md section 12).  Algebra: src/quiver_trace.py on top of src/trace_algebra.py.

State.  A density matrix rho supported on U(n)^3 gauge singlets with N_A = N_B = N_C = m (degree k = 3m), and the
functional phi(X) = Tr(rho X) on gauge-invariant operators (canonical products of traces of closed walks).  For real
couplings C every operator built here is a real matrix in the occupation basis and H, Q have real eigenbases, so
phi can be taken real: phi(m) = phi(m^dag) in R.  This halves every cone relative to the complex engine
(src/trace_bootstrap.py) and is valid for the ground state of each sector and for any BPS state.

Constraints, each an exact operator statement at the given integer n:
  normalisation  phi(1) = 1;  reality  phi(m^dag) = phi(m);
  sector         phi((N_e - m) X) = 0 for the three edge numbers N_e (gauge-invariant X commute with N_e);
  singlet        phi(Chat_v X) = 0 for the three node Casimirs (n x the SU(n) quadratic Casimir; >= 0, zero only on
                 node-v singlets; the U(1) parts are the sector rows); optional Gauss rows phi(Y Tr[W G_v]) = 0 with
                 G_v the node-v u(n) generator matrix and W an open word v -> v;
  EOM            phi([H, X]) = 0 (eigenstates of H and their mixtures);
  finite n       node-wise antisymmetrisers over n+1 indices (quiver_trace.finite_n_relation);
  BPS (opt.)     phi(X Q) = phi(Q X) = phi(X Qbar) = phi(Qbar X) = 0, phi(X H) = phi(H X) = 0;
  cones          adjoint(-projected) Gram matrices phi(Tr[W_a^dag W_b]) (- phi(Tr W_a^dag Tr W_b)/n if the words
                 are v -> v) over open words with common endpoints and edge charges (Schur: blocks of different
                 endpoints or charges are exactly orthogonal); singlet Gram matrices phi(T_a^dag T_b) over traced
                 closed walks of common edge charges, with the identity in the neutral block.
Modes: energy (min phi(H): a rigorous lower bound on the singlet E_0(k) up to solver tolerance), margin (max t with
every cone >= t 1; t* < 0 means infeasible) and feasibility (pure feasibility with a Farkas certificate).
"""
import itertools, time
import numpy as np
import scipy.sparse as sp
import trace_algebra as ta
from trace_algebra import NPoly, Expr
import quiver_trace as qt


# ------------------------------------------------------------------------------------------------ small helpers
def mono_charges(mono):
    q = [0, 0, 0]
    for w in mono:
        for (e, f), t in w:
            q[e] += 1 if t == 'P' else -1
    return tuple(q)


def ev(expr, n, tol=1e-11):
    """Expr -> {monomial: float} at numerical n.  Coefficients must be real (real couplings)."""
    out = {}
    for m, c in expr.items():
        v = complex(c.at(n))
        if abs(v) > tol:
            if abs(v.imag) > 1e-9 * max(1.0, abs(v)):
                raise ValueError(f'complex coefficient {v} at {ta.mono_label(m)}: the real engine needs real couplings')
            out[m] = v.real
    return out


def expr_key(e):
    """Identity of an Expr up to an overall scalar (for de-duplication of labels)."""
    items = sorted(e.items(), key=lambda kv: ta.mono_key(kv[0]))
    c0 = items[0][1]
    k0 = max(c0.items(), key=lambda kv: abs(kv[1]))[1]
    return tuple((m, tuple(sorted((k, np.round(v / k0, 9)) for k, v in c.items()))) for m, c in items)


def open_blocks(p, L, Lmin=1):
    """Open words of length Lmin..L grouped by (start node, end node, edge charges)."""
    blocks = {}
    for l in range(Lmin, L + 1):
        for w in qt.open_words(p, l):
            blocks.setdefault((qt.row_node(w[0]), qt.col_node(w[-1]), qt.word_charges(w)), []).append(w)
    return blocks


def trace_classes(p, L):
    """Distinct (up to scalar) traces of closed walks of length 2..L, grouped by edge charges: {q: [Expr]}."""
    out, seen = {}, set()
    for l in range(2, L + 1):
        for w in qt.closed_words(p, l):
            e = ta.trace(w)
            if not e:
                continue
            key = expr_key(e)
            if key in seen:
                continue
            seen.add(key)
            out.setdefault(qt.word_charges(w), []).append(e)
    return out


def gauss_expr(p, v, W, left=False):
    """Tr[W G_v] (W to the left of G_v in operator order) or, with left=True, Tr[G_v W]; G_v = X_v - n p 1 is the
    u(n)_v generator matrix (quiver_trace.generator_words), W an open word v -> v (() = identity).  Annihilates
    U(n)_v-invariant states from the right (left=False) or, as phi(Tr[G_v W] Y), from the left."""
    W = tuple(W)
    out = Expr()
    for c, gw in qt.generator_words(p, v):
        word = (gw + W) if left else (W + gw)
        for m2, c2 in ta.trace(word).items():
            out.add(m2, c2 * c)
    # - n p Tr[W]   (Tr[1] = n)
    if W:
        for m2, c2 in ta.trace(W).items():
            out.add(m2, c2 * NPoly.N(1, -p))
    else:
        out.add((), NPoly.N(2, -p))
    return out.clean()


# ------------------------------------------------------------------------------------------------ worker state
_G = {}


def _init(C, n, m, opts):
    p = C.shape[0]
    _G.update(n=n, m=m, opts=opts, p=p)
    _G['H'] = qt.hamiltonian(C)
    _G['Ne'] = [qt.edge_number(p, e) + ta.scalar(-m) for e in range(3)]
    _G['Cas'] = [qt.casimir(p, v) for v in range(3)]
    if opts.get('bps'):
        _G['Q'] = qt.supercharge(C)
        _G['Qb'] = ta.dagger(_G['Q'])


def _x_rows(X):
    """Rows generated by a neutral gauge-invariant operator X (an Expr): EOM, sector, Casimir, (BPS: XH, HX)."""
    n, o = _G['n'], _G['opts']
    L = X.total_length()
    rows = []
    if o.get('L_eom_gram') is None or L <= o['L_eom_gram']:
        rows += [ta.mul(Ne, X) for Ne in _G['Ne']]
        rows.append(ta.commutator(_G['H'], X))
    if o.get('casimir') and (o.get('L_cas') is None or L <= o['L_cas']):
        rows += [ta.mul(Cv, X) for Cv in _G['Cas']]
    if o.get('bps') and o.get('bps_H', True) and (o.get('L_bpsH') is None or L <= o['L_bpsH']):
        rows += [ta.mul(X, _G['H']), ta.mul(_G['H'], X)]
    return [ev(r, n) for r in rows]


def _adj_entry(args):
    wa, wb, proj, u, v = args
    n, o = _G['n'], _G['opts']
    wad = ta.word_dagger(wa)
    A0 = ta.trace(wad + wb)
    TT = ta.canonical((wad, wb)) if u == v else None
    A = (A0 - TT * NPoly.N(-1)).clean() if (proj and u == v) else A0.clean()
    rows = _x_rows(A) if A else []
    if o.get('cone_casimir') and (o.get('L_ccas') is None or len(wa) + len(wb) <= o['L_ccas']):
        # (W_b)_ij psi lies in n_u (x) nbar_v for a singlet psi: Chat_x = n^2 - 1 at an endpoint node x (u != v),
        # 2n^2 on the traceless part and 0 on the trace part at x = u = v, and 0 at nodes that are not endpoints
        for x in range(3):
            S = ta.sandwich(wa, wb, _G['Cas'][x])
            if u != v:
                row = S - A0 * NPoly({2: 1.0, 0: -1.0}) if x in (u, v) else S
            else:
                row = S - (A0 - TT * NPoly.N(-1)) * NPoly.N(2, 2.0) if x == u else S
            rows.append(ev(row.clean(), n))
    return ev(A, n), rows


def _sing_entry(args):
    Ea, Eb = args
    n, o = _G['n'], _G['opts']
    S = ta.mul(ta.dagger(Ea), Eb)
    rows = _x_rows(S) if S else []
    if o.get('cone_casimir') and (o.get('L_ccas') is None or S.total_length() <= o['L_ccas']):
        Ead = ta.dagger(Ea)                         # T_b psi is a singlet: phi(T_a^dag Chat_x T_b) = 0
        for x in range(3):
            rows.append(ev(ta.mul(ta.mul(Ead, _G['Cas'][x]), Eb), n))
    return ev(S, n), rows


def _word_rows(X):
    return _x_rows(X)


def _bps_rows(args):
    X, sign = args                      # sign -1: X has charge -(1,1,1), pair with Q; +1: pair with Qbar
    n = _G['n']
    S = _G['Q'] if sign < 0 else _G['Qb']
    return [ev(ta.mul(X, S), n), ev(ta.mul(S, X), n)]


def _gauss_rows(args):
    Y, v, W = args                      # phi(Y Tr[W G_v]) = 0 and phi(Tr[G_v W] Y) = 0 (G_v kills singlets on either side)
    n, p = _G['n'], _G['p']
    return [ev(ta.mul(Y, gauss_expr(p, v, W)), n), ev(ta.mul(gauss_expr(p, v, W, left=True), Y), n)]


def _fin_rows(blocks):
    return ev(qt.finite_n_relation(blocks, _G['n']), _G['n'])


def _dagger(m):
    return m, ev(ta.canonical(tuple(ta.word_dagger(w) for w in reversed(m))), _G['n'])


def rotate_mono(mono, s=1):
    """Image of a monomial under the quiver rotation tau^s: edge e -> e+s (nodes shift alike), flavours fixed."""
    return tuple(tuple((((e + s) % 3, f), t) for (e, f), t in w) for w in mono)


def _rotate(m):
    return m, ev(ta.canonical(rotate_mono(m)), _G['n'])


# ------------------------------------------------------------------------------------------------ the SDP
class QuiverSDP:
    def __init__(self, C, n, m, L_adj=2, L_sing=3, L_eom=None, L_eom_gram='auto', casimir=True, L_cas='auto', gauss=None,
                 bps=False, L_bps=None, bps_H=True, L_bpsH='auto', finite_n_len=0, finite_n_sorted=True,
                 adjoint_projected=True, z3=False, cone_casimir=False, L_ccas='auto', workers=8, verbose=False):
        """C: real couplings (p x p x p); n: rank; m: edge occupation (k = 3m).  L_adj: open-word length of the
        adjoint cones; L_sing: closed-walk length of the singlet cones; L_eom: extra neutral traced words used as EOM
        generators beyond the Gram entries (every neutral closed walk of length <= 2 L_adj is already a Gram entry, so the
        default None adds nothing).  L_eom_gram, L_cas, L_bpsH: EOM + sector, Casimir and phi(XH) rows only for
        generators X of total length <= the cap ('auto': 2 L_adj; None: all).  casimir: Casimir rows on.  gauss: None or (L_Y, L_W): Gauss rows with Y a
        traced closed walk of length <= L_Y (or 1) and W an open v -> v word of length <= L_W.  bps: BPS rows
        (charge -+(1,1,1) operators from traced walks of length <= L_bps and cross-charge Gram products).
        finite_n_len: node-wise antisymmetriser relations of total length <= finite_n_len (sorted block tuples only
        if finite_n_sorted).  z3: impose invariance under the quiver rotation tau (edge e -> e+1), a symmetry of H and
        Q when C_abc = C_bca; the tau-average of a valid functional is valid, so the optimum (or feasibility) is
        unchanged, and only the node-0 adjoint cones, Gauss and finite-n rows are needed.  cone_casimir: for a singlet
        psi each (W_b)_ij psi lies in a definite irrep (n_u x nbar_v, or adj + 1 when u = v), so the sandwiched
        Casimir sum_ij <W_a psi| Chat_x |W_b psi> equals the Casimir value times the Gram entry (D21.4b); likewise
        phi(T_a^dag Chat_x T_b) = 0 for traces; applied to Gram entries of total length <= L_ccas."""
        import multiprocessing as mp
        t0 = time.time()
        C = np.asarray(C, dtype=float)
        p = C.shape[0]
        self.C, self.n, self.m, self.p = C, n, m, p
        self.z3 = z3
        if z3:
            assert np.allclose(C, C.transpose(1, 2, 0)), 'z3 needs cyclic couplings C_abc = C_bca'
        nodes = (0,) if z3 else (0, 1, 2)
        # Row generation is capped by the total length of the generating operator X: the long singlet-cone entries
        # (e.g. Tr[BBB] Tr[PPP], length 6) multiplied by H, Chat_v (length 4) would create length-10 monomials that
        # appear in nothing else and dominate the problem size.  Any subset of valid rows is valid.
        if L_eom_gram == 'auto':
            L_eom_gram = 2 * L_adj
        if L_cas == 'auto':
            L_cas = 2 * L_adj
        if L_bpsH == 'auto':
            L_bpsH = 2 * L_adj
        if L_ccas == 'auto':
            L_ccas = 2 * L_adj
        opts = dict(casimir=casimir, L_cas=L_cas, bps=bps, bps_H=bps_H, L_bpsH=L_bpsH, L_eom_gram=L_eom_gram,
                    cone_casimir=cone_casimir, L_ccas=L_ccas)
        self.level = dict(L_adj=L_adj, L_sing=L_sing, L_eom=L_eom, L_eom_gram=L_eom_gram, casimir=casimir, L_cas=L_cas,
                          cone_casimir=cone_casimir, L_ccas=L_ccas, z3=z3, gauss=gauss, bps=bps,
                          L_bps=L_bps, bps_H=bps_H, L_bpsH=L_bpsH, finite_n_len=finite_n_len, adjoint_projected=adjoint_projected)
        _init(C, n, m, opts)
        self.H = ev(_G['H'], n)
        self.idx = {(): 0}               # monomial -> integer
        self.monos = [()]
        self.cones, rows = [], []
        Pool = mp.get_context('fork').Pool
        with Pool(workers, initializer=_init, initargs=(C, n, m, opts)) as pool:
            # ---- adjoint cones over open words
            blocks = {key: ws for key, ws in open_blocks(p, L_adj).items() if key[0] in nodes}
            for key, ws in sorted(blocks.items()):
                u, v, q = key
                pairs = [(ws[a], ws[b], adjoint_projected, u, v) for a in range(len(ws)) for b in range(a, len(ws))]
                ab = [(a, b) for a in range(len(ws)) for b in range(a, len(ws))]
                ent = {}
                for (a, b), (A, rr) in zip(ab, pool.imap(_adj_entry, pairs, chunksize=16)):
                    if A:
                        ent[(a, b)] = self._intern(A)
                    rows += rr
                self.cones.append(dict(name=f'adj {u}->{v} q={q}', size=len(ws), entries=ent))
            if verbose:
                print(f'    adjoint cones {[c["size"] for c in self.cones]}  [{time.time() - t0:.0f}s]', flush=True)
            # ---- singlet cones over traced closed walks (+ identity in the neutral block)
            sing = trace_classes(p, L_sing)
            for q, es in sorted(sing.items()):
                labels = ([Expr({(): NPoly.const(1)})] if q == (0, 0, 0) else []) + es
                pairs = [(labels[a], labels[b]) for a in range(len(labels)) for b in range(a, len(labels))]
                ab = [(a, b) for a in range(len(labels)) for b in range(a, len(labels))]
                ent = {}
                for (a, b), (S, rr) in zip(ab, pool.imap(_sing_entry, pairs, chunksize=8)):
                    if S:
                        ent[(a, b)] = self._intern(S)
                    rows += rr
                self.cones.append(dict(name=f'sing q={q}', size=len(labels), entries=ent))
            if verbose:
                print(f'    + singlet cones {[c["size"] for c in self.cones if c["name"].startswith("sing")]}  '
                      f'[{time.time() - t0:.0f}s]', flush=True)
            # ---- extra neutral traced words as EOM generators; bare sector / Casimir rows
            extra = trace_classes(p, L_eom).get((0, 0, 0), []) if L_eom else []
            for rr in pool.imap(_word_rows, extra, chunksize=8):
                rows += rr
            rows += [ev(Ne, n) for Ne in _G['Ne']]
            if casimir:
                rows += [ev(Cv, n) for Cv in _G['Cas']]
            self.n_basic_rows = len(rows)
            # ---- Gauss rows
            self.n_gauss_rows = 0
            if gauss is not None:
                LY, LW = gauss
                Ys = [Expr({(): NPoly.const(1)})] + (sum(trace_classes(p, LY).values(), []) if LY >= 2 else [])
                tasks = []
                for v in nodes:
                    Ws = [()] + [w for l in range(2, LW + 1) for w in qt.open_words(p, l, v, v)]
                    for W in Ws:
                        for Y in Ys:
                            qY = mono_charges(next(iter(Y)))
                            if all(a + b == 0 for a, b in zip(qY, qt.word_charges(W))):
                                tasks.append((Y, v, W))
                g = [r for rr in pool.imap(_gauss_rows, tasks, chunksize=16) for r in rr]
                self.n_gauss_rows = len(g)
                rows += g
            # ---- BPS rows
            self.n_bps_rows = 0
            if bps:
                Lb = L_bps if L_bps is not None else max(L_adj, L_sing)
                Xm, Xp = [], []
                for q, es in trace_classes(p, Lb).items():
                    if q == (-1, -1, -1):
                        Xm += es
                    elif q == (1, 1, 1):
                        Xp += es
                for (u, v, q), ws in blocks.items():
                    for (u2, v2, q2), ws2 in blocks.items():
                        if (u2, v2) != (u, v):
                            continue
                        dq = tuple(b - a for a, b in zip(q, q2))
                        if dq not in ((-1, -1, -1), (1, 1, 1)):
                            continue
                        for wa in ws:
                            wad = ta.word_dagger(wa)
                            for wb in ws2:
                                for e in ((ta.trace(wad + wb), ta.canonical((wad, wb))) if u == v else (ta.trace(wad + wb),)):
                                    if e:
                                        (Xm if dq == (-1, -1, -1) else Xp).append(e)
                tasks = [(X, -1) for X in Xm] + [(X, 1) for X in Xp]
                b = [r for rr in pool.imap(_bps_rows, tasks, chunksize=16) for r in rr]
                self.n_bps_rows = len(b)
                rows += b
            # ---- finite-n relations
            self.n_fin_rows = 0
            if finite_n_len:
                tasks = []
                for v in nodes:
                    vv = {l: qt.open_words(p, l, v, v) for l in range(2, finite_n_len + 1)}
                    vv = {l: ws for l, ws in vv.items() if ws}
                    for total in range(2 * (n + 1), finite_n_len + 1):
                        for comp in itertools.product(sorted(vv), repeat=n + 1):
                            if sum(comp) != total:
                                continue
                            for blks in itertools.product(*[vv[c] for c in comp]):
                                if finite_n_sorted and list(blks) != sorted(blks):
                                    continue
                                if tuple(map(sum, zip(*[qt.word_charges(b) for b in blks]))) != (0, 0, 0):
                                    continue
                                tasks.append(blks)
                fr = [r for r in pool.imap(_fin_rows, tasks, chunksize=32) if r]
                self.n_fin_rows = len(fr)
                rows += fr
            # ---- intern rows, then close the variable set under the dagger
            self.rows = []
            seen = set()
            for r in rows:
                if not r:
                    continue
                ir = self._intern(r)
                key = self._row_key(ir)
                if key in seen:
                    continue
                seen.add(key)
                self.rows.append(ir)
            self._intern(self.H)
            for mo in self.monos:
                assert mono_charges(mo) == (0, 0, 0), ta.mono_label(mo)
            self.dag, self.rot = {}, {}
            todo = list(self.monos)
            while todo:
                new = []
                maps = [(_dagger, self.dag)] + ([(_rotate, self.rot)] if z3 else [])
                for fn, store in maps:
                    for mo, md in pool.imap_unordered(fn, todo, chunksize=64):
                        store[self.idx[mo]] = md
                        for m2 in md:
                            if m2 not in self.idx:
                                self.idx[m2] = len(self.monos); self.monos.append(m2); new.append(m2)
                todo = new
        conv = lambda d: {i: (np.array([self.idx[m2] for m2 in md], dtype=np.int64), np.array(list(md.values())))
                          for i, md in d.items()}
        self.dag, self.rot = conv(self.dag), conv(self.rot)
        self.build_time = time.time() - t0
        if verbose:
            print(f'  QuiverSDP (n,p,m)=({n},{p},{m}) L_adj={L_adj} L_sing={L_sing} L_eom={L_eom}: {len(self.monos)} monomials, '
                  f'cones {sorted((c["size"] for c in self.cones), reverse=True)[:8]}..., {len(self.rows)} rows '
                  f'(basic {self.n_basic_rows}, gauss {self.n_gauss_rows}, bps {self.n_bps_rows}, finite-n {self.n_fin_rows}; '
                  f'before dedupe), build {self.build_time:.0f}s', flush=True)

    # ------------------------------------------------------------------------------------------------------------
    def _intern(self, d):
        ix = np.empty(len(d), dtype=np.int64)
        for t, mo in enumerate(d):
            i = self.idx.get(mo)
            if i is None:
                i = self.idx[mo] = len(self.monos); self.monos.append(mo)
            ix[t] = i
        return ix, np.fromiter(d.values(), dtype=float, count=len(d))

    @staticmethod
    def _row_key(ir):
        ix, v = ir
        o = np.argsort(ix)
        ix, v = ix[o], v[o]
        return tuple(ix.tolist()), tuple(np.round(v / v[0], 9).tolist())

    # ------------------------------------------------------------------------------------------------------------
    def _reduction(self):
        """Variable reduction from reality: phi(m) = sum_j d_j phi(m_j) with m^dag = sum_j d_j m_j.  Single-term
        daggers with |d| = 1 identify variables (x_m = d x_m'; d = -1 with m' = m forces x_m = 0); multi-term
        daggers become equality rows.  Returns (col, sign) per monomial (col -1: identically zero), the number of
        reduced variables, and the extra rows."""
        nm = len(self.monos)
        parent = list(range(nm)); sgn = np.ones(nm)           # x_i = sgn[i] * x_parent[i] (union-find with signs)
        zero = np.zeros(nm, dtype=bool)

        def find(i):
            s = 1.0
            while parent[i] != i:
                s *= sgn[i]; i = parent[i]
            return i, s
        extra = []
        rels = [(i, *self.dag[i]) for i in range(nm)] + [(i, *self.rot[i]) for i in range(nm) if i in self.rot]
        for i, ix, v in rels:
            if len(ix) == 1 and abs(abs(v[0]) - 1) < 1e-12:
                ri, si = find(i); rj, sj = find(int(ix[0]))
                # x_i = v x_j  ->  si x_ri = v sj x_rj
                if ri == rj:
                    if abs(si - v[0] * sj) > 1e-12:
                        zero[ri] = True
                else:
                    a, b = (ri, rj) if ri > rj else (rj, ri)
                    # x_a = s x_b with s from si x_ri = v sj x_rj
                    s = (v[0] * sj / si) if a == ri else (si / (v[0] * sj))
                    parent[a] = b; sgn[a] = s
                    if zero[a]:
                        zero[b] = True
            else:
                extra.append((np.append(i, ix), np.append(1.0, -v)))
        col = -np.ones(nm, dtype=np.int64); sign = np.zeros(nm)
        roots = {}
        for i in range(nm):
            r, s = find(i)
            if zero[r]:
                continue
            if r not in roots:
                roots[r] = len(roots)
            col[i] = roots[r]; sign[i] = s
        return col, sign, len(roots), extra

    def assemble(self, mode='energy', solver='clarabel', prune=True, max_dense=6e7):
        """Real standard form  min c.x  s.t.  A x + s = b,  s in {0}^n_eq x PSD x ...  over the reduced variables.
        mode: 'energy' (c = H), 'margin' (extra variable t: cones - t 1 >= 0, maximise t, t <= 1 by a 1x1 cone),
        'feasibility' (c = 0)."""
        col, sign, nv, extra = self._reduction()
        self._col, self._sign, self._nv = col, sign, nv
        R, Cc, V, b = [], [], [], []
        r = 0

        def add_row(ix, v, rhs=0.0):
            nonlocal r
            c = col[ix]; ok = c >= 0
            if not ok.any():
                if abs(rhs) > 1e-12:
                    raise RuntimeError('inconsistent row')
                return
            vals = v[ok] * sign[ix[ok]]
            u, inv = np.unique(c[ok], return_inverse=True)
            s = np.zeros(len(u)); np.add.at(s, inv, vals)
            keep = np.abs(s) > 1e-10 * max(1.0, np.abs(s).max())
            if not keep.any():
                if abs(rhs) > 1e-12:
                    raise RuntimeError('inconsistent row')
                return
            R.extend([r] * int(keep.sum())); Cc.extend(u[keep].tolist()); V.extend(s[keep].tolist()); b.append(rhs); r += 1
        add_row(np.array([0]), np.array([1.0]), 1.0)                 # phi(1) = 1
        for ix, v in extra:
            add_row(ix, v)
        for ix, v in self.rows:
            add_row(ix, v)
        ncol = nv + (1 if mode == 'margin' else 0)
        A_eq = sp.csr_matrix((V, (R, Cc)), shape=(r, ncol)); b_eq = np.array(b)
        rn = np.sqrt(np.asarray(A_eq.multiply(A_eq).sum(axis=1)).ravel()); rn[rn == 0] = 1
        A_eq = (sp.diags(1 / rn) @ A_eq).tocsr(); b_eq = b_eq / rn
        # PSD cones
        R, Cc, V = [], [], []
        r = 0
        dims = []
        for cone in self.cones:
            d = cone['size']; dims.append(d)
            order = ([(i, j) for j in range(d) for i in range(j + 1)] if solver == 'clarabel'
                     else [(j, i) for j in range(d) for i in range(j, d)])
            for (a, bb) in order:
                e = cone['entries'].get((a, bb))
                scale = 1.0 if a == bb else np.sqrt(2.0)
                if e is not None:
                    ix, v = e
                    c = col[ix]; ok = c >= 0
                    if ok.any():
                        R.extend([r] * int(ok.sum())); Cc.extend(c[ok].tolist()); V.extend((-scale * v[ok] * sign[ix[ok]]).tolist())
                if mode == 'margin' and a == bb:
                    R.append(r); Cc.append(nv); V.append(1.0)        # s = svec(M) - t svec(1)
                r += 1
        A_psd = sp.csr_matrix((V, (R, Cc)), shape=(r, ncol))
        b_psd = np.zeros(r)
        if mode == 'margin':
            A_psd = sp.vstack([A_psd, sp.csr_matrix(([1.0], ([0], [nv])), shape=(1, ncol))])
            b_psd = np.append(b_psd, 1.0); dims = dims + [1]                       # 1 - t >= 0
        c = np.zeros(ncol)
        if mode == 'energy':
            for mo, cf in self.H.items():
                i = self.idx[mo]
                if col[i] >= 0:
                    c[col[i]] += cf * sign[i]
        elif mode == 'margin':
            c[-1] = -1.0
        n_eq0, ncol0 = A_eq.shape[0], ncol
        A_eq, b_eq, A_psd, c, ncol = self._presolve(A_eq, b_eq, A_psd.tocsr(), c)
        pruned = False
        if prune and A_eq.shape[0] * A_eq.shape[1] <= max_dense:
            from scipy.linalg import qr
            _, Rq, piv = qr(A_eq.toarray().T, mode='economic', pivoting=True)
            d = np.abs(np.diag(Rq)); rk = int((d > 1e-10 * d[0]).sum())
            keep = np.sort(piv[:rk]); A_eq = A_eq[keep]; b_eq = b_eq[keep]; pruned = True
        n_eq = A_eq.shape[0]
        A = sp.vstack([A_eq, A_psd]).tocsc(); bvec = np.concatenate([b_eq, b_psd])
        return dict(A=A, b=bvec, c=c, n_eq=n_eq, dims=dims, n=ncol, pruned=pruned, n_eq_before_presolve=n_eq0,
                    n_before_presolve=ncol0)

    @staticmethod
    def _presolve(A_eq, b_eq, A_psd, c):
        """Exact reduction: a variable that occurs in exactly one equality row and in no cone and not in the objective
        can always be chosen to satisfy that row, so the row and the variable are dropped (repeated until stable);
        variables left in no row, cone or objective are dropped too.  Feasible set (projected on the kept variables)
        and optimal value are unchanged; a Farkas certificate of the reduced problem certifies the original."""
        A_eq = A_eq.tocsr(); Acsc = A_eq.tocsc()
        nr, nc = A_eq.shape
        protected = np.zeros(nc, dtype=bool)
        protected[np.unique(A_psd.indices)] = True
        protected[np.abs(c) > 0] = True
        cnt = np.diff(Acsc.indptr).astype(np.int64)
        active = np.ones(nr, dtype=bool)
        stack = [j for j in np.where((cnt == 1) & ~protected)[0]]
        while stack:
            j = stack.pop()
            if cnt[j] != 1 or protected[j]:
                continue
            rows_j = Acsc.indices[Acsc.indptr[j]:Acsc.indptr[j + 1]]
            rr = [i for i in rows_j if active[i]]
            if len(rr) != 1:
                continue
            i = rr[0]
            active[i] = False
            for jj in A_eq.indices[A_eq.indptr[i]:A_eq.indptr[i + 1]]:
                cnt[jj] -= 1
                if cnt[jj] == 1 and not protected[jj]:
                    stack.append(jj)
        keep_c = (cnt > 0) | protected
        A_eq = A_eq[active][:, keep_c]; b_eq = b_eq[active]
        A_psd = A_psd[:, keep_c]; c = c[keep_c]
        # rows emptied of all kept variables (cannot happen for active rows, but guard)
        rn = np.diff(A_eq.tocsr().indptr) > 0
        assert np.all(np.abs(b_eq[~rn]) < 1e-12)
        return A_eq[rn], b_eq[rn], A_psd, c, int(keep_c.sum())

    def _solve(self, data, solver, eps, max_iters, verbose):
        A, b, c = data['A'], data['b'], data['c']
        t0 = time.time()
        if solver == 'clarabel':
            import clarabel
            P = sp.csc_matrix((data['n'], data['n']))
            cones = [clarabel.ZeroConeT(data['n_eq'])] + [clarabel.PSDTriangleConeT(d) for d in data['dims']]
            st = clarabel.DefaultSettings(); st.verbose = verbose
            st.chordal_decomposition_enable = False
            st.tol_gap_abs = eps; st.tol_gap_rel = eps; st.tol_feas = eps; st.max_iter = 400
            st.static_regularization_constant = 1e-7 if data['pruned'] else 1e-5
            sol = clarabel.DefaultSolver(P, c, A, b, cones, st).solve()
            out = dict(status=str(sol.status), x=np.array(sol.x), y=np.array(sol.z), obj=float(sol.obj_val))
        else:
            import scs, platform
            lin = {'linear_solver': 'accelerate'} if platform.system() == 'Darwin' else {}
            sol = scs.SCS(dict(A=A, b=b, c=c), dict(z=data['n_eq'], s=data['dims']), eps_abs=eps, eps_rel=eps,
                          eps_infeas=eps, max_iters=max_iters, verbose=verbose, rho_x=1e-3, **lin).solve()
            out = dict(status=sol['info']['status'], x=np.array(sol['x']), y=np.array(sol['y']),
                       obj=float(sol['info']['pobj']))
        out['solve_time'] = time.time() - t0
        out.update(n=data['n'], n_eq=data['n_eq'], n_rows=A.shape[0], nnz=A.nnz, dims=data['dims'], pruned=data['pruned'])
        return out

    def energy(self, solver='clarabel', eps=1e-8, max_iters=100000, verbose=False, prune=None):
        """min phi(H): a lower bound on the lowest singlet energy at degree 3m (when the SDP is solved)."""
        data = self.assemble('energy', solver, prune=(solver == 'clarabel') if prune is None else prune)
        return self._solve(data, solver, eps, max_iters, verbose)

    def margin(self, solver='clarabel', eps=1e-8, max_iters=100000, verbose=False, prune=None):
        data = self.assemble('margin', solver, prune=(solver == 'clarabel') if prune is None else prune)
        out = self._solve(data, solver, eps, max_iters, verbose)
        out['margin'] = float(out['x'][-1]) if len(out['x']) else float('nan')
        return out

    def feasibility(self, solver='clarabel', eps=1e-8, max_iters=100000, verbose=False, prune=None):
        """Pure feasibility; for an infeasible problem the dual vector y is a Farkas certificate
        (A^T y = 0, y in K*, b.y < 0), checked here."""
        data = self.assemble('feasibility', solver, prune=(solver == 'clarabel') if prune is None else prune)
        out = self._solve(data, solver, eps, max_iters, verbose)
        y = out['y']
        if y is not None and len(y) == data['A'].shape[0]:
            out.update(self.check_certificate(data, y, solver))
        return out

    @staticmethod
    def check_certificate(data, y, solver):
        A, b = data['A'], data['b']
        ny = np.linalg.norm(y) or 1.0
        res = dict(cert_ATy=float(np.linalg.norm(A.T @ y) / ny), cert_by=float(b @ y / ny))
        off = data['n_eq']; worst = 0.0
        for d in data['dims']:
            L = d * (d + 1) // 2
            blk = y[off:off + L]; off += L
            M = np.zeros((d, d))
            pairs = ([(i, j) for j in range(d) for i in range(j + 1)] if solver == 'clarabel'
                     else [(i, j) for j in range(d) for i in range(j, d)])
            for t, (i, j) in enumerate(pairs):
                v = blk[t] / (1.0 if i == j else np.sqrt(2)); M[i, j] = v; M[j, i] = v
            worst = min(worst, float(np.linalg.eigvalsh(M)[0]) / ny)
        res['cert_min_eig'] = worst
        return res

    # ------------------------------------------------------------------------------------------------------------
    def check_functional(self, phi):
        """Evaluate every constraint on a functional phi = {monomial: value} (e.g. from an exact state): returns the
        worst row residual (relative to the row's coefficient scale x the typical |phi|), the worst reality
        residual and the smallest eigenvalue of each cone (relative to its largest)."""
        x = np.array([phi[mo] for mo in self.monos])
        scale = max(1.0, np.abs(x).max())
        worst_row = 0.0
        for ix, v in self.rows:
            worst_row = max(worst_row, abs(v @ x[ix]) / (np.abs(v).max() * scale))
        worst_real = 0.0
        for i, (ix, v) in self.dag.items():
            worst_real = max(worst_real, abs(x[i] - v @ x[ix]) / scale)
        eig = []
        for cone in self.cones:
            d = cone['size']; M = np.zeros((d, d))
            for (a, b), (ix, v) in cone['entries'].items():
                M[a, b] = M[b, a] = v @ x[ix]
            w = np.linalg.eigvalsh(M)
            eig.append(w[0] / max(1.0, abs(w[-1])))
        return dict(worst_row=worst_row, worst_reality=worst_real, min_cone_eig=min(eig), energy=sum(cf * phi[mo] for mo, cf in self.H.items()))
