"""
Exact finite-n trace-word algebra for the U(n)^3 fermionic quiver (research/notes/quiver_project.md section 12;
docs/derivations.md D21).  Built on src/trace_algebra.py, whose swap/contraction machinery is node-agnostic.

Letters.  ((e, f), 'P') = the creation fermion of edge e (0: A, 1: B, 2: C) and flavour f, a matrix with row index
at node e and column index at node e+1 (mod 3); ((e, f), 'B') = its conjugate, B_{kl} = (P_{lk})^dagger, with row
index at node e+1 and column index at node e.  The only non-zero anticommutator {P_ij, B_kl} = delta_il delta_jk
contracts indices of the same node, so trace_algebra's contraction rule applies verbatim, and a closed index loop
gives a factor n (equal ranks).  Gauge-invariant words are closed walks: consecutive letters must match
column-node -> row-node.

Contents: words and node bookkeeping; Q = sum C_abc Tr(A^a B^b C^c); H = {Q, Q^+} by exact multiplication; the
fermion number N_tot = sum Tr(B P) (= sum P^+ P); per-node gauge generators and quadratic Casimirs (zero exactly on
singlets); finite-n relations per node (antisymmetriser over n+1 indices of one node, blocks = open words from that
node back to itself); an explicit Fock-space evaluator for verification.
"""
import itertools
import math
import numpy as np
import scipy.sparse as sp
import trace_algebra as ta
from trace_algebra import NPoly, Expr


# ------------------------------------------------------------------------------------------------ letters and words
def P(e, f):
    return ((e, f), 'P')


def B(e, f):
    return ((e, f), 'B')


def row_node(letter):
    (e, f), t = letter
    return e if t == 'P' else (e + 1) % 3


def col_node(letter):
    (e, f), t = letter
    return (e + 1) % 3 if t == 'P' else e


def letters(p):
    return [P(e, f) for e in range(3) for f in range(p)] + [B(e, f) for e in range(3) for f in range(p)]


def is_walk(word):
    return all(col_node(word[i]) == row_node(word[i + 1]) for i in range(len(word) - 1))


def is_closed(word):
    return len(word) > 0 and is_walk(word) and col_node(word[-1]) == row_node(word[0])


def open_words(p, L, u=None, v=None):
    """Walks of length L (as tuples of letters) from node u to node v (None: any)."""
    out = []

    def rec(w):
        if len(w) == L:
            if v is None or col_node(w[-1]) == v:
                out.append(tuple(w))
            return
        for l in letters(p):
            if not w:
                if u is not None and row_node(l) != u:
                    continue
            elif row_node(l) != col_node(w[-1]):
                continue
            w.append(l)
            rec(w)
            w.pop()
    rec([])
    return out


def closed_words(p, L):
    return [w for w in open_words(p, L) if col_node(w[-1]) == row_node(w[0])]


def word_charges(word):
    """(net creation numbers on edges 0, 1, 2) of a word."""
    q = [0, 0, 0]
    for (e, f), t in word:
        q[e] += 1 if t == 'P' else -1
    return tuple(q)


# ------------------------------------------------------------------------------------------------ operators
def supercharge(C):
    p = C.shape[0]
    out = Expr()
    for a, b, c in itertools.product(range(p), repeat=3):
        if abs(C[a, b, c]) > 1e-15:
            out = out + ta.trace((P(0, a), P(1, b), P(2, c))) * NPoly.const(C[a, b, c])
    return out


def hamiltonian(C):
    Q = supercharge(C)
    Qd = ta.dagger(Q)
    return (ta.mul(Q, Qd) + ta.mul(Qd, Q)).clean()


def number_operator(p):
    """N_tot = sum_{e,f} Tr(P B) = sum_{e,f,ij} a^+_ij a_ij (total fermion number; P is a creation operator)."""
    out = Expr()
    for e in range(3):
        for f in range(p):
            out = out + ta.trace((P(e, f), B(e, f)))
    return out.clean()


def edge_number(p, e):
    out = Expr()
    for f in range(p):
        out = out + ta.trace((P(e, f), B(e, f)))
    return out.clean()


def generator_words(p, v):
    """The gauge generator matrix of node v, G_v = X_v - n p 1 with X_v = sum_f [(P_v B_v) + (B_{v-1} P_{v-1})], as a
    list of (coefficient, open word from v to v).  Out-edge fields a_{xy} (x at node v): E_ij = sum_y a^+_iy a_jy =
    (P B)_ij.  In-edge fields c_{lx} (x at node v, antifundamental): E_ij = -sum_l c^+_lj c_li = (B P)_ij - n delta_ij.
    Check: Tr G_v = N_out - N_in, the U(1) charge of node v.  The constant drops out of the traceless Casimir."""
    out = []
    for f in range(p):
        out.append((1.0, (P(v, f), B(v, f))))
        e_in = (v - 1) % 3
        out.append((1.0, (B(e_in, f), P(e_in, f))))
    return out


def casimir(p, v):
    """n * (SU(n) quadratic Casimir of node v, unnormalised) = n Tr(X_v^2) - (Tr X_v)^2, X_v as in generator_words
    (the shift by -n p 1 cancels in this traceless combination).  Polynomial in n; zero exactly on node-v SU(n)
    singlets and positive otherwise."""
    gw = generator_words(p, v)
    e = Expr()
    for c1, w1 in gw:
        for c2, w2 in gw:
            for m, c in ta.trace(w1 + w2).items():
                e.add(m, c * c1 * c2 * NPoly.N(1))
    tr = Expr()
    for c1, w1 in gw:
        for m, c in ta.trace(w1).items():
            tr.add(m, c * c1)
    sq = ta.mul(tr, tr)
    out = Expr(e)
    for m, c in sq.items():
        out.add(m, c * -1.0)
    return out.clean()


def finite_n_relation(blocks, n):
    """Antisymmetriser over n+1 indices of ONE node: blocks must all be open words from that node back to itself."""
    nodes = {(row_node(b[0]), col_node(b[-1])) for b in blocks}
    assert len(nodes) == 1 and next(iter(nodes))[0] == next(iter(nodes))[1], 'blocks must be v -> v words'
    return ta.finiteN_relations([tuple(b) for b in blocks], n)


# ------------------------------------------------------------------------------------------------ explicit evaluator
class ExplicitQuiver:
    """Explicit Fock-space operators for small (n, p): modes (e, f, x, y), creation P^{(e,f)}_{xy}."""

    def __init__(self, n, p, kmax=None):
        self.n, self.p = n, p
        self.nm = 3 * p * n * n
        self.kmax = kmax if kmax is not None else self.nm
        states = [s for k in range(self.kmax + 1) for s in itertools.combinations(range(self.nm), k)]
        self.states = states
        self.index = {s: i for i, s in enumerate(states)}
        self.dim = len(states)
        self._cre = {}

    def mode(self, e, f, x, y):
        return ((e * self.p + f) * self.n + x) * self.n + y

    def creation(self, m):
        if m not in self._cre:
            rows, cols, vals = [], [], []
            for j, s in enumerate(self.states):
                if m in s or len(s) == self.kmax:
                    continue
                t = tuple(sorted(s + (m,)))
                sign = (-1) ** sum(1 for x in s if x < m)
                rows.append(self.index[t]); cols.append(j); vals.append(sign)
            self._cre[m] = sp.csr_matrix((vals, (rows, cols)), shape=(self.dim, self.dim))
        return self._cre[m]

    def letter_op(self, letter, row, col):
        (e, f), t = letter
        if t == 'P':
            return self.creation(self.mode(e, f, row, col))
        return self.creation(self.mode(e, f, col, row)).T                # B_{row,col} = (P_{col,row})^dagger

    def mono_apply(self, mono, v):
        """Apply Tr[w_1] ... Tr[w_m] (operator order) to v by summing over index assignments."""
        w = v
        for word in reversed(mono):
            if len(word) == 0:
                w = self.n * w
                continue
            acc = np.zeros_like(w)
            L = len(word)
            for idx in itertools.product(range(self.n), repeat=L):      # idx[i] = column index of letter i
                u = w
                for i in reversed(range(L)):
                    u = self.letter_op(word[i], idx[i - 1], idx[i]) @ u   # row of letter i = column of letter i-1
                acc = acc + u
            w = acc
        return w

    def expr_apply(self, expr, v):
        out = np.zeros_like(v, dtype=complex)
        for m, c in expr.items():
            cv = c.at(self.n)
            if abs(cv) > 1e-15:
                out = out + cv * self.mono_apply(m, v)
        return out

    def degree(self):
        return np.array([len(s) for s in self.states])
