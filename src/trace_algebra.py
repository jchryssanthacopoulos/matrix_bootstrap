"""
Exact finite-N algebra of trace words for fermionic U(N) matrix models (docs/trace_bootstrap_plan.md, M1).

Letters: (a, 'P') = Psi^a, (a, 'B') = Psibar^a, all Grassmann-odd, with the only non-vanishing anticommutator
    {Psi^a_ij, Psibar^b_kl} = delta^ab delta_il delta_jk
(this holds in the flavor basis and in the Z_p Fourier basis alike).  A *routed term* is a list of letters in
OPERATOR order together with a routing permutation nxt: the column index of letter i is contracted with the row
index of letter nxt[i].  A trace word Tr[l_1...l_L] is the routed term with nxt[i] = i+1 mod L; a monomial
Tr[w_1]...Tr[w_m] (operator order = listed order, letters within each trace in listed order) is a routed term whose
cycles are contiguous blocks.  The single primitive is the swap of two adjacent letters in operator order,
    x y = - y x + {x, y},
where the anticommutator contracts two index pairs and removes the two letters (rewiring the routing; a closed
index loop with no letters gives a factor N).  Everything else -- canonical cyclic rotation of a trace, reordering
of trace factors, products, daggers, commutators, finite-N (antisymmetriser) relations -- is a composition of swaps.

Coefficients are polynomials in N with complex coefficients (NPoly: dict power -> complex); expressions (Expr) are
dicts monomial -> NPoly, monomial = tuple of trace words, each a tuple of letters, in canonical form.

Every identity produced here is checked numerically against explicit operators in tests/test_trace_algebra.py
(N=2, p=3 and N=2,3, p=1), see the plan, section 3.
"""
import itertools
from collections import defaultdict
import numpy as np

# ----------------------------------------------------------------------------------------------------------------
# letters
# ----------------------------------------------------------------------------------------------------------------

def conj_letter(l):
    a, t = l
    return (a, 'B' if t == 'P' else 'P')


def letter_q(l):
    return 1 if l[1] == 'P' else -1


def word_q(word):
    return sum(letter_q(l) for l in word)


def letter_z(l, p):
    """Z_p charge of a Fourier-basis letter (a,'P') -> +a, (a,'B') -> -a (mod p)."""
    return (l[0] if l[1] == 'P' else -l[0]) % p


def word_z(word, p):
    return sum(letter_z(l, p) for l in word) % p


def word_dagger(word):
    """(Tr w)^dagger = Tr[conj(l_L) ... conj(l_1)]."""
    return tuple(conj_letter(l) for l in reversed(word))


# ----------------------------------------------------------------------------------------------------------------
# polynomials in N and expressions
# ----------------------------------------------------------------------------------------------------------------

class NPoly(dict):
    """Polynomial in N: {power: complex coefficient}."""
    @staticmethod
    def const(c):
        return NPoly({0: complex(c)}) if c != 0 else NPoly()

    @staticmethod
    def N(power=1, c=1.0):
        return NPoly({power: complex(c)})

    def __add__(self, other):
        out = NPoly(self)
        for k, v in other.items():
            out[k] = out.get(k, 0) + v
        return out.clean()

    def __mul__(self, other):
        if not isinstance(other, NPoly):
            return NPoly({k: v * other for k, v in self.items()}).clean()
        out = NPoly()
        for k1, v1 in self.items():
            for k2, v2 in other.items():
                out[k1 + k2] = out.get(k1 + k2, 0) + v1 * v2
        return out.clean()

    __rmul__ = __mul__

    def clean(self, tol=1e-14):
        for k in [k for k, v in self.items() if abs(v) <= tol]:
            del self[k]
        return self

    def at(self, N):
        return sum(v * N ** k for k, v in self.items()) if self else 0j

    def conj(self):
        return NPoly({k: np.conj(v) for k, v in self.items()})

    def __repr__(self):
        if not self:
            return '0'
        return ' + '.join(f"({v:.6g})N^{k}" if k else f"({v:.6g})" for k, v in sorted(self.items()))


class Expr(dict):
    """Linear combination of canonical monomials: {monomial: NPoly}."""
    def add(self, mono, coeff):
        if not coeff:
            return
        c = self.get(mono)
        self[mono] = (c + coeff) if c is not None else NPoly(coeff)
        if not self[mono]:
            del self[mono]

    def __add__(self, other):
        out = Expr(self)
        for m, c in other.items():
            out.add(m, c)
        return out

    def __sub__(self, other):
        return self + other * (-1.0)

    def __mul__(self, scalar):
        if isinstance(scalar, Expr):
            return mul(self, scalar)
        if isinstance(scalar, NPoly):
            return Expr({m: c * scalar for m, c in self.items()}).clean()
        return Expr({m: c * scalar for m, c in self.items()}).clean()

    __rmul__ = __mul__

    def clean(self):
        for m in [m for m, c in self.items() if not c]:
            del self[m]
        return self

    def at(self, N):
        """Evaluate coefficients at numerical N: {monomial: complex}."""
        out = {m: c.at(N) for m, c in self.items()}
        return {m: v for m, v in out.items() if abs(v) > 1e-14}

    def total_length(self):
        return max((sum(len(w) for w in m) for m in self), default=0)

    def __repr__(self):
        return ' + '.join(f"[{c}] {mono_label(m)}" for m, c in sorted(self.items(), key=lambda kv: mono_key(kv[0])))


def mono_key(mono):
    return (sum(len(w) for w in mono), len(mono), mono)


def mono_label(mono):
    if not mono:
        return '1'
    lab = lambda l: ('P' if l[1] == 'P' else 'B') + str(l[0])
    return ' '.join('Tr[' + ' '.join(lab(l) for l in w) + ']' for w in mono)


def mono_q(mono):
    return sum(word_q(w) for w in mono)


def mono_z(mono, p):
    return sum(word_z(w, p) for w in mono) % p


def mono_parity(mono):
    return sum(len(w) for w in mono) % 2


# ----------------------------------------------------------------------------------------------------------------
# routed terms and the canonicalisation engine
# ----------------------------------------------------------------------------------------------------------------

def canonical_rotation_start(word):
    """Index j such that word[j:]+word[:j] is the lexicographically smallest rotation (ties -> smallest j)."""
    L = len(word)
    best, bj = None, 0
    for j in range(L):
        rot = word[j:] + word[:j]
        if best is None or rot < best:
            best, bj = rot, j
    return bj


def _cycles(nxt):
    seen = [False] * len(nxt); cyc = []
    for s in range(len(nxt)):
        if seen[s]:
            continue
        c = []; i = s
        while not seen[i]:
            seen[i] = True; c.append(i); i = nxt[i]
        cyc.append(c)
    return cyc


def _target_order(ops, nxt):
    """Canonical operator order for a routing: each cycle in matrix order starting at its canonical rotation,
    cycles sorted by their canonical word (stable in the original minimal position)."""
    cycles = []
    for c in _cycles(nxt):
        word = tuple(ops[i] for i in c)
        j = canonical_rotation_start(word)
        cpos = c[j:] + c[:j]
        cycles.append((tuple(ops[i] for i in cpos), min(c), cpos))
    cycles.sort(key=lambda t: (t[0], t[1]))
    order = [i for _, _, cpos in cycles for i in cpos]
    mono = tuple(w for w, _, _ in cycles)
    return order, mono


def _contract(ops, nxt, order, px_pos, py_pos):
    """Anticommutator branch: remove the letters x, y at operator positions px_pos, py_pos (adjacent in `order`)
    and rewire the routing.  Index edges: e_i = column index of letter i = row index of nxt[i].  The anticommutator
    {x_{r_x c_x}, y_{r_y c_y}} = delta_{r_x c_y} delta_{c_x r_y} identifies e_{prev[x]} with e_y and e_x with
    e_{prev[y]} (union-find on edges); the surviving letters are re-linked through the merged classes and every
    class containing no surviving letter is a closed index loop, i.e. a factor N.
    Returns (new_ops, new_nxt, new_order, power_of_N)."""
    x, y = order[px_pos], order[py_pos]
    L = len(ops)
    prev = [0] * L
    for i, j in enumerate(nxt):
        prev[j] = i
    parent = list(range(L))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    def union(i, j):
        parent[find(i)] = find(j)
    union(prev[x], y); union(x, prev[y])
    keep = [i for i in range(L) if i != x and i != y]
    row_class = {find(prev[j]): j for j in keep}          # class of the row edge of each surviving letter
    remap = {old: new for new, old in enumerate(keep)}
    n_ops = [ops[i] for i in keep]
    n_nxt = [remap[row_class[find(i)]] for i in keep]     # successor = letter whose row class is my column class
    used = {find(i) for i in keep}
    pw = len({find(i) for i in range(L)}) - len(used)      # classes with no surviving letter: closed loops
    n_order = [remap[i] for i in order if i != x and i != y]
    return n_ops, n_nxt, n_order, pw


def _period(word):
    L = len(word)
    for d in range(1, L):
        if L % d == 0 and word[d:] + word[:d] == word:
            return d
    return L


def _self_symmetry(mono):
    """A rearrangement of the operator order of the canonical monomial that reproduces the same letter sequence
    with an odd fermionic sign, if one exists: ('rot', f, d) -- factor f has period d with (L-1)d odd, so
    T = -T + lower; ('swap', f) -- factors f and f+1 are equal words of odd length, so T T = -T T + merging.
    Returns None if the monomial has no such symmetry."""
    for f, w in enumerate(mono):
        L = len(w); d = _period(w)
        if d < L and ((L - 1) * d) % 2 == 1:
            return ('rot', f, d)
    for f in range(len(mono) - 1):
        if mono[f] == mono[f + 1] and len(mono[f]) % 2 == 1:
            return ('swap', f)
    return None


def _add_canonical(mono, coeff, out):
    """Add coeff * mono to out, first resolving self-symmetric monomials: T = -T + (lower)  =>  T = (lower)/2."""
    sym = _self_symmetry(mono)
    if sym is None:
        out.add(mono, coeff)
        return
    ops, nxt = routed_from_monomial(mono)
    starts = [sum(len(w) for w in mono[:f]) for f in range(len(mono))]
    order = list(range(len(ops)))
    if sym[0] == 'rot':
        _, f, d = sym; s0, L = starts[f], len(mono[f])
        blk = order[s0:s0 + L]; order[s0:s0 + L] = blk[d:] + blk[:d]
    else:
        _, f = sym; s0, L = starts[f], len(mono[f])
        order[s0:s0 + L], order[s0 + L:s0 + 2 * L] = order[s0 + L:s0 + 2 * L], order[s0:s0 + L]
    # the rearranged order is the same operator; reducing it gives  -mono + lower  (the sign is -1 by construction)
    tmp = Expr()
    _reduce(ops, nxt, order, NPoly.const(1), tmp, resolve=False)
    main = tmp.pop(mono, None)
    assert main is not None and abs(main.at(1) + 1) < 1e-12, (mono, sym, main)
    for m2, c2 in tmp.items():                # lower terms, already canonical (branches resolve recursively)
        for m3, c3 in canonical(m2, c2 * coeff * 0.5).items():
            out.add(m3, c3)


def _reduce(ops, nxt, order, coeff, out, resolve=True):
    """Bubble-sort the operator order `order` (a permutation of range(len(ops))) into the canonical order for this
    routing by adjacent swaps, accumulating the resulting canonical monomials in `out`.  Each swap of a conjugate
    pair spawns an anticommutator branch (two letters fewer), reduced recursively.  resolve=True applies the
    self-symmetry reduction of _add_canonical to the final monomial."""
    target, mono = _target_order(ops, nxt)
    order = list(order); sign = 1
    L = len(order)
    for i in range(L):
        if order[i] == target[i]:
            continue
        j = order.index(target[i])
        while j > i:
            x, y = order[j - 1], order[j]
            if ops[x] == conj_letter(ops[y]):
                n_ops, n_nxt, n_order, pw = _contract(ops, nxt, order, j - 1, j)
                _reduce(n_ops, n_nxt, n_order, coeff * NPoly.N(pw, sign), out)
            order[j - 1], order[j] = y, x
            sign = -sign
            j -= 1
    if resolve:
        _add_canonical(mono, coeff * sign, out)
    else:
        out.add(mono, coeff * sign)


def canonicalize_routed(ops, nxt, coeff=None):
    out = Expr()
    _reduce(list(ops), list(nxt), list(range(len(ops))), coeff if coeff is not None else NPoly.const(1), out)
    return out


def routed_from_monomial(mono):
    """Routed term for an ordered product of trace words (operator order = listed order)."""
    ops, nxt, pos = [], [], 0
    for w in mono:
        L = len(w)
        ops += list(w)
        nxt += [pos + (i + 1) % L for i in range(L)]
        pos += L
    return ops, nxt


def canonical(mono, coeff=None):
    """Canonical Expr of the (not necessarily canonical) monomial mono."""
    ops, nxt = routed_from_monomial(mono)
    return canonicalize_routed(ops, nxt, coeff)


def canonical_expr(expr):
    out = Expr()
    for m, c in expr.items():
        for m2, c2 in canonical(m, c).items():
            out.add(m2, c2)
    return out


def trace(word):
    """Expr of Tr[word]."""
    return canonical((tuple(word),))


def mul(e1, e2):
    out = Expr()
    for m1, c1 in e1.items():
        for m2, c2 in e2.items():
            for m, c in canonical(m1 + m2, c1 * c2).items():
                out.add(m, c)
    return out


def dagger(expr):
    out = Expr()
    for m, c in expr.items():
        md = tuple(word_dagger(w) for w in reversed(m))
        for m2, c2 in canonical(md, c.conj()).items():
            out.add(m2, c2)
    return out


def commutator(e1, e2):
    return mul(e1, e2) - mul(e2, e1)


def scalar(c):
    return Expr({(): NPoly.const(c)}) if c != 0 else Expr()


# ----------------------------------------------------------------------------------------------------------------
# model operators as expressions
# ----------------------------------------------------------------------------------------------------------------

def fourier_C(C):
    """Couplings in the Z_p Fourier letter basis: Q = C_abc Tr[Psi^a Psi^b Psi^c] = C'_mnr Tr[Psi^(m)Psi^(n)Psi^(r)]
    with Psi^a = p^{-1/2} sum_m w^{-ma} Psi^(m)  (inverse of Psi^(m) = p^{-1/2} sum_a w^{ma} Psi^a)."""
    p = C.shape[0]; w = np.exp(2j * np.pi / p)
    U = np.array([[w ** (-m * a) for a in range(p)] for m in range(p)]) / np.sqrt(p)   # U[m,a]
    return np.einsum('abc,ma,nb,rc->mnr', C, U, U, U)


def supercharge(C):
    """Q = C_abc Tr[Psi^a Psi^b Psi^c] as an Expr (letters in whatever basis C is given)."""
    p = C.shape[0]; out = Expr()
    for a, b, c in itertools.product(range(p), repeat=3):
        if abs(C[a, b, c]) > 1e-15:
            out = out + trace(((a, 'P'), (b, 'P'), (c, 'P'))) * NPoly.const(C[a, b, c])
    return out


def number_operator(p):
    out = Expr()
    for a in range(p):
        out = out + trace(((a, 'P'), (a, 'B')))
    return out


def hamiltonian(C):
    """H = {Q, Qbar} from D2.11 for a cyclic coupling C (any basis with the delta^ab anticommutator):
    9 sum_a Tr[X_a Xbar_a] - 9N M_cd Tr[Psi^c Psibar^d] + 9 L_cd Tr Psi^c Tr Psibar^d + 3N^3|C|^2 - 3N<C,C^rev>,
    X_a = C_abc Psi^b Psi^c, Xbar_a = conj(C_abc) Psibar^c Psibar^b."""
    p = C.shape[0]; Cb = np.conj(C)
    M = np.einsum('abc,abd->cd', C, Cb); Lm = np.einsum('abc,bad->cd', C, Cb)
    normC = np.einsum('abc,abc->', C, Cb); rev = np.einsum('abc,cba->', C, Cb)
    H = Expr()
    for a in range(p):
        for b, c, d, e in itertools.product(range(p), repeat=4):
            coef = 9 * C[a, b, c] * Cb[a, d, e]
            if abs(coef) > 1e-15:
                H = H + trace(((b, 'P'), (c, 'P'), (e, 'B'), (d, 'B'))) * NPoly.const(coef)
    for c, d in itertools.product(range(p), repeat=2):
        if abs(M[c, d]) > 1e-15:
            H = H + trace(((c, 'P'), (d, 'B'))) * NPoly.N(1, -9 * M[c, d])
        if abs(Lm[c, d]) > 1e-15:
            H = H + canonical((((c, 'P'),), ((d, 'B'),))) * NPoly.const(9 * Lm[c, d])
    H = H + Expr({(): NPoly({3: 3 * normC, 1: -3 * rev}).clean()})
    return H.clean()


def sandwich(wa, wb, S):
    """Sum_ij (wa)_ij^dag S (wb)_ij for open words wa, wb and a scalar operator S (an Expr): the open indices are
    contracted into a single cycle (as in Tr[wa^dag wb]) while S sits between them in OPERATOR order, its own traces
    forming separate cycles.  Used for the ground-state / BPS positivity cone [phi(X^dag H Y)] >= 0."""
    wad = word_dagger(wa)
    out = Expr()
    for mono, coeff in S.items():
        ops = list(wad)
        cyc_starts = []
        for word in mono:
            cyc_starts.append((len(ops), len(word))); ops += list(word)
        nb = len(ops)
        ops += list(wb)
        L = len(ops)
        nxt = [0] * L
        outer = list(range(len(wad))) + list(range(nb, L))     # wa^dag then wb, contracted into one cycle
        for t, pos in enumerate(outer):
            nxt[pos] = outer[(t + 1) % len(outer)]
        for st, ln in cyc_starts:                              # S's traces: own cycles
            for t in range(ln):
                nxt[st + t] = st + (t + 1) % ln
        c = coeff if outer else coeff * NPoly.N(1)             # both words empty: a closed index loop gives N
        for m2, c2 in canonicalize_routed(ops, nxt, c).items():
            out.add(m2, c2)
    return out.clean()


def casimir(p):
    """Quadratic gauge Casimir as a trace polynomial (docs/derivations.md D7):
        C2 = 1/2 Tr[M^2] - (Tr M)^2/(2N),   M = sum_a (Psi^a Psibar^a + Psibar^a Psi^a),   Tr M = p N^2,
    in the normalisation Tr T^a T^b = delta/2 (C2 = j(j+1) at N=2).  Flavor-U(p) invariant, so valid in the Fourier
    letter basis as well.  Verified against the explicit operator at (N,p) = (2,3), (2,1), (3,1)."""
    P = lambda c: (c, 'P'); B = lambda c: (c, 'B')
    e = Expr()
    for a in range(p):
        for b in range(p):
            for w in ((P(a), B(a), P(b), B(b)), (P(a), B(a), B(b), P(b)), (B(a), P(a), P(b), B(b)), (B(a), P(a), B(b), P(b))):
                for m2, c2 in trace(w).items():
                    e.add(m2, c2 * 0.5)
    return e + Expr({(): NPoly({3: -p * p / 2.0})})


def _M_terms(p):
    """M = sum_a (Psi^a Psibar^a + Psibar^a Psi^a) as a list of length-2 matrix words (the gauge-generator matrix;
    its trace is the c-number p N^2)."""
    out = []
    for a in range(p):
        out.append(((a, 'P'), (a, 'B')))
        out.append(((a, 'B'), (a, 'P')))
    return out


def casimir3(p):
    """Cubic gauge Casimir as a trace polynomial: Tr[Mtilde^3] with Mtilde = M - (p N) 1 (traceless part, since
    Tr M = p N^2).  Central like C2 (a Gelfand invariant of gl(N)); flavour-U(p) invariant, so basis-independent."""
    terms = _M_terms(p)
    e = Expr()
    for w1 in terms:                                            # Tr[M^3]
        for w2 in terms:
            for w3 in terms:
                for m, c in trace(w1 + w2 + w3).items():
                    e.add(m, c)
    e2 = Expr()                                                 # -3 p N Tr[M^2]
    for w1 in terms:
        for w2 in terms:
            for m, c in trace(w1 + w2).items():
                e2.add(m, c * NPoly.N(1, -3.0 * p))
    e3 = Expr()                                                 # +3 (p N)^2 Tr[M]
    for w1 in terms:
        for m, c in trace(w1).items():
            e3.add(m, c * NPoly({2: 3.0 * p * p}))
    return (e + e2 + e3 + Expr({(): NPoly({4: -1.0 * p ** 3})})).clean()   # - N (p N)^3 / N^0: Tr[1 (pN)^3] = N (pN)^3


def casimir3_value(N, lam):
    """Eigenvalue of casimir3(p) = Tr[Mtilde^3] on the U(N) irrep lam (sum lam_i = 0):
        c3 = f(l) - f(l0),  f(x) = sum_i [x_i^3 + (3/2 - N) x_i^2],  l_i = lam_i + N - i,  l0_i = N - i.
    Fitted and then verified against the joint (C2, C3) spectrum of the k = 1, 2 sectors at N = 2, 3, 4 (full
    multiset match).  At N = 2 it gives c3 = 2 c2, as it must (su(2) has no independent cubic Casimir)."""
    f = lambda x: sum(t ** 3 + (1.5 - N) * t ** 2 for t in x)
    return f([lam[i] + N - (i + 1) for i in range(N)]) - f([N - (i + 1) for i in range(N)])


def casimir_value(N, lam):
    """C2 eigenvalue of the U(N) irrep lam (sum lam_i = 0) in the same normalisation: 1/2 sum_i lam_i (lam_i + N + 1 - 2i)."""
    return 0.5 * sum(l * (l + N + 1 - 2 * (i + 1)) for i, l in enumerate(lam))


# ----------------------------------------------------------------------------------------------------------------
# finite-N relations
# ----------------------------------------------------------------------------------------------------------------

def finiteN_relations(blocks, N):
    """Relations from the vanishing antisymmetriser over N+1 matrix indices,
        delta^{[i_1}_{j_1} ... delta^{i_{N+1}]}_{j_{N+1}} = 0,
    contracted with N+1 matrix blocks (each a tuple of letters, i.e. an open word; a single letter is a block of
    length 1) in the fixed operator order given:  sum_sigma sgn(sigma) prod_a (M_a)_{i_a i_sigma(a)} = 0.
    Each sigma is a routing of the blocks; the sum of the canonicalised routed terms is an Expr that vanishes
    identically as an operator (the zero Expr means no new relation).  Total length sum(len(block)) >= N+1."""
    n = len(blocks); assert n == N + 1
    blocks = [(b,) if (len(b) == 2 and isinstance(b[1], str)) else tuple(b) for b in blocks]   # bare letters allowed
    starts, ops = [], []
    for b in blocks:
        starts.append(len(ops)); ops += list(b)
    out = Expr()
    for perm in itertools.permutations(range(n)):
        sgn = 1; seen = [False] * n
        for s0 in range(n):
            if seen[s0]:
                continue
            l = 0; i = s0
            while not seen[i]:
                seen[i] = True; i = perm[i]; l += 1
            if l % 2 == 0:
                sgn = -sgn
        nxt = [0] * len(ops)
        for a, b in enumerate(blocks):
            for t in range(len(b) - 1):
                nxt[starts[a] + t] = starts[a] + t + 1
            nxt[starts[a] + len(b) - 1] = starts[perm[a]]
        for m, c in canonicalize_routed(list(ops), nxt, NPoly.const(sgn)).items():
            out.add(m, c)
    return out.clean()


# ----------------------------------------------------------------------------------------------------------------
# numerical evaluation against explicit operators (harness helpers)
# ----------------------------------------------------------------------------------------------------------------

def routed_apply(model, ops, nxt, v):
    """Apply the routed operator (sum over all index assignments consistent with the routing) to vector v."""
    from fermion_matrix_model import letter_op
    N = model['N']
    cycles = _cycles(nxt)
    # free index per routing edge: edge e_i = column index of letter i (= row index of nxt[i])
    out = np.zeros_like(v)
    for assign in itertools.product(range(N), repeat=len(ops)):
        # assign[i] = value of edge e_i; row of letter i = assign[prev[i]]
        prev = [0] * len(nxt)
        for i, j in enumerate(nxt):
            prev[j] = i
        w = v
        for i in reversed(range(len(ops))):            # rightmost operator acts first
            w = letter_op(model, ops[i], assign[prev[i]], assign[i]) @ w
        out = out + w
    return out


def mono_apply(model, mono, v):
    """Apply Tr[w_1]...Tr[w_m] (operator order) to v using explicit word matrices."""
    from fermion_matrix_model import word_matrix
    N = model['N']
    w = v
    for word in reversed(mono):
        if len(word) == 0:
            w = N * w
            continue
        W = word_matrix(model, word)
        w = sum(W[i][i] @ w for i in range(N))
    return w


def expr_apply(model, expr, v):
    N = model['N']
    out = np.zeros_like(v)
    for m, c in expr.items():
        cv = c.at(N)
        if abs(cv) > 1e-15:
            out = out + cv * mono_apply(model, m, v)
    return out
