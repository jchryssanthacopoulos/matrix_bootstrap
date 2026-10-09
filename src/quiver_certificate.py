"""
Exact verification of a Farkas (infeasibility) certificate of the quiver singlet bootstrap (src/quiver_bootstrap.py;
docs/derivations.md D21.8).

The bootstrap problem is  A_eq x = b_eq,  M_c(x) = sum_j x_j F_{c,j} >= 0  for the cones c, over the reduced variables
x_j = phi(rep_j) (one representative monomial per class identified by the reality relations).  A certificate is
(mu, Y_c) with Y_c >= 0 and
    r_j := sum_k mu_k (A_eq)_{kj} - sum_c <Y_c, F_{c,j}>   (= 0 ideally),        mu . b_eq < 0 .
If a singlet BPS state existed, its functional x* would satisfy every row exactly and make every M_c(x*) >= 0, so
    mu . b_eq = r . x* + sum_c <Y_c, M_c(x*)>  >=  - sum_j |r_j| |x*_j|  >=  - sum_j |r_j| n^{L_j},
using |phi(m)| <= ||m|| <= n^{L(m)} for a product of traces of total length L (n^L index assignments, each a product
of fermion operators of norm <= 1).  Hence  mu . b_eq + sum_j |r_j| n^{L_j} < 0  proves that no such state exists.

Everything below is done in exact integer arithmetic:
  * constraint coefficients: the algebra produces rationals whose denominators divide n * 2^a (one 1/n from the
    adjoint projection, powers of 1/2 from the self-symmetry resolution, integer couplings); each float coefficient c
    is replaced by round(c D) / D with D = n * 2^20, and the largest rounding deviation is reported (it must be << 1);
  * the solver's certificate (mapped back through the presolve and the row normalisation) is rounded to the grid
    2^-K, and each dual block is shifted by a small delta * 1 (doubling from 2^-40 relative) until it is positive
    definite, verified exactly by Sylvester's criterion with Bareiss elimination;
  * r_j, mu . b_eq and the bound are then exact integers (times the common factor 2^K D).
What remains unverified by this script is the correctness of the constraint generation itself (the trace algebra and
the row families), which is tested against explicit Fock-space operators in tests/test_quiver_trace.py and
tests/test_quiver_bootstrap.py.
"""
import math
import numpy as np


def _svec_pairs(d, solver):
    """(a, b) with a <= b in the order of the solver's packed PSD rows (as in QuiverSDP.assemble)."""
    if solver == 'clarabel':
        return [(i, j) for j in range(d) for i in range(j + 1)]
    return [(j, i) for j in range(d) for i in range(j, d)]


def _positive_definite(Mint):
    """Exact test that an integer symmetric matrix (list of lists) is positive definite: Sylvester's criterion with
    the leading principal minors computed by Bareiss fraction-free elimination (exact integer divisions)."""
    d = len(Mint)
    A = [row[:] for row in Mint]
    prev = 1
    for k in range(d):
        if A[k][k] <= 0:                 # = the (k+1)-th leading principal minor
            return False
        for i in range(k + 1, d):
            for j in range(k + 1, d):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return True


def verify_farkas_exact(S, data, y, D=None, K_bits=62, verbose=True):
    """S: the QuiverSDP; data: the dict returned by S.assemble(...) for the solve that produced y (mode
    'feasibility'); y: the solver's dual vector (Clarabel sol.z or SCS sol['y']).  Returns a report dict."""
    n = S.n
    D = D if D is not None else n * 2 ** 20
    col, sign = S._col, S._sign
    n_eq = data['n_eq']; solver = data['solver']
    row_map, col_map, rn, srcs, extra = data['row_map'], data['col_map'], data['rn'], data['srcs'], data['extra']
    assert data['mode'] == 'feasibility', 'verify the pure feasibility problem'
    y = np.asarray(y, float)
    # ---- certificate pieces (floats)
    mu = np.zeros(n_eq)
    for k in range(n_eq):
        mu[k] = y[k] / rn[row_map[k]]                 # multiplier of the UNnormalised assembled row
    Ys, off = [], n_eq
    for d_c in data['dims']:
        L = d_c * (d_c + 1) // 2
        blk = y[off:off + L]; off += L
        Y = np.zeros((d_c, d_c))
        for t, (a, b) in enumerate(_svec_pairs(d_c, solver)):
            v = blk[t] / (1.0 if a == b else math.sqrt(2.0)); Y[a, b] = v; Y[b, a] = v
        Ys.append(Y)
    scale = max(np.abs(mu).max(), max(np.abs(Y).max() for Y in Ys))
    K = K_bits - int(math.ceil(math.log2(scale)))
    Mint = [int(round(v * 2.0 ** K)) for v in mu]
    # ---- exact integer rows (times D) for the final equality rows
    maxdev = [0.0]

    def to_int(v):
        q = v * D
        Q = round(q)
        maxdev[0] = max(maxdev[0], abs(q - Q))
        return int(Q)

    def exact_row(src):
        kind, i = src
        if kind == 'norm':
            return {int(col[0]): D}, D
        ix, v = (extra[i] if kind == 'dag' else S.rows[i])
        out = {}
        for m, c in zip(ix, v):
            j = int(col[m])
            if j < 0:
                continue
            out[j] = out.get(j, 0) + int(round(sign[m])) * to_int(c)
        return {j: c for j, c in out.items() if c != 0}, 0
    r = {}
    b_term = 0
    for k in range(n_eq):
        if Mint[k] == 0:
            continue
        row, rhs = exact_row(srcs[row_map[k]])
        for j, a in row.items():
            r[j] = r.get(j, 0) + Mint[k] * a
        b_term += Mint[k] * rhs
    # ---- cones: exact PSD-shifted integer duals and their contribution
    shifts, worst_rel_shift = [], 0.0
    for c_idx, (cone, Y) in enumerate(zip(S.cones, Ys)):
        d_c = cone['size']
        Yint = [[int(round(Y[a, b] * 2.0 ** K)) for b in range(d_c)] for a in range(d_c)]
        ymax = max(1, max(abs(v) for row in Yint for v in row))
        delta = max(1, ymax >> 40)                     # start at ~1e-12 relative, double until positive definite
        while True:
            Ysh = [[Yint[a][b] + (delta if a == b else 0) for b in range(d_c)] for a in range(d_c)]
            if _positive_definite(Ysh):
                break
            delta *= 2
            if delta > ymax:
                raise RuntimeError(f'cone {c_idx}: dual block not positive semidefinite within a reasonable shift')
        shifts.append(delta); worst_rel_shift = max(worst_rel_shift, delta / ymax)
        for (a, b), (ix, v) in cone['entries'].items():
            w = 1 if a == b else 2
            yab = Ysh[a][b] * w
            if yab == 0:
                continue
            for m, cval in zip(ix, v):
                j = int(col[m])
                if j < 0:
                    continue
                r[j] = r.get(j, 0) - yab * int(round(sign[m])) * to_int(cval)
    # ---- the bound: |x*_j| <= n^{L_j}
    Lmax = {}
    for i, mo in enumerate(S.monos):
        j = int(col[i])
        if j >= 0:
            Lmax[j] = max(Lmax.get(j, 0), sum(len(w) for w in mo))
    bound = sum(abs(v) * n ** Lmax[j] for j, v in r.items() if v != 0)
    r_l2 = math.sqrt(sum(float(v) ** 2 for v in r.values()))
    y_l2 = math.sqrt(sum(float(v) ** 2 for v in Mint))
    verified = (b_term + bound) < 0
    rep = dict(verified=bool(verified), D=D, K=K, max_coeff_rounding_deviation=maxdev[0],
               b_term=float(b_term), bound_term=float(bound), ratio_bound_over_b=float(bound) / abs(float(b_term)) if b_term else float('inf'),
               r_l2_over_mu_l2=r_l2 / y_l2 if y_l2 else float('nan'), n_nonzero_residuals=sum(1 for v in r.values() if v),
               max_rel_psd_shift=worst_rel_shift, n_cones=len(S.cones), n_eq=n_eq)
    if verbose:
        print(f"  exact check: coefficient rounding deviation {maxdev[0]:.1e} (of 1), grid 2^-{K}, PSD shifts <= {worst_rel_shift:.1e} "
              f"relative; mu.b = {float(b_term):.4e}, sum|r_j| n^L_j = {float(bound):.4e}, ratio {rep['ratio_bound_over_b']:.3e} "
              f"-> {'VERIFIED' if verified else 'NOT verified'}", flush=True)
    return rep
