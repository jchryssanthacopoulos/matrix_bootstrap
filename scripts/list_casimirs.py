#!/usr/bin/env python
"""Enumerate the U(N) irreps of the p-flavour fermionic matrix Fock space and their quadratic Casimirs.

    python scripts/list_casimirs.py --N 3 --p 3
    python scripts/list_casimirs.py --N 5 --p 3 --no-multiplicities   # skip the Fock-space check (memory)

Which irreps can appear, and why:

* The modes are Psi^a_ij with U(N) weight e_i - e_j, whose entries sum to zero.  Every Fock state is a product of
  modes, so every weight in the Fock space has total zero, and only highest weights with sum(lambda) = 0 occur.
* lambda_1 = (row 1 sum) - (column 1 sum) of the occupation matrix; the diagonal entry n_11 cancels between the
  two, leaving at most the p(N-1) off-diagonal modes of row 1.  Hence |lambda_i| <= p(N-1) for every i.
* That box bound is necessary but NOT sufficient: (12,12,0,-12,-12) at N=5 respects it yet cannot occur, because
  row 1 being full forces n_12 = p, which contradicts column 2 being empty.  The exact criterion is DOMINANCE by
  the maximal weight,

      lambda occurs  <=>  lambda dominant, sum(lambda) = 0, and lambda <= lambda*,

  where lambda <= mu means every partial sum satisfies sum_{i<=m} lambda_i <= sum_{i<=m} mu_i (totals equal).
  Verified exactly against the Weyl alternating sum (bps_index.irrep_multiplicities): 25 of 25 candidates at
  N=3, 213 of 241 at N=4 -- the dominance filter reproduces the occurring set on the nose in both cases.
  (Verified at N=3,4; not proved here.)

Quadratic Casimir (derivations D7), in the normalisation where H = 3N(N^2-1) - 9 C_2 for the single-matrix model:

    c_2(lambda) = (1/2) sum_i lambda_i (lambda_i + N + 1 - 2i).

The largest value is attained at the maximal weight lambda*_i = p(N + 1 - 2i), for which
lambda*_i + N+1-2i = (p+1)(N+1-2i) and therefore

    C_2^max = (1/2) p(p+1) sum_i (N+1-2i)^2 = p(p+1) N(N^2-1)/6,

which for p = 1 reduces to Chen 2025 (2.18), C_2(r_*) = N(N^2-1)/3.
"""
import argparse
import itertools
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
import trace_algebra as ta


def candidate_weights(N, p):
    """Non-increasing integer N-tuples with sum 0 and entries in [-p(N-1), p(N-1)]."""
    R = p * (N - 1)
    out = []
    for lam in itertools.combinations_with_replacement(range(R, -R - 1, -1), N):
        if sum(lam) == 0:
            out.append(tuple(lam))
    return out


def dominated(lam, mu):
    """lam <= mu in the dominance order: all partial sums of lam are <= those of mu, with equal totals."""
    s1 = s2 = 0
    for a, b in zip(lam, mu):
        s1 += a; s2 += b
        if s1 > s2:
            return False
    return s1 == s2


def maximal_weight(N, p):
    """The highest weight occurring in the Fock space: lambda*_i = p(N + 1 - 2i)."""
    return tuple(p * (N + 1 - 2 * i) for i in range(1, N + 1))


def occurring_weights(N, p):
    """The U(N) irreps of the Fock space: dominant, traceless, and dominated by lambda*."""
    star = maximal_weight(N, p)
    return [lam for lam in candidate_weights(N, p) if dominated(lam, star)]


def casimir_max(N, p):
    """Closed form p(p+1)N(N^2-1)/6, attained at lambda*_i = p(N+1-2i)."""
    return p * (p + 1) * N * (N * N - 1) / 6


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=int, required=True)
    ap.add_argument('--p', type=int, default=3)
    ap.add_argument('--no-multiplicities', action='store_true',
                    help='skip the Fock-space multiplicity check (it needs the full weight array)')
    a = ap.parse_args()
    cands = candidate_weights(a.N, a.p)
    occ_dom = set(occurring_weights(a.N, a.p))
    mult = None
    if not a.no_multiplicities:
        import bps_index as bi
        m = bi.irrep_multiplicities(a.N, a.p)
        mult = {}
        for (k, w, l), v in m.items():
            mult[l] = mult.get(l, 0) + v
    if mult is not None:
        assert set(mult) == occ_dom, "dominance filter disagrees with the Weyl multiplicities"
    rows = []
    for lam in cands:
        if lam not in occ_dom:
            continue
        rows.append((ta.casimir_value(a.N, lam), lam, None if mult is None else mult[lam]))
    rows.sort(key=lambda r: (-r[0], r[1]))
    star = maximal_weight(a.N, a.p)
    print(f"N={a.N} p={a.p}: {len(cands)} weights in the box, {len(rows)} occurring (dominated by lambda*)"
          + ("; multiplicities cross-checked against the Weyl sum" if mult is not None
             else "; multiplicities not computed"))
    print(f"maximal weight lambda* = {star}, C_2^max = p(p+1)N(N^2-1)/6 = {casimir_max(a.N, a.p):g} "
          f"(check: {ta.casimir_value(a.N, star):g})\n")
    print(f"{'C_2':>8s}  {'irrep':>22s}  {'multiplicity':>14s}")
    for c2, lam, mm in rows:
        print(f"{c2:>8.1f}  {str(lam):>22s}  {'-' if mm is None else format(mm, ',d'):>14s}")
    vals = sorted({r[0] for r in rows}, reverse=True)
    print(f"\n{len(vals)} distinct Casimir values: {[f'{v:g}' for v in vals]}")


if __name__ == '__main__':
    main()
