#!/usr/bin/env python3
"""
Exact gauge-singlet multiplicities n(k) and refined index of the U(n)^3 fermionic quiver (p flavours per edge),
by dual Cauchy decomposition instead of a weight expansion of the 3pn^2-mode Fock character.

    Lambda(C^p (x) V_0 (x) Vbar_1) = (x)_{f=1..p} (+)_{lambda in n x n box} S_lambda(V_0) (x) S_lambda'(Vbar_1)

so a Fock state is labelled by p-tuples of box partitions P_A, P_B, P_C (edges A: 0->1, B: 1->2, C: 2->0), and the
number of U(n)^3 singlets is
    n(t) = sum_{A,B,C} t^{|A|+|B|+|C|} G[A, C'] G[B, A'] G[C, B'] = Tr[(D_t G T)^3],
with G[P,Q] = int_{U(n)} s_P conj(s_Q) dU (s_P = prod_f s_{lambda_f}), T the transpose map P -> P', D_t the grading.
G is computed by exact quadrature on a shifted uniform grid of the maximal torus (Weyl integration formula; the
integrand is a Laurent polynomial of degree < M in each variable, so the grid sum is exact), then rounded.

Refined index per class c = k mod 3 in the convention of quiver_index.py / run_quiver_irrep_cohomology.py:
I_c = sum_{k = c mod 3} (-1)^((k-c)/3) n(k).

    python scripts/quiver_singlet_index.py --n 2 --p 2 --check      # validates against a direct weight count
    python scripts/quiver_singlet_index.py --n 3 --p 2
"""
import argparse, itertools, json, math, os, sys, time
import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')


def box_partitions(n):
    """Partitions with at most n parts, each <= n (the n x n box), as length-n tuples."""
    out = []
    def rec(prefix, maxpart):
        if len(prefix) == n:
            out.append(tuple(prefix)); return
        for x in range(maxpart, -1, -1):
            rec(prefix + [x], x)
    rec([], n)
    return out


def transpose(lam, n):
    return tuple(sum(1 for x in lam if x > j) for j in range(n))


def schur_on_grid(lam, X):
    """s_lam(x) for each grid point (rows of X, shape (pts, n)) via the bialternant."""
    n = X.shape[1]
    rho = np.arange(n - 1, -1, -1)
    num = np.linalg.det(X[:, None, :] ** (np.array(lam) + rho)[None, :, None])
    den = np.linalg.det(X[:, None, :] ** rho[None, :, None])
    return num / den, den


def singlet_series(n, p, M=None):
    parts = box_partitions(n)
    idx = {l: i for i, l in enumerate(parts)}
    tuples = list(itertools.product(range(len(parts)), repeat=p))
    tidx = {t: i for i, t in enumerate(tuples)}
    deg = np.array([sum(sum(parts[i]) for i in t) for t in tuples])
    tau = np.array([tidx[tuple(idx[transpose(parts[i], n)] for i in t)] for t in tuples])
    # torus grid: integrand exponents bounded by 2pn + (n-1) per variable
    M = M or (2 * p * n + n + 2)
    rng = np.random.default_rng(0)
    shifts = rng.uniform(0, 2 * np.pi / M, size=n)            # generic offsets avoid x_i = x_j
    grid = np.array(list(itertools.product(range(M), repeat=n)), dtype=float) * 2 * np.pi / M + shifts
    X = np.exp(1j * grid)
    S = []
    for lam in parts:
        s, den = schur_on_grid(lam, X)
        S.append(s)
    S = np.array(S)                                          # (parts, pts)
    w = np.abs(den) ** 2 / (math.factorial(n) * M ** n)   # Weyl measure, |a_rho|^2 / n!
    Sp = np.ones((len(tuples), X.shape[0]), dtype=complex)
    for f in range(p):
        Sp *= S[[t[f] for t in tuples]]
    Gc = (Sp * w) @ Sp.conj().T
    G = np.rint(Gc.real).astype(np.int64)
    err = float(np.abs(Gc - G).max())
    Tm = np.zeros((len(tuples), len(tuples)), dtype=np.int64)
    Tm[np.arange(len(tuples)), tau] = 1                      # (G T)[A, C] = G[A, tau C]
    GT = G @ Tm
    dmax = int(deg.max())
    Md = [GT * (deg == d)[:, None] for d in range(dmax + 1)]
    nk = np.zeros(3 * dmax + 1, dtype=object)
    # Tr(M_{d1} M_{d2} M_{d3}) for all degree triples (exact integers)
    pair = {}
    for d2 in range(dmax + 1):
        for d3 in range(dmax + 1):
            pair[(d2, d3)] = Md[d2] @ Md[d3]
    for d1 in range(dmax + 1):
        for d2 in range(dmax + 1):
            for d3 in range(dmax + 1):
                v = int(np.einsum('ij,ji->', Md[d1], pair[(d2, d3)]))
                if v:
                    nk[d1 + d2 + d3] += v
    return [int(x) for x in nk], err, len(tuples)


def refined_index(nk):
    I = [0, 0, 0]
    for k, m in enumerate(nk):
        I[k % 3] += (-1) ** ((k - k % 3) // 3) * m
    return I


def direct_count(n, p):
    """Singlet multiplicities by the Weyl alternating sum over weight-space dimensions (small n, p only)."""
    from collections import defaultdict
    nm = 3 * p * n * n
    W = []
    for e, f, x, y in itertools.product(range(3), range(p), range(n), range(n)):
        w = [0] * (3 * n); w[e * n + x] += 1; w[((e + 1) % 3) * n + y] -= 1; W.append(tuple(w))
    Z = {(0, tuple([0] * (3 * n))): 1}
    for w in W:
        new = defaultdict(int)
        for (k, v), c in Z.items():
            new[(k, v)] += c; new[(k + 1, tuple(a + b for a, b in zip(v, w)))] += c
        Z = new
    rho = tuple(n - 1 - i for i in range(n))
    def sgn(pm):
        s = 1
        for i in range(n):
            for j in range(i + 1, n):
                if pm[i] > pm[j]: s = -s
        return s
    shifts = []
    for p3 in itertools.product(list(itertools.permutations(range(n))), repeat=3):
        eps, sh = 1, [0] * (3 * n)
        for e, pm in enumerate(p3):
            eps *= sgn(pm)
            for i in range(n): sh[e * n + i] = rho[i] - rho[pm[i]]
        shifts.append((eps, tuple(sh)))
    return [sum(eps * Z.get((k, sh), 0) for eps, sh in shifts) for k in range(nm + 1)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--n', type=int, default=3); ap.add_argument('--p', type=int, default=2)
    ap.add_argument('--check', action='store_true', help='compare with the direct weight count (small cases)')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    t0 = time.time()
    nk, err, ntup = singlet_series(a.n, a.p)
    I = refined_index(nk)
    nz = {k: m for k, m in enumerate(nk) if m}
    print(f"(n,p)=({a.n},{a.p}): {ntup} partition tuples, quadrature rounding error {err:.1e}; "
          f"{3 * a.p * a.n ** 2} modes, total singlet states {sum(nk)}  [{time.time() - t0:.1f}s]")
    print(f"  singlet multiplicity n(k): {nz}")
    print(f"  refined singlet index I_0, I_1, I_2 = {I}")
    rec = dict(n=a.n, p=a.p, n_k=nk, index=I, quadrature_error=err)
    if a.check:
        dc = direct_count(a.n, a.p)
        ok = dc == nk
        print(f"  direct weight count agrees: {ok}" + ("" if ok else f"  direct={dc}"))
        rec['direct_check'] = ok
    if a.out:
        json.dump(rec, open(a.out, 'w'), indent=1); print('  saved', a.out)


if __name__ == '__main__':
    main()
