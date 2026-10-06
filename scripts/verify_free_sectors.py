#!/usr/bin/env python3
"""
Numerical checks of derivation D16 (free sectors of Chen's three-matrix model).

  characters --N 4 --kmax 15     irreducible content of Lambda^k sl(N): the largest Casimir c_max(k), whether it
                                 equals N k (D16 equality), the size of the c = N k component (ground degeneracy)
                                 and the variational bound E_c - 38Nk + (100/3)(Nk - c_max).
  kernel --N 4                   exact dim ker(delta) on Lambda^k sl(N) for all k, against the c = N k component.
  state --N 4 --k 4 --kind block energy of prod phi^dag(a_m)|0> as (|Q psi|^2 + |Q^dag psi|^2)/|psi|^2 (independent of
                                 D2.13), plus the D2.13 residual |H psi - E psi| on the state's weight block.
                                 --kind: block (abelian matrix units), cartan, control (E_01, E_10: non-commuting).
  lanczos --N 4 --k 4            lowest eigenvalues of H (D2.13) on the zero-weight block of the k-particle sector,
                                 which contains every irrep and therefore the sector ground energy; also checks the
                                 D2.13 builder against {Q, Q^dag} on a random vector.

Deterministic except the random cross-check vector (--seed).
"""
import argparse
import math
import os
import sys
import time

import numpy as np
import scipy.sparse.linalg as sla

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from fermion_matrix_model import chen_C                                       # noqa: E402
from cohomology import weight_basis, kostant                                  # noqa: E402
import free_sectors as fs                                                     # noqa: E402


def cmd_characters(a):
    N = a.N
    kmax = min(a.kmax, N * N - 1)
    Ec = fs.vacuum_energy(N, chen_C())
    chars = fs.wedge_sl_characters(N, kmax)
    thr = N * N // 4
    print(f"N={N}: Lambda^k sl(N), k=0..{kmax}; floor(N^2/4)={thr}; E_c={Ec:.0f}")
    print(f"{'k':>3} {'c_max':>8} {'N k':>5} {'equal':>6} {'deg(c=Nk)':>10} {'E_free':>10} {'E_var':>12}  "
          f"top irreps (lambda: mult x dim)")
    for k in range(kmax + 1):
        irr = fs.decompose(N, chars[k])
        assert all(m > 0 for m in irr.values()), "negative multiplicity: decomposition failed"
        total = sum(m * fs.weyl_dim(l) for l, m in irr.items())
        assert total == math.comb(N * N - 1, k), (total, math.comb(N * N - 1, k))
        c2max = max(fs.casimir2(l) for l in irr)
        top = {l: m for l, m in irr.items() if fs.casimir2(l) == c2max}
        deg = sum(m * fs.weyl_dim(l) for l, m in irr.items() if fs.casimir2(l) == 2 * N * k)
        Efree = fs.free_energy(N, k)
        Evar = Efree + 100 / 3 * (N * k - c2max / 2)
        tops = ', '.join(f"{l}: {m}x{fs.weyl_dim(l)}" for l, m in sorted(top.items()))
        print(f"{k:>3} {c2max / 2:>8.1f} {N * k:>5} {str(c2max == 2 * N * k):>6} {deg:>10} {Efree:>10} "
              f"{Evar:>12.4f}  {tops}")
        assert (c2max == 2 * N * k) == (k <= thr) or k == 0, f"threshold mismatch at k={k}"
        assert c2max <= 2 * N * k, "Casimir exceeds N k: contradicts the D16 identity"
    print("threshold check: c_max(k) = N k exactly for k <= floor(N^2/4), and < N k above  [passed]")


def cmd_kernel(a):
    N = a.N
    t0 = time.time()
    ker = fs.delta_kernel_dims(N)
    chars = fs.wedge_sl_characters(N, N * N - 1)
    print(f"N={N}: dim ker(delta) on Lambda^k sl(N), exact over F_P, against the c = N k component  "
          f"({time.time() - t0:.1f}s)")
    ok = True
    for k in range(N * N):
        irr = fs.decompose(N, chars[k])
        pred = sum(m * fs.weyl_dim(l) for l, m in irr.items() if fs.casimir2(l) == 2 * N * k)
        flag = 'ok' if pred == ker[k] else 'MISMATCH'
        ok &= pred == ker[k]
        print(f"  k={k:>2}: dim Lambda^k = {math.comb(N * N - 1, k):>6}, dim ker delta = {ker[k]:>5}, "
              f"dim(c=Nk) = {pred:>5}  {flag}")
    print("all agree" if ok else "DISAGREEMENT")


def _state_matrices(N, k, kind):
    if kind == 'block':
        mats = fs.block_matrices(N)
    elif kind == 'cartan':
        mats = fs.cartan_matrices(N)
    elif kind == 'control':
        E01 = np.zeros((N, N)); E01[0, 1] = 1
        E10 = np.zeros((N, N)); E10[1, 0] = 1
        mats = [E01, E10]
    else:
        raise ValueError(kind)
    if k > len(mats):
        raise SystemExit(f"--kind {kind} supplies only {len(mats)} commuting matrices at N={N}")
    return mats[:k]


def cmd_state(a):
    N, k = a.N, a.k
    C = chen_C()
    mats = _state_matrices(N, k, a.kind)
    comm = max((np.abs(x @ y - y @ x).max() for x in mats for y in mats), default=0)
    psi = fs.free_state(N, mats)
    print(f"N={N}, k={k}, kind={a.kind}: {len(psi)} Fock components, max |[a_l,a_m]| = {comm:g}")
    Efree = fs.free_energy(N, k)
    pred = Efree
    if a.kind == 'control':                  # D16 identity: <H_4> = (100/3) |[a,b]|_F^2 for |psi| = 1
        br = mats[0] @ mats[1] - mats[1] @ mats[0]
        pred = Efree + 100 / 3 * np.sum(br ** 2)
    if not a.skip_q:
        t0 = time.time()
        E = fs.energy_from_Q(N, 3, C, psi)
        print(f"  <H> from |Q psi|^2 + |Q^dag psi|^2 : {E:.10f}   ({time.time() - t0:.1f}s)")
    print(f"  D16 prediction                     : {pred:.10f}   (free value 16N^3-(15+38k)N = {Efree})")
    # residual with the D2.13 operator on the state's weight block
    lam = [0] * N
    for st in psi:
        for m in st:
            rem = m % (N * N); i, j = divmod(rem, N)
            lam[i] += 1; lam[j] -= 1
        break
    t0 = time.time()
    basis = weight_basis(N, 3, k, tuple(lam))
    H1, Yd = fs.sector_operators(N, C, basis)
    index = {s: n for n, s in enumerate(basis)}
    v = np.zeros(len(basis))
    for st, amp in psi.items():
        v[index[st]] = amp
    v /= np.linalg.norm(v)
    Hv = fs.vacuum_energy(N, C) * v + H1 @ v + 9 * (Yd.T @ (Yd @ v))
    E2 = v @ Hv
    print(f"  D2.13 on weight block {tuple(lam)} (dim {len(basis)}): <H> = {E2:.10f}, "
          f"|H psi - <H> psi| = {np.linalg.norm(Hv - E2 * v):.2e}, |Ydag psi| = {np.linalg.norm(Yd @ v):.2e}  "
          f"({time.time() - t0:.1f}s)")


def cmd_lanczos(a):
    N, k = a.N, a.k
    C = chen_C()
    Ec = fs.vacuum_energy(N, C)
    lam = tuple([0] * N)
    t0 = time.time()
    basis = weight_basis(N, 3, k, lam)
    print(f"N={N}, k={k}: zero-weight block dim {len(basis)} (sector dim {math.comb(3 * N * N, k)})  "
          f"[{time.time() - t0:.1f}s]", flush=True)
    t0 = time.time()
    H1, Yd = fs.sector_operators(N, C, basis)
    print(f"  operators: H1 nnz {H1.nnz}, Ydag {Yd.shape[0]}x{Yd.shape[1]} nnz {Yd.nnz}  [{time.time() - t0:.1f}s]",
          flush=True)
    herm = abs(H1 - H1.T).max()
    assert herm < 1e-10, f"H1 not symmetric: {herm}"
    # cross-check the D2.13 builder against {Q, Q^dag} on a random vector with small support
    rng = np.random.default_rng(a.seed)
    supp = rng.choice(len(basis), size=min(a.check_support, len(basis)), replace=False)
    v = np.zeros(len(basis)); v[supp] = rng.standard_normal(len(supp)); v /= np.linalg.norm(v)
    E_d213 = v @ (Ec * v + H1 @ v + 9 * (Yd.T @ (Yd @ v)))
    t0 = time.time()
    E_q = fs.energy_from_Q(N, 3, C, {basis[n]: v[n] for n in supp})
    print(f"  builder check on a random {len(supp)}-state vector: D2.13 {E_d213:.10f} vs {{Q,Q^dag}} {E_q:.10f}, "
          f"diff {abs(E_d213 - E_q):.1e}  [{time.time() - t0:.1f}s]", flush=True)
    op = sla.LinearOperator((len(basis), len(basis)), dtype=float,
                            matvec=lambda x: Ec * x + H1 @ x + 9 * (Yd.T @ (Yd @ x)))
    t0 = time.time()
    vals = np.sort(sla.eigsh(op, k=a.n_eig, which='SA', tol=a.tol, return_eigenvectors=False))
    print(f"  lowest {a.n_eig} eigenvalues  [{time.time() - t0:.1f}s]:")
    groups = []
    for x in vals:
        if groups and abs(x - groups[-1][0]) < 1e-6 * max(1, abs(x)):
            groups[-1][1] += 1
        else:
            groups.append([x, 1])
    for x, m in groups:
        print(f"    {x:.8f}  x{m}")
    print(f"  D16 free value 16N^3-(15+38k)N = {fs.free_energy(N, k)}  "
          f"(k {'<=' if k <= N * N // 4 else '>'} floor(N^2/4) = {N * N // 4})")
    if a.zero_mult:
        for lamtop in a.zero_mult:
            lt = tuple(int(x) for x in lamtop.split(','))
            print(f"  zero-weight multiplicity of V{lt}: {kostant(N, lt, lam)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('characters'); s.add_argument('--N', type=int, required=True)
    s.add_argument('--kmax', type=int, default=99)
    s = sub.add_parser('kernel'); s.add_argument('--N', type=int, required=True)
    s = sub.add_parser('state'); s.add_argument('--N', type=int, required=True); s.add_argument('--k', type=int,
                                                                                             required=True)
    s.add_argument('--kind', default='block', choices=['block', 'cartan', 'control'])
    s.add_argument('--skip-q', action='store_true', help='skip the (slow at N>=6) {Q,Q^dag} evaluation')
    s = sub.add_parser('lanczos'); s.add_argument('--N', type=int, required=True)
    s.add_argument('--k', type=int, required=True); s.add_argument('--n-eig', type=int, default=12)
    s.add_argument('--tol', type=float, default=1e-10); s.add_argument('--seed', type=int, default=1)
    s.add_argument('--check-support', type=int, default=60)
    s.add_argument('--zero-mult', nargs='*', help='highest weights "2,2,-2,-2" whose zero-weight mult to print')
    a = ap.parse_args()
    {'characters': cmd_characters, 'kernel': cmd_kernel, 'state': cmd_state, 'lanczos': cmd_lanczos}[a.cmd](a)


if __name__ == '__main__':
    main()
