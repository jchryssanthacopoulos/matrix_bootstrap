#!/usr/bin/env python3
"""
Sector ground states of Chen's three-matrix model (p = q = 3, U(N)) by Lanczos on the zero-weight block of the
k-particle sector, with irrep labels from the gauge Casimir.

    python scripts/run_sector_ed.py --N 3 --k 8 --builder python --n-eig 12
    python scripts/run_sector_ed.py --N 3 --k 11 --builder fast --n-eig 12

The zero-weight block contains every irrep of the sector, so its lowest eigenvalue is the sector ground energy; a
level of irrep lambda appears K_{lambda,0} times per copy (zero-weight multiplicity).  Output: printed summary and
results/data/sector_ed_N{N}_k{k}.json (parameters, block size, eigenvalues with <C_2> and candidate irreps, the
D2.13-vs-{Q,Q^dag} cross-check, timings, peak RSS, git revision).  Deterministic given --seed.
"""
import argparse
import json
import math
import os
import resource
import subprocess
import sys
import time

import numpy as np
import scipy.sparse.linalg as sla

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
from fermion_matrix_model import chen_C                                       # noqa: E402
from cohomology import weight_basis                                           # noqa: E402
import free_sectors as fs                                                     # noqa: E402
import sector_ed as se                                                        # noqa: E402


def peak_gb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1e9 if sys.platform == 'darwin' else r / 1e6


def git_rev():
    try:
        rev = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'], cwd=ROOT, text=True).strip()
        dirty = subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True).strip()
        return rev + ('+uncommitted' if dirty else '')
    except Exception:
        return 'unknown'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--N', type=int, required=True)
    ap.add_argument('--k', type=int, required=True)
    ap.add_argument('--builder', choices=['python', 'fast'], default='fast')
    ap.add_argument('--n-eig', type=int, default=12)
    ap.add_argument('--tol', type=float, default=1e-10)
    ap.add_argument('--ncv', type=int, default=0, help='Lanczos basis size (default max(2 n_eig + 1, 40))')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--check-support', type=int, default=40)
    ap.add_argument('--out', default=None)
    ap.add_argument('--no-save', action='store_true')
    a = ap.parse_args()
    N, k, p = a.N, a.k, 3
    C = chen_C()
    Ec = fs.vacuum_energy(N, C)
    rng = np.random.default_rng(a.seed)
    rec = dict(N=N, k=k, p=p, builder=a.builder, n_eig=a.n_eig, tol=a.tol, seed=a.seed, git=git_rev(),
               date=time.strftime('%Y-%m-%d %H:%M'), sector_dim=math.comb(p * N * N, k))
    zero = (0,) * N

    t0 = time.time()
    basis = weight_basis(N, p, k, zero)
    masks = se.masks_from_tuples(basis)
    if a.builder == 'fast':
        order = np.argsort(masks)
        masks = masks[order]
        del basis, order
    dim = len(masks)
    rec['block_dim'] = dim
    rec['t_basis'] = round(time.time() - t0, 1)
    print(f"N={N} k={k}: zero-weight block {dim} of sector {rec['sector_dim']}  [{rec['t_basis']}s, "
          f"{peak_gb():.2f} GB]", flush=True)

    t0 = time.time()
    if a.builder == 'python':
        H1, Yd = fs.sector_operators(N, C, basis)
        del basis
    else:
        H1, Yd = se.build_fast(N, C, masks)
    rec['t_build'] = round(time.time() - t0, 1)
    rec['nnz'] = dict(H1=int(H1.nnz), Ydag=int(Yd.nnz), Ydag_rows=int(Yd.shape[0]))
    print(f"  {a.builder} builder: H1 nnz {H1.nnz}, Ydag {Yd.shape[0]}x{Yd.shape[1]} nnz {Yd.nnz}  "
          f"[{rec['t_build']}s, {peak_gb():.2f} GB]", flush=True)
    herm = abs(H1 - H1.T).max() if H1.nnz else 0.0
    assert herm < 1e-10, f"H1 not symmetric: {herm}"

    def matvec(x):
        return Ec * x + H1 @ x + 9 * (Yd.T @ (Yd @ x))

    # independent cross-check: <v|H|v> from D2.13 against |Q v|^2 + |Q^dag v|^2 on a small-support random vector
    supp = rng.choice(dim, size=min(a.check_support, dim), replace=False)
    v = np.zeros(dim); v[supp] = rng.standard_normal(len(supp)); v /= np.linalg.norm(v)
    E_d213 = float(v @ matvec(v))
    t0 = time.time()
    E_q = fs.energy_from_Q(N, p, C, {se.tuple_from_mask(masks[n]): v[n] for n in supp})
    rec['crosscheck'] = dict(support=int(len(supp)), D213=E_d213, Q=E_q, diff=abs(E_d213 - E_q))
    print(f"  cross-check on a random {len(supp)}-state vector: D2.13 {E_d213:.10f} vs {{Q,Q^dag}} {E_q:.10f} "
          f"(diff {abs(E_d213 - E_q):.1e})  [{time.time() - t0:.1f}s]", flush=True)
    assert abs(E_d213 - E_q) < 1e-8 * max(1.0, abs(E_q)), "D2.13 builder disagrees with {Q, Q^dag}"

    t0 = time.time()
    if dim <= 3000:
        Hd = np.column_stack([matvec(e) for e in np.eye(dim)])
        vals, vecs = np.linalg.eigh(Hd)
        vals, vecs = vals[:a.n_eig], vecs[:, :a.n_eig]
        rec['solver'] = 'dense'
    else:
        op = sla.LinearOperator((dim, dim), dtype=float, matvec=matvec)
        ncv = a.ncv or max(2 * a.n_eig + 1, 40)
        vals, vecs = sla.eigsh(op, k=a.n_eig, which='SA', tol=a.tol, ncv=ncv, v0=rng.standard_normal(dim))
        order = np.argsort(vals); vals, vecs = vals[order], vecs[:, order]
        rec['solver'] = f'eigsh(ncv={ncv})'
    rec['t_eig'] = round(time.time() - t0, 1)
    resid = [float(np.linalg.norm(matvec(vecs[:, n]) - vals[n] * vecs[:, n])) for n in range(len(vals))]
    print(f"  {rec['solver']}: {len(vals)} lowest eigenvalues, max residual {max(resid):.1e}  "
          f"[{rec['t_eig']}s, {peak_gb():.2f} GB]", flush=True)

    t0 = time.time()
    cas = se.casimir_expectations(N, p, masks, vecs)
    irreps = se.sector_irreps(N, p, k)
    rec['t_label'] = round(time.time() - t0, 1)
    groups = []
    for x, c in zip(vals, cas):
        if groups and abs(x - groups[-1]['E']) < 1e-6 * max(1.0, abs(x)):
            groups[-1]['mult'] += 1; groups[-1]['C2'].append(float(c))
        else:
            groups.append(dict(E=float(x), mult=1, C2=[float(c)]))
    print(f"  levels (E, multiplicity in block, <C_2>, irreps with that Casimir: lambda[K0]):")
    for g in groups:
        c = float(np.mean(g['C2']))
        g['C2_mean'], g['C2_spread'] = c, float(np.ptp(g['C2']))
        g['candidates'] = [dict(lam=list(l), mult=m, dim=d, K0=k0) for l, m, cc, d, k0 in irreps if abs(cc - c) < 1e-4]
        cands = ', '.join(f"{tuple(x['lam'])}[{x['K0']}]" for x in g['candidates']) or 'mixed/none'
        print(f"    {g['E']:>14.8f}  x{g['mult']:<3d} <C2> = {c:9.5f} (spread {g['C2_spread']:.1e})  {cands}")
    rec['levels'] = groups
    rec['eigenvalues'] = [float(x) for x in vals]
    rec['residuals'] = resid
    rec['peak_rss_gb'] = round(peak_gb(), 2)
    rec['free_value'] = fs.free_energy(N, k)
    print(f"  D16 free value {rec['free_value']} (k {'<=' if k <= N * N // 4 else '>'} floor(N^2/4)); "
          f"peak RSS {rec['peak_rss_gb']} GB")
    if not a.no_save:
        out = a.out or os.path.join(ROOT, 'results', 'data', f'sector_ed_N{N}_k{k}.json')
        with open(out, 'w') as f:
            json.dump(rec, f, indent=1)
        print(f"  saved {os.path.relpath(out, ROOT)}")


if __name__ == '__main__':
    main()
