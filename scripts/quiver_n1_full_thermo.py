#!/usr/bin/env python3
"""
Exact thermodynamics of the FULL (gauge-unprojected) n = 1 quiver with p flavours, i.e. tripartite N = 2 SYK with
3p complex fermions and Q = sum C_fgh a_f^+ b_g^+ c_h^+, for comparison with the large-N Schwinger-Dyson solution
(docs/derivations.md D22: large p -> N = 2 SYK, qh = 3, N = 3p, J = p^2 <C^2>).

H = {Q, Q^+} is block diagonal in the occupations (x, y, z) of the three species; each block is diagonalised densely.
Block (x,y,z) -> (x+1,y+1,z+1):  a_f^+ b_g^+ c_h^+ |S_a,S_b,S_c> = (-1)^y s(f,S_a) s(g,S_b) s(h,S_c) |...> (c^+ passes
x+y letters, b^+ passes x).  Couplings: real Gaussian, variance 1/p^2 (so J = 1 on average); --coupling int uses the integer ensemble {1..5}.

Recorded besides the thermodynamics (D22):
  * exact infinite-temperature moments <H>/N, Var(H)/N, kappa_3/N (to compare with J/12, the disorder average
    (J^2/16)(1 + 1/p + 1/p^2), and the large-N value J^2/16);
  * zero modes per block, and per flavour sector (x-y, y-z) against that sector's Euler characteristic (index),
    to test the BPS window |J_R| < 3/2 with J_R = x + y + z - 3p/2.

    python scripts/quiver_n1_full_thermo.py --p 6 --seed 1 --out results/data/quiver_n1_full_thermo_p6.json
    python scripts/quiver_n1_full_thermo.py --p 5 --seed 2 --jsonl results/data/quiver_n1_full_thermo.jsonl
"""
import argparse, json, math, os, resource, subprocess, sys, time
import numpy as np
import scipy.sparse as sp

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
from quiver_n1 import creation_matrices                                         # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--p', type=int, default=6)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--coupling', choices=['gauss', 'int'], default='gauss', help='gauss: N(0, 1/p^2); int: uniform {1..5}')
    ap.add_argument('--zero-tol', type=float, default=1e-12, help='eigenvalues below this count as BPS (numerical zeros are ~1e-14; genuine levels in the gapless q = 0 sectors of odd N reach 1e-10)')
    ap.add_argument('--out', default=None, help='write one JSON record')
    ap.add_argument('--jsonl', default=None, help='append one JSON record to this file')
    ap.add_argument('--eigs-out', default=None, help='npz file for the full spectrum (by block)')
    ap.add_argument('--moments-only', action='store_true', help='store only moments and BPS totals (no rows, blocks, sectors)')
    a = ap.parse_args()
    p = a.p
    t0 = time.time()
    rng = np.random.default_rng(a.seed)
    if a.coupling == 'gauss':
        C = rng.standard_normal((p, p, p)) / p
    else:   # the integer ensemble of the exact and bootstrap computations (not zero mean; see D22)
        C = rng.integers(1, 6, (p, p, p)).astype(float)
    J = p ** 2 * float((C ** 2).mean())                  # realised J = p^2 <C^2> (close to 1)
    E = {m: creation_matrices(p, m) for m in range(p)}   # E[m][f]: Lambda^m -> Lambda^{m+1}

    def Qblock(x, y, z):
        """Q: (x,y,z) -> (x+1,y+1,z+1) as a sparse matrix (row-major S_a, S_b, S_c)."""
        Q = None
        for f in range(p):
            for g in range(p):
                for h in range(p):
                    c = C[f, g, h]
                    t = sp.kron(E[x][f], sp.kron(E[y][g], E[z][h], format='csr'), format='csr') * c
                    Q = t if Q is None else Q + t
        return ((-1) ** y) * Q.tocsr()

    eigs, blocks, eig_store = [], [], {}
    for x in range(p + 1):
        for y in range(p + 1):
            for z in range(p + 1):
                D = math.comb(p, x) * math.comb(p, y) * math.comb(p, z)
                H = np.zeros((D, D))
                if x < p and y < p and z < p:
                    Q = Qblock(x, y, z); H += (Q.T @ Q).toarray()
                if x > 0 and y > 0 and z > 0:
                    Q = Qblock(x - 1, y - 1, z - 1); H += (Q @ Q.T).toarray()
                w = np.linalg.eigvalsh(H)
                zero = w < a.zero_tol
                blocks.append(dict(x=x, y=y, z=z, dim=D, n_zero=int(zero.sum()),
                                   max_zero=float(w[zero].max()) if zero.any() else None,
                                   min_nonzero=float(w[~zero].min()) if (~zero).any() else None))
                eigs.append(w)
                if a.eigs_out:
                    eig_store[f'b{x}_{y}_{z}'] = w
    w = np.concatenate(eigs)
    N = 3 * p
    mean, var = float(w.mean()), float(w.var())
    k3 = float(((w - mean) ** 3).mean())
    n_zero = int(sum(b['n_zero'] for b in blocks))
    print(f'p={p} seed={a.seed}: {len(w)} eigenvalues (2^{N}), realised J = {J:.4f}, <H>/N = {mean/N:.5f} (J/12 = {J/12:.5f}), '
          f'Var(H)/N = {var/N:.5f} (avg formula (J^2/16)(1+1/p+1/p^2) = {J**2/16*(1+1/p+1/p**2):.5f}; large N {J**2/16:.5f}), '
          f'#BPS = {n_zero}  [{time.time()-t0:.0f}s]', flush=True)
    worst_zero = max((b['max_zero'] for b in blocks if b['max_zero'] is not None), default=None)
    lowest_nonzero = min(b['min_nonzero'] for b in blocks if b['min_nonzero'] is not None)
    print(f'  separation: largest eigenvalue counted as zero {worst_zero}, smallest nonzero {lowest_nonzero:.3e}', flush=True)

    # flavour sectors (s, t) = (x - y, y - z): chains (z+s+t, z+t, z), graded by J_R = x + y + z - 3p/2 in steps of 3
    by = {(b['x'], b['y'], b['z']): b for b in blocks}
    sectors = []
    for s in range(-p, p + 1):
        for t in range(-p, p + 1):
            chain = [by[(z + s + t, z + t, z)] for z in range(p + 1) if (z + s + t, z + t, z) in by]
            if not chain:
                continue
            chi = sum((-1) ** (b['x'] + b['y'] + b['z']) * b['dim'] for b in chain)
            bps = [(b['x'] + b['y'] + b['z'] - 1.5 * p, b['n_zero']) for b in chain if b['n_zero']]
            sectors.append(dict(s=s, t=t, index=int(chi), n_bps=int(sum(nb for _, nb in bps)),
                                bps_JR=[float(jr) for jr, _ in bps], bps_counts=[int(nb) for _, nb in bps]))
    by_JR = {}
    for b in blocks:
        jr = b['x'] + b['y'] + b['z'] - 1.5 * p
        by_JR[jr] = by_JR.get(jr, 0) + b['n_zero']
    by_JR = {f'{k:+.1f}': v for k, v in sorted(by_JR.items()) if v}
    eq_index = sum(1 for s in sectors if s['n_bps'] == abs(s['index']))
    one_block = sum(1 for s in sectors if len(s['bps_JR']) <= 1)
    in_window = sum(nb for s in sectors for jr, nb in zip(s['bps_JR'], s['bps_counts']) if abs(jr) < 1.5)
    print(f'  BPS by J_R: {by_JR};  sectors with #BPS = |index|: {eq_index}/{len(sectors)};  BPS in a single block: '
          f'{one_block}/{len(sectors)};  BPS with |J_R| < 3/2: {in_window}/{n_zero}', flush=True)

    bJ = np.concatenate([np.linspace(0.05, 2, 40), np.geomspace(2, 200, 60)])
    rows = []
    for b in bJ:
        beta = b / J
        xx = -beta * w                                     # SUSY: E >= 0, BPS states at 0
        lz = np.log(np.exp(xx - xx.max()).sum()) + xx.max()
        Eavg = float((w * np.exp(xx - xx.max())).sum() / np.exp(xx - xx.max()).sum())
        rows.append(dict(betaJ=float(b), logZ_per_N=float(lz / N), E_per_N=float(Eavg / N), S_per_N=float((lz + beta * Eavg) / N)))
    try:
        rev = subprocess.check_output(['git', '-C', ROOT, 'rev-parse', '--short', 'HEAD'], text=True).strip()
    except Exception:
        rev = None
    rec = dict(params=vars(a), p=p, seed=a.seed, N=N, J=J, n_states=int(len(w)), mean_H_per_N=mean / N,
               var_H_per_N=var / N, kappa3_per_N=k3 / N, var_formula_per_N=J ** 2 / 16 * (1 + 1 / p + 1 / p ** 2),
               n_zero=n_zero, bps_by_JR=by_JR, sectors_bps_eq_index=eq_index, sectors_single_block=one_block,
               n_sectors=len(sectors), bps_in_window=in_window, max_zero=worst_zero, min_nonzero=lowest_nonzero,
               blocks=blocks, sectors=sectors, rows=rows, time_s=round(time.time() - t0, 1),
               peak_rss_gb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 2), git=rev,
               date=time.strftime('%Y-%m-%d %H:%M'))
    if a.moments_only:
        for k in ('blocks', 'sectors', 'rows'):
            rec.pop(k)
    if a.out:
        json.dump(rec, open(a.out, 'w'), indent=1)
        print('saved', a.out)
    if a.jsonl:
        with open(a.jsonl, 'a') as f:
            f.write(json.dumps(rec) + '\n')
        print('appended', a.jsonl)
    if a.eigs_out:
        np.savez_compressed(a.eigs_out, **eig_store)


if __name__ == '__main__':
    main()
