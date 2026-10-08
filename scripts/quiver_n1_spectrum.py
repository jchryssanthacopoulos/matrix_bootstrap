#!/usr/bin/env python3
"""
Large-p test of the quiver family at n = 1 (research/notes/quiver_project.md section 10; src/quiver_n1.py): the
three-species N = 2 SYK model Q = sum C_fgh a_f^+ b_g^+ c_h^+ in its charge-balanced (gauge-singlet) sector, blocks
B_m with N_a = N_b = N_c = m, dim C(p,m)^3.

For each p it records
  * BPS counts per block (dense blocks: exact zero modes; large blocks: the lowest eigenvalue, to show there are none),
    against the index sum_m (-1)^m C(p,m)^3 (concentration test);
  * the lowest multiplet energy E_0 of each Q-pair (m, m+1), labelled by q = 3m + 3/2 - 3p/2, from classified
    low-lying eigenvectors of H_m (lower member: Q_{m-1}^+ psi = 0; upper member: Q_m psi = 0);
  * full multiplet spectra of small pairs (non-zero eigenvalues of Q_m^+ Q_m) for level statistics and the
    Turiaci-Witten counting test.
Couplings: Gaussian, <|C_fgh|^2> = J/p^2 with J = 1 (real by default), so energies are O(p) and the
super-Schwarzian scale of N = 2 SYK is O(1/p) in these units.

    scripts/run_guarded.sh 10 log python3 scripts/quiver_n1_spectrum.py --p 2 3 4 5 6 7 8 --seed 1
    scripts/run_guarded.sh 10 log python3 scripts/quiver_n1_spectrum.py --p 10 --blocks 3 4 --nev 16
"""
import argparse, json, math, os, resource, sys, time
import numpy as np
import scipy.sparse.linalg as sla

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
from quiver_n1 import QuiverN1, dixon_index          # noqa: E402


def peak_gb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1e9 if sys.platform == 'darwin' else r / 1e6


def couplings(p, seed, kind):
    rng = np.random.default_rng(seed)
    if kind == 'gauss':
        return rng.standard_normal((p, p, p)) / p
    if kind == 'complex':
        return (rng.standard_normal((p, p, p)) + 1j * rng.standard_normal((p, p, p))) / (p * np.sqrt(2))
    if kind == 'int':
        return rng.integers(1, 6, (p, p, p)).astype(float)
    raise ValueError(kind)


def r_stats(levels, trim=0.1):
    E = np.sort(np.asarray(levels))
    lo, hi = int(trim * len(E)), int((1 - trim) * len(E))
    E = E[lo:hi]
    s = np.diff(E)
    s = s[s > 1e-12 * max(1.0, abs(E).max())]
    r = np.minimum(s[1:], s[:-1]) / np.maximum(s[1:], s[:-1])
    return dict(mean_r=float(r.mean()), sem_r=float(r.std(ddof=1) / np.sqrt(len(r))), n_ratios=int(len(r)))


class Blocks:
    """H_m and Q_m as linear operators, sparse where affordable, matrix-free otherwise."""

    def __init__(self, model, sparse_max):
        self.q = model
        self.sparse_max = sparse_max
        self.cache = {}

    def Q(self, m):
        if m < 0 or m >= self.q.p:
            return None
        if m not in self.cache:
            if self.q.dim(m) * (self.q.p - m) ** 3 <= self.sparse_max:
                self.cache[m] = self.q.Q_sparse(m)
            else:
                self.cache[m] = 'free'
        return self.cache[m]

    def apply_Q(self, m, v):
        Qm = self.Q(m)
        return Qm @ v if not isinstance(Qm, str) else self.q.apply_Q(m, v)

    def apply_Qdag(self, m, w):
        Qm = self.Q(m)
        return Qm.conj().T @ w if not isinstance(Qm, str) else self.q.apply_Qdag(m, w)

    def H(self, m):
        n = self.q.dim(m)
        dt = np.result_type(self.q.C, float)

        def mv(v):
            v = np.asarray(v).ravel()
            out = np.zeros(n, dtype=np.result_type(v, dt))
            if m < self.q.p:
                out += self.apply_Qdag(m, self.apply_Q(m, v))
            if m > 0:
                out += self.apply_Q(m - 1, self.apply_Qdag(m - 1, v))
            return out
        return sla.LinearOperator((n, n), matvec=mv, dtype=dt)

    def dense_H(self, m):
        n = self.q.dim(m)
        H = np.zeros((n, n), dtype=np.result_type(self.q.C, float))
        if m < self.q.p:
            Q = self.Q(m).toarray()
            H += Q.conj().T @ Q
        if m > 0:
            Qp = self.Q(m - 1).toarray()
            H += Qp @ Qp.conj().T
        return H


def classify(blocks, m, vecs, vals, tol=1e-6):
    out = []
    for i, E in enumerate(vals):
        v = vecs[:, i]
        lo = np.linalg.norm(blocks.apply_Qdag(m - 1, v)) if m > 0 else 0.0     # small -> lower member of (m, m+1)
        up = np.linalg.norm(blocks.apply_Q(m, v)) if m < blocks.q.p else 0.0    # small -> upper member of (m-1, m)
        s = np.sqrt(max(E, 1e-300))
        kind = 'bps' if E < tol else ('lower' if lo < 1e-5 * s and up > 1e-3 * s else
                                      ('upper' if up < 1e-5 * s and lo > 1e-3 * s else 'mixed'))
        out.append(dict(E=float(E), kind=kind, res_lower=float(lo), res_upper=float(up)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--p', type=int, nargs='+', default=[2, 3, 4, 5, 6])
    ap.add_argument('--blocks', type=int, nargs='*', default=None, help='restrict to these m (default: all)')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--coupling', choices=['gauss', 'complex', 'int'], default='gauss')
    ap.add_argument('--nev', type=int, default=12)
    ap.add_argument('--dense-max', type=int, default=9300, help='full dense diagonalisation up to this block dimension')
    ap.add_argument('--pair-dense-max', type=int, default=23000, help='dense Q^+Q pair spectra up to this dimension')
    ap.add_argument('--sparse-max', type=float, default=3e7, help='build Q_m explicitly if nnz below this')
    ap.add_argument('--tol', type=float, default=1e-10, help='Lanczos tolerance')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    t0 = time.time()
    rec = dict(params=vars(a), results={})
    if a.selftest:
        for p in (3, 4, 5):
            model = QuiverN1(p, couplings(p, a.seed, a.coupling))
            for m in range(p):
                Qs = model.Q_sparse(m)
                v = np.random.default_rng(0).standard_normal(model.dim(m))
                w = np.random.default_rng(1).standard_normal(model.dim(m + 1))
                d1 = np.abs(Qs @ v - model.apply_Q(m, v)).max()
                d2 = np.abs(Qs.T @ w - model.apply_Qdag(m, w)).max()
                nil = np.abs((model.Q_sparse(m + 1) @ Qs)).max() if m + 1 < p else 0.0
                print(f"selftest p={p} m={m}: matrix-free vs sparse {d1:.1e}, {d2:.1e}; Q^2 = 0 to {nil:.1e}")
        return
    for p in a.p:
        model = QuiverN1(p, couplings(p, a.seed, a.coupling))
        blocks = Blocks(model, a.sparse_max)
        I0 = dixon_index(p)
        res = dict(p=p, index=I0, dims={}, bps={}, lowest={}, pairs={}, pair_spectra={})
        ms = a.blocks if a.blocks is not None else list(range(p + 1))
        print(f"p = {p}: index {I0}, block dims {[model.dim(m) for m in range(p + 1)]}", flush=True)
        for m in ms:
            n = model.dim(m)
            res['dims'][m] = n
            t1 = time.time()
            if n <= a.dense_max:
                H = blocks.dense_H(m)
                w, V = np.linalg.eigh(H)
                tol = 1e-8 * max(1.0, w.max())
                nb = int((w < tol).sum())
                low = classify(blocks, m, V[:, nb:nb + a.nev], w[nb:nb + a.nev])
                res['bps'][m] = nb
            else:
                k = min(a.nev, n - 2)
                w, V = sla.eigsh(blocks.H(m), k=k, which='SA', tol=a.tol, ncv=max(3 * k, 20 if k <= 4 else 40))
                order = np.argsort(w); w, V = w[order], V[:, order]
                tol = 1e-8 * max(1.0, abs(w).max())
                nb = int((w < tol).sum())
                res['bps'][m] = nb if nb < k else f'>= {nb}'
                low = classify(blocks, m, V[:, nb:], w[nb:])
            res['lowest'][m] = low
            for L in low:
                if L['kind'] == 'lower':
                    key = m
                elif L['kind'] == 'upper':
                    key = m - 1
                else:
                    continue
                q = 3 * key + 1.5 - 1.5 * p
                cur = res['pairs'].get(key)
                if cur is None or L['E'] < cur['E0']:
                    res['pairs'][key] = dict(q=q, E0=L['E'], from_block=m)
            print(f"  m = {m}: dim {n}, BPS {res['bps'][m]}, lowest {[round(L['E'], 6) for L in low[:6]]} "
                  f"({[L['kind'][0] for L in low[:6]]})  [{time.time() - t1:.1f}s, {peak_gb():.2f} GB]", flush=True)
        # full multiplet spectra of pairs (m, m+1) where the smaller side is dense-affordable
        for m in ms:
            if m >= p:
                continue
            nsm = model.dim(m)
            if nsm <= a.pair_dense_max and model.dim(m + 1) * nsm < 4e9:
                Q = blocks.Q(m)
                if isinstance(Q, str):
                    continue
                G = (Q.conj().T @ Q).toarray()
                ev = np.linalg.eigvalsh(G)
                lev = np.sort(ev[ev > 1e-8 * max(1.0, ev.max())])
                if len(lev) >= 30:
                    res['pair_spectra'][m] = dict(q=3 * m + 1.5 - 1.5 * p, n_levels=int(len(lev)), E0=float(lev[0]),
                                                  levels=lev.tolist(), stats=r_stats(lev))
                    st = res['pair_spectra'][m]['stats']
                    print(f"  pair ({m},{m + 1}) q = {3 * m + 1.5 - 1.5 * p}: {len(lev)} multiplets, E0 = {lev[0]:.6g}, "
                          f"<r> = {st['mean_r']:.4f} +- {st['sem_r']:.4f}", flush=True)
                del G
        nz = {m: v for m, v in res['bps'].items() if v not in (0,)}
        res['concentration_check'] = dict(nonzero_bps_blocks=nz, index=I0)
        print(f"  BPS by block: {res['bps']}; index {I0}; pair edges "
              f"{ {res['pairs'][k]['q']: round(res['pairs'][k]['E0'], 6) for k in sorted(res['pairs'])} }", flush=True)
        rec['results'][str(p)] = res
        rec['time_s'] = round(time.time() - t0, 1)
        rec['peak_rss_gb'] = round(peak_gb(), 2)
        out = a.out or os.path.join(ROOT, 'results', 'data', f'quiver_n1_spectrum_{a.coupling}_seed{a.seed}.json')
        json.dump(rec, open(out, 'w'), indent=1, default=lambda x: x if not isinstance(x, np.generic) else x.item())
        print('saved', os.path.relpath(out, ROOT), f'[{rec["time_s"]} s, peak {rec["peak_rss_gb"]} GB]', flush=True)


if __name__ == '__main__':
    main()
