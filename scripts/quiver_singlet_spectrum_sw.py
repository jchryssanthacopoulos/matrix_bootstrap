#!/usr/bin/env python3
"""
Singlet spectra of the U(n)^3 quiver in the Schur-Weyl singlet basis (src/quiver_singlet_basis.py; research/notes/
quiver_project.md section 11).  Targets the next-rank laboratory (n,p) = (3,2), validated at n = 2 against every
stored pair spectrum.

For each degree k = 3m up to half filling it records
  * the singlet dimension (checked against the exact count n(k));
  * BPS states: dense degrees -> exact zero modes; Lanczos degrees -> the lowest eigenvalue (positive = none);
    together with particle-hole symmetry and the index this decides concentration at half filling;
  * the lowest singlet levels, classified as lower (Q^+_{m-1} psi = 0, pair (k, k+3)) or upper (Q_m psi = 0,
    pair (k-3, k)) members, giving the pair edges E_0(q), q = k + 3/2 - 3pn^2/2;
  * full pair spectra where Q_m is dense-affordable (level statistics);
  * optionally (--n4-test, p = 2): the hidden N = 4 relations for Q~ = Q((eps x eps x eps) conj(C)) on singlets.

    scripts/run_guarded.sh 8 log python3 scripts/quiver_singlet_spectrum_sw.py --n 3 --p 2 --seed 3 --n4-test
"""
import argparse, json, os, resource, sys, time
import numpy as np
import scipy.sparse.linalg as sla

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
from quiver_singlet_basis import SingletComplex          # noqa: E402
from quiver_index import singlet_series_recursive         # noqa: E402


def peak_gb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1e9 if sys.platform == 'darwin' else r / 1e6


def r_stats(levels, trim=0.1):
    E = np.sort(np.asarray(levels))
    E = E[int(trim * len(E)):int((1 - trim) * len(E))]
    s = np.diff(E)
    s = s[s > 1e-12 * max(1.0, abs(E).max())]
    r = np.minimum(s[1:], s[:-1]) / np.maximum(s[1:], s[:-1])
    return dict(mean_r=float(r.mean()), sem_r=float(r.std(ddof=1) / np.sqrt(len(r))), n_ratios=int(len(r)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--n', type=int, default=3)
    ap.add_argument('--p', type=int, default=2)
    ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('--mmax', type=int, default=None, help='largest m (default: the last degree below half filling)')
    ap.add_argument('--nev', type=int, default=12)
    ap.add_argument('--dense-max', type=float, default=6e7, help='dense Q_m if dim(m+1)*dim(m) below this')
    ap.add_argument('--n4-test', action='store_true')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    t0 = time.time()
    n, p = a.n, a.p
    C = np.random.default_rng(a.seed).integers(1, 6, size=(p, p, p)).astype(float)
    sc = SingletComplex(n, p, C)
    ref, _ = singlet_series_recursive(n, p)
    half = 3 * p * n * n / 2
    mmax = a.mmax if a.mmax is not None else int(np.ceil(half / 3)) - 1
    rec = dict(params=vars(a), couplings=C.tolist(), half_filling=half, degrees={}, pairs={}, pair_spectra={})
    dense_Q = {}

    def Qd(m):
        if m not in dense_Q:
            dense_Q[m] = sc.Q_dense(m) if sc.dim(m) * sc.dim(m + 1) <= a.dense_max else None
        return dense_Q[m]

    def Hop(m):
        N = sc.dim(m)

        def mv(v):
            v = np.asarray(v).ravel()
            out = sc.apply_Qdag(m, sc.apply_Q(m, v)) if sc.dim(m + 1) else np.zeros(N)
            if m > 0:
                out = out + sc.apply_Q(m - 1, sc.apply_Qdag(m - 1, v))
            return out
        return sla.LinearOperator((N, N), matvec=mv, dtype=float)

    for m in range(0, mmax + 1):
        k = 3 * m
        N = sc.dim(m)
        t1 = time.time()
        deg = dict(k=k, dim=N, exact_dim=int(ref[k]), dim_ok=bool(N == ref[k]))
        Qm, Qp = Qd(m), (Qd(m - 1) if m > 0 else None)
        if Qm is not None and (m == 0 or Qp is not None) and N <= 8000:
            H = Qm.T @ Qm + (Qp @ Qp.T if Qp is not None else 0)
            w, V = np.linalg.eigh(H)
            tol = 1e-8 * max(1.0, w.max())
            nb = int((w < tol).sum())
            vals, vecs = w[nb:nb + a.nev], V[:, nb:nb + a.nev]
            deg['method'] = 'dense'
        else:
            kk = min(a.nev, N - 2)
            w, V = sla.eigsh(Hop(m), k=kk, which='SA', tol=1e-10, ncv=max(3 * kk, 40))
            o = np.argsort(w); w, V = w[o], V[:, o]
            tol = 1e-8 * max(1.0, abs(w).max())
            nb = int((w < tol).sum())
            vals, vecs = w[nb:], V[:, nb:]
            deg['method'] = 'lanczos'
        deg['bps'] = nb
        deg['lowest_eigenvalue'] = float(w[0])
        low = []
        for i, E in enumerate(vals):
            v = vecs[:, i]
            lo = np.linalg.norm(sc.apply_Qdag(m - 1, v)) if m > 0 else 0.0
            up = np.linalg.norm(sc.apply_Q(m, v)) if sc.dim(m + 1) else 0.0
            s = np.sqrt(E)
            kind = 'lower' if (lo < 1e-5 * s and up > 1e-3 * s) else ('upper' if (up < 1e-5 * s and lo > 1e-3 * s) else 'mixed')
            low.append(dict(E=float(E), kind=kind))
            key = k if kind == 'lower' else (k - 3 if kind == 'upper' else None)
            if key is not None:
                q = key + 1.5 - half
                if key not in rec['pairs'] or E < rec['pairs'][key]['E0']:
                    rec['pairs'][key] = dict(q=q, E0=float(E), from_degree=k)
        deg['lowest'] = low
        rec['degrees'][k] = deg
        print(f"k={k}: dim {N} ({'ok' if deg['dim_ok'] else 'MISMATCH'} vs {ref[k]}), {deg['method']}, BPS {nb}, lowest "
              f"{[round(L['E'], 6) for L in low[:6]]} ({[L['kind'][0] for L in low[:6]]})  [{time.time() - t1:.1f}s, "
              f"{peak_gb():.2f} GB]", flush=True)
        if Qm is not None and Qm.shape[0] and m + 1 <= mmax + 1:
            ev = np.linalg.eigvalsh(Qm.T @ Qm)
            lev = np.sort(ev[ev > 1e-8 * max(1.0, ev.max())])
            entry = dict(q=k + 1.5 - half, n_levels=int(len(lev)), E0=float(lev[0]) if len(lev) else None,
                         levels=lev.tolist())
            if len(lev) >= 40:
                entry['stats'] = r_stats(lev)
            rec['pair_spectra'][k] = entry
            print(f"   pair ({k},{k + 3}) q={k + 1.5 - half}: {len(lev)} multiplets, E0 {entry['E0']}"
                  + (f", <r> {entry['stats']['mean_r']:.4f} +- {entry['stats']['sem_r']:.4f}" if 'stats' in entry else ''),
                  flush=True)
    if a.n4_test and p == 2:
        eps = np.array([[0.0, 1.0], [-1.0, 0.0]])
        Ct = np.einsum('ai,bj,ck,ijk->abc', eps, eps, eps, np.conj(C))
        st = SingletComplex(n, p, Ct)
        rng = np.random.default_rng(7)
        res = {}
        for m in range(1, mmax + 1):
            x = rng.standard_normal(sc.dim(m))
            Hx = sc.apply_Qdag(m, sc.apply_Q(m, x)) + sc.apply_Q(m - 1, sc.apply_Qdag(m - 1, x))
            Htx = st.apply_Qdag(m, st.apply_Q(m, x)) + st.apply_Q(m - 1, st.apply_Qdag(m - 1, x))
            mix = sc.apply_Qdag(m, st.apply_Q(m, x)) + st.apply_Q(m - 1, sc.apply_Qdag(m - 1, x))   # {Q^+, Q~}
            anti = sc.apply_Q(m + 1, st.apply_Q(m, x)) + st.apply_Q(m + 1, sc.apply_Q(m, x)) if sc.dim(m + 2) else np.zeros(1)
            res[3 * m] = dict(H_minus_Ht=float(np.linalg.norm(Hx - Htx) / np.linalg.norm(Hx)),
                              Qdag_Qt=float(np.linalg.norm(mix) / np.linalg.norm(Hx)),
                              Q_Qt=float(np.linalg.norm(anti) / max(np.linalg.norm(Hx), 1e-300)))
            print(f"N=4 test k={3 * m}: |{{Q~,Q~+}} - H|/|H| = {res[3 * m]['H_minus_Ht']:.1e}, |{{Q+,Q~}}|/|H| = "
                  f"{res[3 * m]['Qdag_Qt']:.1e}, |{{Q,Q~}}|/|H| = {res[3 * m]['Q_Qt']:.1e}", flush=True)
        rec['n4_test'] = res
    rec['time_s'] = round(time.time() - t0, 1)
    rec['peak_rss_gb'] = round(peak_gb(), 2)
    out = a.out or os.path.join(ROOT, 'results', 'data', f'quiver_singlet_sw_n{n}_p{p}_seed{a.seed}.json')
    json.dump(rec, open(out, 'w'), indent=1)
    print('pair edges:', {v['q']: round(v['E0'], 6) for kk, v in sorted(rec['pairs'].items())})
    print('saved', os.path.relpath(out, ROOT), f"[{rec['time_s']} s, peak {rec['peak_rss_gb']} GB]")


if __name__ == '__main__':
    main()
