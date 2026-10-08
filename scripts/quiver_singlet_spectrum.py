#!/usr/bin/env python3
"""
Step 1 of the quiver programme (research/notes/quiver_project.md): the singlet near-BPS spectrum and chaos at n = 2.

All computations live in the W-symmetric zero-weight sectors of src/quiver_n2.py.  There the Casimir penalty
P = sum_v J_v^2 commutes with H, vanishes exactly on gauge singlets and is >= 2 on every non-singlet.

Modes
  pairs  Exact singlet bases U_k = ker P_sym (dense), checked against the dual-Cauchy multiplicities n(k), for every
         degree whose W-symmetric dimension is <= --dense-max.  For each such k the nonzero squared singular values of
         Q: S_k -> S_{k+3} restricted to singlets are the energies of the multiplets of the pair (k, k+3)
         (Turiaci-Witten: each pair is an independent random-matrix ensemble), with multiplet charge
         q = k + 3/2 - 3pn^2/2.  Level statistics: mean ratio of consecutive spacings <r>, against size-matched
         references (real/complex Gaussian square matrices for singular values, and Poisson), sampled here.
  edge   Lowest eigenpairs of K = H_sym + mu P at one degree (Lanczos), each classified as singlet (<P> = 0), lower
         member of (k, k+3) (Q^dag psi = 0) or upper member of (k-3, k) (Q psi = 0).  Eigenvalues below 2 mu are exact
         singlet energies.
  lmrs   BPS chaos (the projected-operator statistic of Lin-Maldacena-Rozenberg-Shan as used by Chang-Chen-Sia-Yang
         section 3.4): eigenvalues of P_BPS O P_BPS on the singlet BPS space at the window centre, for simple
         gauge-invariant operators O, against size-matched GOE and Poisson references.

    python scripts/quiver_singlet_spectrum.py pairs --n 2 --p 2
    scripts/run_guarded.sh 12 log python3 scripts/quiver_singlet_spectrum.py edge --n 2 --p 3 --k 15 --mu 1 --nev 12
    python scripts/quiver_singlet_spectrum.py lmrs --n 2 --p 2
"""
import argparse, json, os, resource, sys, time
import numpy as np
import scipy.sparse.linalg as sla

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
from quiver_n2 import QuiverN2, sector, Q_sym, raising_sym, op_sym           # noqa: E402
from quiver_index import singlet_series                                       # noqa: E402


def peak_gb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1e9 if sys.platform == 'darwin' else r / 1e6


# ----------------------------------------------------------------------------------------------- statistics

def r_values(levels, trim=0.1):
    """Ratios of consecutive spacings r_i = min(s_i, s_{i+1}) / max(...) over the central (1 - 2 trim) of the levels."""
    E = np.sort(np.asarray(levels, float))
    lo, hi = np.quantile(E, [trim, 1 - trim]) if trim > 0 else (E[0], E[-1])
    E = E[(E >= lo) & (E <= hi)]
    s = np.diff(E)
    a, b = s[:-1], s[1:]
    mx = np.maximum(a, b)
    ok = mx > 0
    return np.minimum(a, b)[ok] / mx[ok], int((~ok).sum())


def r_stats(levels, trim=0.1):
    r, ndeg = r_values(levels, trim)
    E = np.sort(levels)
    s = np.diff(E)
    tiny = int((s < 1e-9 * max(np.mean(s), 1e-300)).sum()) if len(s) else 0
    return dict(n_levels=int(len(levels)), n_ratios=int(len(r)), mean_r=float(np.mean(r)) if len(r) else float('nan'),
                sem_r=float(np.std(r) / np.sqrt(len(r))) if len(r) > 1 else float('nan'), near_degenerate=tiny)


def reference_r(n, kind, trim=0.1, samples=None, seed=0):
    """Size-matched <r> distribution: 'loe'/'lue' = squared singular values of n x n real/complex Gaussian matrices,
    'goe' = eigenvalues of n x n real symmetric Gaussian matrices, 'poisson' = n sorted uniforms."""
    rng = np.random.default_rng(seed)
    samples = samples or max(20, int(40000 / max(n, 1)))
    out = []
    for _ in range(samples):
        if kind == 'poisson':
            lev = np.sort(rng.uniform(size=n))
        elif kind == 'goe':
            X = rng.standard_normal((n, n)); lev = np.linalg.eigvalsh((X + X.T) / 2)
        elif kind == 'loe':
            X = rng.standard_normal((n, n)); lev = np.linalg.eigvalsh(X.T @ X)
        elif kind == 'lue':
            X = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)); lev = np.linalg.eigvalsh(X.conj().T @ X)
        r, _ = r_values(lev, trim)
        out.append(np.mean(r))
    return dict(kind=kind, n=n, samples=samples, mean=float(np.mean(out)), std=float(np.std(out)))


def unfolded_spacings(levels, trim=0.1, deg=7):
    """Spacings after unfolding with a polynomial fit of the cumulative count (central part only)."""
    E = np.sort(np.asarray(levels, float))
    N = np.arange(1, len(E) + 1)
    x = (E - E.mean()) / E.std()
    coef = np.polyfit(x, N, deg)
    Nf = np.polyval(coef, x)
    lo, hi = int(trim * len(E)), int((1 - trim) * len(E))
    s = np.diff(Nf[lo:hi])
    return s / s.mean()


# ----------------------------------------------------------------------------------------------- singlet bases

def singlet_basis(q, sec, k, tol=1e-6):
    """Orthonormal basis of the singlets inside the W-symmetric zero-weight sector (columns), and P's spectrum gap."""
    sb = sec['sb']
    dim = len(sb['reps'])
    if dim == 0:
        return np.zeros((0, 0)), (np.nan, np.nan)
    P = np.zeros((dim, dim))
    for v in range(3):
        W = raising_sym(q, sec['masks'], sb, v, k)
        P += (W.T @ W).toarray()
    vals, vecs = np.linalg.eigh((P + P.T) / 2)
    sel = vals < tol
    nz = int(sel.sum())
    gap = (float(vals[nz - 1]) if nz else float('nan'), float(vals[nz]) if nz < dim else float('nan'))
    return vecs[:, sel], gap


def charge(k, n, p):
    return k + 1.5 - 1.5 * p * n * n


# ----------------------------------------------------------------------------------------------- modes

def mode_pairs(a):
    t0 = time.time()
    q = QuiverN2(a.n, a.p, a.seed)
    nk, _, _ = singlet_series(a.n, a.p)
    secs, U, rec = {}, {}, dict(params=vars(a), couplings=q.C.tolist(), degrees={}, pairs={})
    ks = list(range(0, q.nm + 1, 3))
    for k in ks:
        secs[k] = sector(q, k)
        dim = len(secs[k]['sb']['reps'])
        if dim <= a.dense_max:
            U[k], gap = singlet_basis(q, secs[k], k)
            ok = U[k].shape[1] == nk[k]
            rec['degrees'][k] = dict(zero_weight=int(len(secs[k]['masks'])), w_sym=dim, singlets=int(U[k].shape[1]),
                                     n_k_index=int(nk[k]), P_gap=gap)
            print(f"k={k:2d}: zero-weight {len(secs[k]['masks'])}, W-sym {dim}, singlets {U[k].shape[1]} "
                  f"(dual Cauchy n(k) = {nk[k]}: {'ok' if ok else 'MISMATCH'}); P gap {gap[0]:.1e} / {gap[1]:.3f}  "
                  f"[{time.time()-t0:.0f}s, {peak_gb():.2f} GB]", flush=True)
            assert ok, "singlet count disagrees with the dual-Cauchy multiplicity"
        else:
            rec['degrees'][k] = dict(zero_weight=int(len(secs[k]['masks'])), w_sym=dim, singlets=None, n_k_index=int(nk[k]))
            print(f"k={k:2d}: W-sym {dim} > dense-max: no dense singlet basis", flush=True)
            # keep only the symmetric data needed as a Q target
    ranks = {}
    for k in ks:
        if k not in U or k + 3 > q.nm or U[k].shape[1] == 0:
            continue
        Qk = Q_sym(q, secs[k]['masks'], secs[k]['sb'], secs[k + 3]['sb'])
        S = Qk @ U[k]
        if k + 3 in U and U[k + 3].shape[1]:
            leak = float(np.linalg.norm(S - U[k + 3] @ (U[k + 3].T @ S)) / max(np.linalg.norm(S), 1e-300))
        else:
            leak = float('nan')
        sig2 = np.linalg.eigvalsh(S.T @ S)
        thr = 1e-9 * max(sig2.max(), 1.0)
        lev = np.sort(sig2[sig2 > thr])
        ranks[k] = len(lev)
        qch = charge(k, a.n, a.p)
        st = r_stats(lev, a.trim) if len(lev) >= 8 else None
        rec['pairs'][f"{k},{k+3}"] = dict(q=qch, n_levels=int(len(lev)), kernel=int((sig2 <= thr).sum()),
                                          E0=float(lev[0]) if len(lev) else None, levels=lev.tolist(),
                                          singlet_leak=leak, stats=st)
        msg = (f"pair ({k},{k+3}) q={qch:+.1f}: {len(lev)} multiplets, E0 = {lev[0]:.6g}" if len(lev) else
               f"pair ({k},{k+3}): no multiplets")
        if st:
            msg += f", <r> = {st['mean_r']:.4f} +- {st['sem_r']:.4f} ({st['n_ratios']} ratios, {st['near_degenerate']} near-degenerate)"
        print(msg + (f"; Q U_k stays in singlets to {leak:.1e}" if leak == leak else ""), flush=True)
    # BPS counts where both neighbouring ranks are known
    for k in ks:
        if k in U and (k in ranks or k + 3 > q.nm) and (k - 3 in ranks or k - 3 < 0):
            h = U[k].shape[1] - ranks.get(k, 0) - ranks.get(k - 3, 0)
            rec['degrees'][k]['bps_singlets'] = int(h)
            if h:
                print(f"  singlet BPS states at k={k}: {h}")
    # size-matched references for the pairs with enough levels
    refs = {}
    for key, pr in rec['pairs'].items():
        if pr['stats'] and pr['n_levels'] >= 40:
            n = pr['n_levels']
            refs[key] = {kind: reference_r(n, kind, a.trim) for kind in ('loe', 'lue', 'poisson')}
            pr['unfolded_spacings'] = unfolded_spacings(pr['levels'], a.trim).tolist()
            print(f"  references for {n} levels: " + ", ".join(
                f"{kd} {v['mean']:.4f}+-{v['std']:.4f}" for kd, v in refs[key].items()))
    rec['references'] = refs
    rec.update(time_s=round(time.time() - t0, 1), peak_rss_gb=round(peak_gb(), 2))
    out = a.out or os.path.join(ROOT, 'results', 'data', f'quiver_singlet_pairs_n{a.n}_p{a.p}.json')
    json.dump(rec, open(out, 'w'), indent=1)
    print(f"saved {os.path.relpath(out, ROOT)}  [{rec['time_s']}s, peak {rec['peak_rss_gb']} GB]")


def mode_edge(a):
    t0 = time.time()
    q = QuiverN2(a.n, a.p, a.seed)
    k = a.k
    S = {kk: sector(q, kk) for kk in (k - 3, k, k + 3)}
    Qu = Q_sym(q, S[k]['masks'], S[k]['sb'], S[k + 3]['sb']) if S[k + 3] else None
    Qd = Q_sym(q, S[k - 3]['masks'], S[k - 3]['sb'], S[k]['sb']) if S[k - 3] else None
    Ws = [raising_sym(q, S[k]['masks'], S[k]['sb'], v, k) for v in range(3)]
    dim = len(S[k]['sb']['reps'])
    print(f"(n,p)=({a.n},{a.p}) k={k} mu={a.mu}: W-sym dim {dim}; operators built  [{time.time()-t0:.0f}s, "
          f"{peak_gb():.2f} GB]", flush=True)
    def parts(x):
        up = Qu @ x if Qu is not None else np.zeros(1)
        dn = Qd.T @ x if Qd is not None else np.zeros(1)
        pw = [W @ x for W in Ws]
        return up, dn, pw
    def Kmv(x):
        y = np.zeros_like(x)
        if Qu is not None:
            y += Qu.T @ (Qu @ x)
        if Qd is not None:
            y += Qd @ (Qd.T @ x)
        for W in Ws:
            y += a.mu * (W.T @ (W @ x))
        return y
    op = sla.LinearOperator((dim, dim), matvec=Kmv, dtype=float)
    vals, vecs = sla.eigsh(op, k=a.nev, which='SA', tol=1e-10, ncv=max(4 * a.nev, 40),
                           v0=np.random.default_rng(a.v0_seed).standard_normal(dim))
    order = np.argsort(vals); vals, vecs = vals[order], vecs[:, order]
    rows = []
    for i in range(len(vals)):
        x = vecs[:, i]
        up, dn, pw = parts(x)
        eu, ed, c = float(up @ up), float(dn @ dn), float(sum(w @ w for w in pw))
        res = float(np.linalg.norm(Kmv(x) - vals[i] * x))
        E = eu + ed
        kind = 'singlet' if c < 1e-8 else f'non-singlet (sum J^2 = {c:.3f})'
        memb = ('lower' if ed < 1e-8 * max(E, 1e-12) else 'upper' if eu < 1e-8 * max(E, 1e-12) else 'mixed')
        pair = (f"({k},{k+3})" if memb == 'lower' else f"({k-3},{k})" if memb == 'upper' else '?')
        qch = charge(k if memb == 'lower' else k - 3, a.n, a.p) if memb in ('lower', 'upper') else None
        rows.append(dict(K=float(vals[i]), E=E, casimir=c, kind=kind, member=memb, pair=pair, q=qch, residual=res))
        print(f"  K = {vals[i]:.9f}  E = {E:.9f}  {kind:28s} {memb:5s} member of {pair}"
              + (f" (q = {qch:+.1f})" if qch is not None else "") + f"  res {res:.1e}")
    rec = dict(params=vars(a), dim=dim, levels=rows, time_s=round(time.time() - t0, 1), peak_rss_gb=round(peak_gb(), 2))
    out = a.out or os.path.join(ROOT, 'results', 'data', f'quiver_singlet_edge_n{a.n}_p{a.p}_k{k}.json')
    json.dump(rec, open(out, 'w'), indent=1)
    print(f"saved {os.path.relpath(out, ROOT)}  [{rec['time_s']}s, peak {rec['peak_rss_gb']} GB]")


def lmrs_operators(q):
    """Simple gauge-invariant operators preserving k: (name, terms) with terms = [(coef, [(mode, dagger), ...])]."""
    n, p = q.n, q.p
    A, B = 0, 1
    ops = []
    # Tr(X^a Xbar^b) = sum_ij c^dag_(X,a,i,j) c_(X,b,i,j)
    def bil(e, a_, b_):
        return [(1.0, [(q.mode(e, a_, i, j), True), (q.mode(e, b_, i, j), False)]) for i in range(n) for j in range(n)]
    ops.append(('N_A0 = Tr(A^0 Abar^0)', bil(A, 0, 0)))
    ops.append(('Tr(A^0 Abar^1) + h.c.', bil(A, 0, 1)))
    ops.append(('Tr(A^0 Abar^1) + Tr(B^1 Bbar^0) + h.c.', bil(A, 0, 1) + bil(B, 1, 0)))
    # Tr(A^a B^b Bbar^c Abar^d) = sum A_ij B_jk (Bbar)_kl (Abar)_li,  (Bbar)_kl = c_(B,l,k), (Abar)_li = c_(A,i,l)
    def quart(a_, b_, c_, d_):
        out = []
        for i in range(n):
            for j in range(n):
                for kk in range(n):
                    for l in range(n):
                        out.append((1.0, [(q.mode(A, a_, i, j), True), (q.mode(B, b_, j, kk), True),
                                          (q.mode(B, c_, l, kk), False), (q.mode(A, d_, i, l), False)]))
        return out
    ops.append(('Tr(A^0 B^0 Bbar^0 Abar^0) + h.c.', quart(0, 0, 0, 0)))
    ops.append(('Tr(A^0 B^1 Bbar^0 Abar^1) + h.c.', quart(0, 1, 0, 1)))
    return ops


def mode_lmrs(a):
    t0 = time.time()
    q = QuiverN2(a.n, a.p, a.seed)
    kc = 3 * a.p * a.n * a.n // 2
    S = {kk: sector(q, kk) for kk in (kc - 3, kc, kc + 3)}
    U, gap = singlet_basis(q, S[kc], kc)
    Qu = Q_sym(q, S[kc]['masks'], S[kc]['sb'], S[kc + 3]['sb'])
    Qd = Q_sym(q, S[kc - 3]['masks'], S[kc - 3]['sb'], S[kc]['sb'])
    X1, X2 = Qu @ U, Qd.T @ U
    Hs = X1.T @ X1 + X2.T @ X2
    ev, evec = np.linalg.eigh((Hs + Hs.T) / 2)
    nb = int((ev < 1e-9 * max(ev.max(), 1)).sum())
    V = U @ evec[:, :nb]                                         # BPS singlets in W-sym coordinates
    print(f"(n,p)=({a.n},{a.p}) centre k={kc}: {U.shape[1]} singlets, {nb} singlet BPS states; next singlet level "
          f"{ev[nb]:.6f}  [{time.time()-t0:.0f}s]", flush=True)
    rec = dict(params=vars(a), k=kc, n_singlets=int(U.shape[1]), n_bps=nb, operators=[])
    refs = {kind: reference_r(nb, kind, a.trim) for kind in ('goe', 'poisson')}
    print("  size-matched references: " + ", ".join(f"{kd} {v['mean']:.4f}+-{v['std']:.4f}" for kd, v in refs.items()))
    for name, terms in lmrs_operators(q):
        O = op_sym(S[kc]['masks'], S[kc]['sb'], S[kc]['sb'], terms)
        Ohat = V.T @ (O @ V)
        Ohat = (Ohat + Ohat.T) / 2
        lam = np.linalg.eigvalsh(Ohat)
        st = r_stats(lam, a.trim)
        # how much of O V stays in the BPS space (O is not BPS-preserving in general)
        OV = O @ V
        frac = float(np.linalg.norm(V.T @ OV) ** 2 / max(np.linalg.norm(OV) ** 2, 1e-300))
        distinct = int(len(np.unique(np.round(lam, 8))))
        print(f"  {name:40s}: {distinct} distinct eigenvalues of {nb}; <r> = {st['mean_r']:.4f} +- {st['sem_r']:.4f} "
              f"({st['near_degenerate']} near-degenerate spacings); BPS weight of O|BPS> {frac:.3f}")
        rec['operators'].append(dict(name=name, eigenvalues=lam.tolist(), stats=st, distinct=distinct, bps_fraction=frac))
    rec['references'] = refs
    rec.update(time_s=round(time.time() - t0, 1), peak_rss_gb=round(peak_gb(), 2))
    out = a.out or os.path.join(ROOT, 'results', 'data', f'quiver_singlet_lmrs_n{a.n}_p{a.p}.json')
    json.dump(rec, open(out, 'w'), indent=1)
    print(f"saved {os.path.relpath(out, ROOT)}  [{rec['time_s']}s, peak {rec['peak_rss_gb']} GB]")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mode', choices=['pairs', 'edge', 'lmrs'])
    ap.add_argument('--n', type=int, default=2); ap.add_argument('--p', type=int, default=2)
    ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('--k', type=int, default=None)
    ap.add_argument('--mu', type=float, default=1.0); ap.add_argument('--nev', type=int, default=8)
    ap.add_argument('--v0-seed', type=int, default=1)
    ap.add_argument('--dense-max', type=int, default=4000)
    ap.add_argument('--trim', type=float, default=0.1)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    {'pairs': mode_pairs, 'edge': mode_edge, 'lmrs': mode_lmrs}[a.mode](a)


if __name__ == '__main__':
    main()
