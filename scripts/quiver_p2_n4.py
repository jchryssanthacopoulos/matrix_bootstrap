#!/usr/bin/env python3
"""
Hidden N = 4 supersymmetry of the p = 2 quiver on gauge singlets, and level statistics resolved by it
(research/notes/quiver_project.md section 6).

Observed 2026-10-07 at (n, p) = (2, 2), for every coupling tensor tested (seed-3 integers, random real, random
complex): on the singlet sectors the second cubic supercharge Qt = Q(eps.eps.eps Cbar), eps = [[0,1],[-1,0]], obeys
Qt^2 = 0 and {Q, Qt} = 0 (automatic: creation operators only), {Qt, Qt^dag} = H and {Q^dag, Qt} = 0 to machine
precision.  On the full (non-singlet) zero-weight sectors the last two fail at the 30% level, so the extended algebra
closes only up to gauge transformations.  On singlets the flavour map F (eps on every edge) is a unitary symmetry of
H with F Q F^-1 = Qt and F^2 = (-1)^k; with complex conjugation K (real couplings), Theta = F K is antiunitary with
Theta^2 = (-1)^k.

This script
  (i)   verifies these relations on every singlet sector (dense, n = 2, p = 2);
  (ii)  builds the N = 4 bottom spaces B_k = ker Q^dag ∩ ker Qt^dag (minus BPS) and checks the long-multiplet count
        n(k) = b(k) + 2 b(k-3) + b(k-6) + h(k);
  (iii) resolves each B_k by F (eigenvalues +-1 for even k, +-i for odd k; the +-i blocks are exchanged by K and so
        isospectral) and computes <r> per block against size-matched GOE (even k) or GUE (odd k) references;
  (iv)  recomputes the BPS projected-operator (LMRS) statistic at the centre, with F taken into account.

    python scripts/quiver_p2_n4.py                 # seed-3 integer couplings
    python scripts/quiver_p2_n4.py --couplings complex --seed 7
"""
import argparse, json, os, resource, sys, time
import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'scripts')); sys.path.insert(0, os.path.join(ROOT, 'src'))   # src first: scripts/quiver_index.py is a different module
from quiver_n2 import QuiverN2, sector, op_sym, flavour_eps_map, perm_sym          # noqa: E402
from quiver_index import singlet_series                                             # noqa: E402
from quiver_singlet_spectrum import singlet_basis, r_stats, lmrs_operators              # noqa: E402

EPS = np.array([[0., 1.], [-1., 0.]])


def peak_gb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1e9 if sys.platform == 'darwin' else r / 1e6


def cubic_terms(q, Ct):
    n, p = q.n, q.p
    out = []
    for a in range(p):
        for b in range(p):
            for c in range(p):
                if abs(Ct[a, b, c]) < 1e-15:
                    continue
                for i in range(n):
                    for j in range(n):
                        for kk in range(n):
                            out.append((Ct[a, b, c], [(q.mode(0, a, i, j), True), (q.mode(1, b, j, kk), True),
                                                      (q.mode(2, c, kk, i), True)]))
    return out


def null_space(M, tol=1e-9):
    if M.shape[0] == 0:
        return np.eye(M.shape[1], dtype=M.dtype)
    G = M.conj().T @ M
    w, V = np.linalg.eigh((G + G.conj().T) / 2)
    return V[:, w < tol * max(1.0, w.max())]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('--couplings', choices=['integer', 'real', 'complex'], default='integer',
                    help="integer: the seed-3 convention of all earlier quiver runs; real/complex: Gaussian, from --seed")
    ap.add_argument('--trim', type=float, default=0.1)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    t0 = time.time()
    n, p = 2, 2
    q = QuiverN2(n, p, a.seed if a.couplings == 'integer' else 3)
    rng = np.random.default_rng(a.seed)
    if a.couplings == 'integer':
        C = q.C.astype(float)
    elif a.couplings == 'real':
        C = rng.standard_normal((2, 2, 2))
    else:
        C = rng.standard_normal((2, 2, 2)) + 1j * rng.standard_normal((2, 2, 2))
    Ct = np.einsum('ax,by,cz,xyz->abc', EPS, EPS, EPS, C.conj())
    cplx = np.iscomplexobj(C)
    nk, _, _ = singlet_series(n, p)
    ks = list(range(0, q.nm + 1, 3))
    S = {k: sector(q, k) for k in ks}
    U = {k: singlet_basis(q, S[k], k)[0] for k in ks}
    assert all(U[k].shape[1] == nk[k] for k in ks)
    tQ, tQt = cubic_terms(q, C), cubic_terms(q, Ct)
    Qs, Qts = {}, {}
    for k in ks[:-1]:
        Qs[k] = U[k + 3].T @ (op_sym(S[k]['masks'], S[k]['sb'], S[k + 3]['sb'], tQ) @ U[k])
        Qts[k] = U[k + 3].T @ (op_sym(S[k]['masks'], S[k]['sb'], S[k + 3]['sb'], tQt) @ U[k])
    sig, sgn, inv = flavour_eps_map(q)
    Fs = {k: U[k].T @ (perm_sym(S[k]['masks'], S[k]['sb'], S[k]['sb'], sig, sgn, inv, q.nm) @ U[k]) for k in ks}
    def H(k):
        h = np.zeros((nk[k], nk[k]), dtype=complex if cplx else float)
        if k in Qs:
            h += Qs[k].conj().T @ Qs[k]
        if k - 3 in Qs:
            h += Qs[k - 3] @ Qs[k - 3].conj().T
        return h
    def Ht(k):
        h = np.zeros((nk[k], nk[k]), dtype=complex if cplx else float)
        if k in Qts:
            h += Qts[k].conj().T @ Qts[k]
        if k - 3 in Qts:
            h += Qts[k - 3] @ Qts[k - 3].conj().T
        return h
    Hs = {k: H(k) for k in ks}
    rec = dict(params=vars(a), couplings=np.round(C, 10).tolist() if not cplx else [str(x) for x in C.ravel()],
               checks={}, bottoms={}, lmrs={})
    print(f"(n,p)=(2,2), {a.couplings} couplings: checks of the hidden N = 4 algebra on singlets")
    for k in ks:
        Hk = Hs[k]; nh = max(np.linalg.norm(Hk), 1e-300)
        c = dict(Htilde_minus_H=float(np.linalg.norm(Ht(k) - Hk) / nh),
                 F_H_commutator=float(np.linalg.norm(Fs[k] @ Hk - Hk @ Fs[k]) / nh),
                 F_squared=float(np.linalg.norm(Fs[k] @ Fs[k] - (-1) ** k * np.eye(nk[k]))),
                 F_unitary=float(np.linalg.norm(Fs[k].T @ Fs[k] - np.eye(nk[k]))))
        X = np.zeros_like(Hk)
        if k in Qs:
            X = X + Qs[k].conj().T @ Qts[k]
        if k - 3 in Qs:
            X = X + Qts[k - 3] @ Qs[k - 3].conj().T
        c['Qdag_Qt_anticommutator'] = float(np.linalg.norm(X) / nh)
        if k in Qs:
            c['FQF_minus_Qt'] = float(np.linalg.norm(Fs[k + 3] @ Qs[k] @ Fs[k].T - (Qts[k] if not cplx else
                                      Fs[k + 3] @ Qs[k] @ Fs[k].T)) / max(np.linalg.norm(Qs[k]), 1e-300)) if not cplx else None
        rec['checks'][k] = c
        print(f"  k={k:2d}: |Ht-H| {c['Htilde_minus_H']:.1e}, |{{Q^dag,Qt}}| {c['Qdag_Qt_anticommutator']:.1e}, "
              f"|[F,H]| {c['F_H_commutator']:.1e}, |F^2-(-1)^k| {c['F_squared']:.1e}"
              + (f", |FQF^-1 - Qt| {c['FQF_minus_Qt']:.1e}" if c.get('FQF_minus_Qt') is not None else ""))
    # (ii)-(iii) bottom spaces
    print("N = 4 bottom spaces B_k = ker Q^dag ∩ ker Qt^dag (non-BPS), resolved by F:")
    bcount, hcount = {}, {}
    for k in ks:
        if k - 3 in Qs:
            Z = null_space(np.vstack([Qs[k - 3].conj().T, Qts[k - 3].conj().T]))
        else:
            Z = np.eye(nk[k])
        HB = Z.conj().T @ Hs[k] @ Z
        e, W = np.linalg.eigh((HB + HB.conj().T) / 2)
        thr = 1e-9 * max(1.0, abs(e).max() if len(e) else 1.0)
        hcount[k] = int((e < thr).sum())
        N = Z @ W[:, e >= thr]
        bcount[k] = N.shape[1]
        rec['bottoms'][k] = dict(dim=int(N.shape[1]), bps=hcount[k])
        if N.shape[1] == 0:
            continue
        Fb = N.conj().T @ Fs[k] @ N
        leak = float(np.linalg.norm(Fs[k] @ N - N @ Fb))
        info = dict(F_leak=leak)
        blocks = {}
        if not cplx:
            # real couplings: F is a unitary symmetry of H on singlets -> resolve by F
            G = Fb if k % 2 == 0 else 1j * Fb
            g, Wg = np.linalg.eigh((G + G.conj().T) / 2)
            for sgn_, lab in ((1, '+'), (-1, '-')):
                sel = np.abs(g - sgn_) < 1e-6
                if sel.any():
                    Hb = (N @ Wg[:, sel]).conj().T @ Hs[k] @ (N @ Wg[:, sel])
                    blocks[lab] = np.sort(np.linalg.eigvalsh((Hb + Hb.conj().T) / 2))
            lab_F = ('F = +1', 'F = -1') if k % 2 == 0 else ('F = +i', 'F = -i')
            info['block_dims'] = {kk: int(len(v)) for kk, v in blocks.items()}
            if k % 2 == 1 and '+' in blocks and '-' in blocks and len(blocks['+']) == len(blocks['-']):
                info['odd_k_blocks_isospectral'] = float(np.abs(blocks['+'] - blocks['-']).max())
            msg = f"  k={k:2d}: bottom dim {N.shape[1]} ({lab_F[0]}: {len(blocks.get('+', []))}, {lab_F[1]}: {len(blocks.get('-', []))})"
            if 'odd_k_blocks_isospectral' in info:
                msg += f", +-i blocks isospectral to {info['odd_k_blocks_isospectral']:.1e}"
        else:
            # complex couplings: F is not a symmetry (only Theta = K F, antiunitary, Theta^2 = (-1)^k)
            HB2 = N.conj().T @ Hs[k] @ N
            lev = np.sort(np.linalg.eigvalsh((HB2 + HB2.conj().T) / 2))
            if k % 2 == 1:
                dbl = float(np.abs(lev[0::2] - lev[1::2]).max()) if len(lev) % 2 == 0 else float('nan')
                info['kramers_doubling'] = dbl
                lev = lev[0::2]
            blocks['+'] = lev
            msg = f"  k={k:2d}: bottom dim {N.shape[1]} (complex couplings; F leak {leak:.1e} = F not a symmetry)" + \
                  (f", Kramers-doubled to {info['kramers_doubling']:.1e}" if k % 2 == 1 else "")
        for lab, lev in blocks.items():
            st = r_stats(lev, a.trim) if len(lev) >= 12 else None
            info[f'levels_{lab}'] = lev.tolist(); info[f'stats_{lab}'] = st
            if st and (k % 2 == 0 or lab == '+'):
                msg += f"; <r>[{lab}] = {st['mean_r']:.3f}+-{st['sem_r']:.3f} (n={len(lev)}, near-deg {st['near_degenerate']})"
        print(msg + ("" if cplx else f"; F leak {leak:.1e}"))
        rec['bottoms'][k].update(info)
    ok = True
    for k in ks:
        pred = bcount.get(k, 0) + 2 * bcount.get(k - 3, 0) + bcount.get(k - 6, 0) + hcount.get(k, 0)
        ok &= pred == nk[k]
    print(f"  long-multiplet count n(k) = b(k) + 2 b(k-3) + b(k-6) + h(k) holds in every degree: {ok}")
    rec['multiplet_count_identity'] = bool(ok)
    # size-matched references for the largest blocks
    refs = {}
    for k in ks:
        for lab in ('+', '-'):
            lev = rec['bottoms'][k].get(f'levels_{lab}')
            if lev and len(lev) >= 20 and (k % 2 == 0 or lab == '+'):
                kinds = ('goe', 'poisson') if k % 2 == 0 else (('gse', 'poisson') if cplx else ('gue', 'poisson'))
                refs[f"{k}{lab}"] = {kd: ref_r(len(lev), kd, a.trim) for kd in kinds}
                print(f"  references for k={k} block {lab} ({len(lev)} levels): " +
                      ", ".join(f"{kd} {v['mean']:.3f}+-{v['std']:.3f}" for kd, v in refs[f'{k}{lab}'].items()))
    rec['references'] = refs
    # (iv) LMRS at the centre with F
    kc = 12
    Hc = Hs[kc]
    e, W = np.linalg.eigh((Hc + Hc.conj().T) / 2)
    V = W[:, e < 1e-9 * max(1.0, e.max())]
    FB = V.conj().T @ Fs[kc] @ V
    gB = np.linalg.eigvalsh((FB + FB.conj().T) / 2)
    print(f"BPS singlets at k={kc}: {V.shape[1]}; F eigenvalues on them: +1 x {int((gB > 0).sum())}, -1 x {int((gB < 0).sum())}")
    rec['lmrs']['n_bps'] = int(V.shape[1]); rec['lmrs']['F_plus'] = int((gB > 0).sum()); rec['lmrs']['F_minus'] = int((gB < 0).sum())
    rec['lmrs']['operators'] = []
    for name, terms in lmrs_operators(q):
        O = U[kc].T @ (op_sym(S[kc]['masks'], S[kc]['sb'], S[kc]['sb'], terms) @ U[kc])
        Oh = V.conj().T @ O @ V
        Oh = (Oh + Oh.conj().T) / 2
        comm = float(np.linalg.norm(Oh @ FB - FB @ Oh) / max(np.linalg.norm(Oh), 1e-300))
        acomm = float(np.linalg.norm(Oh @ FB + FB @ Oh) / max(np.linalg.norm(Oh), 1e-300))
        lam = np.linalg.eigvalsh(Oh)
        if cplx or (comm > 1e-8 and acomm > 1e-8):
            treat, seqs = 'generic', [lam]
        elif comm <= 1e-8:
            gF, WF = np.linalg.eigh((FB + FB.conj().T) / 2)
            seqs = []
            for sg in (1, -1):
                Wsel = WF[:, np.abs(gF - sg) < 1e-6]
                Ob = Wsel.conj().T @ Oh @ Wsel
                seqs.append(np.linalg.eigvalsh((Ob + Ob.conj().T) / 2))
            treat = 'F-commuting: per F block'
        else:
            treat, seqs = 'F-anticommuting: spectrum symmetric, positive half', [np.sort(lam[lam > 1e-12])]
        parts = []
        for sq in seqs:
            st = r_stats(sq, a.trim) if len(sq) >= 12 else None
            ref = ref_r(len(sq), 'goe', a.trim) if st else None
            parts.append(dict(n=int(len(sq)), stats=st, ref_goe=ref))
        txt = "; ".join(f"n={pp['n']}: <r> = {pp['stats']['mean_r']:.3f}+-{pp['stats']['sem_r']:.3f} (GOE {pp['ref_goe']['mean']:.3f}+-{pp['ref_goe']['std']:.3f})"
                        for pp in parts if pp['stats'])
        print(f"  {name:40s}: [{treat}] {txt}")
        st = parts[0]['stats']
        rec['lmrs']['operators'].append(dict(name=name, F_commutator=comm, F_anticommutator=acomm, treatment=treat,
                                             eigenvalues=lam.tolist(), parts=parts))
    rec.update(time_s=round(time.time() - t0, 1), peak_rss_gb=round(peak_gb(), 2))
    out = a.out or os.path.join(ROOT, 'results', 'data', f'quiver_p2_n4_{a.couplings}_seed{a.seed}.json')
    json.dump(rec, open(out, 'w'), indent=1, default=float)
    print(f"saved {os.path.relpath(out, ROOT)}  [{rec['time_s']}s, peak {rec['peak_rss_gb']} GB]")


def ref_r(n, kind, trim):
    """Size-matched <r> for eigenvalues of n x n GOE / GUE matrices, or Poisson."""
    rng = np.random.default_rng(0)
    samples = max(20, int(40000 / max(n, 1)))
    from quiver_singlet_spectrum import r_values
    out = []
    for _ in range(samples):
        if kind == 'poisson':
            lev = np.sort(rng.uniform(size=n))
        elif kind == 'goe':
            X = rng.standard_normal((n, n)); lev = np.linalg.eigvalsh((X + X.T) / 2)
        elif kind == 'gse':
            A_ = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)); A_ = (A_ + A_.conj().T) / 2
            B_ = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)); B_ = (B_ - B_.T) / 2
            M_ = np.block([[A_, B_], [-B_.conj(), A_.conj()]])
            lev = np.sort(np.linalg.eigvalsh(M_))[0::2]            # quaternion self-dual: Kramers pairs, keep one
        else:
            X = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)); lev = np.linalg.eigvalsh((X + X.conj().T) / 2)
        r, _ = r_values(lev, trim)
        out.append(np.mean(r))
    return dict(kind=kind, n=n, samples=samples, mean=float(np.mean(out)), std=float(np.std(out)))


if __name__ == '__main__':
    main()
