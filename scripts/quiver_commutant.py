#!/usr/bin/env python3
"""
Conserved-charge search on the gauge-singlet sector of the n = 2 quiver (research/notes/quiver_project.md section 6).

Span of gauge-invariant, k-preserving operators with at most four fermion operators:
  identity; bilinears Tr(X^a Xbar^b) (X in {A, B, C}); single-trace quartic loops (two arrows forward, two backward
  around the triangle, every start node and order, all flavours); double traces Tr(X^a Xbar^b) Tr(Y^c Ybar^d).
Restricted to the singlet sector at one degree k (dense singlet basis), find the dimension of the span and of its
commutant with H: {O in span : [H, O] = 0}.  Generic expectation: only the identity and H itself (H = {Q, Q^dag} lies
in this span).  Extra solutions are conserved charges (hidden symmetries or integrability).

    scripts/run_guarded.sh 3 log python3 scripts/quiver_commutant.py --p 2 --k 9
"""
import argparse, itertools, json, os, sys, time
import numpy as np
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'scripts')); sys.path.insert(0, os.path.join(ROOT, 'src'))
from quiver_n2 import QuiverN2, sector, op_sym, Q_sym                       # noqa: E402
from quiver_singlet_spectrum import singlet_basis                             # noqa: E402

LET = {'A': (0, True, 0, 1), 'B': (1, True, 1, 2), 'C': (2, True, 2, 0),     # edge, creation, from-node, to-node
       'a': (0, False, 1, 0), 'b': (1, False, 2, 1), 'c': (2, False, 0, 2)}  # barred letters go backwards


def letter_op(q, L, f, r, s):
    """Matrix element (r, s) of a letter as (mode, dagger): X_rs = c^dag_(X,f,r,s); (Xbar)_rs = c_(X,f,s,r)."""
    e, cre, _, _ = LET[L]
    return (q.mode(e, f, r, s), True) if cre else (q.mode(e, f, s, r), False)


def trace_word(q, word, flav):
    """Tr(L_1 ... L_m) as terms [(1.0, [(mode, dag), ...])] summing over internal indices."""
    n = q.n; m = len(word); terms = []
    for idx in itertools.product(range(n), repeat=m):
        ops = [letter_op(q, word[t], flav[t], idx[t], idx[(t + 1) % m]) for t in range(m)]
        terms.append((1.0, ops))
    return terms


def quartic_words():
    out = []
    for start in range(3):
        for pat in itertools.permutations('FFBB'):
            node, w = start, []
            for st in pat:
                if st == 'F':
                    L = 'ABC'[node]; w.append(L); node = LET[L][3]
                else:
                    L = {0: 'c', 1: 'a', 2: 'b'}[node]; w.append(L); node = LET[L][3]
            assert node == start
            out.append(''.join(w))
    return sorted(set(out))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--p', type=int, default=2); ap.add_argument('--k', type=int, default=9)
    ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    t0 = time.time()
    q = QuiverN2(2, a.p, a.seed)
    k = a.k
    S = {kk: sector(q, kk) for kk in (k - 3, k, k + 3)}
    U = singlet_basis(q, S[k], k)[0]
    d = U.shape[1]
    Qu = U_up = None
    Qup = Q_sym(q, S[k]['masks'], S[k]['sb'], S[k + 3]['sb']) @ U
    Qdn = Q_sym(q, S[k - 3]['masks'], S[k - 3]['sb'], S[k]['sb'])
    Hk = Qup.T @ Qup + (U.T @ Qdn) @ (U.T @ Qdn).T
    def restrict(terms):
        return U.T @ (op_sym(S[k]['masks'], S[k]['sb'], S[k]['sb'], terms) @ U)
    ops, labels = [np.eye(d)], ['1']
    bil = {}
    for X in 'ABC':
        e = LET[X][0]
        for f1 in range(a.p):
            for f2 in range(a.p):
                terms = [(1.0, [(q.mode(e, f1, r, s), True), (q.mode(e, f2, r, s), False)]) for r in range(2) for s in range(2)]
                bil[(X, f1, f2)] = restrict(terms); ops.append(bil[(X, f1, f2)]); labels.append(f'Tr({X}{f1} {X.lower()}bar{f2})')
    for w in quartic_words():
        for flav in itertools.product(range(a.p), repeat=4):
            ops.append(restrict(trace_word(q, w, flav))); labels.append(f'Tr[{w}]{flav}')
    keys = list(bil)
    for i in range(len(keys)):
        for j in range(i, len(keys)):
            ops.append(bil[keys[i]] @ bil[keys[j]]); labels.append(f'{keys[i]}*{keys[j]}')
    print(f"(n,p)=(2,{a.p}) k={k}: singlet dim {d}; {len(ops)} operators "
          f"({len(quartic_words())} quartic loop shapes)  [{time.time()-t0:.0f}s]", flush=True)
    M = np.array([o.ravel() for o in ops]).T                 # (d^2, N)
    # basis of the span (operators as vectors)
    G = M.T @ M; w, V = np.linalg.eigh(G)
    keep = w > 1e-9 * w.max()
    Bv = V[:, keep] / np.sqrt(w[keep])                       # columns: orthonormal combinations
    span = int(keep.sum())
    Ob = (M @ Bv)                                             # orthonormal operator basis, (d^2, span)
    # commutator map on the span
    C = np.array([(Hk @ Ob[:, i].reshape(d, d) - Ob[:, i].reshape(d, d) @ Hk).ravel() for i in range(span)]).T
    gc, Vc = np.linalg.eigh(C.T @ C)
    null = int((gc < 1e-9 * gc.max()).sum())
    # is H itself in the span?
    h = Hk.ravel(); hres = np.linalg.norm(h - Ob @ (Ob.T @ h)) / np.linalg.norm(h)
    print(f"  span dimension {span}; commutant dimension {null}; H in span: residual {hres:.1e}")
    print(f"  smallest commutator Gram eigenvalues {np.round(gc[:min(8, len(gc))] / gc.max(), 12).tolist()}")
    rec = dict(params=vars(a), singlet_dim=d, n_ops=len(ops), span=span, commutant=null, H_in_span_residual=float(hres),
               smallest_gram=[float(x / gc.max()) for x in gc[:12]], time_s=round(time.time() - t0, 1))
    if null > 2:
        # characterise the extra conserved operators: spectra on the singlet sector
        Xs = [(Ob @ Vc[:, i]).reshape(d, d) for i in range(null)]
        for i, X in enumerate(Xs):
            Xh = (X + X.T) / 2
            ev = np.linalg.eigvalsh(Xh)
            print(f"  conserved operator {i}: |X - X^T|/|X| = {np.linalg.norm(X - X.T)/np.linalg.norm(X):.1e}, "
                  f"distinct eigenvalues {len(np.unique(np.round(ev, 6)))} of {d}")
    out = a.out or os.path.join(ROOT, 'results', 'data', f'quiver_commutant_n2_p{a.p}_k{k}.json')
    json.dump(rec, open(out, 'w'), indent=1)
    print(f"saved {os.path.relpath(out, ROOT)}  [{rec['time_s']}s]")


if __name__ == '__main__':
    main()
