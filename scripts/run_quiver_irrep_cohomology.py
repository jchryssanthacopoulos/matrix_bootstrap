"""Exact BPS cohomology per U(n)^3 irrep for the p-flavour quiver Q = sum_{abc} C_abc Tr(A^a B^b C^c).

Only the weight spaces lambda + rho - w rho (w in W) needed for each target irrep lambda are materialised, found by
vectorised enumeration of all 2^(3 p n^2) occupation masks.  Irrep multiplicities per degree come from the
Weyl-denominator formula; ranks are exact over two primes.  Targets: the singlet and the irreps whose refined index
is largest (the most macroscopic complexes), plus optional extra highest weights.

    python scripts/run_quiver_irrep_cohomology.py --n 2 --p 2 --top 6 --seed 3
"""
import argparse, itertools, os, sys
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cohomology import rank_mod_p_blocked, SMALL_PRIMES  # noqa: E402
from multitrace import apply_omega, _add_term  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=2); ap.add_argument("--p", type=int, default=2)
    ap.add_argument("--top", type=int, default=6); ap.add_argument("--seed", type=int, default=3)
    a = ap.parse_args(); n, P = a.n, a.p
    nm = 3 * P * n * n
    # modes: edge e (0:A,1:B,2:C), flavour f, fundamental x, antifundamental y
    mode = lambda e, f, x, y: ((e * P + f) * n + x) * n + y
    W = np.zeros((nm, 3 * n), dtype=np.int8)
    for e in range(3):
        for f in range(P):
            for x in range(n):
                for y in range(n):
                    W[mode(e, f, x, y), e * n + x] += 1
                    W[mode(e, f, x, y), ((e + 1) % 3) * n + y] -= 1
    rng = np.random.default_rng(a.seed)
    Cc = rng.integers(1, 6, size=(P, P, P))
    om = {}
    for fa, fb, fc in itertools.product(range(P), repeat=3):
        for i, j, k in itertools.product(range(n), repeat=3):
            _add_term(om, (mode(0, fa, i, j), mode(1, fb, j, k), mode(2, fc, k, i)), int(Cc[fa, fb, fc]))

    # --- per-irrep index to pick targets (Fock character)
    Z = {(0, tuple([0] * (3 * n))): 1}
    for i in range(nm):
        w = tuple(int(x) for x in W[i]); new = defaultdict(int)
        for (k, v), c in Z.items():
            new[(k, v)] += c; new[(k + 1, tuple(p_ + q_ for p_, q_ in zip(v, w)))] += c
        Z = new
    rho = tuple(n - 1 - i for i in range(n)); perms = list(itertools.permutations(range(n)))
    def sgn(p):
        s = 1
        for i in range(n):
            for j in range(i + 1, n):
                if p[i] > p[j]: s = -s
        return s
    shifts = []
    for p3 in itertools.product(perms, repeat=3):
        eps = 1; sh = [0] * (3 * n)
        for e, p in enumerate(p3):
            eps *= sgn(p)
            for i in range(n): sh[e * n + i] = rho[i] - rho[p[i]]
        shifts.append((eps, tuple(sh)))
    dom = {v for (k, v) in Z if all(all(v[e*n+i] >= v[e*n+i+1] for i in range(n - 1)) for e in range(3))}
    index = {}
    for lam in dom:
        I = [0, 0, 0]
        for k in range(nm + 1):
            m = sum(eps * Z.get((k, tuple(l + s for l, s in zip(lam, sh))), 0) for eps, sh in shifts)
            I[k % 3] += (-1) ** ((k - k % 3) // 3) * m
        if any(I): index[lam] = I
    sing = tuple([0] * (3 * n))
    targets = [sing] + [l for l, _ in sorted(index.items(), key=lambda kv: -max(abs(x) for x in kv[1])) if l != sing][:a.top]
    need = {tuple(l + s for l, s in zip(lam, sh)) for lam in targets for _, sh in shifts}

    # --- enumerate masks, keep only needed weights
    masks = np.arange(1 << nm, dtype=np.uint32)
    wv = np.zeros((1 << nm, 3 * n), dtype=np.int8); kk = np.zeros(1 << nm, dtype=np.int8)
    for i in range(nm):
        bit = ((masks >> i) & 1).astype(np.int8)
        kk += bit
        wv += bit[:, None] * W[i][None, :]
    keep = np.zeros(1 << nm, dtype=bool)
    needarr = np.array(sorted(need), dtype=np.int8)
    for t in needarr:
        keep |= np.all(wv == t[None, :], axis=1)
    sel = np.nonzero(keep)[0]
    blocks = defaultdict(list)
    for mk in sel:
        s = tuple(i for i in range(nm) if (int(mk) >> i) & 1)
        blocks[(len(s), tuple(int(x) for x in wv[mk]))].append(s)
    del masks, wv, kk, keep
    print("(n,p)=(%d,%d): %d modes; materialised %d states in %d blocks, largest %d" %
          (n, P, nm, sum(len(b) for b in blocks.values()), len(blocks), max(len(b) for b in blocks.values())))

    big = max((len(blocks.get((k + 3, w), [])) * len(b) for (k, w), b in blocks.items()), default=0)
    print("largest Q-matrix: %.2f GB as float64" % (big * 8 / 1e9)); sys.stdout.flush()
    rk = {}
    for (k, w), b in blocks.items():
        tgt = blocks.get((k + 3, w))
        if not tgt: rk[(k, w)] = 0; continue
        idx = {s: r for r, s in enumerate(tgt)}; ent = {}
        for c, s in enumerate(b):
            for s2, v in apply_omega(om, s).items():
                if s2 in idx: ent[(idx[s2], c)] = ent.get((idx[s2], c), 0) + v
        rs = []
        for PR in SMALL_PRIMES:                         # in-place blocked LU on float64: one matrix, no temporaries
            M = np.zeros((len(tgt), len(b)), dtype=np.float64)
            for (r_, c_), v in ent.items(): M[r_, c_] = v % PR
            rs.append(rank_mod_p_blocked(M, PR, copy=False)); del M
        assert rs[0] == rs[1]; rk[(k, w)] = rs[0]
    H = {key: len(b) - rk[key] - rk.get((key[0] - 3, key[1]), 0) for key, b in blocks.items()}
    print("%-26s %-16s %-10s %s" % ("irrep A|B|C", "index I_0,I_1,I_2", "window", "BPS multiplicity by degree"))
    for lam in targets:
        hk = {}
        for k in range(nm + 1):
            m = sum(eps * H.get((k, tuple(l + s for l, s in zip(lam, sh))), 0) for eps, sh in shifts)
            if m: hk[k] = m
        cls = defaultdict(list)
        for k in hk: cls[k % 3].append(k)
        conc = all(len(v) == 1 for v in cls.values())
        lab = "|".join(",".join(str(x) for x in lam[e*n:(e+1)*n]) for e in range(3))
        ks = sorted(hk)
        print("%-26s %-16s %-10s %s  %s" % (lab, str(index.get(lam, [0,0,0])), ("%d..%d" % (ks[0], ks[-1])) if ks else "-",
              [hk[k] for k in ks], "concentrated" if conc and ks else ("EMPTY" if not ks else "NOT concentrated")))


if __name__ == "__main__":
    main()
