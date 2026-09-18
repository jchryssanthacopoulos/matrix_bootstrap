"""Refined Witten index per (N_Psi mod 3, U(N) irrep, Z_p flavor charge) for the fermionic matrix models (src/bps_index.py).

    cd src && PYTHONPATH=. ../.venv/bin/python ../scripts/compute_bps_index.py --N 3 --p 3 --out ../results/data/refined_index_N3_p3.json
Prints class totals (with the closed form), per-class statistics, and writes {"c|w|lambda": I} to --out.  --check_ed (N=2, p=3)
compares with the BPS multiplets per (k, spin, flavor) from exact diagonalisation (index saturation test).
"""
import argparse, json, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
import bps_index as bi


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--N', type=int, required=True); ap.add_argument('--p', type=int, default=3)
    ap.add_argument('--out', default=None); ap.add_argument('--check_ed', action='store_true')
    a = ap.parse_args()
    mult = bi.irrep_multiplicities(a.N, a.p, verbose=True)
    idx = bi.refined_index(mult, a.p)
    print(f"class totals sum_R I dim R: {[sum(v*bi.irrep_dim(a.N, lam) for (c, w, lam), v in idx.items() if c == cc) for cc in range(3)]}"
          f"   closed form {bi.class_totals(a.N, a.p)}")
    for c in range(3):
        cells = {key: v for key, v in idx.items() if key[0] == c}
        if not cells:
            print(f"  class {c}: no complex with non-zero index"); continue
        pos = sum(1 for v in cells.values() if v > 0); neg = len(cells) - pos
        print(f"  class {c}: {len(cells)} complexes ({pos} +, {neg} -), max |I| {max(abs(v) for v in cells.values())}, "
              f"singlet I (w=0) {cells.get((c, 0, tuple([0]*a.N)), 0)}")
    if a.out:
        json.dump({f"{c}|{w}|{','.join(map(str, lam))}": v for (c, w, lam), v in idx.items()}, open(a.out, 'w'), indent=0)
        print("wrote", a.out)
    if a.check_ed:
        from fermion_matrix_model import build_model, chen_C
        from symmetry_reduction import IsotypicReducer
        m = build_model(a.N, a.p, chen_C(a.p) if a.p > 1 else np.ones((1, 1, 1)), fourier=(a.p > 1))
        bad = 0; n_bps = 0
        for k in range(a.p * a.N * a.N + 1):
            ix = m['sectors'][k]; Hk = m['H'][ix][:, ix].toarray(); lam, U = np.linalg.eigh(Hk); ker = U[:, lam < 1e-8]
            if not ker.shape[1]:
                continue
            R = IsotypicReducer(m, k)
            for b in R.blocks:
                j = b['label'][0]; fe = b['label'][1:]
                wq = int(round(np.angle(complex(fe[0], fe[1])) / (2 * np.pi / a.p))) % a.p if len(fe) == 2 else 0
                P = sum(V @ V.conj().T for V in b['V'])
                nb = np.trace(P @ ker @ ker.conj().T).real / b['dim']
                if nb > 0.5:
                    I = idx.get((k % 3, wq, (j, -j)), 0); n_bps += 1
                    if abs(abs(I) - round(nb)) > 1e-6:
                        bad += 1; print(f"  NOT saturated: k={k} j={j} w={wq}: BPS multiplets {nb:.2f}, index {I}")
        print(f"ED check: {n_bps} BPS cells, {bad} unsaturated")


if __name__ == '__main__':
    main()
