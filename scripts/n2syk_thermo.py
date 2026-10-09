#!/usr/bin/env python3
"""
Thermodynamics of the large-N supersymmetric SYK saddle (src/n2syk_sd.py), used for the large-p limit of the quiver
(docs/derivations.md D22).  Per Majorana fermion of the N = 1 model; the N = 2 model at zero R-chemical potential has
twice these values per complex fermion.

  E(T)/N = J/(qh 2^qh) - (J/2) int_0^beta Gt_b G_psi^(qh-1)   (true energy, including the c-number part of H = Q^2)
  log Z/N = (1/2) log 2 - int_0^beta E dbeta',   S = log Z + beta E,   S_0 = (1/2) log 2 - int_0^infty E dbeta
  Schwarzian coupling (FGMS (5.38)-(5.39) convention):  E/N -> 2 pi^2 alpha_s T^2 / J  as T -> 0.
  Fit: E = a2 T^2 + ... + a5 T^5 on beta J >= 50; the spread over other windows/orders is quoted as systematic.
Grids are refined until beta J / M <= --dtau and the energy is Richardson-extrapolated in M (midpoint error ~ dtau^2).

    python scripts/n2syk_thermo.py --qh 3 --out results/data/n2syk_thermo_qh3.json
    python scripts/n2syk_thermo.py --refit results/data/n2syk_thermo_qh3.json      # re-analyse saved rows
"""
import argparse, json, os, sys, time
import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
from n2syk_sd import solve, energy as true_energy, thermodynamics            # noqa: E402


def energy(bJ, qh, dtau, tol):
    M = 2 ** max(10, int(np.ceil(np.log2(bJ / dtau))))
    vals = []
    for MM in (M, 2 * M):
        r = solve(bJ, qh=qh, M=MM, tol=tol)
        vals.append(true_energy(r))
    E_rich = vals[1] + (vals[1] - vals[0]) / 3.0
    return E_rich, vals, M


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--qh', type=int, default=3)
    ap.add_argument('--bmax', type=float, default=500.0)
    ap.add_argument('--dtau', type=float, default=2e-3, help='largest beta J / M')
    ap.add_argument('--tol', type=float, default=1e-13)
    ap.add_argument('--out', default=None)
    ap.add_argument('--refit', default=None, help='re-analyse the rows of an existing output file (no new solves) and update it')
    a = ap.parse_args()
    t0 = time.time()
    if a.refit:
        d = json.load(open(a.refit))
        b = np.array([r['betaJ'] for r in d['rows']]); E = np.array([r['E'] for r in d['rows']])
        qh = d['params']['qh']
        summary = analyse(b, E, qh, d['params']['bmax'])
        print_summary(summary, qh)
        for k in ('fit', 'alpha_s_N1', 'alpha_s_N2', 'S0_N1', 'S0_exact_N1', 'S0_N1_trapezoid', 'integral', 'tail'):
            d.pop(k, None)
        d.update(summary)
        json.dump(d, open(a.refit, 'w'), indent=1)
        return
    qh = a.qh
    betas = np.unique(np.concatenate([np.linspace(0.0, 2.0, 11)[1:], np.geomspace(2.0, a.bmax, 40)]))
    rows = []
    for bJ in betas:
        E, vals, M = energy(bJ, qh, a.dtau, a.tol)
        rows.append(dict(betaJ=float(bJ), E=float(E), E_M=float(vals[0]), E_2M=float(vals[1]), M=int(M)))
        print(f'betaJ={bJ:8.3f}  M=2^{int(np.log2(M))}  E/N={E:.10e}  (M:{vals[0]:.10e}, 2M:{vals[1]:.10e})  [{time.time()-t0:.0f}s]', flush=True)
    b = np.array([r['betaJ'] for r in rows]); E = np.array([r['E'] for r in rows])
    summary = analyse(b, E, qh, a.bmax)
    print_summary(summary, qh)
    print(f'[{time.time()-t0:.0f}s]', flush=True)
    if a.out:
        json.dump(dict(params=vars(a), rows=rows, **summary, time_s=round(time.time() - t0, 1)), open(a.out, 'w'), indent=1)


def fit_alpha(b, E, bmin, order):
    """Least-squares E = sum_{k=2}^{order} a_k T^k on beta J >= bmin; returns (alpha_s per Majorana, coefficients)."""
    sel = b >= bmin
    T = 1.0 / b[sel]
    A = np.vstack([T ** k for k in range(2, order + 1)]).T
    coef, *_ = np.linalg.lstsq(A, E[sel], rcond=None)
    return coef[0] / (2 * np.pi ** 2), coef


def analyse(b, E, qh, bmax):
    """alpha_s from the T^2 coefficient (central fit: beta J >= 50, powers T^2..T^5; systematic = spread over
    beta J >= 20, 30, 50, 80 with powers up to T^4, T^5), and S_0 = (1/2) log 2 - int_0^inf E dbeta."""
    alpha, coef = fit_alpha(b, E, 50, 5)
    alts = [fit_alpha(b, E, bmin, order)[0] for bmin in (20, 30, 50, 80) for order in (4, 5)]
    syst = float(max(abs(x - alpha) for x in alts))
    bb = np.concatenate([[0.0], b]); EE = np.concatenate([[1.0 / (2 ** qh * qh)], E])
    integral = np.trapezoid(EE, bb)                       # kept for reference (quadrature error ~ 1e-3)
    tail = 2 * np.pi ** 2 * alpha / b[-1]                 # int_bmax^inf 2 pi^2 alpha_s / beta^2
    S0 = float(thermodynamics(b, E, qh, alpha)(1e12)[2][0])   # spline quadrature + Schwarzian tail
    S0_trap = 0.5 * np.log(2) - integral - tail
    S0_exact = 0.5 * np.log(2 * np.cos(np.pi / (2 * qh)))
    return dict(alpha_s_N1=float(alpha), alpha_s_N1_syst=syst, alpha_s_N2=float(2 * alpha), fit=[float(c) for c in coef],
                fit_window='beta J >= 50, E = a2 T^2 + ... + a5 T^5', S0_N1=float(S0), S0_exact_N1=float(S0_exact),
                S0_N1_trapezoid=float(S0_trap), integral=float(integral), tail=float(tail))


def print_summary(r, qh):
    print(f"qh={qh}: alpha_s (N=1, per Majorana) = {r['alpha_s_N1']:.7f} +- {r['alpha_s_N1_syst']:.1e} (fit-window spread); "
          f"T^2 coefficient a2 = {r['fit'][0]:.7f}", flush=True)
    print(f"S_0/N = {r['S0_N1']:.7f} (exact {r['S0_exact_N1']:.7f}; spline quadrature; trapezoid {r['S0_N1_trapezoid']:.6f}); "
          f"tail {r['tail']:.2e}", flush=True)
    print(f"N = 2 per complex fermion: alpha_s = {r['alpha_s_N2']:.7f}, S_0 = {2*r['S0_N1']:.7f} (exact {2*r['S0_exact_N1']:.7f})", flush=True)


if __name__ == '__main__':
    main()
