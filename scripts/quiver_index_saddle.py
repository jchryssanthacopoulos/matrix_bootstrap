#!/usr/bin/env python3
"""
Large-n saddle of the singlet index (t = -1) and singlet count (t = +1) of the U(n)^3 fermionic quiver
(research/notes/quiver_project.md section 8; docs/derivations.md D18; code in src/quiver_saddle.py).

Exact representation (D18):  I_0(n,p) = (-i)^{3pn^2}/(n!)^3 int prod dtheta/(2pi) e^{L(theta)} sgn(theta)^p, a
three-species log gas on the circle (like species repel with weight 2, unlike with weight p).  The count has the cross
kernel ln|2 cos(alpha/2)| (attraction) and no phase.  Hence ln|I_0| = F*(p) n^2 + o(n^2), ln(count) = G*(p) n^2 + o(n^2).

This script
  1. solves the continuum problem (src/quiver_saddle.Continuum) for F*(p) (index: three arcs rotated by 2pi/3) and
     G*(p) (count: three coincident species), with convergence and optimality checks, and compares with the analytic
     limits: F* -> (p/2 - 1) 39 zeta(3)/(2 pi^2) as p -> 2+, F* ~ (3p/2) ln 3 - (3/2) ln(2p/3) - 9/4 as p -> oo;
     G* -> (p - 1) 21 zeta(3)/(2 pi^2) as p -> 1+, G* ~ 3p ln 2 - (3/2) ln(p/2) - 9/4 as p -> oo;
  2. maximises the discrete L for n up to --nmax (symmetry-reduced Newton, started at continuum quantiles) and
     extrapolates L_max - 3 n ln n = F n^2 + b n + c ln n + d + e/n as a cross-check of F*; for n <= --nfull it checks
     the symmetric maximum against unconstrained L-BFGS from random starts;
  3. forms the Laplace (Gaussian-fluctuation) estimate of ln|I_0| and ln(count) for n <= --nlap and compares it, and
     F* n^2, with every exact value in results/data/quiver_singlet_index_dixon.json;
  4. checks the exact consequences of D18 against the exact data (p = 2 closed form, sign (-1)^{pn/2} for even p,
     zeros for odd p and odd n, the Hoelder bound |I_0(n,p)|^2 <= |I_0(n,p-1)| |I_0(n,p+1)| for odd p).

    scripts/run_guarded.sh 3 log python3 scripts/quiver_index_saddle.py --nmax 2048
"""
import argparse, json, math, os, sys, time
import numpy as np
from scipy.special import zeta

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
from quiver_saddle import (Continuum, newton_reduced, reduced_LgH, full_LgH, expand, lbfgs_full, laplace_lnI,  # noqa: E402
                           laplace_corrected, initial_reduced, sign_factor, CRYSTAL_DELTA)

C_ARC = 39 * zeta(3) / (2 * np.pi ** 2)          # index slope at threshold (uniform arcs of length 2pi/3)
C_HALF = 21 * zeta(3) / (2 * np.pi ** 2)         # count slope at threshold (uniform half circle)


def limits(p, kind):
    if kind == 'index':
        return (p / 2 - 1) * C_ARC, 1.5 * p * np.log(3) - 1.5 * np.log(2 * p / 3) - 2.25
    return (p - 1) * C_HALF, 3 * p * np.log(2) - 1.5 * np.log(p / 2) - 2.25


def continuum_scan(kind, ps):
    out = {}
    for p in ps:
        t0 = time.time()
        res = {}
        for K, J in ((40, 400), (60, 700)):
            c = Continuum(p, kind, K, J)
            a, F = c.optimise()
            res[(K, J)] = (a, F)
        hi = (2 * np.pi / 3 - c.a) if kind == 'index' else (np.pi - c.a)
        Ph = c.Phi(np.linspace(-0.999 * c.a, 0.999 * c.a, 61))
        Pg = c.Phi(np.linspace(1.0001 * c.a, 0.995 * hi, 80))
        G = c.G(np.linspace(-0.9999, 0.9999, 4001))
        thr, large = limits(p, kind)
        rec = dict(p=p, a=c.a, F=c.F, lam=c.lam, F_K40=res[(40, 400)][1], dF_resolution=abs(res[(40, 400)][1] - c.F),
                   Phi_spread=float(Ph.max() - Ph.min()), gap_margin=float(Pg.max() - c.lam), min_G=float(G.min()),
                   entropy_term=c.entropy_term(), threshold_form=float(thr), large_p_form=float(large),
                   d=c.d.tolist(), time_s=round(time.time() - t0, 2))
        out[float(p)] = (rec, c)
        print(f"  {kind:5s} p = {p:5}: a* = {c.a:.10f} (arc gap {2*np.pi/3 - 2*c.a if kind == 'index' else np.pi - 2*c.a:.4f}), "
              f"F* = {c.F:.12f} (K40 diff {rec['dF_resolution']:.0e}); Phi spread {rec['Phi_spread']:.0e}, gap margin "
              f"{rec['gap_margin']:.1e}, min G {rec['min_G']:.1e}; threshold form {thr:.5f}, large-p form {large:.5f}", flush=True)
    return out


def discrete_scan(kind, p, ns, cont, nfull, nlap, exact, rng):
    rows = []
    for n in ns:
        t0 = time.time()
        if cont is not None:
            th0 = cont.quantiles(n)
        else:                                    # index p = 2: the 3n-point crystal with arc labelling
            th0 = initial_reduced(n, np.pi / 3, kind)
        L, th, gmax, it = newton_reduced(th0, p, kind)
        row = dict(n=n, L_max=float(L), grad=gmax, newton_iters=it, span=float(th.max() - th.min()))
        if n <= nfull:                           # unconstrained check of the symmetric maximum
            best = -np.inf
            for s in range(6):
                Lf, thf, gf = lbfgs_full(rng.uniform(0, 2 * np.pi, 3 * n), n, p, kind)
                best = max(best, Lf)
            row['L_full_random_best'] = float(best)
            row['symmetric_is_max'] = bool(best <= L + 1e-7 * max(1.0, abs(L)))
        if n <= nlap:
            thF = expand(th, kind)
            Lf, gF, H = full_LgH(thF, n, p, kind)
            lap, logdet, zero = laplace_lnI(n, p, kind, Lf, H)
            row.update(laplace=float(lap), laplace_corr=float(laplace_corrected(n, lap)), logdet=logdet,
                       zero_mode=float(zero), L_full=float(Lf), full_grad=float(np.abs(gF).max()))
            if kind == 'index' and p == 2:          # exact: (3n)!/(n!)^3
                row['exact_closed_form'] = math.lgamma(3 * n + 1) - 3 * math.lgamma(n + 1)
            if kind == 'index' and p % 2 == 1:
                row['sign_at_max'] = sign_factor(thF, n)
        key = f'{n},{p}'
        if key in exact:
            v = exact[key]
            if kind == 'index':
                I = v['index']
                row['exact'] = math.log(abs(I)) if I else None
                row['exact_sign'] = int(np.sign(I))
            elif v.get('singlet_total'):
                row['exact'] = math.log(float(v['singlet_total']))
        row['time_s'] = round(time.time() - t0, 2)
        rows.append(row)
        ex = row.get('exact')
        lp = row.get('laplace_corr')
        print(f"    n = {n:5d}: L_max = {L:16.6f}  (L_max - 3n ln n)/n^2 = {(L - 3 * n * np.log(n)) / n ** 2:10.6f}  |grad| {gmax:.0e} "
              f"[{it} it]" + (f"  sym=max: {row['symmetric_is_max']}" if 'symmetric_is_max' in row else '')
              + (f"  Laplace-corr {lp:10.4f}" if lp is not None else '') + (f"  exact {ex:9.4f}" if ex is not None else '')
              + f"  [{row['time_s']}s]", flush=True)
    return rows


def fit(rows, F_cont, nmin):
    big = [r for r in rows if r['n'] >= nmin]
    if len(big) < 6:
        return None
    N = np.array([r['n'] for r in big], float)
    y = np.array([r['L_max'] - 3 * r['n'] * np.log(r['n']) for r in big])
    A = np.column_stack([N ** 2, N, np.log(N), np.ones_like(N), 1 / N])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    A1 = np.column_stack([N, np.log(N), np.ones_like(N), 1 / N])
    coef1, *_ = np.linalg.lstsq(A1, y - (F_cont * N ** 2 if F_cont is not None else 0), rcond=None)
    Af = np.column_stack([N ** 2, N * np.log(N), N, np.log(N), np.ones_like(N), 1 / N])
    coeff, *_ = np.linalg.lstsq(Af, y + 3 * N * np.log(N), rcond=None)
    return dict(nmin=nmin, F=float(coef[0]), b=float(coef[1]), c=float(coef[2]), d=float(coef[3]), e=float(coef[4]),
                F_minus_continuum=float(coef[0] - F_cont) if F_cont is not None else None,
                with_F_fixed=dict(b=float(coef1[0]), c=float(coef1[1]), d=float(coef1[2]), e=float(coef1[3])),
                free_nlogn=dict(F=float(coeff[0]), nlogn=float(coeff[1])))


def fit_laplace(rows, F_cont, nmin=16):
    pts = [r for r in rows if r['n'] >= nmin and r.get('laplace_corr') is not None and np.isfinite(r['laplace_corr'])]
    if len(pts) < 5 or F_cont is None:
        return None
    N = np.array([r['n'] for r in pts], float)
    y = np.array([r['laplace_corr'] for r in pts]) - F_cont * N ** 2
    A = np.column_stack([N, np.log(N), np.ones_like(N), 1 / N])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return dict(nmin=nmin, nmax=int(N.max()), b=float(coef[0]), c=float(coef[1]), d=float(coef[2]), e=float(coef[3]),
                max_abs_resid=float(np.abs(A @ coef - y).max()))


def structure_checks(exact):
    out = []
    for key, v in sorted(exact.items(), key=lambda kv: tuple(int(x) for x in kv[0].split(','))[::-1]):
        n, p = (int(x) for x in key.split(','))
        I = v['index']
        if p == 2:
            pred = (-1) ** n * math.factorial(3 * n) // math.factorial(n) ** 3
            out.append(dict(n=n, p=p, check='p=2 closed form', ok=(I == pred), value=I, predicted=pred))
        if p % 2 == 0:
            out.append(dict(n=n, p=p, check='sign (-1)^{pn/2}, nonzero', ok=(I != 0 and np.sign(I) == (-1) ** (p * n // 2)), value=I))
        elif n % 2 == 1:
            out.append(dict(n=n, p=p, check='odd p, odd n: zero', ok=(I == 0), value=I))
        elif p == 1:
            out.append(dict(n=n, p=p, check='p = 1 < 2 (uniform saddle, no growth): observed I_0 = 0', ok=(I == 0), value=I))
        else:
            lo, hi = exact.get(f'{n},{p-1}'), exact.get(f'{n},{p+1}')
            if lo and hi and lo['index'] and hi['index'] and p > 1:
                out.append(dict(n=n, p=p, check='Hoelder |I(p)|^2 <= |I(p-1) I(p+1)|',
                                ok=(I * I <= abs(lo['index'] * hi['index'])), value=I,
                                ratio=float(abs(lo['index'] * hi['index']) / (I * I)) if I else None))
            out.append(dict(n=n, p=p, check='odd p, even n: sign + (saddle prediction, not a theorem)', ok=(I > 0), value=I))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--p-index', type=int, nargs='+', default=[2, 3, 4, 5, 6])
    ap.add_argument('--p-count', type=int, nargs='+', default=[2, 3, 4, 5, 6])
    ap.add_argument('--p-cont-index', type=float, nargs='+',
                    default=[2.5, 2.75, 3, 3.5, 4, 5, 6, 8, 10, 12, 16, 24, 32, 64, 128])
    ap.add_argument('--p-cont-count', type=float, nargs='+',
                    default=[1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10, 12, 16, 24, 32, 64, 128])
    ap.add_argument('--nmax', type=int, default=2048)
    ap.add_argument('--nfull', type=int, default=12, help='unconstrained random-start check up to this n')
    ap.add_argument('--nlap', type=int, default=256, help='full-Hessian Laplace estimate up to this n')
    ap.add_argument('--fit-nmin', type=int, default=32)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out', default=os.path.join(ROOT, 'results', 'data', 'quiver_index_saddle.json'))
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    exact = json.load(open(os.path.join(ROOT, 'results', 'data', 'quiver_singlet_index_dixon.json')))
    for q in range(1, 9):                         # n = 1 in closed form: Dixon (index) and Franel numbers (count)
        dix = sum((-1) ** m * math.comb(q, m) ** 3 for m in range(q + 1))
        fra = sum(math.comb(q, m) ** 3 for m in range(q + 1))
        e = exact.setdefault(f'1,{q}', dict(index=dix, source='closed form (n = 1)'))
        assert e['index'] == dix
        e.setdefault('singlet_total', fra)
    ns = [n for n in (1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048)
          if n <= a.nmax]
    t0 = time.time()
    rec = dict(params=vars(a), C_arc=float(C_ARC), C_half=float(C_HALF), crystal_delta=CRYSTAL_DELTA, continuum={},
               discrete={}, fits={})

    print('1. continuum solutions')
    cont = {}
    for kind, ps in (('index', sorted(set(a.p_cont_index) | {float(p) for p in a.p_index if p > 2})),
                     ('count', sorted(set(a.p_cont_count) | {float(p) for p in a.p_count if p > 1}))):
        out = continuum_scan(kind, ps)
        rec['continuum'][kind] = {str(float(p)): r for p, (r, c) in out.items()}
        cont[kind] = {p: c for p, (r, c) in out.items()}

    print('2-3. discrete maxima, extrapolation, Laplace estimates')
    for kind, ps in (('index', a.p_index), ('count', a.p_count)):
        for p in ps:
            c = cont[kind].get(float(p))
            Fc = c.F if c is not None else (0.0 if kind == 'index' and p == 2 else None)
            print(f"  {kind} p = {p}: continuum F* = {Fc}")
            rows = discrete_scan(kind, p, ns, c, a.nfull, a.nlap, exact, rng)
            rec['discrete'][f'{kind},{p}'] = rows
            if kind == 'index' and p == 2:          # crystal: L_max = 3n ln(3n); corrected Laplace vs the closed form
                dev = [r['laplace_corr'] - r['exact_closed_form'] for r in rows if 'exact_closed_form' in r]
                print(f"    p = 2 check: corrected Laplace - ln((3n)!/(n!)^3) = " +
                      ", ".join(f"{d:+.5f}" for d in dev[:6]) + f", ..., {dev[-1]:+.2e} (n = {rows[len(dev)-1]['n']}); "
                      f"expected -1/(36n)")
                continue
            ft = fit(rows, Fc, a.fit_nmin)
            if ft:
                ent = c.entropy_term() if c is not None else None
                ft['b_heuristic_3_int_g_ln_2pi_g'] = 3 * ent if ent is not None else None
                ft['laplace_subleading'] = fit_laplace(rows, Fc)
                rec['fits'][f'{kind},{p}'] = ft
                print(f"    fit n >= {a.fit_nmin}: F = {ft['F']:.8f} (continuum {Fc:.8f}, diff {ft['F_minus_continuum']:.1e}); "
                      f"with F fixed: b = {ft['with_F_fixed']['b']:.5f} (heuristic 3 int g ln 2pi g = {ft['b_heuristic_3_int_g_ln_2pi_g']:.5f}), "
                      f"c = {ft['with_F_fixed']['c']:.4f}; free n ln n coefficient {ft['free_nlogn']['nlogn']:.4f} (expected 3)", flush=True)
                ls = ft['laplace_subleading']
                if ls:
                    print(f"    corrected Laplace - F* n^2 = {ls['b']:.4f} n + {ls['c']:.4f} ln n + {ls['d']:.4f} + {ls['e']:.3f}/n "
                          f"(n = {ls['nmin']}..{ls['nmax']}, max resid {ls['max_abs_resid']:.1e})", flush=True)

    print('4. exact consequences of D18 against the exact data')
    rec['structure_checks'] = structure_checks(exact)
    for s in rec['structure_checks']:
        print(f"  ({s['n']},{s['p']}) {s['check']}: {'ok' if s['ok'] else 'FAILS'}" + (f" (ratio {s['ratio']:.2f})" if s.get('ratio') else ''))

    print('5. exact vs saddle')
    table = []
    for kind in ('index', 'count'):
        for key, rows in rec['discrete'].items():
            if not key.startswith(kind):
                continue
            p = int(key.split(',')[1])
            Fc = rec['continuum'][kind].get(str(float(p)), {}).get('F', 0.0 if (kind == 'index' and p == 2) else None)
            for r in rows:
                if r.get('exact') is None and not (kind == 'index' and r.get('exact_sign') == 0):
                    continue
                if 'exact_sign' in r and r['exact_sign'] == 0:
                    table.append(dict(kind=kind, n=r['n'], p=p, exact='0', laplace=r.get('laplace')))
                    continue
                t = dict(kind=kind, n=r['n'], p=p, exact=r['exact'], laplace=r.get('laplace'),
                         laplace_corr=r.get('laplace_corr'), leading=Fc * r['n'] ** 2, residual=r['exact'] - Fc * r['n'] ** 2)
                table.append(t)
                print(f"  {kind:5s} (n,p)=({r['n']},{p}): exact {r['exact']:8.4f}  Laplace {r['laplace']:8.4f} "
                      f"(diff {r['laplace'] - r['exact']:+.3f})  corrected {r['laplace_corr']:8.4f} (diff {r['laplace_corr'] - r['exact']:+.3f})"
                      f"  F* n^2 {Fc * r['n'] ** 2:8.4f} (exact - F* n^2 = {t['residual']:+.3f})")
    rec['exact_vs_saddle'] = table
    rec['time_s'] = round(time.time() - t0, 1)
    json.dump(rec, open(a.out, 'w'), indent=1, default=float)
    print('saved', os.path.relpath(a.out, ROOT), f'[{rec["time_s"]} s]')


if __name__ == '__main__':
    main()
