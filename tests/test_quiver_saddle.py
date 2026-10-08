#!/usr/bin/env python3
"""
Checks for src/quiver_saddle.py (docs/derivations.md D18).  Run as a script: python3 tests/test_quiver_saddle.py

  1. reduced energy/gradient = full energy/gradient on symmetric configurations (index and count);
  2. reduced and full Hessians against finite differences;
  3. p = 2: the corrected Laplace estimate reproduces ln((3n)!/(n!)^3) to O(1/n);
  4. continuum F*(3), G*(3) against the recorded values; Phi flat on the support;
  5. the eigenvalue representation against direct Monte Carlo-free quadrature at n = 1 (I_0(1,p) = Dixon).
"""
import math, os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
from quiver_saddle import (reduced_LgH, full_LgH, expand, newton_reduced, laplace_lnI, laplace_corrected,  # noqa: E402
                           initial_reduced, Continuum, f0)


def check(name, ok, detail=''):
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}")
    if not ok:
        sys.exit(1)


def test_reduced_vs_full(rng):
    for kind in ('index', 'count'):
        for p in (2, 3, 5):
            n = 7
            th = np.sort(rng.uniform(-0.8, 0.8, n))
            Lr, gr = reduced_LgH(th, p, kind, hess=False)
            Lf, gf = full_LgH(expand(th, kind), n, p, kind, hess=False)
            gsum = gf[:n] + gf[n:2 * n] + gf[2 * n:]
            check(f'reduced = full ({kind}, p={p})', abs(Lr - Lf) < 1e-9 * abs(Lf) and np.allclose(gr, gsum, atol=1e-9),
                  f'dL {abs(Lr - Lf):.1e}, dg {np.abs(gr - gsum).max():.1e}')


def test_hessians(rng):
    eps = 1e-6
    for kind in ('index', 'count'):
        p, n = 3, 5
        th = initial_reduced(n, 0.7) + 0.02 * rng.standard_normal(n)        # well separated (spacing ~0.28)
        _, _, H = reduced_LgH(th, p, kind)
        Hfd = np.column_stack([(reduced_LgH(th + eps * e, p, kind, hess=False)[1]
                                - reduced_LgH(th - eps * e, p, kind, hess=False)[1]) / (2 * eps) for e in np.eye(n)])
        rel = np.abs(H - Hfd).max() / np.abs(H).max()
        check(f'reduced Hessian vs finite differences ({kind})', rel < 1e-7, f'max relative dev {rel:.1e}')
        thf = expand(th, kind) + 1e-3 * rng.standard_normal(3 * n)
        _, _, Hf = full_LgH(thf, n, p, kind)
        Hffd = np.column_stack([(full_LgH(thf + eps * e, n, p, kind, hess=False)[1]
                                 - full_LgH(thf - eps * e, n, p, kind, hess=False)[1]) / (2 * eps) for e in np.eye(3 * n)])
        relf = np.abs(Hf - Hffd).max() / np.abs(Hf).max()
        check(f'full Hessian vs finite differences ({kind})', relf < 1e-7, f'max relative dev {relf:.1e}')


def test_p2_closed_form():
    for n in (1, 2, 3, 5, 8):
        L, th, g, _ = newton_reduced(initial_reduced(n, np.pi / 3), 2, 'index')
        Lf, gf, H = full_LgH(expand(th, 'index'), n, 2, 'index')
        lap, _, _ = laplace_lnI(n, 2, 'index', Lf, H)
        exact = math.lgamma(3 * n + 1) - 3 * math.lgamma(n + 1)
        dev = laplace_corrected(n, lap) - exact
        check(f'p=2 corrected Laplace vs (3n)!/(n!)^3, n={n}', abs(dev + 1 / (36 * n)) < 2e-3 / n,
              f'dev {dev:+.5f}, expected ~ {-1 / (36 * n):+.5f}')


def test_continuum():
    for kind, p, ref in (('index', 3, 1.237284226192), ('count', 3, 3.042661737425), ('index', 4, 2.554963345090)):
        c = Continuum(p, kind, 40, 400)
        a, F = c.optimise()
        Ph = c.Phi(np.linspace(-0.99 * a, 0.99 * a, 21))
        check(f'continuum {kind} p={p}', abs(F - ref) < 1e-9 and Ph.max() - Ph.min() < 1e-10,
              f'F* = {F:.12f} (ref {ref}), Phi spread {Ph.max() - Ph.min():.1e}')


def test_eigenvalue_form_n1():
    """n = 1: I_0 = (-i)^{3p} int d^3theta/(2pi)^3 e^L sgn^p, by a product-trapezoid rule (exact for trig polynomials)."""
    M = 64
    t = 2 * np.pi * np.arange(M) / M
    T0, T1, T2 = np.meshgrid(t, t + 1e-9, t + 2e-9, indexing='ij')         # offsets avoid the integrable zeros
    for p in range(1, 7):
        val = 1.0 + 0j
        for a, b in ((T0, T1), (T1, T2), (T2, T0)):
            val = val * (1 - np.exp(1j * (a - b))) ** p                         # the original integrand
        direct = val.mean()
        L = p * (f0(T0 - T1) + f0(T1 - T2) + f0(T2 - T0))
        sg = np.sign(np.sin((T0 - T1) / 2) * np.sin((T1 - T2) / 2) * np.sin((T2 - T0) / 2)) ** p
        rep = (-1j) ** (3 * p) * np.mean(np.exp(L) * sg)
        dixon = sum((-1) ** m * math.comb(p, m) ** 3 for m in range(p + 1))
        check(f'eigenvalue form at n=1, p={p}', abs(direct - dixon) < 1e-6 and abs(rep - dixon) < 1e-6,
              f'direct {direct.real:+.6f}, representation {rep.real:+.6f}, Dixon {dixon}')


if __name__ == '__main__':
    rng = np.random.default_rng(1)
    test_reduced_vs_full(rng)
    test_hessians(rng)
    test_p2_closed_form()
    test_continuum()
    test_eigenvalue_form_n1()
    print('all quiver_saddle checks passed')
