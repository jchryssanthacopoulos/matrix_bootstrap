"""
Large-n saddle point of the gauge-singlet index (t = -1) and singlet count (t = +1) of the U(n)^3 fermionic quiver
(research/notes/quiver_project.md section 8; docs/derivations.md D18).

Exact eigenvalue representation (D18).  Let x_{vi} = e^{i theta_{vi}} be the eigenvalues of U_v (v = 0, 1, 2), with
edges v -> v+1 (mod 3) carrying p flavours, and f(alpha) = ln|2 sin(alpha/2)|.  Since
1 - e^{i phi} = 2 sin(phi/2) e^{i(phi - pi)/2} and the phases telescope around the 3-cycle,
    I_0(n,p) = (-i)^{3 p n^2} / (n!)^3  int prod_{v,i} d theta_{vi}/(2 pi)  exp L(theta) * sgn(theta)^p,
    L = sum_{same-node pairs} 2 f(theta_a - theta_b) + p sum_{different-node pairs} f(theta_a - theta_b),
    sgn(theta) = prod_{v} prod_{i,j} sign sin((theta_{vi} - theta_{v+1,j})/2)      (only matters for odd p).
The singlet count has the same form with f -> h(alpha) = ln|2 cos(alpha/2)| on different-node pairs and no phase.

Contents
  * kernels f, h and derivatives;
  * discrete maximisation of L (a three-species log gas on the circle: like species repel with weight 2, unlike with
    weight p (index) or attract with weight p (count)): symmetry-reduced Newton (index: three arcs rotated by 2 pi/3;
    count: coincident species) and unconstrained L-BFGS for checking the reduction;
  * the full 3n x 3n Hessian and the Laplace (Gaussian fluctuation) estimate of ln|I_0| about the maximum;
  * the continuum limit: arc densities expanded as sqrt(1-t^2) U_k(t); at fixed arc half-width a the functional is a
    quadratic form maximised exactly, then a is optimised.  ln|I_0| = F*(p) n^2 + o(n^2).
"""
import math
import numpy as np
from scipy.optimize import minimize, minimize_scalar

TWO_PI_3 = 2 * np.pi / 3


# ----------------------------------------------------------------------------------------------------------- kernels
def f0(x):
    return np.log(np.abs(2 * np.sin(x / 2)))


def f1(x):
    return 0.5 / np.tan(x / 2)


def f2(x):
    return -0.25 / np.sin(x / 2) ** 2


def h0(x):
    return np.log(np.abs(2 * np.cos(x / 2)))


def h1(x):
    return -0.5 * np.tan(x / 2)


def h2(x):
    return -0.25 / np.cos(x / 2) ** 2


def cross_kernel(kind):
    """(x, x', x'') of the even cross kernel of the reduced problem: index -> (f(a - 2pi/3) + f(a + 2pi/3))/2, count -> h."""
    if kind == 'index':
        return (lambda a: 0.5 * (f0(a - TWO_PI_3) + f0(a + TWO_PI_3)),
                lambda a: 0.5 * (f1(a - TWO_PI_3) + f1(a + TWO_PI_3)),
                lambda a: 0.5 * (f2(a - TWO_PI_3) + f2(a + TWO_PI_3)))
    if kind == 'count':
        return h0, h1, h2
    raise ValueError(kind)


# -------------------------------------------------------------------------------------------- full discrete problem
def species(n):
    return np.repeat(np.arange(3), n)


def full_LgH(th, n, p, kind='index', hess=True):
    """L, gradient and Hessian of the full 3n-angle problem (species v = entries v*n ... v*n+n-1)."""
    sp = species(n)
    same = sp[:, None] == sp[None, :]
    D = th[:, None] - th[None, :]
    N = 3 * n
    eye = np.eye(N, dtype=bool)
    Ds = np.where(eye, 1.0, D)                                  # avoid the diagonal singularity of f
    with np.errstate(divide='ignore', invalid='ignore'):       # masked entries (coincident cross pairs) may be singular
        return _full_LgH(D, Ds, same, eye, p, kind, hess)


def _full_LgH(D, Ds, same, eye, p, kind, hess):
    if kind == 'index':
        K0 = np.where(same, 2 * f0(Ds), p * f0(Ds))
        K1 = np.where(same, 2 * f1(Ds), p * f1(Ds))
        K2 = np.where(same, 2 * f2(Ds), p * f2(Ds)) if hess else None
    else:
        K0 = np.where(same, 2 * f0(Ds), p * h0(D))
        K1 = np.where(same, 2 * f1(Ds), p * h1(D))
        K2 = np.where(same, 2 * f2(Ds), p * h2(D)) if hess else None
    K0[eye] = 0.0
    K1[eye] = 0.0
    L = 0.5 * K0.sum()
    g = K1.sum(axis=1)
    if not hess:
        return L, g
    K2[eye] = 0.0
    H = -K2
    H[eye] = K2.sum(axis=1)
    return L, g, H


def sign_factor(th, n):
    """sgn(theta) = prod_v prod_ij sign sin((theta_vi - theta_{v+1,j})/2) (the odd-p sign of the index integrand)."""
    s = 1
    for v in range(3):
        a = th[v * n:(v + 1) * n]
        b = th[((v + 1) % 3) * n:((v + 1) % 3 + 1) * n]
        s *= int(np.prod(np.sign(np.sin((a[:, None] - b[None, :]) / 2))))
    return s


def expand(th0, kind):
    """Full configuration from the reduced one: index -> arcs rotated by 2 pi v/3; count -> coincident species."""
    if kind == 'index':
        return np.concatenate([th0 + TWO_PI_3 * v for v in range(3)])
    return np.concatenate([th0, th0, th0])


# ----------------------------------------------------------------------------------------- reduced discrete problem
def reduced_LgH(th, p, kind='index', hess=True):
    """L(expand(th)) = 3 sum_{i != j} f(th_i - th_j) + 3 p sum_{i,j} x(th_i - th_j), with x the even cross kernel."""
    x0, x1, x2 = cross_kernel(kind)
    n = len(th)
    D = th[:, None] - th[None, :]
    eye = np.eye(n, dtype=bool)
    Ds = np.where(eye, 1.0, D)
    F0 = f0(Ds); F0[eye] = 0.0
    X0 = x0(D)
    L = 3 * F0.sum() + 3 * p * X0.sum()
    F1 = f1(Ds); F1[eye] = 0.0
    X1 = x1(D); X1[eye] = 0.0
    g = 6 * F1.sum(axis=1) + 6 * p * X1.sum(axis=1)
    if not hess:
        return L, g
    K2 = 6 * f2(Ds) + 6 * p * x2(D)
    K2[eye] = 0.0
    H = -K2
    H[eye] = K2.sum(axis=1)
    return L, g, H


def newton_reduced(th, p, kind='index', tol=1e-11, maxit=200):
    """Damped Newton ascent of the reduced L from an ordered starting configuration; returns (L, th, |grad|, iters)."""
    th = np.sort(np.array(th, float))
    n = len(th)
    L, g, H = reduced_LgH(th, p, kind)
    it = 0
    for it in range(maxit):
        if np.abs(g).max() < tol * max(1.0, n):
            break
        A = -(H - np.ones((n, n)) / n)                         # positive definite near a maximum (rotation mode lifted)
        try:
            c = np.linalg.cholesky(A)
            d = np.linalg.solve(c.T, np.linalg.solve(c, g))
        except np.linalg.LinAlgError:
            d = g / max(np.abs(np.diag(A)).max(), 1.0)          # fall back to a scaled gradient step
        s = 1.0
        while s > 1e-12:
            tn = th + s * d
            Ln, gn = reduced_LgH(tn, p, kind, hess=False)
            if np.isfinite(Ln) and Ln >= L - 1e-13 * abs(L) and np.all(np.diff(tn) > 0):   # no crossings
                break
            s *= 0.5
        th = tn
        L, g, H = reduced_LgH(th, p, kind)
    return L, th, float(np.abs(g).max()), it


def lbfgs_full(th, n, p, kind='index'):
    res = minimize(lambda t: tuple(-x for x in full_LgH(t, n, p, kind, hess=False)), th, jac=True, method='L-BFGS-B',
                   options=dict(maxiter=50000, gtol=1e-11, ftol=1e-16))
    L, g = full_LgH(res.x, n, p, kind, hess=False)
    return L, res.x, float(np.abs(g).max())


def initial_reduced(n, a, kind='index'):
    """Equally spaced (midpoint) starting points on [-a, a]."""
    return -a + 2 * a * (np.arange(n) + 0.5) / n


# ------------------------------------------------------------------------------------------------------- Laplace
CRYSTAL_DELTA = 1 - 0.5 * math.log(2 * math.pi)     # = ln(e / sqrt(2 pi)): Gaussian-about-crystal excess per beta = 2 point

def laplace_lnI(n, p, kind, L_max, H_full):
    """
    Gaussian estimate of ln|I_0| (index) or ln(count) about the maximum family:
        ln c + (1/2) ln(3n) + ((1 - 3n)/2) ln(2 pi) + L_max - (1/2) ln det'(-H),
    with the (n!)^3 cancelled by permutations inside the species and c the number of remaining maximum families
    (rotation orbits): index p > 2 -> 2 cyclic orientations of the arcs (they cancel for odd p and odd n: returns -inf);
    index p = 2 -> all labellings of the 3n-point crystal, (3n-1)!/(n!)^3; count -> 1.
    """
    ev = np.linalg.eigvalsh(-H_full)
    zero = ev[0]
    logdet = float(np.sum(np.log(ev[1:])))
    if kind == 'count':
        logc = 0.0
    elif p == 2:
        logc = math.lgamma(3 * n) - 3 * math.lgamma(n + 1)
    elif p % 2 == 1 and n % 2 == 1:
        return -np.inf, logdet, zero
    else:
        logc = math.log(2)
    val = logc + 0.5 * math.log(3 * n) + 0.5 * (1 - 3 * n) * math.log(2 * np.pi) + L_max - 0.5 * logdet
    return val, logdet, zero


def laplace_corrected(n, lap):
    """
    Laplace estimate minus the universal beta = 2 crystal excess 3n * ln(e/sqrt(2 pi)).  For one species on the full
    circle (CUE) the Gaussian approximation about the equally spaced crystal gives sqrt(N) N^N (2 pi)^{(1-N)/2} against
    the exact N!, an excess of N ln(e/sqrt(2 pi)) - 1/(12N) + O(N^-3) (Stirling); at p = 2 the quiver integrand is
    exactly CUE(3n).  Applying the same per-point excess to every species for p > 2 is a heuristic (local beta = 2 log
    gas in each arc; the inter-species interaction is smooth at the scale of the spacing).
    """
    return lap - 3 * n * CRYSTAL_DELTA


# ------------------------------------------------------------------------------------------------ continuum problem
def _r0(alpha):
    """r(alpha) = ln|2 sin(alpha/2)/alpha| (smooth for |alpha| < 2 pi, r(0) = 0)."""
    alpha = np.asarray(alpha, float)
    small = np.abs(alpha) < 1e-3
    with np.errstate(divide='ignore', invalid='ignore'):
        r = np.log(np.abs(2 * np.sin(alpha / 2) / alpha))
    return np.where(small, -alpha ** 2 / 24 - alpha ** 4 / 2880, r)


def _r1(alpha):
    """r'(alpha) = cot(alpha/2)/2 - 1/alpha."""
    alpha = np.asarray(alpha, float)
    small = np.abs(alpha) < 1e-3
    with np.errstate(divide='ignore', invalid='ignore'):
        r = 0.5 / np.tan(alpha / 2) - 1 / alpha
    return np.where(small, -alpha / 12 - alpha ** 3 / 720, r)


class Continuum:
    """
    Continuum (n -> infinity) maximum of F[rho] = lim L/n^2 within the symmetric ansatz
      index: rho_v(theta) = g(theta - 2 pi v/3), g even on [-a, a] with a < pi/3,
             F = 3 <g, f g> + 3 p <g, x g>,  x(alpha) = (f(alpha - 2pi/3) + f(alpha + 2pi/3))/2;
      count: rho_v = g for all v, a < pi/2,  F = 3 <g, f g> + 3 p <g, h g>.
    Method.  Write g(a t) = G(t)/a, G(t) = sum_{k even <= 2K} d_k T_k(t) / (pi sqrt(1-t^2)), d_0 = 1 (int G = 1).
    The Euler-Lagrange condition Phi'(theta) = 0 on (-a, a), with
        Phi(theta) = int g(theta') [f + p x](theta - theta') d theta'   (= delta F / (6 delta g)),
    splits into the Cauchy part (exact: PV int T_k / ((t - t') sqrt(1-t'^2)) = -pi U_{k-1}(t)) and a smooth part
    (Gauss-Chebyshev quadrature, J nodes); its projection on sqrt(1-t^2) U_{l-1}, l = 2..2K, is linear in d.  For each
    a this gives the hard-wall solution; its edge coefficient e(a) = sum_k d_k multiplies 1/sqrt(1-t^2) at t = +-1
    (e > 0: support too small, e < 0: negative density).  The physical maximum has a soft edge, e(a*) = 0.  Then
    F* = 3 lambda with lambda = Phi on the support (int g Phi = <g,fg> + p <g,xg> = F/3).
    """

    def __init__(self, p, kind='index', K=40, J=400):
        self.p, self.kind, self.K, self.J = p, kind, K, J
        self.ks = np.arange(0, 2 * K + 1, 2)
        jj = np.arange(1, J + 1)
        self.tau = np.cos((2 * jj - 1) * np.pi / (2 * J))                      # first-kind nodes (1/sqrt weight)
        self.T_tau = np.cos(np.outer(self.ks, np.arccos(self.tau)))            # T_k(tau_j)
        phi = jj * np.pi / (J + 1)
        self.t2 = np.cos(phi)                                                   # second-kind nodes (sqrt weight)
        self.w2 = np.pi / (J + 1) * np.sin(phi) ** 2
        self.Ul = np.sin(np.outer(self.ks[1:], phi)) / np.sin(phi)              # U_{l-1}(t2_mu), l = 2, 4, ..., 2K
        x0, x1, _ = cross_kernel(kind)
        self.x0, self.x1 = x0, x1
        self.amax = np.pi / 3 if kind == 'index' else np.pi / 2

    def S0(self, alpha):
        return _r0(alpha) + self.p * self.x0(alpha)

    def S1(self, alpha):
        return _r1(alpha) + self.p * self.x1(alpha)

    def solve_a(self, a):
        """Hard-wall Euler-Lagrange solution on [-a, a]: returns (d, edge coefficient e(a))."""
        S1m = self.S1(a * (self.t2[:, None] - self.tau[None, :]))              # (J, J)
        N = (self.Ul * self.w2) @ S1m @ self.T_tau.T / self.J                   # (K, K+1)
        A = np.pi / (2 * a) * np.eye(self.K) - N[:, 1:]
        dr = np.linalg.solve(A, N[:, 0])
        d = np.concatenate([[1.0], dr])
        return d, float(d.sum())

    def optimise(self, ngrid=200, xtol=1e-14):
        from scipy.optimize import brentq
        grid = np.linspace(0.02, 0.9995, ngrid) * self.amax
        es = np.array([self.solve_a(a)[1] for a in grid])
        idx = np.where(np.sign(es[:-1]) != np.sign(es[1:]))[0]
        if len(idx) == 0:
            raise RuntimeError(f'no soft-edge root for p={self.p} ({self.kind}); e(a) in [{es.min():.3g}, {es.max():.3g}]')
        roots = [brentq(lambda a: self.solve_a(a)[1], grid[i], grid[i + 1], xtol=xtol) for i in idx]
        self.roots = roots
        self.a = roots[0]
        self.d, self.edge = self.solve_a(self.a)
        self.lam = float(self.Phi(0.0)[0])
        self.F = 3 * self.lam
        return self.a, self.F

    def numerator(self, t):
        return np.cos(np.outer(np.arccos(np.clip(t, -1, 1)), self.ks)) @ self.d

    def G(self, t):
        t = np.asarray(t, float)
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.where(np.abs(t) < 1, self.numerator(t) / (np.pi * np.sqrt(1 - t ** 2)), 0.0)

    def density(self, theta):
        """Species-0 density g(theta) (zero outside [-a, a])."""
        return self.G(np.asarray(theta, float) / self.a) / self.a

    def cdf(self, theta):
        """int_{-a}^{theta} g: with t = cos(phi), int_t^1 T_k / (pi sqrt(1-t^2)) = phi/pi (k = 0), sin(k phi)/(k pi)."""
        scalar = np.ndim(theta) == 0
        t = np.clip(np.atleast_1d(np.asarray(theta, float)) / self.a, -1, 1)
        ph = np.arccos(t)
        ks = self.ks[1:]
        upper = ph / np.pi + (np.sin(np.outer(ph, ks)) / (ks * np.pi)) @ self.d[1:]   # mass in [t, 1]
        return float(1 - upper[0]) if scalar else 1 - upper

    def quantiles(self, n):
        """Positions theta_i with cdf = (i - 1/2)/n (starting points for the discrete problem)."""
        from scipy.optimize import brentq
        return np.array([brentq(lambda th: self.cdf(th) - (i + 0.5) / n, -self.a, self.a, xtol=1e-15) for i in range(n)])

    def Phi(self, theta):
        """Phi(theta) = int g(theta') [f + p x](theta - theta') d theta'; equals F/3 on [-a, a] at the optimum."""
        theta = np.atleast_1d(np.asarray(theta, float))
        out = np.empty_like(theta)
        ks = self.ks[1:]
        for i, th in enumerate(theta):
            x = th / self.a
            if abs(x) <= 1:
                logpart = -np.log(2) - np.sum(self.d[1:] / ks * np.cos(ks * np.arccos(x)))
            else:
                z = abs(x) - np.sqrt(x * x - 1)
                logpart = -np.log(2 * z) - np.sum(self.d[1:] / ks * z ** ks)
            smooth = np.sum((self.d @ self.T_tau) * self.S0(th - self.a * self.tau)) / self.J
            out[i] = np.log(self.a) + logpart + smooth
        return out

    def entropy_term(self):
        """int g ln(2 pi g) over one species."""
        ph = np.linspace(0, np.pi, 20001)[1:-1]
        t = np.cos(ph)
        G = self.G(t)
        integrand = np.where(G > 0, G * np.log(2 * np.pi * G / self.a), 0.0) * np.sin(ph)      # dt = sin(phi) dphi
        return float(np.trapezoid(integrand, ph))
