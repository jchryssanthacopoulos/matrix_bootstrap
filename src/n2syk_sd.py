"""
Large-N Schwinger-Dyson equations of supersymmetric SYK with a q-hat-fermion supercharge, used for the large-p limit
of the fermionic U(n)^3 quiver (docs/derivations.md D22; research/notes/quiver_project.md section 14).

Equations (Fu-Gaiotto-Maldacena-Sachdev 2016, eqs. (2.11)-(2.12), N = 1 with Majorana fermions; at zero R-chemical
potential the N = 2 complex equations are identical, and the N = 2 free energy per complex fermion is twice the N = 1
free energy per Majorana):
    Sigma_psi(tau) = (qh - 1) J G_b(tau) G_psi(tau)^(qh-2),      Sigma_b(tau) = J G_psi(tau)^(qh-1),
    G_psi(iw) = 1 / (-iw - Sigma_psi(iw))   (fermionic w),         G_b(inu) = -1 / (1 + Sigma_b(inu))  (bosonic nu).
G_b = -delta(tau) + Gt_b(tau); the delta part drops out of Sigma_psi at zero chemical potential (G_psi(0) = 0 by
particle-hole symmetry).  Imaginary time is discretised on the midpoint grid tau_k = (k + 1/2) beta / M, the free
fermion tail 1/(-iw) is transformed analytically.

Energy (true <H>/N per Majorana, including the c-number part of H = Q^2):
    E/N = J / (qh 2^qh) - (J/2) int_0^beta Gt_b(tau) G_psi(tau)^(qh-1) dtau.
The first term is the normal-ordering constant sum_I C_I^2 / 2^qh of Q^2 (= <H>_infinity / N); the second is
J d/dJ of the bilocal action from the smooth part of G_b.  FGMS (2.37) gives J d_J log Z through d_tau G_psi(0+),
which drops the constant (their large-qh free energy (2.38) has a ground-state energy -J/(4 qh^2) to subtract); the
quantity returned by solve() as 'E' is that slope form and is NOT the true energy; use energy() or
scripts/n2syk_thermo.py.  Checks (D22): E -> J/24 and dE/dbeta -> -J^2/32 at beta -> 0 (qh = 3), E -> 0 as T -> 0,
and S_0 = (1/2) log 2 - int_0^infty E dbeta = (1/2) log(2 cos(pi/2qh)) to 3e-6.
"""
import numpy as np


def _grid(beta, M):
    dt = beta / M
    tau = (np.arange(M) + 0.5) * dt
    n = np.fft.fftfreq(M, d=1.0 / M).astype(int)          # 0, 1, ..., M/2-1, -M/2, ..., -1  (FFT order)
    wf = np.pi * (2 * n + 1) / beta                      # fermionic Matsubara frequencies in FFT order
    wb = 2 * np.pi * n / beta                            # bosonic
    return dt, tau, n, wf, wb


def f_to_w(F, beta, M, n):
    """Fermionic: F(iw_n) = dt sum_k F_k exp(i pi (2n+1)(k+1/2)/M)."""
    dt = beta / M
    k = np.arange(M)
    S = M * np.fft.ifft(F * np.exp(1j * np.pi * k / M))   # sum_k f_k e^{2 pi i n k / M}
    return dt * np.exp(1j * np.pi * (2 * n + 1) / (2 * M)) * S


def f_to_tau(Fw, beta, M, n):
    """Inverse fermionic: F_k = (1/beta) sum_n Fw_n exp(-i pi (2n+1)(k+1/2)/M)."""
    k = np.arange(M)
    S = np.fft.fft(Fw * np.exp(-1j * np.pi * n / M))     # sum_n g_n e^{-2 pi i n k / M}
    return (np.exp(-1j * np.pi * (k + 0.5) / M) * S / beta)


def b_to_w(F, beta, M, n):
    dt = beta / M
    S = M * np.fft.ifft(F)
    return dt * np.exp(1j * np.pi * n / M) * S


def b_to_tau(Fw, beta, M, n):
    S = np.fft.fft(Fw * np.exp(-1j * np.pi * n / M))
    return S / beta


def solve(betaJ, qh=3, M=2 ** 14, J=1.0, tol=1e-12, maxit=20000, init=None, mix=0.5, verbose=False):
    """Solve the SD equations at inverse temperature beta = betaJ / J.  Returns dict with tau, G_psi, Gt_b, energy."""
    beta = betaJ / J
    dt, tau, n, wf, wb = _grid(beta, M)
    if init is not None and len(init['G']) == M:
        G = init['G'].copy(); Gb = init['Gb'].copy()
    else:
        G = 0.5 * np.ones(M)                              # free fermion, 0 < tau < beta
        Gb = np.zeros(M)
    x = mix
    diff_old = np.inf
    for it in range(maxit):
        Sp = (qh - 1) * J * Gb * G ** (qh - 2)
        Sb = J * G ** (qh - 1)
        Spw = f_to_w(Sp, beta, M, n)
        Sbw = b_to_w(Sb, beta, M, n)
        Gw = 1.0 / (-1j * wf - Spw)
        Rw = Gw - 1.0 / (-1j * wf)                        # subtract the free tail 1/(-iw) <-> 1/2
        Gnew = 0.5 + f_to_tau(Rw, beta, M, n).real
        Gbw = Sbw / (1.0 + Sbw)                           # Gt_b(inu) = G_b + 1 = Sigma_b / (1 + Sigma_b)
        Gbnew = b_to_tau(Gbw, beta, M, n).real
        diff = np.abs(Gnew - G).max() + np.abs(Gbnew - Gb).max()
        if diff > diff_old:
            x = max(x / 2, 0.01)
        diff_old = diff
        G = (1 - x) * G + x * Gnew
        Gb = (1 - x) * Gb + x * Gbnew
        if diff < tol:
            break
    # energy from the slope at tau = 0+, extrapolated from the first grid points (quadratic fit)
    k = np.arange(4)
    c = np.polyfit(tau[k], G[k], 2)
    dG0 = c[1]                                            # derivative at tau = 0
    E = -dG0 / (2 * (qh - 1))                             # FGMS (2.37) slope form; NOT the true energy (see docstring)
    return dict(tau=tau, G=G, Gb=Gb, beta=beta, J=J, qh=qh, M=M, E=E, iterations=it + 1, converged=diff < tol,
                residual=diff)


def conformal_G(tau, beta, J=1.0, qh=3):
    """FGMS (2.29)-(2.30): G = b_psi [pi / (beta sin(pi tau / beta))]^(1/qh), b_psi = [tan(pi/2qh) / (2 pi J)]^(1/qh)."""
    b = (np.tan(np.pi / (2 * qh)) / (2 * np.pi * J)) ** (1.0 / qh)
    return b * (np.pi / (beta * np.sin(np.pi * tau / beta))) ** (1.0 / qh)


def energy(r):
    """True energy per Majorana (docstring formula) from a solve() result, midpoint rule."""
    dt = r['beta'] / r['M']
    return r['J'] / (2 ** r['qh'] * r['qh']) - 0.5 * r['J'] * np.sum(r['Gb'] * r['G'] ** (r['qh'] - 1)) * dt


def thermodynamics(betaJ, E, qh, alpha_s=None):
    """Spline thermodynamics from the energy E(beta J) per Majorana (J = 1) on a grid.

    Cubic spline of E (1 + beta) in u = log(1 + beta) (smooth at both ends), with E(0) = 1/(qh 2^qh) prepended;
    log Z = (1/2) log 2 - int_0^beta E, S = log Z + beta E.  Beyond the grid the Schwarzian tail
    E = 2 pi^2 alpha_s / beta^2 is used if alpha_s is given.  Returns a function bJ -> (logZ, E, S) (arrays)."""
    from scipy.interpolate import CubicSpline
    b = np.concatenate([[0.0], np.asarray(betaJ, float)])
    e = np.concatenate([[1.0 / (2 ** qh * qh)], np.asarray(E, float)])
    u = np.log1p(b)
    F = CubicSpline(u, e * (1 + b))
    bmax = b[-1]

    def f(bJ):
        bJ = np.atleast_1d(np.asarray(bJ, float))
        logZ, EE = np.empty_like(bJ), np.empty_like(bJ)
        for i, x in enumerate(bJ):
            if x <= bmax:
                ux = np.log1p(x)
                logZ[i] = 0.5 * np.log(2) - F.integrate(0.0, ux)
                EE[i] = F(ux) / (1 + x)
            else:
                if alpha_s is None:
                    raise ValueError('beta J beyond the grid; pass alpha_s for the Schwarzian tail')
                c = 2 * np.pi ** 2 * alpha_s
                logZ[i] = 0.5 * np.log(2) - F.integrate(0.0, np.log1p(bmax)) - c * (1 / bmax - 1 / x)
                EE[i] = c / x ** 2
        return logZ, EE, logZ + bJ * EE
    return f
