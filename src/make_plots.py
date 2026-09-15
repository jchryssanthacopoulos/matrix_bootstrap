"""Example figures for the matrix-SYK bootstrap program."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from math import comb
from msyk_model import build_msyk

# ---- data ----
# single matrix E_0(N_Psi) (exact, verified N=2,3,4)
sm = {
    2: [18, 0, 0, 0, 18],
    3: [72, 45, 18, 0, 0, 0, 0, 18, 45, 72],
    4: [180, 144, 108, 72, 36, 18, 0, 0, 0, 0, 0, 18, 36, 72, 108, 144, 180],
}
# 3-matrix SYK N=2 profile (dense per-charge)
M = build_msyk(N=2, p=3); Hd = M['H'].toarray(); Hd = (Hd + Hd.conj().T) / 2
pc = np.array([bin(i).count("1") for i in range(M['dim'])])
syk2 = []
for k in range(13):
    ix = np.where(pc == k)[0]
    syk2.append(max(float(np.linalg.eigvalsh(Hd[np.ix_(ix, ix)])[0].real), 0.0))
# BPS counts
def bps_counts(profile):
    return [1 if e < 1e-4 else 0 for e in profile]

fig, ax = plt.subplots(2, 2, figsize=(12, 9))

# Panel A: E_0(N_Psi) intensive, single matrix vs SYK
axA = ax[0, 0]
for N, col in zip((2, 3, 4), ('#1f77b4', '#2ca02c', '#9467bd')):
    xs = np.arange(len(sm[N])) / N**2
    axA.plot(xs, np.array(sm[N]) / (3 * N * (N**2 - 1)), 'o-', color=col,
             label=f'single matrix N={N}')
xs = np.arange(13) / 12
axA.plot(xs, np.array(syk2) / max(syk2), 's--', color='#d62728', label='3-matrix SYK N=2')
axA.set_xlabel(r'$N_\Psi / (\#\mathrm{fermions})$'); axA.set_ylabel(r'$E_0 / E_{\max}$')
axA.set_title(r'(A) $R$-charge-resolved ground energy: flat-zero = BPS window')
axA.legend(fontsize=8); axA.grid(alpha=0.3)

# Panel B: BPS window WIDTH vs N
axB = ax[0, 1]
Ns = np.arange(2, 8)
axB.plot(Ns, Ns + 1, 'o-', color='#1f77b4', label='single matrix (exact: $N+1$)')
axB.plot([2], [3], 's', ms=12, color='#d62728', label='3-matrix SYK N=2 (=3)')
axB.axhline(3, ls=':', color='#d62728', alpha=0.5)
axB.annotate('concentration = this stays flat\nwhile the blue line grows',
             xy=(4, 3), xytext=(4.2, 5.5), fontsize=9,
             arrowprops=dict(arrowstyle='->', color='gray'))
axB.set_xlabel('N'); axB.set_ylabel('BPS window width (# charge sectors)')
axB.set_title('(B) The concentration question, visually')
axB.legend(fontsize=8); axB.grid(alpha=0.3)

# Panel C: BPS distributions (normalized), single matrix binomial vs SYK
axC = ax[1, 0]
for N, col in zip((2, 3, 4), ('#1f77b4', '#2ca02c', '#9467bd')):
    w0 = N * (N - 1) // 2
    cnt = np.array([comb(N, k - w0) if 0 <= k - w0 <= N else 0 for k in range(N * N + 1)],
                   float)
    xs = (np.arange(N * N + 1) - N * N / 2)
    axC.plot(xs, cnt / cnt.max(), 'o-', color=col, label=f'single matrix N={N} (binomial)')
syk_bps = np.array([243, 486, 243], float)
axC.plot([-1, 0, 1], syk_bps / syk_bps.max(), 's--', color='#d62728',
         label='3-matrix SYK N=2 (243:486:243)')
axC.set_xlabel(r'$N_\Psi - (\#\mathrm{fermions})/2$'); axC.set_ylabel('#BPS (normalized)')
axC.set_title('(C) BPS distribution: single matrix spreads (binomial), SYK is pinned')
axC.legend(fontsize=8); axC.grid(alpha=0.3)

# Panel D: spectra -- solvable (integer, degenerate) vs interacting (irrational)
axD = ax[1, 1]
Ms = build_msyk(N=3, p=1); Hs = Ms['H'].toarray(); Hs = (Hs + Hs.conj().T) / 2
ev_sm = np.linalg.eigvalsh(Hs)
ev_syk = np.linalg.eigvalsh(Hd)
axD.eventplot([ev_sm], colors='#1f77b4', lineoffsets=1, linelengths=0.8)
axD.eventplot([ev_syk], colors='#d62728', lineoffsets=0, linelengths=0.8)
axD.set_yticks([0, 1]); axD.set_yticklabels(['3-matrix SYK\nN=2 (168 levels)',
                                             'single matrix\nN=3 (4 levels)'])
axD.set_xlabel('E'); axD.set_title('(D) Spectra: solvable (4 integer levels) vs interacting')
axD.grid(alpha=0.3, axis='x')

plt.tight_layout()
plt.savefig("matrix_syk_figures.pdf", bbox_inches='tight')
plt.savefig("matrix_syk_figures.png", dpi=140, bbox_inches='tight')
print("saved matrix_syk_figures.pdf / .png")
