"""
R-charge-resolved ground-state energy  E_0(N_Psi)  of the single-matrix model
Q=Tr[Psi^3], as a bootstrap validation of the target observable.

Since H = 3N(N^2-1) - 9 C2, the minimum energy at fixed fermion number is
    E_0(k) = 3N(N^2-1) - 9 * C2max(k),   C2max(k) = max Casimir in ^k adj.
The bootstrap MAXIMIZES phi(C2) at sharp charge phi(N_Psi)=k -> upper bound on C2max
-> rigorous LOWER bound on E_0(k). Compared against exact diagonalization.
"""
import numpy as np
from mm_model import build_model
from mm_bootstrap import bootstrap, enumerate_words, independent_subset

N = 2
M = build_model(N)
Q, Qbar = M['Q'], M['Qbar']
H, C2, Npsi = M['H'], M['C2'], M['Npsi']
const = 3 * N * (N * N - 1)
dim = M['dim']

# ---- exact E_0(k) per charge sector via ED ----
pc = np.array([bin(i).count("1") for i in range(dim)])
Hd = H.copy()
print(f"single matrix N={N}: exact BPS window = [{N*(N-1)//2}, {N*(N+1)//2}]  (E_0=0 there)")
print(f"{'N_Psi':>5} {'dim':>4} {'E0_exact':>9} {'C2max_ex':>9} {'E0_boot(L=1+Q)':>15} {'E0_boot(L=2+Q)':>15}")
for k in range(N * N + 1):
    ix = np.where(pc == k)[0]
    blk = Hd[np.ix_(ix, ix)]
    e = np.linalg.eigvalsh(blk); e0 = e.min()
    c2max_ex = (const - e0) / 9.0
    # bootstrap: max phi(C2) at sharp charge k -> lower bound on E0
    out = []
    for L in (1, 2):
        r = bootstrap(M, L=L, Lc=3, objective='maxC2',
                      extra_ops=[Q, Qbar], charge=k)
        c2b = r['value']
        e0b = const - 9 * c2b
        out.append(e0b)
    print(f"{k:>5} {len(ix):>4} {e0:>9.3f} {c2max_ex:>9.3f} {out[0]:>15.3f} {out[1]:>15.3f}")
print("\n(bootstrap gives a LOWER bound on E_0(k); it should rise toward E0_exact as L grows,")
print(" and hit 0 exactly inside the BPS window.)")
