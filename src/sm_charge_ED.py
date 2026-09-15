"""
Single-matrix model Q=Tr[Psi^3] via sparse ED at N=2,3,4.
Compute E_0(N_Psi) per R-charge sector; locate the BPS window.
Expected (exact): E_0=0 iff N_Psi in [N(N-1)/2, N(N+1)/2]  (window width N).
This is the large-N baseline: the window GROWS with N -> NOT R-charge concentrated,
in contrast to the p=3 matrix SYK model whose window stays O(1).
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from msyk_model import build_msyk

def E0_by_charge(N):
    M = build_msyk(N=N, p=1)          # single fermion matrix
    H = M['H'].tocsr(); dim = M['dim']; n = M['n']
    const = 3 * N * (N * N - 1)
    pc = np.array([bin(i).count("1") for i in range(dim)])
    res = {}
    for k in range(n + 1):
        ix = np.where(pc == k)[0]
        d = len(ix)
        Hb = H[ix][:, ix]
        Hb = (Hb + Hb.getH()) * 0.5
        if d == 1:
            e0 = float(Hb.toarray()[0, 0].real)
        elif d <= 300:
            e0 = float(np.linalg.eigvalsh(Hb.toarray())[0].real)
        else:
            # smallest algebraic eigenvalue via shift-invert (H is PSD)
            try:
                e0 = float(sla.eigsh(Hb, k=1, sigma=-1e-6, which='LM',
                                     return_eigenvectors=False)[0].real)
            except Exception:
                e0 = float(sla.eigsh(Hb, k=1, which='SA',
                                     return_eigenvectors=False)[0].real)
        res[k] = (d, max(e0, 0.0), const)
    return res

for N in (2, 3, 4):
    r = E0_by_charge(N)
    lo, hi = N * (N - 1) // 2, N * (N + 1) // 2
    window = [k for k, (_, e0, _) in r.items() if e0 < 1e-5]
    print(f"\n===== single matrix N={N}  (n={N*N} fermions) =====")
    print(f"   predicted BPS window [N(N-1)/2, N(N+1)/2] = [{lo},{hi}]  (width {N}, {hi-lo+1} sectors)")
    print(f"   {'N_Psi':>5} {'dim':>6} {'E_0':>10} {'C2max=(const-E0)/9':>18}")
    for k, (d, e0, const) in r.items():
        star = "  <-- BPS (E_0=0)" if e0 < 1e-5 else ""
        print(f"   {k:>5} {d:>6} {e0:>10.4f} {(const-e0)/9:>18.4f}{star}")
    print(f"   -> measured BPS window = {window}  (matches [{lo},{hi}]: {window==list(range(lo,hi+1))})")
