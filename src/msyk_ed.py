"""
Exact diagonalization of the p=3, N=2 matrix SYK model (dim 4096).
Quantifies (1) BPS spectrum & R-charge concentration, (2) chaos via level statistics.
"""
import numpy as np
import scipy.sparse as sp
from msyk_model import build_msyk

np.set_printoptions(suppress=True)
M = build_msyk(N=2, p=3)
dim, n = M['dim'], M['n']
H = M['H']
Hd = H.toarray()
Hd = (Hd + Hd.conj().T) / 2

# occupation (R-charge) of each computational basis state = popcount of index
popcount = np.array([bin(i).count("1") for i in range(dim)])

print(f"===== matrix SYK  p=3, N=2   (n={n} complex fermions, dim={dim}) =====")

# ---------- (0) global spectrum ----------
evals = np.linalg.eigvalsh(Hd)
E0 = evals.min()
vals, cts = np.unique(np.round(evals, 4), return_counts=True)
print(f"\n[spectrum]  E_min = {E0:.6f}   E_max = {evals.max():.4f}")
nbps = int(np.sum(np.abs(evals) < 1e-6))
print(f"            #(E=0 BPS) states = {nbps}   ({100*nbps/dim:.1f}% of Hilbert space)")
print(f"            # DISTINCT energies = {len(vals)}  (out of {dim} states)")
print(f"            low levels are irrational, e.g. {', '.join(f'{v:.4f}' for v in vals[1:6])} ...")
print(f"            -> unlike the single-matrix model (only integer levels 0,18,45,72): genuinely interacting")

# ---------- (1) BPS states resolved by R-charge  -> concentration ----------
print("\n[R-charge concentration]  BPS (E=0) states per fermion number N_Psi:")
bps_per_k = {}
for k in range(n + 1):
    idx = np.where(popcount == k)[0]
    if len(idx) == 0:
        continue
    block = Hd[np.ix_(idx, idx)]
    ev = np.linalg.eigvalsh(block)
    z = int(np.sum(np.abs(ev) < 1e-6))
    bps_per_k[k] = (len(idx), z)
tot = sum(z for _, z in bps_per_k.values())
mean = sum(k * z for k, (_, z) in bps_per_k.items()) / tot
var = sum((k - mean) ** 2 * z for k, (_, z) in bps_per_k.items()) / tot
print(f"   {'N_Psi':>5} {'dim sector':>10} {'#BPS':>6} {'frac':>7}")
for k, (d, z) in bps_per_k.items():
    bar = "#" * int(40 * z / max(1, max(zz for _, zz in bps_per_k.values())))
    print(f"   {k:>5} {d:>10} {z:>6} {z/tot:>7.3f}  {bar}")
print(f"   total BPS = {tot};  <N_Psi> = {mean:.3f} (center = n/2 = {n/2});  std = {np.sqrt(var):.3f}")
print(f"   BPS occupy N_Psi in [{min(k for k,(_,z) in bps_per_k.items() if z>0)}, "
      f"{max(k for k,(_,z) in bps_per_k.items() if z>0)}] out of full range [0,{n}]")

# ---------- (2) chaos: level statistics in a symmetry-resolved sector ----------
# Fix N_Psi = 6 (central). Resolve gauge SU(2) via  K = a*C2gauge + b*Jz  (commute with H).
print("\n[chaos]  level-spacing ratio <r> in the N_Psi=6 sector, resolved by gauge SU(2)")
k0 = n // 2
idx = np.where(popcount == k0)[0]
Hb = Hd[np.ix_(idx, idx)]
C2g = M['C2gauge'].toarray()[np.ix_(idx, idx)]
Jz = M['Ja'][2].toarray()[np.ix_(idx, idx)]     # Cartan generator of gauge su(2)
K = 0.7123 * C2g + 0.31719 * Jz                 # generic commuting resolver
K = (K + K.conj().T) / 2
# simultaneous eigenbasis: diagonalize K, then block H by K-eigenvalue
wK, VK = np.linalg.eigh(K)
# group by K eigenvalue
groups, cur, curv = [], [wK[0]], [0]
for a in range(1, len(wK)):
    if abs(wK[a] - cur[-1]) < 1e-6:
        cur.append(wK[a]); curv.append(a)
    else:
        groups.append(curv); cur, curv = [wK[a]], [a]
groups.append(curv)

def r_ratios(energies, tol=1e-7):
    e = np.sort(np.array(energies).real)
    e = e[np.concatenate(([True], np.diff(e) > tol))]   # drop residual exact degeneracies
    if len(e) < 3:
        return []
    s = np.diff(e)
    s = s[s > tol]
    if len(s) < 2:
        return []
    return [min(s[i], s[i+1]) / max(s[i], s[i+1]) for i in range(len(s) - 1)]

allr, nsec, nlev = [], 0, 0
for gv in groups:
    sub = VK[:, gv].conj().T @ Hb @ VK[:, gv]
    esub = np.linalg.eigvalsh((sub + sub.conj().T) / 2)
    rr = r_ratios(esub)
    if rr:
        allr += rr; nsec += 1; nlev += len(esub)
allr = np.array(allr)
ev6 = np.linalg.eigvalsh(Hb)
ev6d = ev6[np.concatenate(([True], np.diff(ev6) > 1e-7))]
print(f"   N_Psi=6 sector: dim={len(idx)}, but only {len(ev6d)} DISTINCT energies (huge symmetry degeneracy)")
print(f"   resolved sub-sectors used = {nsec}, pooled ratios = {len(allr)},  <r> = {allr.mean():.3f}")
print(f"   (Poisson 0.386, GOE 0.5996, GUE 0.5307)")
print("   CAVEAT: with only ~50 distinct levels split among gauge SO(3) x flavor S3 x SUSY multiplets,")
print("   N=2 has too few independent levels for a clean RMT diagnosis. Spectral chaos is a LARGE-N")
print("   question -> exactly the regime inaccessible to ED and targeted by the bootstrap.")
