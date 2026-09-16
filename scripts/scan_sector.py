#!/usr/bin/env python3
"""
Allowed-region ("island / peninsula / archipelago") scans for one R-charge sector of a fermionic matrix model.

For a grid of energies E, imposes phi(H) = E on the sector SDP (docs/derivations.md D3) and
  * tests feasibility  (optionally with the eigenstate constraints phi(O H) = E phi(O)  -> archipelago),
  * minimises / maximises chosen observables at fixed E (-> allowed region in the (E, observable) plane).
Also records the exact (E, <observable>) pairs of the sector's eigenstates for comparison.

Example:
  cd src && PYTHONPATH=. ../.venv/bin/python ../scripts/scan_sector.py --N 3 --p 1 --k 2 --L_adj 2 --L_sing 2 \
      --Emin -5 --Emax 80 --nE 86 --eigenstate --out ../results/data/scan_p1_N3_k2_L2.json
"""
import argparse, json, os, sys, time, warnings
import numpy as np, scipy.sparse as sp
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
from fermion_matrix_model import build_model, chen_C
from sector_bootstrap import SectorSDP, _restrict


def observables_for(model, k):
    """Block-restricted observables: gauge Casimir, QbarQ, QQbar, flavor-symmetric number N_s."""
    ix = model['sectors'][k]; n = model['n']; d = len(ix); p = model['p']; N = model['N']
    obs = {'C2gauge': _restrict(model['C2'], ix, ix)}
    if k + 3 <= n:
        SQ = _restrict(model['Q'], model['sectors'][k + 3], ix); obs['QbarQ'] = (SQ.getH() @ SQ).tocsr()
    if k - 3 >= 0:
        SB = _restrict(model['Qbar'], model['sectors'][k - 3], ix); obs['QQbar'] = (SB.getH() @ SB).tocsr()
    if p > 1:
        Fs = sp.csr_matrix((model['dim'], model['dim']), dtype=complex)
        for c in range(p):
            for dd in range(p):
                for i in range(N):
                    for j in range(N):
                        Fs = Fs + model['Psi'](c, i, j) @ model['Psibar'](dd, j, i)
        obs['Ns'] = _restrict(Fs / p, ix, ix)
    return obs


def exact_pairs(S, obs):
    """Joint eigen-decomposition of H_k with each observable that commutes with it (C2gauge, Ns need not)."""
    Hd = S.Hk.toarray(); w, v = np.linalg.eigh(Hd)
    out = {'E': w.tolist()}
    for name, X in obs.items():
        Xd = X.toarray()
        vals = []
        # within each degenerate eigenspace diagonalise X (exact if [H,X]=0, otherwise gives the range of <X>)
        i = 0
        while i < len(w):
            j = i
            while j + 1 < len(w) and abs(w[j + 1] - w[i]) < 1e-7: j += 1
            B = v[:, i:j + 1]; Xs = B.conj().T @ Xd @ B
            vals += np.linalg.eigvalsh((Xs + Xs.conj().T) / 2).real.tolist()
            i = j + 1
        out[name] = vals
        out[name + '_commutes'] = bool(abs(Hd @ Xd - Xd @ Hd).max() < 1e-9)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=int, required=True); ap.add_argument('--p', type=int, required=True)
    ap.add_argument('--k', type=int, required=True)
    ap.add_argument('--L_adj', type=int, default=2); ap.add_argument('--L_sing', type=int, default=2); ap.add_argument('--L_eom', type=int, default=3)
    ap.add_argument('--gs', action='store_true'); ap.add_argument('--eigenstate', action='store_true')
    ap.add_argument('--Emin', type=float, required=True); ap.add_argument('--Emax', type=float, required=True); ap.add_argument('--nE', type=int, default=41)
    ap.add_argument('--observables', default='C2gauge,QbarQ,Ns'); ap.add_argument('--solver', default='SCS')
    ap.add_argument('--budget_gb', type=float, default=6.0); ap.add_argument('--out', required=True)
    a = ap.parse_args()
    C = chen_C(a.p) if a.p == 3 else np.ones((a.p, a.p, a.p))
    model = build_model(a.N, a.p, C)
    obs_all = observables_for(model, a.k)
    obs = {n_: obs_all[n_] for n_ in a.observables.split(',') if n_ in obs_all}
    t0 = time.time()
    S = SectorSDP(model, a.k, a.L_adj, a.L_sing, a.L_eom, add_Q=True, gs=a.gs, observables=obs,
                  eigenstate=a.eigenstate, verbose=True, budget_gb=a.budget_gb)
    res = dict(params=vars(a), r=S.r, d=S.d, n_eom_independent=int(S.eom_rows.shape[0]), build_time_s=S.build_time,
               exact=exact_pairs(S, obs), E_bound=S.solve('H', 'min', solver=a.solver)['value'],
               E_max_bound=S.solve('H', 'max', solver=a.solver)['value'], scan=[])
    print(f"sector k={a.k}: E_0 bound {res['E_bound']:.4f} (exact {S.exact:.4f}); max-E bound {res['E_max_bound']:.4f}", flush=True)
    for E in np.linspace(a.Emin, a.Emax, a.nE):
        row = dict(E=float(E))
        f = S.solve(None, 'min', fixE=E, solver=a.solver)
        row['feasible'] = f['status'] in ('optimal', 'optimal_inaccurate')
        row['feas_status'] = f['status']
        if row['feasible']:
            for name in obs:
                lo = S.solve(name, 'min', fixE=E, solver=a.solver); hi = S.solve(name, 'max', fixE=E, solver=a.solver)
                row[name] = [lo['value'], hi['value']]; row[name + '_status'] = [lo['status'], hi['status']]
        res['scan'].append(row)
        print(f"  E={E:8.3f}: feasible={row['feasible']} " + " ".join(f"{n_}=[{row[n_][0]:.3f},{row[n_][1]:.3f}]" for n_ in obs if n_ in row and row[n_][0] == row[n_][0]), flush=True)
    res['total_time_s'] = time.time() - t0
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    json.dump(res, open(a.out, 'w'), indent=1)
    print("wrote", a.out)


if __name__ == '__main__':
    main()
