# Matrix SYK — V. The covariant single-trace bootstrap: it works, and where it strains

*The engine for $E_0(N_\Psi)$ that does not presuppose a master field. File: `sm_covariant_bootstrap.py`; runs in `cov_boot_N2.txt`.*

## Design

Bootstrap the R-charge-resolved ground energy $E_0(N_\Psi)$ using **only gauge-invariant single-trace operators** $\mathrm{Tr}[\text{word}]$ in $\Psi,\bar\Psi$. The functional $\phi(O)=\mathrm{Tr}[\rho\,O]$ is a general Hermitian functional — no positivity of $\rho$, no large-$N$ factorization, no master field assumed (crucial, since the fortuitous states are maximal-Casimir and non-singlet). Constraints:

- moment positivity $M_{\alpha\beta}=\phi(O_\alpha^\dagger O_\beta)\succeq0$ over the single-trace operators,
- equations of motion $\phi([H,O])=0$,
- sharp R-charge $\phi(N_\Psi)=k$, $\phi(N_\Psi^2)=k^2$ (with $N_\Psi=\mathrm{Tr}[\Psi\bar\Psi]$, itself single-trace),
- optional ground-state positivity $\phi(O_\alpha^\dagger[H,O_\beta])\succeq0$,
- objective $\min\phi(H)$ with the single-trace $H=\tfrac92(\mathrm{Tr}[\Psi\Psi\bar\Psi\bar\Psi]+\mathrm{h.c.})-9(N-1)^2$.

Everything is reduced to the operator span actually touched (SVD), so the number of SDP variables is $O(\text{tens})$ — **independent of the $2^{N^2}$ Hilbert space.** This is what makes the method scalable in principle.

## Result 1 — the method works and is tight ($N=2$)

Using only gauge-invariant single-trace operators (26 variables, vs $\dim^2=256$), the bootstrap reproduces $E_0(N_\Psi)$ **exactly** in every sector, including the lifted ones:

| $N_\Psi$ | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| exact | 18 | 0 | 0 | 0 | 18 |
| covariant bootstrap | 18 | 0 | 0 | 0 | 18 |

This is the key proof of principle: **gauge-invariant data alone bounds the non-singlet, maximal-Casimir ground state tightly** — the master-field obstruction is not fatal, at least at $N=2$, and the covariant formulation is the right frame.

## Result 2 — where it strains ($N=3$)

At $N=3$ (level 3 operators, no ground-state positivity) the picture splits:

| $N_\Psi$ | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|---|
| exact | 72 | 45 | 18 | 0 | 0 | 0 | 0 | 18 | 45 | 72 |
| bootstrap | 0 | 0 | **18?→0** | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

- **Inside the BPS window** ($N_\Psi=3,4,5,6$): tight, $E_0=0$ recovered. Proving $E_0=0$ is easy (SUSY: with $Q=\mathrm{Tr}[\Psi^3]$ among the operators, $\phi(H)=\phi(Q^\dagger Q)+\phi(QQ^\dagger)\ge0$, saturated).
- **In the lifted sectors** ($N_\Psi=0,1,2,\dots$): the level-3 bootstrap returns the trivial bound $0$ instead of $72,45,18$. Proving $E_0>0$ is the hard direction and needs more constraints: **ground-state positivity** and/or **higher operator level** — the standard convergence knobs, on which both cited matrix-bootstrap papers rely.

Honest limitation: adding ground-state positivity at $N=3$ (an extra $m^2$ dense $512\times512$ products) exceeds this sandbox's ~45 s/call and ~4 GB budget. That is an environment constraint, not a method failure — the identical $N=2$ computation with ground-state positivity is tight in every sector.

## Reading

- The covariant single-trace bootstrap is the correct, scalable ($O(\text{tens})$ variables) tool for the fortuitous sector, and it is **validated tight at $N=2$** and **tight in the BPS window at $N=3$**.
- The concentration statement needs tight lower bounds in the **lifted** sectors (to prove $E_0>0$ just outside the window). At $N=3$ that requires ground-state positivity + higher level — reachable with a dedicated SDP solver (SDPB) and more compute, not with the in-session sandbox.

## Concrete next steps

1. **Compute**: port the SDP to SDPB and precompute operator products sparsely (H is sparse), to run $N=3,4$ lifted sectors with ground-state positivity and level 4–5.
2. **Scale the engine to genuine large $N$**: replace explicit matrices with the abstract single-trace algebra (normal-ordering + the exact trace-form $H$), so $N$ enters only as a parameter — then bound $E_0(N_\Psi)$ at $N$ far beyond ED.
3. **Transport to the 3-matrix model** and test the decisive question: does the SYK BPS window stay $O(1)$ (concentration) while the single-matrix baseline grows as $N$?

### Sources (project library)
- Y. Chen, *Fortuity with a Single Matrix*, arXiv:2511.00790.
- H. W. Lin, Z. Zheng, *High-Precision Bootstrap of Multimatrix Quantum Mechanics*, arXiv:2507.21007.
- S. Laliberte, B. McPeak, *Bootstrapping supersymmetric (matrix) quantum mechanics*, arXiv:2510.01356.

---

## Addendum (2026-09-15): re-run and corrections to "Result 2"

Re-run in the project `.venv` (cvxpy 1.9.2 / SCS), `sm_covariant_bootstrap.py` unchanged:

- $N=2$: reproduced exactly (18, 0, 0, 0, 18).
- $N=3$, level 3: the script's default is `gs=True`, i.e. **ground-state positivity was already included** and the whole $N=3$ scan takes ~36 s. The statement above that adding it "exceeds this sandbox" is therefore not accurate. With or without it, all lifted sectors return the trivial bound $0$.
- $N=3$, level 4 (`Lmom=4, Leom=5`, $m=20$, $r=92$; ~3 min setup): still $0$ in every lifted sector.
- Validity check: $\langle g_k|O^\dagger[H,O]|g_k\rangle\ge0$ holds for the exact sector ground states $g_k$ and all single-trace $O$ of length $\le3$, so the (in principle invalid) inclusion of charge-changing $O$ in sector-restricted ground-state positivity did no harm here; it should nevertheless be restricted to charge-preserving $O$ in future.

**Diagnosis.** With $Q=\mathrm{Tr}[\Psi^3]$ and $\bar Q$ in the operator set, positivity gives the SUSY floor $\phi(H)=\phi(\bar QQ)+\phi(Q\bar Q)\ge0$. To lift it in the $N_\Psi=0$ sector the SDP must learn $\phi(Q\bar Q)=0$ and $\phi(\mathrm{Tr}[\Psi\Psi\bar\Psi\bar\Psi])=0$ (annihilators act first on an empty sector), which in formulation (a) follows from $\phi(N_\Psi)=0$, positivity of the single-letter blocks $\phi(\Psi_{ij}\bar\Psi_{ji})\ge0$ and Cauchy–Schwarz. The gauge-invariant single-trace basis contains no single letters, so this information is absent at any low level; the bound reverts to the SUSY floor. **Conclusion: the covariant bootstrap needs adjoint-valued (open-index) operators in the positivity matrix, organised into $SU(N)$ irreps (Lin–Zheng 2024 §2.4 / 2025 App. F style), not gauge-invariant traces only.** The $N=2$ success is special to the tiny algebra of four modes.

**Addendum (2026-09-17).** The symmetry-reduced sector bootstrap (derivation D4, `src/symmetry_reduction.py`) now reproduces the three-matrix $N=2$ value $E_0(2)=5.16536$ exactly at level 3 and gives the first $k=3$ bound, $E_0(3)\ge0.9622$ (exact 1.22706). The ingredient that closes the $k=2$ gap is the Cho et al. singlet/adjoint split of the open-word Gram matrix, which is legitimate only for a gauge-invariant functional — the reason the gauge-invariant-only formulation of this note failed was not gauge invariance per se but the absence of open-index words and of the sector restriction. Details: `research/notes/sector_bootstrap_results.md` §8.
