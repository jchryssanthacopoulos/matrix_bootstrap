# To do

*Updated 2026-09-15.*

## Setup
- [ ] Create a Python environment (numpy, scipy, cvxpy or an SDP solver) so that `src/` can be run; re-verify the ED numbers recorded in `research/notes/` (single-matrix $E_0(N_\Psi)$ for $N\le4$; 3-matrix $N=2$ BPS counts 243:486:243).
- [ ] Decide on an SDP backend (SCS via cvxpy for prototyping; SDPB/SDPA-GMP for production — Lin–Zheng 2025 App. A warns that these problems are not strictly feasible).

## Theory (docs/derivations.md)
- [ ] Write out the refined definition of R-charge concentration for the matrix models (irreducible complexes labelled by $N_\Psi$ mod 3, $SU(N)$ irrep, flavor) and derive what the single-matrix window (Chen 2.22) implies per complex.
- [ ] Derive the refined index $\mathrm{Tr}[(-1)^F q^{N_\Psi}]$ per $SU(N)$ irrep for the 3-matrix model as a unitary matrix integral with character insertion (Chen §3, Anninos–Silva §4); evaluate at $N=2$ and compare with ED.
- [ ] Derive the exact single-trace form of $H=\{Q,\bar Q\}$ for the 3-matrix model (analogue of the single-matrix identity recorded in research note IV) and check which parts are Casimirs.
- [ ] Large-$N$ moments $\mathrm{tr}H,\mathrm{tr}H^2$ per charge sector via planar diagrams (Gaitan et al. §4 method).
- [ ] State precisely when ground-state positivity may be imposed sector-by-sector (only for sector-preserving operators).

## Numerics (after theory)
- [ ] Reproduce the single-matrix covariant bootstrap at $N=2,3$; push lifted sectors at $N=3$ to tightness (higher level + ground-state positivity).
- [ ] Port to the 3-matrix model at $N=2$ with $S_3$-flavor blocking; validate against ED.
- [ ] Attempt $N=3$ lifted sectors for the 3-matrix model.
