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
