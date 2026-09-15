# Chang, Lin (2024) — Holographic covering and the fortuity of black holes

**File:** `papers/chang_lin_2024.pdf`

## Full citation

Chi-Ming Chang, Ying-Hsuan Lin, *Holographic covering and the fortuity of black holes*, arXiv:2402.10129v2 [hep-th] (17 Apr 2024).

## Main question

How should BPS states of a holographic CFT at different ranks $N$ be related, and can one classify them by their large-$N$ behaviour? Proposes the *monotone* vs *fortuitous* taxonomy and conjectures monotone ↔ smooth horizonless geometries, fortuitous ↔ typical black-hole microstates.

## Physical system

Sequences of SCFTs labelled by $N$ with supercharge $Q_N$; case studies: $SU(N)$ $\mathcal N=4$ SYM (covering = formal multi-traces, $I_N$ = trace relations) and symmetric product orbifolds $\mathrm{Sym}^N(T^4)$, $\mathrm{Sym}^N(K3)$ (LLM and Lunin–Mathur geometries).

## Important definitions

- **Covering** (Def. 1): a vector space $\tilde{\mathcal H}$ and subspaces $I_N$ with $\mathcal H_N\simeq\tilde{\mathcal H}/I_N$, short exact sequence $0\to I_N\to\tilde{\mathcal H}\to\mathcal H_N\to0$ (2.2).
- **Strictness** (Def. 2): $I_N\supsetneq I_{N+1}$ and every $v\in\tilde{\mathcal H}$ leaves $I_{N'}$ for some $N'$ (2.3).
- **Long exact sequence** (2.5) $\cdots\to H^n(I_N)\to H^n(\tilde{\mathcal H})\xrightarrow{\pi^*}H^n(\mathcal H_N)\xrightarrow{f}H^{n+1}(I_N)\to\cdots$ when $\pi_N\tilde Q=Q_N\pi_N$.
- **Monotone / fortuitous** (Def. 4): $\mathcal H^{\rm mon}_N\simeq\mathrm{im}\,\pi^*\simeq H^*(\tilde{\mathcal H})/\mathrm{im}\,i^*$; $\mathcal H^{\rm fts}_N\simeq\mathrm{im}\,f\simeq H^*(\mathcal H_N)/\mathrm{im}\,\pi^*$ (2.6). Equivalently: monotone lift $\tilde O$ has $\tilde Q\tilde O=0$ (2.7); fortuitous lift has $\tilde Q\tilde O\ne0$ but $\tilde Q\tilde O\in I_N$ (2.8) — it is BPS only because of a trace relation.
- Sequences (2.9)–(2.11): monotone states lift to all $N'>N$; fortuitous ones exist only for $N\le N_{\max}$ and acquire anomalous dimension outside (Fig. 1).

## Main assumptions

Strict covering exists; $Q$-cohomology is one-loop exact / constant on the conformal manifold; Hodge isomorphism cohomology ↔ BPS states (2.4).

## Main analytical results and conjectures

- **Conjecture 1**: $H^*(\tilde{\mathcal H})\simeq$ Fock space of perturbative BPS particles in AdS (3.1).
- **Conjecture 2**: $\mathcal H^{\rm mon}_N$ = quantisation of the moduli space of smooth horizonless SUSY geometries at fixed flux. Verified for LLM (half-BPS, deformation/canonical quantisation reproduces the finite-$N$ Hilbert space, App. A) and Lunin–Mathur/superstrata.
- **Theorem 1** (3.2): $\dim\mathcal H^{\rm mon}_{N,q}\le\dim H^*(\tilde{\mathcal H})_q$, $\dim\mathcal H^{\rm fts}_{N,q}\le\dim H^*(I_N)_q$.
- **Conjecture 3** (from (3.3)–(3.5)): at black-hole charges $q\sim N^{d/2}$, fortuitous states outnumber monotone ones exponentially; in $\mathcal N=4$ SYM $\#\text{mon}\lesssim e^{N^{3/2}}$ vs $\#\text{BPS}\sim e^{N^2}$ (4.1)–(4.3).
- Discussion: the bulk dual of a single fortuitous state is unknown (branes?); the $1/N$ vs $e^{-N^\#}$ organisation; open problem of lifting the covering to the operator algebra.

## Numerical methods

None (structural/representation-theoretic; earlier explicit SU(2) cohomology computations cited).

## Relevant figures/results

Fig. 1 (anomalous dimension of monotone vs fortuitous lifts as $N$ varies).

## Limitations

Definitions are at the vector-space level (no norm); the covering for disordered models (SYK) does not exist — this is exactly why Chang–Chen–Sia–Yang re-define monotony via successive uplifts rather than an infinite-$N$ cohomology.

## Relationship to our project

The conceptual origin of "fortuity". For the matrix models the covering is concrete: $\tilde{\mathcal H}$ = formal polynomials in the matrix entries / trace expressions, $I_N$ = relations that hold at rank $N$ (Chen's Young diagrams with more than $N-1$ rows, trace relations). Chen's fortuitous states are of type (2.8): $Q|\lambda\rangle=0$ only because the staircase diagram cannot be extended. The paper also supplies the *dictionary* one wants the 3-matrix model to reproduce: an exponential fortuitous sector (Conjecture 3) versus a polynomially small monotone sector; in the single-matrix model the BPS multiplicity of $r_*$ is $2^N$, not $e^{cN^2}$ (Chen (2.19)–(2.21)), which is the sharp sense in which it is "skeleton" fortuity. Our question Q1c (multiplicity growth of BPS states with $N$) is Conjecture 3 for matrix SYK.
