"""Numerical harness for src/trace_algebra.py: every identity is checked as an explicit operator equation.
Run:  cd src && PYTHONPATH=. ../.venv/bin/python ../tests/test_trace_algebra.py  (or pytest)."""
import sys, os, random, itertools
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from fermion_matrix_model import build_model, chen_C, letters
import trace_algebra as ta

TOL = 1e-9


def models():
    return [('N=2 p=3 fourier', build_model(2, 3, chen_C(3), fourier=True), ta.fourier_C(chen_C(3))),
            ('N=2 p=3 flavor', build_model(2, 3, chen_C(3)), chen_C(3)),
            ('N=2 p=1', build_model(2, 1, np.ones((1, 1, 1))), np.ones((1, 1, 1))),
            ('N=3 p=1', build_model(3, 1, np.ones((1, 1, 1))), np.ones((1, 1, 1)))]


def rvec(model, seed):
    rng = np.random.default_rng(seed)
    return rng.normal(size=model['dim']) + 1j * rng.normal(size=model['dim'])


def rel(a, b):
    return np.linalg.norm(a - b) / max(np.linalg.norm(a), np.linalg.norm(b), 1e-12)


def test_routed_canonicalisation(n_random=60, maxL=5):
    """Random routed terms (random letters, random routing) equal their canonical trace-form expansion."""
    for name, m, C in models():
        lets = letters(m); v = rvec(m, 1); random.seed(1); worst = 0
        for _ in range(n_random):
            L = random.randint(1, maxL)
            ops = [random.choice(lets) for _ in range(L)]
            nxt = list(range(L)); random.shuffle(nxt)
            lhs = ta.routed_apply(m, ops, nxt, v)
            e = ta.canonicalize_routed(ops, nxt)
            rhs = ta.expr_apply(m, e, v)
            err = np.linalg.norm(lhs - rhs) / max(np.linalg.norm(v), 1e-12)
            worst = max(worst, err)
            assert err < TOL, (name, ops, nxt, err, e)
        print(f"  routed canonicalisation {name}: {n_random} terms, worst abs err {worst:.1e}")


def test_trace_and_product(n_random=40):
    for name, m, C in models():
        lets = letters(m); v = rvec(m, 2); random.seed(2)
        for _ in range(n_random):
            w1 = tuple(random.choice(lets) for _ in range(random.randint(1, 4)))
            w2 = tuple(random.choice(lets) for _ in range(random.randint(1, 3)))
            assert rel(ta.expr_apply(m, ta.trace(w1), v), ta.mono_apply(m, (w1,), v)) < TOL or \
                np.linalg.norm(ta.mono_apply(m, (w1,), v)) < 1e-9, (name, w1)
            prod = ta.mul(ta.trace(w1), ta.trace(w2))
            direct = ta.mono_apply(m, (w1, w2), v)
            assert np.linalg.norm(ta.expr_apply(m, prod, v) - direct) < TOL * max(1, np.linalg.norm(direct)), (name, w1, w2)
        print(f"  trace / product {name}: ok")


def test_dagger(n_random=30):
    for name, m, C in models():
        lets = letters(m); u = rvec(m, 3); v = rvec(m, 4); random.seed(3)
        for _ in range(n_random):
            mono = tuple(tuple(random.choice(lets) for _ in range(random.randint(1, 3))) for _ in range(random.randint(1, 2)))
            e = ta.canonical(mono); ed = ta.dagger(e)
            lhs = np.vdot(u, ta.expr_apply(m, ed, v))          # <u|E^dag|v>
            rhs = np.conj(np.vdot(v, ta.expr_apply(m, e, u)))  # conj <v|E|u>
            assert abs(lhs - rhs) < TOL * max(1, abs(lhs)), (name, mono, lhs, rhs)
        print(f"  dagger {name}: ok")


def test_model_operators():
    for name, m, C in models():
        v = rvec(m, 5)
        Q = ta.supercharge(C); H = ta.hamiltonian(C); Np = ta.number_operator(m['p'])
        eQ = rel(ta.expr_apply(m, Q, v), m['Q'] @ v)
        eH = rel(ta.expr_apply(m, H, v), m['H'] @ v)
        eN = rel(ta.expr_apply(m, Np, v), m['Npsi'] @ v)
        # {Q, Qbar} = H from the algebra itself
        Qb = ta.dagger(Q); anti = ta.mul(Q, Qb) + ta.mul(Qb, Q)
        eA = rel(ta.expr_apply(m, anti, v), m['H'] @ v)
        print(f"  model operators {name}: Q {eQ:.1e}  H(D2.11) {eH:.1e}  N_Psi {eN:.1e}  {{Q,Qbar}} vs H {eA:.1e}")
        assert eQ < TOL and eH < TOL and eN < TOL and eA < TOL


def test_commutators(n_random=15):
    for name, m, C in models():
        lets = letters(m); v = rvec(m, 6); random.seed(6); H = ta.hamiltonian(C)
        for _ in range(n_random):
            w = tuple(random.choice(lets) for _ in range(random.randint(1, 3)))
            X = ta.trace(w)
            comm = ta.commutator(H, X)
            lhs = ta.expr_apply(m, comm, v)
            Xv = ta.mono_apply(m, (w,), v)
            rhs = m['H'] @ Xv - ta.mono_apply(m, (w,), m['H'] @ v)
            assert np.linalg.norm(lhs - rhs) < TOL * max(1, np.linalg.norm(rhs)), (name, w)
        print(f"  [H, Tr w] {name}: ok")


def test_finiteN_relations():
    for name, m, C in models():
        if m['N'] > 3:
            continue
        lets = letters(m); v = rvec(m, 7); random.seed(7); nontriv = 0
        for _ in range(20):
            seq = [random.choice(lets) for _ in range(m['N'] + 1)]
            e = ta.finiteN_relations(seq, m['N'])
            if e:
                nontriv += 1
                assert np.linalg.norm(ta.expr_apply(m, e, v)) < TOL * max(1, np.linalg.norm(v)), (name, seq, e)
        print(f"  finite-N relations {name}: {nontriv}/20 non-trivial, all vanish numerically")


def test_self_symmetry():
    """Words/monomials with an odd self-symmetry reduce to lower terms (Tr[PP]=0, Tr[P]Tr[P]=0, ...)."""
    for name, m, C in models():
        lets = letters(m); v = rvec(m, 9); random.seed(9)
        cases = []
        for l in lets:
            cases += [((l, l),), ((l, l, l, l),), ((l,), (l,))]
        for _ in range(15):
            w = tuple(random.choice(lets) for _ in range(random.choice([1, 3])))
            cases.append((w, w))
            w2 = tuple(random.choice(lets) for _ in range(2))
            cases.append((w2 + w2,)); cases.append((w2 + w2 + w2,))
        for mono in cases:
            e = ta.canonical(mono)
            assert ta._self_symmetry(mono) is None or mono not in e, (name, mono, e)
            direct = ta.mono_apply(m, mono, v)
            assert np.linalg.norm(ta.expr_apply(m, e, v) - direct) < TOL * max(1, np.linalg.norm(v)), (name, mono, e)
        print(f"  self-symmetry reductions {name}: {len(cases)} cases ok")


def test_charge_parity():
    m = build_model(2, 3, chen_C(3), fourier=True); lets = letters(m); random.seed(8)
    for _ in range(30):
        mono = tuple(tuple(random.choice(lets) for _ in range(random.randint(1, 4))) for _ in range(random.randint(1, 2)))
        e = ta.canonical(mono); q, z, par = ta.mono_q(mono), ta.mono_z(mono, 3), ta.mono_parity(mono)
        for m2 in e:
            assert ta.mono_q(m2) == q and ta.mono_z(m2, 3) == z and ta.mono_parity(m2) == par
    print("  charge / Z_3 / parity conservation: ok")


if __name__ == '__main__':
    for t in (test_routed_canonicalisation, test_trace_and_product, test_dagger, test_model_operators,
              test_commutators, test_finiteN_relations, test_self_symmetry, test_charge_parity):
        print(t.__name__); t()
    print("all tests passed")
