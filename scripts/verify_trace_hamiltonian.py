"""verify_trace_hamiltonian.py -- numerical check of docs/derivations.md D2 (and D1 for p=1).

Usage:  cd src && PYTHONPATH=. ../.venv/bin/python ../scripts/verify_trace_hamiltonian.py

Verify the general single-trace form of H={Q,Qbar} for Q = C_abc Tr[Psi^a Psi^b Psi^c]
(C cyclic), for several (p,N) and random complex C, and for Chen's C at p=3,N=2.
Formula (derived 2026-09-16):
  H = 9 K_{bc;ed} Tr[Psi^b Psi^c Psibar^e Psibar^d] - 9N M_{cd} Tr[Psi^c Psibar^d]
      + 9 L_{cd} TrPsi^c TrPsibar^d + 3N^3 |C|^2 - 3N <C,C^rev>
  K_{bc;ed}=sum_a C_abc conj(C_dea), M_cd=sum_ab C_abc conj(C_abd), L_cd=sum_ab C_abc conj(C_bad)
  |C|^2=sum C_abc conj(C_abc), <C,C^rev>=sum C_abc conj(C_cba)."""
import sys, os, time, itertools, numpy as np, scipy.sparse as sp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
from msyk_model import jw_creation_sparse

def cyclic(C):  return (C+np.transpose(C,(1,2,0))+np.transpose(C,(2,0,1)))/3
def chenC(p=3):
    Ct=np.zeros((p,p,p))
    for a,b,c in itertools.product(range(p),repeat=3):
        if a<=b<=c: Ct[a,b,c]=1
    return cyclic(Ct)

def build(N,p,C):
    n=p*N*N; dim=2**n
    cd=jw_creation_sparse(n); an=[x.getH().tocsr() for x in cd]
    Cr=lambda f,i,j: cd[f*N*N+i*N+j]            # Psi^f_ij
    Bb=lambda f,i,j: an[f*N*N+j*N+i]            # Psibar^f_ij = (Psi^f_ji)^dagger
    Z=lambda: sp.csr_matrix((dim,dim),dtype=complex)
    Q=Z()
    for a,b,c in itertools.product(range(p),repeat=3):
        if abs(C[a,b,c])<1e-14: continue
        for i,j,k in itertools.product(range(N),repeat=3):
            Q=Q+C[a,b,c]*(Cr(a,i,j)@Cr(b,j,k)@Cr(c,k,i))
    Q=Q.tocsr(); Qb=Q.getH().tocsr(); H=(Q@Qb+Qb@Q).tocsr()
    return dict(n=n,dim=dim,Cr=Cr,Bb=Bb,Q=Q,H=H,Z=Z,cd=cd,an=an)

def formula(N,p,C,m):
    Cr,Bb,Z,dim=m['Cr'],m['Bb'],m['Z'],m['dim']
    K=np.einsum('abc,dea->bced',C,np.conj(C)); M=np.einsum('abc,dab->cd',C,np.conj(C))
    L=np.einsum('abc,dba->cd',C,np.conj(C)); n2=np.einsum('abc,abc',C,np.conj(C)); rev=np.einsum('abc,cba',C,np.conj(C))
    H4=Z()
    for b,c,e,d in itertools.product(range(p),repeat=4):
        if abs(K[b,c,e,d])<1e-14: continue
        V=Z()
        for i,j,k,l in itertools.product(range(N),repeat=4):   # Tr[Psi^b Psi^c Psibar^e Psibar^d] = Psi^b_ij Psi^c_jk Psibar^e_kl Psibar^d_li
            V=V+Cr(b,i,j)@Cr(c,j,k)@Bb(e,k,l)@Bb(d,l,i)
        H4=H4+9*K[b,c,e,d]*V
    F={}
    for c,d in itertools.product(range(p),repeat=2):
        F[c,d]=sum((Cr(c,i,j)@Bb(d,j,i) for i in range(N) for j in range(N)),Z())
    TrP=[sum((Cr(c,i,i) for i in range(N)),Z()) for c in range(p)]
    TrB=[sum((Bb(c,i,i) for i in range(N)),Z()) for c in range(p)]
    H2=Z()
    for c,d in itertools.product(range(p),repeat=2):
        H2=H2-9*N*M[c,d]*F[c,d]+9*L[c,d]*(TrP[c]@TrB[d])
    H0=3*N**3*n2-3*N*rev
    return (H4+H2+H0*sp.identity(dim,dtype=complex,format='csr')).tocsr(), H0

if __name__ == '__main__':
  rng=np.random.default_rng(1)
  cases=[(2,1,'chen'),(3,1,'chen'),(2,2,'rand'),(3,2,'rand'),(2,3,'rand'),(2,3,'chen'),(2,3,'rand2')]
  for N,p,kind in cases:
    t=time.time()
    if kind=='chen': C=chenC(p) if p==3 else cyclic(np.ones((p,p,p)))  # p<3: all-ones cyclic (p=1: Tr Psi^3)
    else:
        C=cyclic(rng.normal(size=(p,p,p))+1j*rng.normal(size=(p,p,p)))
    m=build(N,p,C); Hf,H0=formula(N,p,C,m)
    v=rng.normal(size=m['dim'])+1j*rng.normal(size=m['dim'])
    res=np.linalg.norm(Hf@v-m['H']@v)/np.linalg.norm(m['H']@v)
    e0=np.zeros(m['dim']); e0[-1]=1   # Fock vacuum = last index in this JW convention (creation clears a bit)
    vac=np.real(e0@(m['H']@e0))
    print(f"N={N} p={p} C={kind:5s} dim={m['dim']:7d}: rel.resid {res:.1e}   <0|H|0>={vac:.6f} vs H0={H0.real:.6f}   [{time.time()-t:.0f}s]")
