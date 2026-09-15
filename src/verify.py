import numpy as np
import itertools

def su_generators(N):
    """Generalized Gell-Mann matrices T^a normalized Tr(T^a T^b)=1/2 delta^ab, a=1..N^2-1."""
    gens=[]
    # symmetric off-diagonal
    for i in range(N):
        for j in range(i+1,N):
            S=np.zeros((N,N),complex); S[i,j]=1; S[j,i]=1
            gens.append(S/2)  # Tr(S^2)=2 -> (S/2): Tr=1/2? Tr((S/2)^2)=Tr(S^2)/4=2/4=1/2 ✓
            A=np.zeros((N,N),complex); A[i,j]=-1j; A[j,i]=1j
            gens.append(A/2)
    # diagonal (Cartan)
    for k in range(1,N):
        d=np.zeros(N,complex)
        for m in range(k): d[m]=1
        d[k]=-k
        D=np.diag(d)*np.sqrt(1.0/(k*(k+1)))  # normalize
        # Tr(D^2)=  (k*1 + k^2) * 1/(k(k+1)) = k(k+1)/(k(k+1))=1 -> need 1/2
        D=D/np.sqrt(2)
        gens.append(D)
    return gens

def struct_consts(gens):
    N2=len(gens)
    f=np.zeros((N2,N2,N2))
    for a in range(N2):
        for b in range(N2):
            comm=gens[a]@gens[b]-gens[b]@gens[a]
            for c in range(N2):
                # [T^a,T^b]=i f^{abc} T^c ; f^{abc}=-2i Tr([T^a,T^b] T^c)
                f[a,b,c]=(-2j*np.trace(comm@gens[c])).real
    return f

def fermions(n):
    """n fermionic creation ops c^dag (psi) as 2^n x 2^n via Jordan-Wigner. {c_i,c_j^dag}=delta."""
    I=np.eye(2); Z=np.array([[1,0],[0,-1]]); sp=np.array([[0,1],[0,0]])  # sp|down> ; use c^dag = sp with JW
    # define c_a^dag
    cdag=[]
    for a in range(n):
        op=np.array([[1.0]])
        for k in range(n):
            if k<a: op=np.kron(op,Z)
            elif k==a: op=np.kron(op,sp)
            else: op=np.kron(op,I)
        cdag.append(op.astype(complex))
    return cdag

def build(N):
    gens=su_generators(N)
    n=N*N-1
    f=struct_consts(gens)
    cdag=fermions(n)              # psi^a = c_a^dag (creation)
    c=[x.conj().T for x in cdag]  # bar psi^a = annihilation
    # check anticommutators
    dim=2**n
    def acomm(A,B): return A@B+B@A
    ok=True
    for a in range(n):
        for b in range(n):
            if not np.allclose(acomm(cdag[a],c[b]), (1 if a==b else 0)*np.eye(dim)): ok=False
            if not np.allclose(acomm(cdag[a],cdag[b]),0): ok=False
    # Q = (i/sqrt2) f^{abc} psi^a psi^b psi^c
    Q=np.zeros((dim,dim),complex)
    for a in range(n):
        for b in range(n):
            for cc in range(n):
                if abs(f[a,b,cc])>1e-12:
                    Q+= f[a,b,cc]*(cdag[a]@cdag[b]@cdag[cc])
    Q=(1j/np.sqrt(2))*Q
    Qbar=Q.conj().T
    H=Q@Qbar+Qbar@Q
    # C2 = sum_a J^a J^a, J^a=-i f^{abc} psi^b barpsi^c
    J=[]
    for a in range(n):
        Ja=np.zeros((dim,dim),complex)
        for b in range(n):
            for cc in range(n):
                if abs(f[a,b,cc])>1e-12:
                    Ja+= f[a,b,cc]*(cdag[b]@c[cc])
        Ja=-1j*Ja
        J.append(Ja)
    C2=sum(Ja@Ja for Ja in J)
    # fermion number in adjoint sector
    Nf=sum(cdag[a]@c[a] for a in range(n))
    return dict(N=N,n=n,anticomm_ok=ok,Q=Q,Qbar=Qbar,H=H,C2=C2,Nf=Nf)

for N in (2,3):
    d=build(N); H=d['H']; C2=d['C2']; Q=d['Q']; n=d['n']
    Q2=Q@Q
    ident = H - (3*N*(N*N-1)*np.eye(2**n) - 9*C2)
    ev=np.round(np.linalg.eigvalsh(H).real,6)
    vals,counts=np.unique(np.round(ev,4),return_counts=True)
    print(f"=== N={N}, n=N^2-1={n} fermions, Hilbert dim={2**n} (su(N) sector) ===")
    print(" anticommutators OK:",d['anticomm_ok'])
    print(" Q^2=0:",np.allclose(Q2,0))
    print(" {Q,Qbar}=H by construction; H Hermitian:",np.allclose(H,H.conj().T))
    print(" H == 3N(N^2-1) - 9 C2  (operator identity):",np.allclose(ident,0))
    print(" spectrum (E: degeneracy):",dict(zip(vals,counts)))
    # C2 eigenvalues at E=0 (BPS)
    w,V=np.linalg.eigh(H)
    bps = np.abs(w)<1e-6
    c2vals=np.round(np.einsum('ij,jk,ki->i',V.conj().T,C2,V).real,4)
    print(" #BPS(E=0) states:",int(bps.sum()),"  expected dim(r*)=3^(N(N-1)/2)=",3**(N*(N-1)//2))
    print(" C2 on BPS states (should be max = N(N^2-1)/3 =",N*(N*N-1)/3,"):",np.unique(c2vals[bps]))

print("\n=== Cross-check: E = 3N(N^2-1) - 9 C2 on every eigenspace ===")
for N in (2,3):
    d=build(N); H=d['H']; C2=d['C2']
    w,V=np.linalg.eigh(H)
    c2=np.einsum('ij,jk,ki->i',V.conj().T,C2,V).real
    const=3*N*(N*N-1)
    pred=const-9*c2
    print(f" N={N}: max|E - (3N(N^2-1)-9C2)| over all states = {np.max(np.abs(w-pred)):.2e}")
    # list distinct (E, C2)
    pairs=sorted(set((round(float(e),3),round(float(cc),3)) for e,cc in zip(w,c2)))
    print(f"        distinct (E, C2):", pairs)
