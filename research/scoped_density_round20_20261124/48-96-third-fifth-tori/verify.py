from fractions import Fraction as F
from math import comb
N=25

def mon(j=0): return [int(i==j) for i in range(N)]
def mul(a,b): return [sum(a[j]*b[i-j] for j in range(i+1)) for i in range(N)]
def pw(a,b):
    out=mon()
    for _ in range(b):out=mul(out,a)
    return out
def prod_eta(a,b,p):
    # Strip leading q powers from eta(tau)^a eta(p*tau)^b.
    out=mon()
    for d in range(1,N):
        k=a+b*int(d%p==0)
        w=[0]*N
        for j in range((N-1)//d+1):
            w[j*d]=((-1)**j*comb(k,j) if j<=k else 0) if k>=0 else comb(-k+j-1,j)
        out=mul(out,w)
    return out
def normalize(basis,const):
    out=[0]*N;cs=[]
    for i,b in enumerate(basis):
        assert b[i]==1 and not any(b[:i])
        c=(const if not i else 0)-out[i];cs.append(c)
        out=[x+c*y for x,y in zip(out,b)]
    assert out[:len(cs)]==[const]+[0]*(len(cs)-1)
    return out,cs
def v5(x):
    assert x>0
    e=0
    while x%5==0:e+=1;x//=5
    return e

A2=[0]*N
for x in range(-8,9):
    for y in range(-8,9):
        h=x*x+x*y+y*y
        if h<N:A2[h]+=1
B3=prod_eta(3,-1,3)
C3=mul(mon(1),prod_eta(-3,9,3))
UB=prod_eta(5,-1,5)
VB=mul(mon(1),prod_eta(-1,5,5))
E4=[1]+[240*sum(d**3 for d in range(1,n+1) if n%d==0) for n in range(1,N)]
Delta=mul(mon(1),prod_eta(24,0,5))
expected_third = {
 24:[(F(2),36,-270),(F(14,9),2,-27),(F(16,9),9,-108)],
 48:[(F(4),1836,-19512),(F(32,9),105,-5904),(F(34,9),450,-11250)],
 72:[(F(6),105084,-1278423),(F(50,9),6060,-480411),(F(52,9),25632,-802872)],
 96:[(F(8),6334524,-82872720),(F(68,9),366750,-34230924),(F(70,9),1541727,-54117315)],
}
expected_fifth={24:(60,7500),48:(3060,104062500),72:(175200,1623632812500),96:(10563300,26685278320312500)}
for s in (1,2,3,4):
    n=24*s;D=F(n,12)
    theta,_=normalize([mul(pw(E4,3*(s-i)),pw(Delta,i)) for i in range(s+1)],1)
    count_min=theta[s+1]
    thirds=[]
    for r in (0,1,2):
        k=12*s-r;J,a=divmod(k,3)
        out,cs=normalize([mul(mul(pw(A2,a),pw(B3,3*(J-i)+r)),pw(C3,i)) for i in range(s+1)],1)
        j=J-s
        assert cs[-1]>0
        e=F(3*j+r-a-9*s,2)
        assert e.denominator==1
        thirds.append((F(2,3)*(j+F(r,3)), cs[-1]*F(3)**int(e),out[s+1]))
    assert thirds==expected_third[n]
    print(n,"thirds (depth, nearest count at endpoint, least-shell twist)",thirds)
    # Whole-deep linear subspace Fourier divisibility.
    discrepancy=count_min-thirds[0][2]
    t=0
    while discrepancy%3==0:t+=1;discrepancy//=3
    assert t=={24:9,48:11,72:13,96:11}[n]
    print(n,"wholly deep linear span rank cap",t)
    # Scalar level 5 sum for isotropic and product for all.
    _,co1=normalize([mul(pw(UB,6*s-i),pw(VB,i)) for i in range(s+1)],2)
    _,co2=normalize([mul(pw(UB,12*s-i),pw(VB,i)) for i in range(s+1)],1)
    iso=co1[-1]
    assert co1[-1]>0 and co2[-1]>0
    general=co2[-1]*5**(3*s)
    assert (iso,general)==expected_fifth[n]
    quartic=F(3*iso*D**4,n*(n+2))
    assert quartic.denominator>1
    print(n,"fifth isotropic paired cap and count",D,iso,"quartic contradiction",quartic)
    if s>=2:
        # Need moments of order s+1, using h<=s+1 vanishing.
        gamma={2:F(n+3,n),3:F((n+3)*(n+5),n*(n+2)),4:F((n+5)*(n+7),n*(n+2))}[s]
        pval=2*(s+1)-v5(gamma.numerator)+v5(gamma.denominator)
        assert v5(general)<2*pval
        print(n,"fifth anisotropic refined sum cap",F(2,5)*(11*s-1),"product valuation",v5(general),"needs",2*pval)
    # Exact all-radii ternary comparison with endpoint mu=D+2.
    kmax={24:9,48:12,72:15,96:17}[n]
    left=F(9,8)*D/(D+2); right=F(18,D+2)
    assert left**n*right**kmax<=1 and left**n*right**(kmax+1)>1
    print(n,"all radius small-palette or lattice vertical dimension cap",kmax)
    if n in (48,72):
        # Unlimited palette, vertical dimension exactly one.  Gram argument.
        for deficit,lw,hi,residues in [
            (F(2,3),D-F(2,3)+F(1,10**4),D if n==48 else D+F(1,12),[0,1,2]),
            (F(2,9),D,D+F(1,12),[0]),
        ]:
            a=min(-thirds[r][2] for r in residues)
            M=1+count_min//a
            assert all(2*M*M*(u/(D+2))**(n+1)/(u-D+deficit)<1 for u in (lw,hi))
            print(n,"unlimited-palette vertical one interval",str(lw),str(hi),"group cap",M)
