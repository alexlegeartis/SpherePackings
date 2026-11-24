import sys
from pathlib import Path
here=Path(__file__).resolve().parent
prior=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else here.parent.parent/'scoped_density_round19_20261124'/'96-golden-three-shell'
sys.path.insert(0,str(prior))
from affine3 import itinerary,relations,classes,F,product,json,jointrows3

from even9 import staircase
from joint3 import parity
from even import mul,powd,J,K,D,product
from anchored_replay import COEF,q
from span import coeffs,s as sg
from math import prod
import sympy as S
from jacobi_replay import theta

S16={0:680642990590923360000,4:509650967828070432000,
     8:209200292309552544000,12:43731241296017184000,
     16:4017806237491956000,20:123081209203680000,
     24:726542605152000,28:231390432000}
S16|={-s:v for s,v in S16.items()}
N16=sum(S16.values())
assert N16==theta(0)[8]
def gram(p,c,s):
    m=(16+s)//2; n=(48+s)//2
    return ((10,0,p,c),(0,10,c,p+c),(p,c,16,m),(c,p+c,m,n))
good={10:{2,6,10},12:{0,4,8,12,16},14:{2,6,10,14,18,22},16:set(S16)}
def valid(p,c,s):
    G=gram(p,c,s)
    a=160-p*p-c*c;b=10*G[3][3]-c*c-(p+c)**2;d=10*G[2][3]-p*c-c*(p+c)
    if min(a,b)<=0 or a*b<=d*d:return False
    for t in coeffs(G,16):
        n=q(G,t)
        if n and n<=16 and (n not in good or abs(sg(G,t)) not in good[n]):
            return False
    return True
def parts16():
    # Leave same-shell affine grouping to the shared exact quotient.
    return [((p,c,s),(-p,-c,s)) if (p,c)!=(0,0) else ((p,c,s),)
            for p,c,s in product(range(-8,9),range(-8,9),sorted(S16))
            if (p,c)>=(0,0) and valid(p,c,s)]
def odd4():
    JJ={(0,1):1,(0,-1):-1}
    for n in range(1,9):
        JJ=mul(JJ,mul(mul({(0,0):1,(n,2):-1},{(0,0):1,(n,-2):-1}),
                        {(n*k,0):k+1 for k in range(9//n+1)}))
    def one(a):
        P=mul(JJ,mul(powd(J,a),powd(K,3-a)))
        return [{p:w for (i,p),w in P.items() if i==n and p>0} for n in range(4)]
    D5=powd(D,5)
    return [[{(p,c):sum(D5.get((n+5-i-j,0),0)*A[i].get(p,0)*B[j].get(c,0)
                     for i in range(n+1) for j in range(n+1-i))
              for p,c in product(range(1,11),repeat=2)} for n in range(4)]
             for a,b in ((2,3),(3,2)) for A,B in [(one(a),one(b))]]
def rows16(P):
    C,E=staircase()
    assert len(C)==9 and all(len(e)==9 for e in E)
    J=[]
    for od,C,Es in [(0,C[5:9],[e[5:9] for e in E]),(1,[{}]*4,odd4())]:
        D=sorted(set(product(range(1 if od else 0,11),repeat=2))
                 |{q for z in [C,*Es] for coefficient in z for q,v in coefficient.items() if v})
        Q=S.Matrix([[e[0].get(k,0) for e in Es] for k in D])
        sel=Q.T.rref()[1];QI=Q[list(sel),:].inv()
        PR=[parity(P[0],D[i],od) for i in sel]
        PV=[C[0].get(D[i],0) for i in sel]
        for k in D:
            rr=[F(x) for x in S.Matrix([[e[3].get(k,0) for e in Es]])*QI]
            wr=parity(P[-1],k,od)
            if not any(rr) and not any(wr):
                assert not C[3].get(k,0)
                continue
            R=[-sum(t*r[i] for t,r in zip(rr,PR))for i in range(len(P[0]))]
            J.append((R+[0]*(len(P[1])+len(P[2]))+wr,
                      C[3].get(k,0)-sum(t*v for t,v in zip(rr,PV)),('q8',od,k)))
    return J
def poprows(P):
    M=[]
    for z,v in S16.items():
        M += [(lambda p,c,s,z=z:int(s==z),F(v),('size16',z)),
              (lambda p,c,s,z=z:(5*p*p+(2*c-p)**2)*int(s==z),F(800-10*z,48)*v,('cov116',z)),
              (lambda p,c,s,z=z:p*(2*c-p)*int(s==z),F(5*z-80,48)*v,('cov216',z))]
    for i in range(11):
        for j in range(11-i):
            if not (i+j) or (i+j)%2:continue
            m=(i+j)//2
            t=0 if i%2 else F(N16*160**m*prod(range(1,i,2))*prod(range(1,j,2)),prod(96+2*h for h in range(m)))
            M.append((lambda p,c,s,i=i,j=j:p**i*c**j,t,('moment16',i,j)))
    return [([0]*sum(map(len,P[:-1]))+[sum(f(*a) for a in bb) for bb in P[-1]],t,l) for f,t,l in M]
def model():
    R,V,M,A,B,C=jointrows3()
    P=[A,B,C,parts16()]
    print('parts',*[len(p) for p in P],flush=True)
    E=relations(P,(10,12,14,16))
    Q=classes(sum(map(len,P)),E)
    print('relations',len(E),'orbits',len(set(Q)-{-1}),'zeros',Q.count(-1),flush=True)
    return [(r+[0]*len(P[-1]),v,m) for r,v,m in zip(R,V,M)]+rows16(P)+poprows(P),P,E,Q
if __name__=='__main__':
    from span import admissible
    def read(path):
        return {int(k):F(v) for k,v in json.loads(path.read_text()).items()}
    M,P,E,Q=model()
    Y=read(here/'joint4-counts.json')
    X=read(prior/'joint3-counts.json')
    upto=sum(map(len,P[:-1]))
    assert X=={j:y for j,y in Y.items() if j<upto}
    assert len(Y)==556 and all(x>0 for x in Y.values()) and len(E)==3257
    for r,b,m in M:assert sum(r[j]*v for j,v in Y.items())==b,m
    for i,j in E:assert Y.get(i,0)==Y.get(j,0)
    n0=0
    for d,B in zip((10,12,14,16),P):
        for j,bb in enumerate(B):
            if n0+j in Y:
                f=valid if d==16 else lambda p,c,s:admissible(d,p,c,s)
                assert all(f(*z) for z in bb)
        n0+=len(B)
    print('PASS: '+str(len(M))+' rows, '+str(len(E))+' affine identities, unchanged first three shells and '+str(len(Y)-len(X))+' new positive orbits.')
