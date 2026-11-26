"""Replay using Python 3 standard library. Run from any directory."""
from fractions import Fraction as F
from functools import reduce
from pathlib import Path
from math import gcd, comb
import json, hashlib

here=Path(__file__).resolve().parent
Gfile=here.parent/'adaptive_layers_20261115'/'sources'/'leterrier54_gram.txt'
data=json.loads((here/'data.json').read_text())
def need(b):
    if not b: raise ArithmeticError("certificate check failed")
need(hashlib.sha256(Gfile.read_bytes()).hexdigest()=="74a193d77f6faa0e94c79fa2795e57bcd7778d31b5a3840e276fc27acacde0f3")
G=[list(map(int,x.split())) for x in Gfile.read_text().splitlines() if x.strip()]
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def mv(m,x): return [dot(r,x) for r in m]
def bil(a,m,b): return dot(a,mv(m,b))
def inverse_and_det(A):
    n=len(A); T=[[F(x) for x in r]+[F(i==j) for j in range(n)] for i,r in enumerate(A)]; det=F(1)
    for i in range(n):
        j=next(j for j in range(i,n) if T[j][i])
        if j!=i: T[i],T[j]=T[j],T[i]; det=-det
        p=T[i][i];det*=p;T[i]=[x/p for x in T[i]]
        for j in range(n):
            if j!=i and T[j][i]:
                p=T[j][i];T[j]=[a-p*b for a,b in zip(T[j],T[i])]
    return [r[n:] for r in T],det
N=len(G);need(N==54 and all(len(r)==54 for r in G))
Gi,det=inverse_and_det(G);need(det==3**25)
H=[[3*x for x in r] for r in Gi]
need(all(H[i][j]==H[j][i] and H[i][j].denominator==1 and G[i][j]==G[j][i] for i in range(N) for j in range(N)))
need(all(H[i][i]%2==G[i][i]%2==0 for i in range(N)))
a=data['normal'];w=mv(H,a)
need(reduce(gcd,a)==1 and reduce(gcd,map(int,w))==1 and dot(a,w)==10)
U=[]
for j,x in enumerate(data['witnesses_labels_0_through_5']):
    t=dot(a,x);u=10*bil(x,G,x)-3*t*t
    need(t%10==j and u>0);U.append(u)
U+=U[-2:0:-1];need(U==[100,97,88,73,72,65,72,73,88,97])
best2=F(0)
for mask in range(1,1024,2):
    S=[i for i in range(10) if mask>>i&1]
    u=min(U[(x-y)%10] for x in S for y in S)
    r2=len(S)**2*F(u,100)**53
    if len(S)>1: best2=max(best2,r2)
    need(r2<=1)
need(best2==4*F(97,100)**53 and best2<1)
# Series truncated mod q^6, d = D/27.
def mul(x,y): return [sum((x[i] if i<len(x) else 0)*(y[k-i] if k-i<len(y) else 0) for i in range(k+1)) for k in range(6)]
def powp(x,k): return reduce(mul,[x]*k,[1])
def fac(n,p):
    z=[0]*6;z[0]=1
    for j in range(1,6//n+1):
        if j*n>5: break
        z[j*n]=(-1)**j*comb(p,j) if p>=0 and j<=p else (comb(-p+j-1,j) if p<0 else 0)
    return z
C=[1];d=[0,1]
for i in range(1,6):
    C=mul(C,mul(fac(i,9),fac(3*i,-3)))
    d=mul(d,mul(fac(i,-3),fac(3*i,9)))
T=[mul(powp(C,9-j),powp(d,j)) for j in range(10)]
def coefficient(s,j):
    b=[]
    for i in range(j+1): b.append(s[i]-sum(x*T[t][i] for t,x in enumerate(b)))
    return b[-1]
def row(j):
    return [coefficient([int(i==k) for k in range(j+1)],j) for i in range(j+1)]
need(row(4)==[193914,16016,1026,45,1] and row(5)==[398034,53463,6058,540,33,1])
need(15309000+398034==81*193914 and 9*193914-398034==1347192)
need(15309000*(1-F(10,54)*F(10,3))==5859000)
K=data['normal_dsw'];B=data['dual_scaled_48'];Ki,dk=inverse_and_det(K);need(dk==3072 and all(48*Ki[i][j]==B[i][j] for i in range(5) for j in range(5)))
for dim,goal,count,d0 in [(44,F(7,24),5,768),(45,F(1,4),4,192)]:
    axes=data['release'+str(dim)];cl=data['cluster'+str(dim)];nn=len(axes)
    M=[[bil(a,Ki,b) for b in axes] for a in axes];inv,dT=inverse_and_det(M)
    need(dT*dk==d0 and len(cl)==count)
    if nn==1: need(reduce(gcd,axes[0])==1)
    else: need(reduce(gcd,(axes[0][i]*axes[1][j]-axes[1][i]*axes[0][j] for i in range(5) for j in range(i)))==1)
    qs=[]
    for i,x in enumerate(cl):
        for y in cl[:i]:
            v=[a-b for a,b in zip(x,y)];t=[bil(v,Ki,a) for a in axes]
            qs.append(bil(v,Ki,v)-bil(t,inv,t))
    need(min(qs)>0 and max(qs)==goal)
print("PASS: scaled dual minimum and shell coefficient identity; full ten-fiber check; 44D and 45D candidate Gram and diameter checks.")
