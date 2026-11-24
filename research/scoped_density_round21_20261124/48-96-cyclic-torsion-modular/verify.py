from fractions import Fraction as F
from math import comb
N=6
def mon(j=0): return [int(i==j) for i in range(N)]
def mul(a,b): return [sum(a[j]*b[i-j] for j in range(i+1)) for i in range(N)]
def pw(a,b):
    out=mon()
    for _ in range(b):out=mul(out,a)
    return out
def ep(a,b,p):
    out=mon()
    for d in range(1,N):
        k=a+b*int(d%p==0)
        w=[0]*N
        for j in range((N-1)//d+1):
            w[j*d]=((-1)**j*comb(k,j) if j<=k else 0) if k>=0 else comb(-k+j-1,j)
        out=mul(out,w)
    return out
def norm(basis,C):
    out=[0]*N;cs=[]
    for i,b in enumerate(basis):
        assert b[i]==1 and not any(b[:i])
        c=(C if not i else 0)-out[i];cs.append(c)
        out=[x+c*y for x,y in zip(out,b)]
    assert out[:len(cs)]==[C]+[0]*(len(cs)-1)
    assert cs[-1]>0
    return cs
def val(q,p):
    q=F(q);a=q.numerator;b=q.denominator
    assert a!=0
    e=0
    while a%p==0:e+=1;a//=p
    while b%p==0:e-=1;b//=p
    return e
for p in (3,5,7,13):
 r=(p-1)//2
 Dstep=F(p*p-1,24)
 eta_step=24//(p-1)
 for s in (1,2,3,4):
    def pbasis(J,i):
       return mul(mon(i),ep(p*J-eta_step*i,-J+eta_step*i,p))
    J=24*s//(p-1) # half-weight exponent total is (p-1)*J; form weight (p-1)*J/2 =12s
    R=norm([pbasis(J,i) for i in range(s+1)],r)
    P=norm([pbasis(J*r,i) for i in range(s+1)],1)
    assert J*Dstep==(p+1)*s
    U=J*r*Dstep-s
    e=F(J*r-eta_step*s,2); assert e.denominator==1
    product=P[-1]*p**int(e)
    if p==3:
       assert R[-1]%r==0
       print(p,24*s,"isotropic each",R[-1]//r)
       continue
    h=max(h for h in range(1,9) if F(p*p-1,24)*(12*s+h)<U+r*(s+1))
    n=24*s
    gamma={
      1:F(1),
      2:F(n+1,n),
      3:F(n+3,n),
      4:F((n+3)*(n+5),n*(n+2)),
      5:F((n+5)*(n+7),n*(n+2)),
    }[h]
    need=r*(2*h-val(gamma,p))
    actual=val(product,p)
    assert (need>actual)==(p>=7 or p==5 and s>=2)
    assert R[-1]%r==0
    e0=F(J-eta_step*s,2);assert e0==0
    print('P=',p,'N=',n,'iso combined',R[-1],'isotropic each',R[-1]//r,'anisotropic sum cap',F(2*(U-int(need>actual)),p),'raw product vP',actual,'needed',need,'harm design',h,flush=True)
for s in (1,2,3,4):
 n=24*s
 p=13
 n_orb=6
 for bsize in (2,3):
  # Degree-t products over cosets of the unique t-subgroup of F13*/±1.
  J=2*s*bsize
  e=J//2-s
  r=n_orb//bsize
  bs=lambda i:mul(mon(i),ep(13*J-2*i,-J+2*i,13))
  cs=norm([bs(i)for i in range(s+1)],r)
  rawcap=14*s*bsize-s
  assert cs[-1]%r==0
  count=cs[-1]//r*13**e
  h=max(h for h in range(1,9) if 7*(12*s+h)<r*rawcap+n_orb*(s+1))
  gamma={1:F(1),2:F(n+1,n),3:F(n+3,n)}[h]
  need=bsize*(2*h-val(gamma,13))
  assert val(count,13)<need
  print(n,'p13 block',bsize,'some block sum cap',F(2*(rawcap-1),13),'count at equality',count,'v13',val(count,13),'needs',need)
# Uniform residue obstruction at every positive genus level.
# E(q)=(q/Delta(q))^s; residue uses E_(s-1) a1+ ... + E_0 a_s.
for s in (1,2,3,4):
    E=ep(-24*s,0,13)
    B={2:3,3:4,4:6} # upper bounds on |a_j| by Deligne
    assert E[s-1]>sum(E[s-j]*B[j] for j in range(2,s+1))
    print("residue",s,"a1 coefficient",E[s-1],"rest bound",sum(E[s-j]*B[j] for j in range(2,s+1)))

# Level 6 mixed cusps: exact finite envelope.
M=((6,-3,-2,1),(-3,6,1,-2),(-2,1,6,-3),(1,-2,-3,6))
def eta6(v):
    arr=[sum(x*y for x,y in zip(row,v)) for row in M]
    assert all(x.denominator==1 for x in arr)
    out=mon()
    for d in range(1,N):
        k=int(sum(x for x,b in zip(arr,(1,2,3,6)) if d%b==0))
        w=[0]*N
        for j in range((N-1)//d+1):
            w[j*d]=((-1)**j*comb(k,j) if j<=k else 0) if k>=0 else comb(-k+j-1,j)
        out=mul(out,w)
    return mul(mon(int(v[-1])),out)
C=(29,1477,84539,5096373)
for s in (1,2,3,4):
    verified=0
    for r in range(6):
        rs=[F(r%k,k) for k in (6,3,2)]+[F(0)]
        for u3 in range(3*s+1-(r%3!=0)):
            for u2 in range(2*s+(r%2==0)):
                v3=F(u3)+rs[1];v2=F(u2)+rs[2]
                W=12*s-v3-v2
                z=norm([eta6((W-i,v3,v2,F(i))) for i in range(s+1)],1)
                assert z[-1]>0
                verified+=1
    z=norm([eta6((F(7*s-i),F(3*s),F(2*s),F(i))) for i in range(s+1)],1)
    assert z[-1]==C[s-1]
    v=(6*s,3*s,2*s,s)
    assert [sum(x*y for x,y in zip(row,v)) for row in M]==[24*s,0,0,0]
    print(24*s,"sixth verified residue/depth rows",verified,"equal nearest cardinality",z[-1])
