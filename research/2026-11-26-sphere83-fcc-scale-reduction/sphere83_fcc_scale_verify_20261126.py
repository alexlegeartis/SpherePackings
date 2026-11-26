import sys,json,base64,zlib,math,itertools
from collections import Counter
from fractions import Fraction as F
from itertools import combinations,product
def det(M):
 a,b,c=M
 return a[0]*(b[1]*c[2]-c[1]*b[2])+a[1]*(b[2]*c[0]-c[2]*b[0])+a[2]*(b[0]*c[1]-c[0]*b[1])
P=[(1,0,0,1),(0,2,0,1),(0,0,2,1),(1,0,0,4),(0,1,0,4),(0,0,1,4)]+[(*v,4) for v in product(range(3),range(5),range(5)) if any(v)]
Q=set()
for t in combinations(P,3):
 m=[list(r[:3])for r in t];d=det(m)
 if not d:continue
 a,b,c=[F(det([row[:i]+[r[3]]+row[i+1:] for row,r in zip(m,t)]),d)for i in range(3)]
 if 1<=a<=4 and F(1,2)<=b<=c<=4:Q.add((a,b,c))
return_Q=sorted(Q)
def X(q):
 return [2*((q>>(i+1))&1)+(q>>5)+4*(i==0)*(q&1) for i in range(4)]
X=[X(x) for x in range(64)]
def d(x):
 r=[(i+2)%4-2 for i in x];return (sum(i*i for i in r)+(sum((a-b)//4 for a,b in zip(x,r))%2)*min(16-8*abs(i)for i in r))//4
def cls(x):return d(X[x])
D=[cls(i) for i in range(64)]
def ch(x,y):
 s=sum(t*u for t,u in zip(X[x],X[y]));assert not(s%2);return [1,0,-1,0][(s//2)%4]
def block(a,bb,cc):
 az=[a(i)for i in range(bb)];n=len(set(az));Fs=[];size=[]
 for k in range(n):
  ks=[j for j in range(bb)if az[j]==k];size.append(len(ks)); row=[]
  for b in range(n):
   rr={sum(cc(i,j) for j in ks)for i in range(bb)if az[i]==b};assert len(rr)==1;row.append(F(rr.pop(),len(ks)))
  Fs.append(row)
 return Fs,size
I,sz=block(cls,64,ch)
outer=lambda x:0 if not x else 2 if x==7 else 1
dual=lambda x:0 if not x else (1 if x.bit_count()%2==0 else 2)
O=[];so=[]
for r in range(3):
 ks=[k for k in range(8)if dual(k)==r];so.append(len(ks));row=[]
 for b in range(3):
  z={sum((-1)**((k&i).bit_count()%2)for k in ks)for i in range(8)if outer(i)==b};assert len(z)==1;row.append(F(z.pop(),len(ks)))
 O.append(row)
lab=list(itertools.product(range(3),range(5),range(5)))
mat=[[O[k][x]*I[h][i]*I[l][j]for x,i,j in lab]for k,h,l in lab][1:]
mult=[so[k]*sz[h]*sz[l]for k,h,l in lab][1:]
scl=math.lcm(*(v.denominator for row in mat for v in row))
P=[[int(x*scl)for x in row]for row in mat]
assert [sum(w*r[u] for w,r in zip(mult,P))for u in range(75)]==[32767*scl]+[-scl]*74
r=json.loads(zlib.decompress(base64.b64decode(open(sys.argv[1]).read())));assert r['L']==scl and r['scale']==10**8
assert Q==set((F(A,d),F(B,d),F(C,d))for d,A,B,C,*_ in r['rows']) and len(Q)==len(r['rows'])==592
surv=[]
for d,A,B,C,t,cap,no,Y,oldbad in r['rows']:
 assert len(Y)==74 and min(Y)>=0
 ok=[u for u,(x,i,j)in enumerate(lab) if A*x+B*i+C*j>=4*d]
 scale=r['scale'];bad=max([0]+[sum(z*P[k][u]for k,z in enumerate(Y))+scl*scale for u in ok]);assert bad==oldbad
 N=scale*scl+scl*sum(Y)+32767*bad;Den=scale*scl
 assert t==math.isqrt(288**2*A**3*B**4*C**4//(2*d**11))+1 and cap==N//Den and no==(cap<t)
 if cap>=t:surv.append(((A,B,C,d),t,N,Den))
assert sorted((x[1],x[0])for x in surv)==[(51,(2,1,2,2)),(204,(1,1,1,1))]
for val,t,N,D in surv:
 print(val,t,'ordinary upper',F(N,D),'odd upper',F(N,D)-1)
 if t==51:assert N/D<52 and (N-D)//D<51
print('Verified all 592 rational vertices; only the unscaled size 204 can win.')
