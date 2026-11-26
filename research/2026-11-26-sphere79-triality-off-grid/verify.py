import sys,itertools,collections,math,json
def xyz(q):return [2*((q>>(i+1))&1)+(q>>5)+4*((q&1) and i==0) for i in range(4)]
# Bool times 4 -> 4 or 0
X=[xyz(q) for q in range(64)]
def norm(q):
 r=[(u+2)%4-2 for u in q];par=sum((u-v)//4 for u,v in zip(q,r))%2
 return (sum(x*x for x in r)+par*min(16-8*abs(x) for x in r))//4
D=[[norm([a-b for a,b in zip(x,y)]) for y in X]for x in X]
def enc(x):return next(j for j,y in enumerate(X) if norm([a-b for a,b in zip(x,y)])==0)
ADD=[[enc([a+b for a,b in zip(x,y)]) for y in X]for x in X]
NEG=[enc([-a for a in x]) for x in X]
def add(a,b):return ((a&~63)^(b&~63))+ADD[a%64][b%64]
def neg(a):return a//64*64+NEG[a%64]
def ds(a,b):return ((a^b)>>6).bit_count()+D[a%64][b%64]
def check(path):
 a=json.load(open(path))['motif']; A=set(a); assert len(A)==len(a)==144 and all(isinstance(v,int) and 0<=v<8192 for v in A)
 assert min(ds(i,j) for i in A for j in A if j!=i)>=4
 u=next(iter(A));H=sorted(v for y in A for v in [add(y,neg(u))] if all(add(x,v) in A for x in A))
 C=collections.Counter(ds(x,y) for x in A for y in A if x<y)
 assert H==[0] and C[4]==4308; print('Verified 144 sites, 4308 quotient pairs at wrapped squared length four, full period covolume 4096.');return H
# Append below the standalone verifier definitions; call backbone(path).
def backbone(path):
 a=json.load(open(path))['motif']; A=set(a)
 def ov(t): return sum(add(x,t) in A for x in A)
 # Ambient vectors (v,z/2) with squared norm 4. The outer signs do not
 # change their class mod 2Z^7, so keep support and squared norm only.
 outer={}
 for v in itertools.product(range(-2,3),repeat=7):
  n=sum(a*a for a in v)
  if n<=4: outer.setdefault(n,set()).add(sum(((u%2)<<i) for i,u in enumerate(v)))
 counts={}; big=set(); mids=set()
 for z in itertools.product(range(-4,5),repeat=4):
  if any((w-z[0])%2 for w in z):continue
  rr=16-sum(w*w for w in z)
  if rr<0 or rr%4:continue
  q=enc(z)
  for p in outer[rr//4]:
   k=64*p+q
   if k not in counts:counts[k]=ov(k)
   d=counts[k]
   if d>136:big.add((p,z))
   if d==136:mids.add((p,z))
 assert big=={(0,(0,0,0,0))}
 # In the latter case the only nonzero representatives of norm 2
 # are the 14 pure translations +/-2 coordinate units.
 assert len(mids)==8 and {p for p,z in mids}=={0}
 K={enc(z) for p,z in mids}
 assert len(K)==1; k=K.pop(); assert add(k,k)==0
 Z=[z for p,z in sorted(mids) if next(w for w in z if w)!=-2]
 assert len(Z)==4 and all(sum(u*v for u,v in zip(x,y))==16*int(x==y) for x in Z for y in Z)
 # The four perpendicular inner physical vectors z/2 together with seven
 # standard outer 2e_i form a cubic period with covolume 2048. They all
 # translate the motif by class k or zero. The intersection with its
 # k-translate is k-invariant and has 136 sites per cell.
 B={v for v in A if add(v,k) in A}
 assert len(B)==136 and all(add(v,k) in B for v in B)
 print('Largest subpacking invariant under any rotated cubic 2Z^11 has intensity 68/2048.')


def escape(path):
 A=set(json.load(open(path))['motif']);z=[[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]]
 # z/2 is orthogonal; physical inner coordinates X[q]/2 transform by z/2.
 assert all(sum(x*y for x,y in zip(a,b))==4*int(a==b) for a in z for b in z)
 assert all(sum(u*v for u,v in zip(r,X[a%64]))%4==0 for a in A for r in z)
 for a,j,N in [(4167,0,57),(6169,5,58)]:
  sh=[x for x in A if x!=a and ds(a,x)==4]
  assert len(sh)==N and all(((a^x)>>(6+j))&1==0 for x in sh)
  assert add(a,30) not in A
 # The displacement 2e_j to a different center label on this axis is
 # impossible since 2e_j is already a period. A nonzero component for
 # a sharp contact would therefore require a different bit in that axis.
 assert (((4167>>6)&1)-((5675>>6)&1))==1
 # Endpoint moves to the two centers of a deeper guest cubic lattice.
 B=A-{4167,6169};U=4167^(1<<6);V=6169^(1<<11)
 assert U==4103 and V==4121 and add(U,30)==V and add(V,30)==U
 for x in (U,V): assert collections.Counter(ds(x,y) for y in B)=={5:129,9:13}
 print('Off-grid displacement and deeper movable guest cubic coset verified.')
if __name__=='__main__':
 check(sys.argv[1]);backbone(sys.argv[1]);escape(sys.argv[1])
