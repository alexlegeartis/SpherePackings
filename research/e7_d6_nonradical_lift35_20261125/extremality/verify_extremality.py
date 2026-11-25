import hashlib,base64,zlib,json,sys,numpy as np,scipy.linalg as la
from pathlib import Path
here=Path(__file__).resolve().parent
pack=json.loads((here/'extreme35.json').read_text());p=json.loads((here/'../../e7_d6_lattice34_20261125/cert.json').read_text())
if 'parent' in p and isinstance(p['parent'],dict):p=p['parent']
V=zlib.decompress(base64.b64decode(pack['zlib64']));assert hashlib.sha256(V).hexdigest()==pack['raw_sha256']
V=np.frombuffer(V,dtype='i1').reshape(630,35).astype('i8')
g=np.array(p['E7_gram']);alph=p['alpha'][-4:];inv=np.argsort(p['phase_isometry']);tags=p['phase_syndromes'];u=p['offset_syndrome'];M=[21,60,41,40];B=[(3 if i==0 else 1 if i==5 else 0) for i in range(6)]
def gf(a,b):
 v=0
 while b:
  if b&1:v^=a
  a*=2
  if a&64:a^=67
  b//=2
 return v
def sy(j,t):
 x=t%64
 return x+64*gf(alph[j],x)+4096*(t//64)
def xm(bits,vs):
 a=0
 for j,v in enumerate(vs):
  if bits>>j&1:a^=v
 return a
def belongs(v):
 w=int(v[-1])
 for m in range(2):
  yy=[int(x)-w*b+m for x,b in zip(v[:6],B)]
  if any(k%2 for k in yy):continue
  yy=[k//2 for k in yy]
  if sum(yy)%2:continue
  s=((w%2)*896)^(m*u)^xm(sum((k%2)<<i for i,k in enumerate(yy[:5]))+32*((sum(yy)//2)%2),tags)
  ok=True
  for j in range(4):
   k=[int(x)-w*(M[j]>>i&1)-m*(j>0)*(82>>i&1) for i,x in enumerate(v[6+7*j:13+7*j])]
   if any(k0%2 for k0 in k):ok=False;break
   mask=sum(((z//2)&1)<<i for i,z in enumerate(k)); s^=sy(j,int(inv[mask]))
  if ok and s==0:return True
 return False
assert all(belongs(v) for v in V)
G=la.block_diag(np.eye(6,dtype='i8'),*([g]*4),np.array([[6]],'i8'))
assert np.all(np.einsum('vi,ij,vj->v',V,G,V)==32)
bmat=np.rint(6*la.inv(G)).astype('i8');assert np.all(bmat@G==6*np.eye(35,dtype=int))
i,j=np.tril_indices(35);A=(V[:,i]*V[:,j]).T.copy();N=630;b=bmat[i,j]
T=2**26;W=2**36
B0=np.rint(T*la.inv(A)).astype('i8');w0=np.rint(W*la.solve(A,b)).astype('i8')
# These two bounds prove that every fixed-width multiplication below is exact.
amax=int(abs(A).max());assert int(abs(B0).max())*amax*N+T<2**63
assert int(abs(w0).max())*amax*N+int(abs(b).max())*W<2**63
E=T*np.eye(N,dtype='i8')-B0@A;r0=W*b-A@w0
assert int(abs(E).max())*N<2**63 and int(abs(B0).max())*N<2**63
e=int(abs(E).sum(1).max());k=int(abs(B0).sum(1).max());r=int(abs(r0).max());m=int(w0.min())
assert e<T and (T-e)*m>k*r
print('630 independent contacts, positive eutaxy:',T,e,W,m,k,r)
