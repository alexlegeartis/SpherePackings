from joint import *
def advance(v,k):
 q,p,c,s=v
 assert k in (1,-1)
 if k==1:return (3*q+s)/2,c,p+c,(5*q+3*s)/2
 return (3*q-s)/2,c-p,p,(-5*q+3*s)/2
def shift(z,u,v):
 q,p,c,s=z
 return (q+10*(u*u+v*v)+2*(p*u+c*v),p+10*u,c+10*v,
          s-10*u*u+40*u*v+10*v*v+2*((2*c-p)*u+(2*p+c)*v))
def unroll(z):
 q,p,c,s=z
 g=q-(p*p+c*c)/10
 if not g:return # fixed exceptional orbit in the anchor plane
 assert g>0
 yield z
 for k in (-1,1):
  v=z
  for j in range(1,100):
   v=advance(v,k)
   q,p,c,s=v
   h=q-(p*p+c*c)/10
   a=s-(-p*p+4*p*c+c*c)/10
   if h>12 and k*a>=0: break
   yield v
  else:raise ValueError(z)
def extra_rows():
 R,V,M,A,B=jointrows()
 K=A+B; H={q:i for i,b in enumerate(K)for p,c,s in b if (q:=(10 if i<len(A)else 12,p,c,s))}
 n0=len(R)
 Z=set()
 for i,b in enumerate(K):
  for p,c,s in b:
   y=(F(10 if i<len(A)else 12),F(p),F(c),F(s))
   for z in unroll(y):
    p0,c0=map(int,z[1:3])
    assert z[1:3]==(p0,c0)
    for u,v in product(range((-10-p0+9)//10,(10-p0)//10+1),range((-10-c0+9)//10,(10-c0)//10+1)):
     w=shift(z,u,v)
     if w[0] not in (10,12):continue
     assert all(t.denominator==1 for t in w)
     j=H.get(tuple(map(int,w)),-1)
     if i!=j:Z.add(tuple(sorted((i,j))))
 for i,j in sorted(Z):
  out=[F(0)]*len(K);out[j]=1
  if i>=0:out[i]=-1
  R.append(out);V.append(0);M.append(('affine',i,j))
 return R,V,M,A,B
if __name__=='__main__':
 R,V,M,A,B=extra_rows()
 import json
 Y={int(k):F(v)for k,v in json.loads(Path('joint-counts.json').read_text()).items()}
 fail=[m for r,b,m in zip(R,V,M) if sum(r[i]*v for i,v in Y.items())!=b]
 print('extra rows',len(R)-652,'failed rows',len(fail), fail[:15])
