"""Verify exact atlas: python3 atlas.py (requires Python with sympy, numpy, Node and C++17)."""
import json, os, gzip, subprocess, tempfile, sys
from pathlib import Path
import numpy as np
import sympy as S
here=Path(__file__).resolve().parent
def matrix(fn): return S.Matrix(json.loads((here/fn).read_text()))
def run(args,wd=None): subprocess.run(list(map(str,args)),cwd=wd,check=True)
with tempfile.TemporaryDirectory() as tmp:
 w=Path(tmp)
 run(["node",here/"glue.js",w/"q.dat"])
 G=matrix("gram.json");T=matrix("transform.json"); assert abs(T.det())==1
 J=30*G.inv();assert all(z.is_Integer for z in J)
 M=T*J*T.T; V=T*J; N=30*M.inv(); assert all(z.is_Integer for z in N)
 L,D=M.LDLdecomposition(hermitian=False)
 # Fixed-point LOWER BOUND to exact Q in every branch, |v_i|<=H from Cauchy Schwarz.
 from math import isqrt
 n=39;cut=100;scale=16384;H=[isqrt((111*int(N[i,i]))//30) for i in range(n)]
 C=[[int(L[j,i]*scale) if j>i else 0 for j in range(n)] for i in range(n)]
 low=[int(D[i,i]*scale) for i in range(n)]
 err=[sum(H[i+1:]) for i in range(n)]
 assert min(low)>0 and all(abs(L[j,i]*scale-C[i][j])<1 for i in range(n) for j in range(i+1,n))
 # Each accepted or first rejected nearest-out integer has displacement <= sqrt(B/d)+2S+E in scaled units. Enforce int limits before compiling.
 for i in range(n):
  t=isqrt((cut*scale**3)//low[i])+2*scale+max(err)+2
  assert t*t*low[i]<2**63
  assert sum(abs(C[i][j])*H[j] for j in range(i+1,n))+t<2**30
 assert sum(abs(int(M[i,j]))*H[i]*H[j] for i in range(n) for j in range(n))<2**30
 values=[list(M),list(V),[C[i][j] for j in range(n) for i in range(n)],low,H,err]
 (w/"exactdual.dat").write_text("\n".join(" ".join(map(str,row)) for row in values))
 run([os.environ.get("CXX","c++"),"-O3","-std=c++17",here/"dual.cpp","-o",w/"dual"])
 run([w/"dual",cut],w)
 G=np.asarray(G,dtype=np.int64);gens=np.asarray(matrix("gens.json"),dtype=np.int64);ids=np.arange(21870);pc=ids[:,None]%10*gens[7]
 for i in range(7):pc+=(ids//(10*3**i)%3)[:,None]*gens[i]
 pc%=30;lookup={bytes(x.astype('u1')):k for k,x in enumerate(pc)};qq=np.fromfile(w/'q.dat',dtype='u1').astype(int)
 z=np.frombuffer(gzip.decompress((here/'witness.gz').read_bytes()),dtype='i1').reshape((-1,39))
 W=30*z.astype(np.int64)+pc[1:];R=(W@G*W).sum(axis=1);assert np.all(R%30==0)
 ub=np.r_[0,R//30];raw=(w/'low.dat').read_bytes();assert len(raw)%43==0
 for k in range(0,len(raw),43):
  h=raw[k:k+39];c=int.from_bytes(raw[k+39:k+43],sys.byteorder);i=lookup[h];j=lookup[bytes((30-u)%30 for u in h)];ub[i]=min(ub[i],c);ub[j]=min(ub[j],c)
 # For a conjugate gluing class the lower bound is 180-q, allowed squared norms rise by 60 in these units. If it is above cut, a witness attains it. Otherwise the exhaustive shell settles attainment/nonattainment.
 u,q=ub[1:],qq[1:]; assert np.all((u+q==180) | ((u+q==240)&(180-q<=cut)))
 print('true lifts exactly 2 farther in squared norm:',{int(a):int(((q==a)&(u+q==240)).sum()) for a in sorted(set(q)) if np.any((q==a)&(u+q==240))})
 assert int((u+q==240).sum())==1475 and int(max(u))==171 and (u==171).sum()==2
 print('verified all 21,869 nonzero true class minima.')

