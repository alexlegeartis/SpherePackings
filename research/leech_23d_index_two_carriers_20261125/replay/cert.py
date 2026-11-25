import sys
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog
import sympy as sp
d={i:[] for i in range(64)};A=[]
for s in open('parts.txt'):
 if s.startswith('part'):
  _,c,o,e=s.split();d[int(c)].append((int(o),int(e)))
 else:A.append(list(map(int,s.split())))
def mul(c,n):
 v=0
 for j in range(n%16):v=A[v][c]
 return v
C=[0,32,5];rows=[];names=[]
def ins(a,h):
 if a not in rows:rows.append(a);names.append(h)
for i in range(3):
 for j in range(i,3):
  for u in range(-12,13):
   for v in range(-5,6):
    if i==j and not(u or v):continue
    x,y=F(3*u+j-i,3),F(3*v+j-i,3);c=A[A[C[j]][mul(C[i],-1)]][A[mul(4,u)][mul(33,v)]]
    for o,e in d[c]:
     val=(4*x*x+4*x*y+37*y*y+o+e)/16;assert val>=4
     if val!=4:continue
     row=[x*x,y*y,2*x*y,0,0,0,0,F(o,16),F(e,16),-1]
     for h,sg in [(i,-1),(j,1)]:
      if h:row[3+2*(h-1)]+=sg*(8*x+4*y)/16;row[4+2*(h-1)]+=sg*(4*x+74*y)/16
     ins(tuple(row),(i,j,u,v,o,e))
ins((0,0,0,0,0,0,0,4,0,-1),'pure-octad');ins((0,0,0,0,0,0,0,0,4,-1),'pure-exterior')
f=(F(-37,18),F(-2,9),F(2,9),0,0,0,0,F(-7,2),-7,F(23,8));R=np.asarray(rows,dtype=float);ff=np.asarray(f,dtype=float)
s=linprog(-ff,A_ub=-R,b_ub=np.zeros(len(R)),bounds=[(-1,1)]*len(f),method='highs')
print('primal',len(rows),-s.fun,s.x)
s=linprog(np.zeros(len(R)),A_eq=R.T,b_eq=-ff,bounds=[(0,None)]*len(R),method='highs')
if s.success:
 I=np.nonzero(s.x>1e-8)[0];M=sp.Matrix([[sp.Rational(rows[i][j]) for i in I] for j in range(len(f))]);y=sp.Matrix([-sp.Rational(x) for x in f]);z=M.gauss_jordan_solve(y)[0];assert M*z==y and all(t>=0 for t in z)
 print('exact dual',list(zip((names[i] for i in I),z)))
