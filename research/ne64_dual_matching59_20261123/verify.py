#!/usr/bin/env python3
"""Exact 59D Ne64 six-translate lower bound using two matched normal edges."""
from pathlib import Path
from itertools import combinations
import sympy as s
from sympy.matrices.normalforms import smith_normal_form
from math import lcm
from fractions import Fraction
import json,sys
assert __debug__
F=Path(__file__).resolve().parent
def read(f):return s.Matrix([list(map(int,r.split()))for r in (F/f).read_text().splitlines()])
def write(p,X):
 assert all(z.q==1 for z in X)
 p.write_text('\n'.join(' '.join(map(str,X.row(i)))for i in range(X.rows))+'\n')
G=s.zeros(64)
lines=(F/'ne64.triangle').read_text().splitlines();assert len(lines)==64
for i,l in enumerate(lines):
 r=list(map(int,l.split()));assert len(r)==i+1
 for j,z in enumerate(r):G[i,j]=G[j,i]=z
# Source: Nebe-Sloane Catalogue entry CQ64b. Inherit exactly its
# published minimum 12 and 3-modularity, which imply min(3 G^-1)=12.
L,P=G.LDLdecomposition(hermitian=False)
assert L*P*L.T==G and all(P[i,i]>0 for i in range(64))
assert s.prod(P[i,i]for i in range(64))==3**32
T=read('change.txt');assert T.shape==(64,64) and abs(T.det())==1
C=T*G*T.T
J=3*C.inv();assert all(z.q==1 for z in J) and J.det()==3**32
Ji=C/3
A=read('normal.txt');assert A.shape==(5,64)
K=A*Ji*A.T
raw=14*s.eye(5)-2*s.ones(5)
for i,j in ((1,2),(3,4)):raw[i,j]=raw[j,i]=-4
assert 3*K==raw and K.det()==s.Rational(32768,81)
I=[0,1,3,5,25];free=[i for i in range(64)if i not in I]
P=A[:,I];assert abs(P.det())==1
S=smith_normal_form(A,domain=s.ZZ);assert all(abs(S[i,i])==1 for i in range(5))
R=s.zeros(64,5);U=P.inv()
for j,v in enumerate(I):R[v,:]=U[j,:]
assert A*R==s.eye(5)
N=s.zeros(59,64);B=-A[:,free].T*U.T
for j,i in enumerate(I):N[:,i]=B[:,j]
for j,i in enumerate(free):N[j,i]=1
assert N*A.T==s.zeros(59,5) and N[:,free]==s.eye(59)
Q=N*J*N.T
det=Q.det();assert det==3**32*K.det()==32768*3**28
labels=s.zeros(1,5).col_join(s.eye(5))
lifts=labels*R.T;assert lifts*A.T==labels
BD=K.inv(); proj=lifts-labels*BD*A*Ji
c=proj[:,free]
assert c*N==proj and proj*A.T==s.zeros(6,5)
diam=[]
for i,j in combinations(range(6),2):
 d=labels[i,:]-labels[j,:];b=(d*BD*d.T)[0];assert 0<b<12;diam.append(b)
beta=max(diam);assert beta==s.Rational(7,16)
m=12-beta
contact=read('contact.txt');assert contact.shape==(1,64)
assert (contact*J*contact.T)[0]==12
ed=labels[0,:]-labels[4,:]
assert contact*A.T==ed and (ed*BD*ed.T)[0]==beta
contact_perp=contact-ed*BD*A*Ji
assert (contact_perp*J*contact_perp.T)[0]==m
den=lcm(*(int(x.q)for x in c))
assert all(any(v.q!=1 for v in (2*c[1,:]-c[j,:]))for j in range(6))
incumbent=s.Rational(3**33*401**59,2**296*25*13**61)
delta2=36*(m/4)**59/det
assert delta2>incumbent
volume=s.pi**s.Rational(59,2)/s.gamma(s.Rational(61,2))
val={'status':'PASS','dimension':59,'cosets':6,'normal_det':str(K.det()),'beta':str(beta),
     'period_det':str(det),'shift_denominator':den,'minimum_squared_exact':str(m),'periodic_not_lattice':True,
     'center_density_squared_floor':str(delta2),'center_density_floor':str(s.N(s.sqrt(delta2),35)),
     'volume_fraction_floor':str(s.N(s.sqrt(delta2)*volume,35)),
     'improvement_percent_floor':str(s.N((s.sqrt(delta2/incumbent)-1)*100,30))}
if len(sys.argv)>1:
 out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True)
 write(out/'ambient_gram.txt',J);write(out/'period_basis.txt',N);write(out/'period_gram.txt',Q);write(out/'shift_numerators.txt',c*den);write(out/'lifts.txt',lifts)
 (out/'results.json').write_text(json.dumps(val,indent=2)+'\n')
print(json.dumps(val,indent=2))
