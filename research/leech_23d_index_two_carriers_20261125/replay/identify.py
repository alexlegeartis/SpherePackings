import numpy as np
Z=np.loadtxt('group.txt',dtype=int);w=Z[:,0];A=Z[:,1:];neg=np.argmin(A,axis=1)
def mul(c,k):
 v=0
 for i in range(k%16):v=A[v,c]
 return v
def canon(s):
 x,y,c=s;u,x=divmod(x,6);v,y=divmod(y,6);return x,y,int(A[c,A[mul(1,-u),mul(6,-v)]])
def add(s,t):return canon((s[0]+t[0],s[1]+t[1],A[s[2],t[2]]))
S={canon((x,y,c)) for x,y,c in [(0,0,0),(0,0,208),(3,3,96),(3,3,176),(3,1,124),(3,1,172),(6,4,3),(6,4,211)]};H=sorted(h for h in S if {add(s,h) for s in S}==S);d=sorted(S-set(H))[0]
print('H',H,'index',len(S)//len(H),'d',d)
g=[(6,0,1,0),(0,6,6,0)]+[(*h,0) for h in H]+[(*d,1)]
def plus(g,h):return g[0]+h[0],g[1]+h[1],A[g[2],h[2]],g[3]+h[3]
def no(v):x,y,c,j=v;return 2*x*x+3*y*y+9*w[c]+24*j*j
print('generator norms x144',*[no(a) for a in g]); print('bad pairs',[(i,j,no(plus(a,b))%288) for i,a in enumerate(g) for j,b in enumerate(g) if no(plus(a,b))%288]);assert all(no(a)%288==0 and all(no(plus(a,b))%288==0 for b in g) for a in g)
sixd=canon((6*d[0],6*d[1],mul(d[2],6)));print('6d',sixd,'in H',sixd in H);assert sixd in H and len(H)==4
for j in range(4):
 dh=canon((j*d[0],j*d[1],mul(d[2],j)));mi=999999
 for h in H:
  x,y,c=add(h,dh)
  for i in range(-8,9):
   for k in range(-8,9):
    v=(x+6*i,y+6*k,A[c,A[mul(1,i),mul(6,k)]],j);z=no(v)
    if z==0:z=576
    mi=min(mi,z)
 print('height',j,'minimum full 24D norm',mi,'/144');assert mi>=576
