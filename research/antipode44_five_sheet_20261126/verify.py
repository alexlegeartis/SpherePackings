"""Exact small-rank and all-superbase checks for the 44D five-sheet result."""
from fractions import Fraction as F
from math import prod
from itertools import combinations

def require(x):
    if not x: raise ArithmeticError("44D certificate check failed")

ref=F(25,605)*F(157,110)**44
rad=lambda b: ((6-b)/4)**44
require(F(25,324)*rad(F(2,3))<ref)
bound=lambda b: F(25,2)*b**3*rad(b)/ref
require(bound(F(18,47))<F(831,100) and bound(F(1,2))<8)
# All r<4 bounds in squared center density, including M in the constant.
require(F(1,324)*rad(F(0))<ref)
for r,c2 in [(1,F(4,108)),(2,F(16*3,4*27)),(2,F(25*4,9*27)),(3,F(16,2*6)),(3,F(25,4*6))]:
    b=F(6*r,44+r); require(c2*b**r*rad(b)<ref)

def d3(A):
    (a,b,c),(d,e,f),(g,h,i)=A
    return a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)
def rem(A,i,j):
    return [[v for l,v in enumerate(row) if l!=j] for k,row in enumerate(A) if k!=i]
def co(A,i,j): return (-1)**(i+j)*d3(rem(A,i,j))
def lim(x,y): return (x+y-6)//2
def choices(x,y): return range(-lim(x,y),lim(x,y)+1)
counts=[];eq=0
for n6 in range(6):
    diag=[6]*n6+[8]*(5-n6);a,b,c,d,e=diag;num=0
    for p in choices(a,b):
     for q in choices(a,c):
      for r in choices(a,d):
       if abs(a+p+q+r)>lim(a,e): continue
       for s in choices(b,c):
        for t in choices(b,d):
         z=(e-a-b-c-d)//2-p-q-r-s-t
         if abs(z)>lim(c,d) or abs(b+p+s+t)>lim(b,e) or abs(c+q+s+z)>lim(c,e) or abs(d+r+t+z)>lim(d,e): continue
         A=[[a,p,q,r],[p,b,s,t],[q,s,c,z],[r,t,z,d]]
         if a*b<=p*p or co(A,3,3)<=0: continue
         C=[[co(A,i,j) for j in range(4)] for i in range(4)]
         det=sum(A[0][j]*C[0][j] for j in range(4))
         if det<=0: continue
         num+=1
         mx=max([C[i][i] for i in range(4)]+[C[i][i]+C[j][j]-2*C[i][j] for i,j in combinations(range(4),2)])
         if mx>=6*det: continue
         lhs=(6*det-mx)**44*605**45;rhs=3454**44*det**45
         require(lhs<=rhs)
         if lhs==rhs:
             eq+=1;require(n6==5 and det==605 and mx==176)
             Gram=[row+[-sum(row)] for row in A]
             Gram.append([-sum(row) for row in A]+[e])
             require(all(sorted(Gram[i][:i]+Gram[i][i+1:])==[-2,-2,-1,-1] for i in range(5)))
    counts.append(num)
require(counts==[10477,6588,3981,2305,1269,654] and sum(counts)==25274 and eq==12)
print("PASS: all ranks of M<=5, 25,274 superbase matrices, 12 labelled equalities, no higher certificate.")
