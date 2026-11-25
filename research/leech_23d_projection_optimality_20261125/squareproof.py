from itertools import combinations
from fractions import Fraction as F
from math import comb
out=[list(map(int,l.split())) for l in open('square.counts')]
a=out[0];V=[5]+[1]*23;W=[0]*24
for j in [0,1,5,6,7,9,11,23]:W[j]=2
W[0]-=4;W[1]-=4;W[5]-=4;W[2]+=4
dot=lambda x,y:sum(a*b for a,b in zip(x,y))
assert (dot(a,a),dot(a,V),dot(a,W))==(48,24,24)
for i,u in enumerate(range(-3,4)):
 for j,v in enumerate(range(-3,4)):
  assert (out[i+1][j]>0)==(abs(u)!=3 or abs(v)!=3)
def D(a,b):
 x=(a%6-b%6)%6;y=(a//6-b//6)%6
 x=min(x,6-x);y=min(y,6-y)
 return 24 if (x,y)==(0,0) else (6 if (x,y)==(3,3) else x*x+y*y)
def add(a,b):return ((a%6+b%6)%6)+6*((a//6+b//6)%6)
def cap(c):
 pu=[i for i in range(1,36) if D(i,0)<c];best=1;out=[]
 for k in range(1,len(pu)+1):
  if k>best:break
  for ls in combinations(pu,k):
   if all(D(a,b)<c for a,b in combinations(ls,2)):
    best=k+1;out.append(set(ls)|{0})
 return best,[l for l in out if len(l)==best]
assert cap(2)[0]==2
sz,S=cap(5);assert sz==5
cross={0,1,5,6,30};crosses=[{add(k,j) for j in cross} for k in range(36)]
assert all(A in crosses for A in S)
b=F(13,8)
assert F(25,4)/(5-b)*((24-b)/23)**23<1
assert F(54,4)*(F(5,6)**22)/(F(23,24)**23)<1 # squared ratio of 3*(5/6)^11 to CHLWW
assert (F(47,46)**23)<2 # (1+.5/23)^23 / (4*.5)<1
print('full distance table, exceptional vector, 5-point cap, rational endpoints verified')

