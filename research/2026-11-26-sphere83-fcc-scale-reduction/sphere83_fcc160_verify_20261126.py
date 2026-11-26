"""Sphere83: self-contained sharp product FCC witness. Python >= 3.8."""
from itertools import combinations
from collections import Counter
import json
def X(q):
 return [2*((q>>(k+1))&1)+(q>>5)+4*(k==0)*(q&1)for k in range(4)]
def dd(u,v):
 t=[a-b for a,b in zip(X(u),X(v))];r=[(x+2)%4-2 for x in t]
 return (sum(x*x for x in r)+(sum((a-b)//4 for a,b in zip(t,r))%2)*min(16-8*abs(x)for x in r))//4
D=[[dd(i,j)for j in range(64)]for i in range(64)]
def out(x):return 0 if x==0 else 2 if x==7 else 1
def distance(x,y):return out((x^y)>>12)+D[(x>>6)&63][(y>>6)&63]+D[x&63][y&63]
def color(q):return 2*(q>>5)+(bin((q>>1)&15).count('1')%2)
def masks():
 def sy(a,b):return bin(a&(((b&5)<<1)|((b&10)>>1))).count('1')%2
 return {tuple(sorted((p,p^a,p^b,p^a^b))) for a in range(1,16)for b in range(a+1,16)if not sy(a,b)for p in range(16)}
def check(P):
 assert len(P)==len(set(P));ds=Counter(distance(a,b)for a,b in combinations(P,2));assert min(ds)==4
 cols=Counter(4*color(x>>6&63)+color(x&63)for x in P)
 return len(P),ds[4],len(P)-max(sum(cols[x]for x in plane)for plane in masks())
def product():
 C=set()
 for w in ('0100000011','0011111101','1100101100','0001010111'):
  z=int(w,2)
  for j in range(10):C.add(z);z=((z<<1)&1023)|(z>>9)
 assert len(C)==40 and min(bin(a^b).count('1')for a,b in combinations(C,2))==4
 P=[]
 for word in sorted(C):
  x=word&255;c=(word>>8)&1;p=(bin(x).count('1')&1)^(7*c)
  for j in range(4):P.append(4096*p+64*(2*(x&15)+(j&1))+2*(x>>4)+(j>>1))
 return sorted(P)
P=product();assert len(P)==160
print('Product FCC:',check(P))
