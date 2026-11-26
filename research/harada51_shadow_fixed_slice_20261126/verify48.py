"""Exact 48D fixed-space complement antipode obstruction; ignores all floating searches."""
from pathlib import Path
import itertools,json
from fractions import Fraction as F
from verify import A,dot,isnum,require
P=Path(__file__).parent
N=lambda i:divmod(i,25)
def ident(x,y):return 25*(x%25)+y%25
def difference(i,j):
 a,b=N(i);c,d=N(j);return ident(a-c,b-d)
D=json.loads((P/'certificate48.json').read_text())
U={}
require(len(D['witnesses'])==307)
for k,w in D['witnesses'].items():
 k=int(k);require(len(w)==51 and isnum(w) and dot(w,w)%18==0)
 a=sum(w[:25]);t=w[25];b=sum(w[26:])
 require(ident(a,b)==k and k>0)
 require((t-a)%3==0 and (b-t+2*a)%9==0 and (a+t+b)%2==0)
 h=25*dot(w,w)-a*a-b*b-25*t*t
 require(0<h<1350 and h%18==0 and k not in U)
 U[k]=U[ident(-a,-b)]=h
require(len(U)==614)
# Mathematical notation: u labels have actual squared distance witnessed at most U[u]/225.
def forbidden(i,j,lo):return i!=j and U.get(difference(i,j),1350)<=lo
old=0;worst=F(0);cut=None
for lo,hi,data in D['partition']:
 require(lo==old and lo<hi<=1350 and data)
 old=hi
 if isinstance(data[0],list):
  # Each color is an independent set in the allowed graph; must cover its entire
  # 0-neighborhood, including unwitnessed classes.
  W=[i for z in data for i in z]
  require(all(z for z in data))
  require(sorted(W)==[i for i in range(1,625) if not forbidden(i,0,lo)])
  for z in data:
   require(all(forbidden(i,j,lo) for i,j in itertools.combinations(z,2)))
  M=1+len(data)
 else:
  require(len(data)==len(set(data)) and all(0<=i<625 for i in data))
  require(all(forbidden(i,j,lo) for i,j in itertools.combinations(data,2)))
  M=625//len(data)
 x=F(M,25)*F(hi,1350)**24
 require(x<F(13,100))
 if x>worst:worst=x;cut=(lo,hi,M)
require(old==1350)
print("48D fixed-space complement: upper ratio <",float(worst),"at integer interval",cut,"via",len(D['partition']),"exact interval certificates")
print("At 6, the allowed graph has a verified one-colour covering of the anchor neighborhood: at most two translates")


# Independent branch-and-bound improving the last intervals (no trusted output file).
# Anchor the motif at 0, and solve maximum clique in its allowed-neighbor graph.
def exact_max(lo):
 N0=[i for i in range(1,625) if not forbidden(i,0,lo)]
 Adj=[sum(1<<j for j,o in enumerate(N0) if i!=o and not forbidden(i,o,lo)) for i in N0]
 best=0
 def search(mask,size):
  nonlocal best
  if not mask:
   best=max(best,size);return
  c=0;pending=mask;order=[]
  # Greedily partition the current graph into independent sets; every prefix
  # of the resulting order needs at most its final color label colors.
  while pending:
   c+=1;q=pending
   while q:
    z=q&-q;i=z.bit_length()-1;q^=z;pending^=z;q&=~Adj[i];order.append((i,c))
  for i,c in reversed(order):
   if size+c<=best:return
   search(mask&Adj[i],size+1);mask^=1<<i
 search((1<<len(N0))-1,0)
 return 1+best
KM=[0,0,-2,-1,1,3,0,0,0,0,1,0,1,0,-1,0,1,0,-1,1,0,-2,0,-1,0,0,-1,-1,0,1,3,1,0,1,-1,0,0,0,0,-1,-1,0,0,-1,1,0,-1,-2,0,0,2]
require(isnum(KM) and dot(KM,KM)==54 and KM[25]==sum(KM[:25])==sum(KM[26:])==0)
bestcuts={lo:exact_max(lo) for lo in range(1188,1350,18)}
require(bestcuts=={1188:9,1206:7,1224:7,1242:7,1260:5,1278:4,1296:4,1314:2,1332:2})
old0=F(0);gates={}
for lo,hi,data in D['partition']:
 M=1+len(data) if isinstance(data[0],list) else 625//len(data)
 if lo in bestcuts:M=min(M,bestcuts[lo])
 gates[lo]=M
 x=F(M,25)*F(hi,1350)**24
 require(x<F(9,100));old0=max(old0,x)
print("48D refined ratio bound < 9%; exact certified maximum of upper expressions",float(old0))
# 49D: an interval of width sqrt(mu-beta) can host at most its 48D motif cap M.
# Live 49D comparator is (3/2)^(49/2)/sqrt(2), so square the ratio.
u=F(51,10)
lower=F(2,625)*625**2/u*(u/6)**49
require(lower<F(1,3)**2)
worst1=(F(0),None)
for m in range(5100,6000):
 l=F(m,1000);u=F(m+1,1000)
 g,b=min((F(2,625)*M*M*F(1,l-F(k,225))*(u/6)**49,k) for k,M in gates.items() if F(k,225)<l)
 require(g<F(1,3)**2)
 worst1=max(worst1,(g,m))
print("49D arbitrary axial relayering: certified ratio",float(worst1[0])**.5,"or lower-radius ratio",float(lower)**.5,"both < 1/3")


from verify import a as integral_upper, N as grid
from math import isqrt
# Doubles-only planar screen at full squared minimum 6.
I=F(integral_upper[isqrt(2*grid*grid//3)],grid*grid)
Area=F(3,50)*(F(314159,100000)-8*I)
r2=F(12,625)/Area**2
require(r2<F(91,100)**2)
print("50D two fibres at each normal site and squared minimum six: cap",float(r2)**.5,"< 0.91")


from verify import R
require(sum(R)==107)
P=((1,1,8),(0,3,3),(0,0,18))
require(all((t-a)%3==0 and (b-t+2*a)%9==0 and (a+t+b)%2==0 for a,t,b in P))
require({ident(a,8*a+18*c) for a in range(25) for c in range(25)}==set(range(625)))
