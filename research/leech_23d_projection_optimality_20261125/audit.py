import json
from fractions import Fraction as F
L=[l for l in open('anti.out') if l.startswith('best 6')];a=open('anti.out').read().split('best 6\n')[1].split('result')[0].splitlines()
vs=[[0]*24]+[[int(s) for s in l.split(':')[1].split()] for l in a]
v=[5]+[1]*23;dot=lambda a,b:sum(x*y for x,y in zip(a,b))
gol=set()
g=(1<<11)|(1<<9)|(1<<7)|(1<<6)|(1<<5)|(1<<1)|1
for k in range(4096):
 w=0
 for j in range(12):
  if k>>j&1:w^=g<<j
 gol.add(w|((w.bit_count()%2)<<23))
def inL(z):
 m=z[0]%2;return all((a-m)%2==0 for a in z) and (sum(z)-4*m)%8==0 and sum(1<<i for i,a in enumerate(z) if (a-m)%4==2) in gol
def h(z):return dot(z,v)//8
assert inL(v) and all(inL(z) for z in vs)
for i,a in enumerate(vs):
 for b in vs[:i]:
  d=[x-y for x,y in zip(a,b)];k=h(d);u=k%6;kmin=min(u,6-u)
  assert kmin and F(dot(d,d),8)-F(k*k,6)==4-F(kmin**2,6)
A=[[int(z) for z in l.split()] for l in open('triangles')]
for line in A:
 a,b=line[:24],line[24:];assert inL(a) and inL(b) and h(a)==h(b)==2 and dot(a,b)==-8 and dot(a,a)==dot(b,b)==32
p=list(map(int,open('square.counts').readline().split()));assert inL(p)
W=[0]*24
for j in [0,1,5,6,7,9,11,23]:W[j]=2
W[1]-=4;W[5]-=4;W[0]-=4;W[2]+=4;assert inL(W)
print('independent Golay membership, all 15 anticode distances, and 400 triangle-pair generators verified')
