#!/usr/bin/env python3
"""Exact scoped four-Leech 96D obstruction. Python standard library and g++ C++20."""
import gzip,hashlib,json,subprocess,tempfile
from fractions import Fraction as F
from pathlib import Path
assert __debug__
p=Path(__file__).resolve().parent
def pin(b,h):
 assert hashlib.sha256(b).hexdigest()==h
 return b
data=pin(gzip.decompress((p/'input.txt.gz').read_bytes()),'25e5c6ed1f425c27c6c051ee54df420e5a11537ec48561db1f48fc50329c91ad')
lines=data.decode().splitlines()
assert len(lines)==4120 and int(lines[24])==4095
Fmat=[int(x)for x in lines[:24]]
fib={r[0]:r[1:] for line in lines[25:] for r in [list(map(int,line.split()))]}
def space(H):
 U={0};basis=[]
 for a in sorted(H):
  if a not in U:U|={a^x for x in U};basis.append(a)
 return U,basis
A,ba=space(fib)
B,bb=space(set().union(*map(set,fib.values())))
assert A=={0,*fib} and len(A)==len(B)==4096 and len(ba)==len(bb)==12
A=list(sorted(A)); AI={a:i for i,a in enumerate(A)}
def right(a,b):
 v=0
 for i in range(24):
  if a>>i&1:v^=Fmat[i]
 return (v&b).bit_count()%2
assert all(right(a,d)==right(b,d)==0 for d in ba for a,b in zip(ba,ba))
assert all(right(a,b)==0 for a in bb for b in bb)
assert len({tuple(right(a,b) for b in bb) for a in A})==4096
assert all(AI[a^b]==i^AI[b] for i,a in enumerate(A)for b in ba)
assert A[0]==0
q4=set()
for a,L in fib.items():
 assert len(set(L))==24 and all(b in B and right(a,b)==0 for b in L)
 q4.update(a^b for b in L)
assert len(q4)==98280
def writeints(X):
 return ''.join(f'{q}\n'for q in sorted(X)).encode()
with tempfile.TemporaryDirectory() as tmp:
 t=Path(tmp)
 for name in ('affine_planes','walsh'):
  subprocess.run(['g++','-O3','-std=c++20',str(p/(name+'.cpp')),'-o',str(t/name)],check=True,cwd=t)
 subprocess.run([str(t/'affine_planes')],input=data,check=True,cwd=t)
 pl=pin((t/'admissible_affine.txt').read_bytes(),'3be19b30a638424f3dccd1e61caafb3f039830c97ffa0800b6d89f6ad4865e37')
 P=[list(map(int,line.split()))for line in pl.splitlines()]
 assert len(P)==975 and all(a^b==c for a,b,c in P)
 C={AI[a]|(AI[b]<<12) for S in P for a in S for b in S if a!=b}
 assert len(C)==5850
 subprocess.run([str(t/'walsh'),'source'],input=writeints(q4),check=True,cwd=t)
 subprocess.run([str(t/'walsh'),'target'],input=writeints(C),check=True,cwd=t)
# PSD bound on a selected vertical subset. A ray for each q4 edge, both orientations.
base=F(98280+24,2**24);assert base==F(3,512)
color_ceil=F(24*(1+15),1)/base;assert color_ceil==2**16
mu6=16*(F(6,4)**48)/(F(6,16)**12)
mu5=2**12*(F(5,4)**48)/(F(5,16)**12)
mu4=2**24*(F(4,4)**48)/(F(4,16)**12)
mu3=2**36*(F(3,4)**48)/(F(3,16)**12)
middle=F(2**36)*color_ceil
relaxed=F(2)**48*F(4,3)**12
bound=max(mu6,mu5,mu4,mu3,middle,relaxed)
assert bound==relaxed
comparison=json.loads(pin((p/'pinned_comparison.json').read_bytes(),
 '3da0691c99ab510a124930a7771972497b0047dc7e62e1f91ab6219e0c78a131'))['rows']
prior=next(F(r['center_density_squared'])for r in comparison if r['dimension']==96 and r['kind']=='periodic existence')
assert bound**2<prior
print('PASS N cap',color_ceil,'max center ceiling',bound,'ceil/prior ~',float(bound)/(float(prior)**.5))
