from pathlib import Path
from fractions import Fraction as Q
from subprocess import run
from tempfile import TemporaryDirectory
from sympy import Matrix
P=Path(__file__).resolve().parent
B=Matrix([list(map(int,r.split()))for r in (P/'odd_basis.txt').read_text().splitlines()])
assert B.shape==(24,24) and abs(B.det())==2**36
cols=[7,1,4,4,2,7,3,3,4,7,2,1]
# Compare the supplied lattice basis to its exact generating stock.
S=[[4*((j==0)+(j==i))for j in range(24)]for i in range(1,24)]+[[8]+[0]*23]
g=[[int((cols[j-i]&1) if 0<=j-i<12 else 0) if j!=23 else 1 for j in range(24)]for i in range(12)]
S += [[2*x for x in a]for a in g]+[[1]*24]
assert all(x.q==1 for x in Matrix(S)*B.inv())
with TemporaryDirectory() as t:
 u=Path(t)/'v'
 run(['g++','-O3','-std=c++17',str(P/'replay.cpp'),'-o',str(u)],check=True)
 r=run([str(u)],cwd=P,text=True,capture_output=True,check=True);print(r.stdout)
 assert '24 2048' in r.stdout.splitlines() and '80 1318912' in r.stdout.splitlines()
 assert max(int(l.split()[0]) for l in r.stdout.splitlines() if l and l[0].isdigit())==80
 assert 'POSITIVE DEFINITE 72 49200 0' in r.stdout
 assert f'POSITIVE DEFINITE 13210188 {(2**24)*4677} 0' in r.stdout
N=2**24;C1=N//2;C2=Q(72*N,49200);C3=Q(4403396,1559)
assert C3 < Q(27,20)**12*Q(6,5)**24
assert C2*C2 < N*Q(6,5)**24
assert C2*C1 < N*Q(27,20)**12*Q(60,47)**24
assert 2+Q(3,2)*Q(27,20)<Q(9,2)
print('PASS exact period basis, Fourier coefficients and density comparisons')
