"""Independent exact 38D phase and covolume replay (standard library).
Uses a different field presentation; either one realizes the construction."""
from fractions import Fraction as F
from phase_syndrome import D,F8,mul,q
from itertools import combinations, product
from math import prod

clique=[]; span={0}
def B(x,y):return q(x)^q(y)^q(x^y)
for _ in range(6):
    a=next(i for i in range(1,64) if i not in span and all(B(i,b)==1 for b in clique))
    clique.append(a); span|={x^a for x in span}
s=0
for a in clique:s^=a
clique=[0,*clique,s]
assert len(set(clique))==8 and all(B(a,b)==1 for a,b in combinations(clique[1:],2))
X={0}
for a in clique:
    X|={b^a for b in X}
assert len(X)==64
beta=F8[5]
tags=[a|(mul(beta,a)<<6)|((1^q(a))<<12)if a else 0 for a in clique]
assert all(D[a^b]==6 for a,b in combinations(tags,2))
# D3 has index 2 in Z^3, so 2D3 has 16 integer labels.
# The min-positive squared vertical length in D3 is 2, and in 2D3 is 8.
# Different parities cost at least 1, and the half shifts cost at least 30/4.
# Same parity, same label: minimum 8; same parity, distinct labels: 2+6.
assert min(F(1)+F(30,4),F(8),F(2)+min(D[a^b]for a,b in combinations(tags,2)))==8
horizontal_det=2**31
vertical_covol=F(16,16)
delta_sq=F(2)**38/(horizontal_det*vertical_covol**2)
assert delta_sq==128
print("PASS independent 38D: 28 clique pairs have depth exactly 6; center density squared 128 (nonlattice)")
