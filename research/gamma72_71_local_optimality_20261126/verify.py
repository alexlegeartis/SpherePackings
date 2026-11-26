"""Rational replay for 71D periodic local rigidity and the rank-two sheet comparison."""
from fractions import Fraction as F
from math import prod
def interp(target, zeros):
    c=[F(1)]
    for j in zeros:
        y=[F(0)]*(len(c)+1)
        for i,b in enumerate(c):
            y[i]-=b*j*j; y[i+1]+=b
        c=y
    d=prod(target**2-j*j for j in zeros)
    return [b/d for b in c]
def value(p,X):return sum(c*X**i for i,c in enumerate(p))
def ratio(p,d):
    return sum(c*F(80**k*prod(range(1,2*k,2)),prod(72+d+2*j for j in range(k))) for k,c in enumerate(p))
p5=interp(5,[0,1,2,3,4])
upper5=interp(5,[0,1,2,3])
lo2=interp(1,[0,2,3,4])
lo4=interp(1,[0,2,5])
for j in range(6):
    assert value(p5,j*j)==(j==5)
    assert value(upper5,j*j)>=(j==5)
    assert value(lo2,j*j)<=(j==1) and value(lo4,j*j)<=(j==1)
assert value(lo2,25)==-210
p=ratio(p5,0);U=F(8,72)*ratio(upper5,2);T=F(11,2)*p
assert p==F(142,77727195) and U==F(8,72)*F(83,2878785)
assert 3*U<T<4*U
A=F(571,500000)
assert 210**2*(3*U**2+(T-3*U)**2-T**2/71)<A**2
B=F(128,72*74)*ratio(lo4,4);C=F(8,72)*ratio(lo2,2)
assert B==F(128,72*74)*F(850,2223) and C==F(8,72)*F(39952,82251)
mu=F(79,10);S=F(1,10000);eta=F(1,100000)
assert C>S and B>A*mu+S
assert 4*(B-A*mu-S)*(C-S)>A**2*mu
assert mu*mu+mu<81 and (9+eta)*eta<S
# At the endpoints of the log-density envelopes and at their transition;
# squares cancel the positive common factors. D(e)^2 = f(e)^2*((8-e)/4)^71/48.
def envelope_sq(e,freq_sq):return freq_sq*((8-e)/4)**71/48
opt=envelope_sq(F(1,8),F(24))
assert opt==4*((8-F(1,8))/4)**71/8
assert opt>envelope_sq(F(0),F(6))
e=F(5,12)
assert opt>envelope_sq(e,F(9)/(F(1,2)-e))
assert opt>envelope_sq(e,F(48**2)/(8-e))
inc=4*((8-F(1,10))/4)**71/10
assert opt<inc
print("scope: 71D primitive two-translate deformations, not supercells")
print("moment majorant",A,"quadratic and linear witnesses",B,C)
print("positive rational discriminant margin",4*(B-A*mu-S)*(C-S)-A**2*mu)
print("rank-two sheet optimum / unrestricted incumbent",float(opt/inc)**.5)
