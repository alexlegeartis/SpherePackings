"""Exact, finite arithmetic replay for the quarter-torus upper bounds."""
from itertools import combinations
from math import ceil, isqrt
from fractions import Fraction as Q
from sympy import Matrix, simplify, rf
from prior_scalar import ring, remaining

# The convention is t = exp(pi i tau) for the half-level block
# and t = exp(2 pi i tau) for the quarter-level block.
N = 8
c, add, mul, power = ring(N)
def scale(z, a):
    return tuple(a*x for x in z)

theta3 = tuple(1 if not k else (2 if isqrt(k)**2 == k else 0) for k in range(N))
theta4 = tuple(1 if not k else (2*(-1)**isqrt(k) if isqrt(k)**2 == k else 0) for k in range(N))
theta2_regular = tuple(2 if any(k==j*(j+1) for j in range(N)) else 0 for k in range(N))
X = (0,) + power(theta2_regular,4)[:-1]
Y = power(theta4,4)
W = add(X,Y)

def min_count(n):
    """Canonical level-one theta series in E4 and Delta."""
    k = n//24
    E4 = tuple(1 if m==0 else 240*sum(i**3 for i in range(1,m+1) if m%i==0) for m in range(N))
    D = (0,1) + (0,)*(N-2)
    for m in range(1,N):
        D=mul(D,power(tuple(1 if j==0 else (-1 if j==m else 0) for j in range(N)),24))
    F=[mul(power(E4,3*(k-i)),power(D,i)) for i in range(k+1)]
    a=Matrix([[f[j] for f in F] for j in range(k+1)]).inv()*Matrix([1]+[0]*k)
    return sum(z*f[k+1] for z,f in zip(a,F))

def dot_upper(n,s,j):
    # Solve the six exact design moments for k=0,...,5 over all 6-site supports.
    # The polytope is bounded by the total shell count, so a maximum is at a vertex.
    m=2*(n//24+1)
    num=min_count(n)
    sites=tuple(range(s//2+1))
    moments=Matrix([num*(m*s)**k*rf(Q(1,2),k)/rf(Q(n,2),k) for k in range(6)])
    bounds=[]
    for support in combinations(sites,6):
        A=Matrix([[x**(2*k) for x in support] for k in range(6)])
        a=A.inv()*moments
        if min(a)>=0:
            bounds.append(a[support.index(j)] if j in support else 0)
    return max(bounds)

def quarter(n,r,H):
    # The source-level form is P Y^i W^j times a homogeneous polynomial of degree d.
    # At M: (X,Y,W) -> (W,-Y,X); the t^alpha factor is omitted.
    k=n//24
    i=n//6
    alpha=Q(4-r,4)
    j=ceil(H/2-alpha)
    d=n//4-1-i-j
    def seq(X,Y,W,P):
        return [mul(P,mul(power(Y,i),mul(power(W,j),mul(power(X,u),power(Y,d-u))))) for u in range(d+1)]
    P=mul(power(theta4,r),power(theta3,4-r))
    fam=seq(X,Y,W,P)
    a,p=Matrix([[f[t] for f in fam] for t in range(k+1)]).gauss_jordan_solve(Matrix([1]+[0]*k))
    PM=mul(power(theta4,r),power(theta2_regular,4-r))
    F=seq(W,scale(Y,-1),X,PM)
    return tuple(simplify(sum(b*f[jj] for b,f in zip(a,F))) for jj in (1,2)),tuple(p)

def half(n,H):
    # At S: (X,Y,W) -> (Y,X,W) after common signs cancel.
    # The t^.5 factor when h is half-integral is omitted.
    odd=int(H.denominator==2)
    e=int(H)
    d=n//4-odd-2*e
    k=n//24
    def seq(X,Y,W,P):
        Z=add(X,scale(Y,2))
        return [mul(P,mul(power(Y,e),mul(power(W,e),mul(power(Z,d-2*u),power(X,2*u))))) for u in range(d//2+1)]
    P=mul(power(theta3,2),power(theta4,2)) if odd else c(1)
    fam=seq(X,Y,W,P)
    a,p=Matrix([[f[2*t] for f in fam] for t in range(k+1)]).gauss_jordan_solve(Matrix([1]+[0]*k))
    PS=mul(power(theta3,2),power(theta2_regular,2)) if odd else c(1)
    F=seq(Y,X,W,PS)
    return tuple(simplify(sum(b*f[jj] for b,f in zip(a,F))) for jj in range(e,e+3)),tuple(p)

assert min_count(48)==52416000 and min_count(72)==6218175600
assert dot_upper(48,10,5)==200
assert dot_upper(72,10,5)==11360
assert dot_upper(72,12,6)==11088

def min_abs_aff(a,b): # intervals in this replay all stay positive, verified below
    assert 0<a<=b
    return a

a,_=quarter(48,1,Q(5,2)); assert abs(a[0])==516
new={(48,1,Q(5,2)):(abs(a[0]),dot_upper(48,10,5)/2)}
(h,j,l),p=half(72,Q(7,2));assert len(p)==2 and simplify(l+4*j+78*h)==6272
for H,A,U in [(Q(5,2),0,dot_upper(72,10,5)/2),(Q(7,2),2,6272)]:
    (a,b),p=quarter(72,1,H);assert len(p)==1 and simplify(b+66*a)==29648
    new[(72,1,H)]=(min_abs_aff(29648-66*A,29648+66*A),U)
(a,b),p=quarter(72,1,Q(9,2));assert not p and b==29648
(h,j,_),p=half(72,Q(9,2));assert len(p)==1 and simplify(j+4*h)==6272
new[(72,1,Q(9,2))]=(b,6272)
(a,b),p=quarter(72,2,Q(3));assert len(p)==1 and simplify(b+68*a)==17376
new[(72,2,Q(3))]=(min_abs_aff(17376-136,17376+136),dot_upper(72,12,6)/2)
(a,b),p=quarter(72,2,Q(4));assert not p and b==17376
(h,j,l),p=half(72,Q(4));assert len(p)==2 and simplify(l+24*j+96*h)==136704
new[(72,2,Q(4))]=(b,Q(136704,24))
(a,b),p=quarter(72,2,Q(5));assert not p and b==17376
(h,j,_),p=half(72,Q(5));assert len(p)==1 and simplify(j+24*h)==136704
new[(72,2,Q(5))]=(b,Q(136704,24))

assert set(new)=={(n,r,H) for n,T in remaining.items() for r,H,_ in T}
for key,(need,available) in new.items():
    assert need>available
    print('SHELL GATE',key,'needs',need,'available at most',available)
print('PASS universal quarter-torus squared depths <=4 in rank48 and <=6 in rank72')
