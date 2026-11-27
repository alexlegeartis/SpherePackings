#!/usr/bin/env python3
"""Exact, standard-library-only projected-shell and density computations."""
from collections import Counter
from fractions import Fraction as F
from functools import cache
from itertools import combinations, product
from math import prod

N6, N8 = 52416000, 39007332000


def exponents(dim, degree):
    if dim == 1:
        yield (degree,)
    else:
        for j in range(degree+1):
            for u in exponents(dim-1, degree-j):
                yield (j, *u)


def moment_factory(diagonal, off):
    @cache
    def wick(e):
        if not sum(e):
            return 1
        if sum(e) % 2:
            return 0
        i = next(j for j, x in enumerate(e) if x)
        a = list(e)
        a[i] -= 1
        z = 0
        for j in range(len(a)):
            if a[j]:
                x = a[j]
                a[j] -= 1
                z += x * (diagonal if j == i else off) * wick(tuple(a))
                a[j] += 1
        return z

    def moment(e, m, cardinality):
        j = sum(e)//2
        return F(cardinality*m**j*wick(e), prod(48+2*k for k in range(j)))
    return moment


def mon(e, a):
    return prod(x**i for i, x in zip(e, a))


def interpolate(xs, at):
    """Unnormalized Lagrange polynomial and its value at 'at'."""
    p = [1]
    for r in xs:
        if r == at:
            continue
        q = [0]*(len(p)+1)
        for i, c in enumerate(p):
            q[i] -= r*c
            q[i+1] += c
        p = q
    return p, sum(x*at**i for i, x in enumerate(p))


def classes(dim, bound, key, norm):
    C = {}
    for a in product(range(-bound, bound+1), repeat=dim):
        q, k = norm(a), key(a)
        if k not in C or q < C[k][0]:
            C[k] = q, [a]
        elif q == C[k][0]:
            C[k][1].append(a)
    return C


def largest_clique(C, modulus, distance, cutoff):
    """Translation accepts a phase at zero without loss; integer graph."""
    zero = (0,)*len(next(iter(C)))
    def edge(a, b):
        if a == b:
            return False
        k = tuple((x-y) % modulus for x, y in zip(a, b))
        return distance(C[k][0]) >= cutoff
    best = 1
    def visit(n, A):
        nonlocal best
        best = max(best, n)
        while len(A)+n > best:
            a = A.pop()
            visit(n+1, [b for b in A if edge(a, b)])
    visit(1, [a for a in C if edge(a, zero)])
    return best


def square_density(dim, count, mu, det):
    return F(count*count, det)*(F(mu)/4)**dim


def verify46():
    root = [(6,-3),(-6,3),(-3,6),(3,-6),(3,3),(-3,-3)]
    moment = moment_factory(6, -3)
    M = {(i,j): moment((i,j),6,N6)-sum(mon((i,j),r) for r in root)
         for i in range(12) for j in range(12-i)}
    S = list(range(-3,4))
    c = {}
    for a in [i for i in S if i]:
        row = [b for b in S if abs(a+b) <= 3]
        P, dp = interpolate(S,a)
        for b in row:
            Q, dq = interpolate(row,b)
            c[a,b] = sum(x*y*M[i,j] for i,x in enumerate(P)
                         for j,y in enumerate(Q))/(dp*dq)
    Z = [M[0,j]-sum(n*b**j for (a,b),n in c.items()) for j in range(7)]
    for b in S:
        Q,dq = interpolate(S,b)
        c[0,b] = sum(x*y for x,y in zip(Q,Z))/dq
    assert len(c)==37 and all(a.denominator==1 and a>0 for a in c.values())
    counts = {0:(1,12309600),1:(6,5365440),3:(6,922185),4:(6,361584),
              7:(6,16767),9:(2,1656)}
    def key(a):
        x,y=a
        return (2*x+y)%9,(x+2*y)%9
    def n(a):
        x,y=a
        return x*x+x*y+y*y
    C = classes(2,4,key,n)
    assert Counter(x[0] for x in C.values())=={k:v[0] for k,v in counts.items()}
    assert all(c[a]==counts[q][1] for q,A in C.values() for a in A)
    # A coordinate outside the box gives n >= 3*25/4 > 9.
    assert largest_clique(C,9,lambda q:6-F(2*q,9),F(52,9)) == 3
    assert C[3,6][0]==C[6,3][0]==9  # This pair + zero is an order-three subgroup.
    best=square_density(46,3,F(52,9),27)
    assert best==F(13**46,3**93)
    for s,mu in [(1,6),(9,F(16,3)),(27,F(46,9))]:
        assert square_density(46,s,mu,27)<best


def solve_param(rows, n):
    """Solve the first n-1 unknowns retaining the last as a free parameter.
    Rows already contain right-hand sides for any additional parameters.
    """
    A = [list(map(int,row)) for row in rows]
    assert all(all(F(x).denominator==1 for x in row) for row in rows)
    R=n-1
    prev=1
    for j in range(R):
        i=next(i for i in range(j,len(A)) if A[i][j])
        A[i],A[j]=A[j],A[i]
        p=A[j][j]
        for i in range(j+1,len(A)):
            u=A[i][j]
            for k in range(j+1,len(A[i])):
                x=p*A[i][k]-u*A[j][k]
                assert x%prev==0
                A[i][k]=x//prev
            A[i][j]=0
        prev=p
    assert all(not any(row[R:]) for row in A[R:])
    m=len(A[0])-n
    V=[[F(0)]*(m+1) for _ in range(n)]
    V[-1][-1]=F(1)
    for j in range(R-1,-1,-1):
        right=A[j][n:]+[0]
        for k in range(j+1,n):
            right=[u-A[j][k]*v for u,v in zip(right,V[k])]
        V[j]=[F(u,A[j][j]) for u in right]
    return V


def verify45():
    def canon(a):
        return next((x for x in a if x),0)
    def tup4(a):
        return (*a,-sum(a))
    def add(a,b,k=1):
        return tuple(x+k*y for x,y in zip(a,b))
    def pairs(a):
        return [x+y for x,y in combinations(a,2)]
    roots=[tuple(6 if i==j else -2 for j in range(4)) for i in range(4)]
    norm8={add(a,b) for a,b in combinations(roots,2)}
    target=(3,1,-1)
    P=[]
    for a in product(range(-3,4),repeat=3):
        p=tup4(a)
        if canon(p)>=0 and abs(p[-1])<=3 and max(map(abs,pairs(p)))<=4:
            P.append(a)
    assert len(P)==101
    P.remove(target)
    P.append(target)
    ex=[e for d in range(0,11,2) for e in exponents(3,d)]
    moment=moment_factory(6,-2)
    R6=[(p,2) for p in roots]
    A6=[]
    for e in ex:
        A6.append([mon(e,p) for p in P]+[moment(e,6,N6)-sum(n*mon(e,p) for p,n in R6)])
    V6=solve_param(A6,len(P))
    assert all(len(v)==2 for v in V6)
    # Deterministic positivity endpoints: these checks prove the parameter range.
    n6=dict(zip(P,V6))
    assert n6[target]==[0,1] and n6[0,0,3]==[736,-2]
    assert all(v[0]>=0 and v[0]+368*v[1]>=0 for v in V6)

    # Exact norm-8 boundary counts inherited from norm-6 labels.
    D={}
    for u,v in n6.items():
        p=tup4(u)
        D[p]=D[tuple(-x for x in p)]=v
    for p in roots:
        D[p]=D[tuple(-x for x in p)]=[2,0]
    B=[]
    for x in product(range(-4,5),repeat=3):
        p=tup4(x)
        if canon(p)<=0 or abs(p[-1])>4:
            continue
        if p in norm8:
            B.append((p,2,0))
            continue
        if max(map(abs,pairs(p)))>5:
            continue
        b=None
        if 4 in map(abs,p):
            i=next(i for i,z in enumerate(p) if abs(z)==4)
            b=add(p,roots[i],-p[i]//4)
        else:
            for i,j in combinations(range(4),2):
                if abs(p[i]+p[j])==5:
                    b=add(p,add(roots[i],roots[j]),-(p[i]+p[j])//5)
                    break
        if b is not None:
            assert any(b)  # All mapped labels are nonzero, so both use pair sums.
            u,v=D.get(b,(0,0))
            B.append((p,u,v))
    assert len(B)==105
    A8=[]
    for e in ex:
        A8.append([mon(e,p) for p in P]+[
            moment(e,8,N8)-sum(u*mon(e,p) for p,u,v in B),
            -sum(v*mon(e,p) for p,u,v in B)])
    V8=solve_param(A8,len(P))
    n8=dict(zip(P,V8))
    assert n8[target]==[0,0,1] and n8[0,0,3]==[8662080,4,-2]

    # Degree-12 harmonic linkage across shells (monomial x1² x2⁴ x3⁶).
    e=(2,4,6)
    m6=[sum(v[j]*mon(e,p) for p,v in n6.items()) for j in range(2)]
    m6[0]+=sum(n*mon(e,p) for p,n in R6)
    m8=[sum(v[j]*mon(e,p) for p,v in n8.items()) for j in range(3)]
    for p,u,v in B:
        m8[0]+=u*mon(e,p)
        m8[1]+=v*mon(e,p)
    const=m8[0]+72*m6[0]-moment(e,8,N8)-72*moment(e,6,N6)
    tau=m8[1]+72*m6[1]
    sig=m8[2]
    assert sig and -const/sig==1887360 and -tau/sig==-40

    def key(a):
        return tuple((x+sum(a))%8 for x in a)
    def q(a):
        return sum(x*x for x in a)+sum(a)**2
    C=classes(3,5,key,q)
    assert Counter(v[0] for v in C.values())=={
        0:1,2:12,4:6,6:24,8:12,10:24,12:4,14:24,16:3,18:12,20:6}
    # Any label outside the box has q>=36, so all shortest representatives are listed.
    for Q,V in C.values():
        for p in V:
            p=p if canon(p)>=0 else tuple(-x for x in p)
            u,v=n6[p]
            if Q==18:
                assert [u,v]==[736,-2]
            elif Q==20:
                assert [u,v]==[0,1]
            else:
                assert u>0 and u+368*v>0
            if Q in (18,20):
                x,y,z=n8[p]
                assert x+1887360*z>0 and x+368*y+(1887360-40*368)*z>0
    def optimistic(q):
        return F(23,4) if q==18 else F(11,2) if q==20 else 6-F(q,8)
    assert largest_clique(C,8,optimistic,F(23,4))==4
    assert largest_clique(C,8,optimistic,F(11,2))==6
    H=[a for a in C if all(v in (0,4) for v in a)]
    assert len(H)==8 and all(optimistic(C[a][0])<=F(9,2) for a in H if any(a))
    best=square_density(45,4,F(23,4),128)
    assert best==F(23**45,2**183)
    for s,mu in [(1,6),(6,F(11,2)),(16,F(21,4)),(128,F(9,2))]:
        assert square_density(45,s,mu,128)<best
    assert (7380720-7379838)//9==98


def verify47():
    best=square_density(47,2,F(47,8),8)
    assert best==F(47**47,2**236)
    for t in (6,8,10):
        moment=moment_factory(t,0)
        xs=[k*k for k in range(t//2+1)]
        M=[moment((2*j,),6,N6)-(2*6**(2*j) if t==6 else 0)
           for j in range(len(xs))]
        for k in xs:
            a,b=interpolate(xs,k)
            count=sum(c*d for c,d in zip(a,M))/b
            assert count.denominator==1 and count>0
        B=F(0)
        for S in range(1,1<<t,2):
            V=[i for i in range(t) if (S>>i)&1]
            gap=max((min(abs(a-b),t-abs(a-b)) for a,b in combinations(V,2)),default=0)
            B=max(B,square_density(47,len(V),6-F(gap*gap,t),t))
        assert B==square_density(47,2,6-F(1,t),t)
        assert (B==best if t==8 else B<best)


verify46()
verify45()
verify47()
print("PASS: exact projected moments, cross-shell closure, and scoped 45D–47D phase optima")
