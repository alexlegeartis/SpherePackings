"""Exact rational polytope checks for the 9D E7/Bravais theorem."""
from fractions import Fraction as F
from itertools import combinations


def dot(x, y):
    return sum(a*b for a,b in zip(x,y))


def norm(x):
    return dot(x,x)


def det(A):
    if not A:
        return 1
    return sum((-1)**j*a*det([r[:j]+r[j+1:] for r in A[1:]])
               for j,a in enumerate(A[0]))


def verts(rows, rhs):
    """An extreme point with k nonzero gaps has k independent active rows."""
    seen = {(F(0),)*7}
    for k in range(1,len(rows)+1):
        for J in combinations(range(7),k):
            for I in combinations(range(len(rows)),k):
                A = [[rows[i][j] for j in J] for i in I]
                d = det(A)
                if not d:
                    continue
                g = [F(0)]*7
                for l,j in enumerate(J):
                    g[j] = F(det([r[:l]+[rhs[i]]+r[l+1:] for r,i in zip(A,I)]),d)
                if min(g)>=0 and all(dot(a,g)<=b for a,b in zip(rows,rhs)):
                    seen.add(tuple(g))
    return [tuple(sum(g[j]*(F(int(i<=j))-F(j+1,8)) for j in range(7))
                  for i in range(8)) for g in sorted(seen)]


def mix(x,y,t):
    return tuple(a+t*(b-a) for a,b in zip(x,y))


def p(x):
    return x[0]+x[1]


def cuts(V, levels):
    yield from V
    for x,y in combinations(V,2):
        a,b=p(x),p(y)
        for c in levels:
            if (c-a)*(c-b)<0:
                yield mix(x,y,(c-a)/(b-a))


def main():
    rows = [[2]*7,[1,2,3,4,3,2,1]]
    Q = verts(rows,[2,2])
    P = verts(rows+[[-1,-2,-1,0,1,2,1]],[2,2,0])
    h=(F(3,4),F(3,4))+(-F(1,4),)*6
    assert len(P)==21 and norm(h)==F(3,2)
    assert h in P and all(norm(y)<=1 for y in P if y!=h)
    crossed=0
    for y in P:
        if y==h:
            continue
        delta=tuple(b-a for a,b in zip(h,y))
        direction=sorted(delta,reverse=True)
        # Rearrangement across coordinates and the sorted-Q vertex maximum.
        gauge=max(dot(direction,z) for z in Q)
        if gauge>=F(1,2):
            crossed+=1
            assert norm(mix(h,y,F(1,2)/gauge))<=1
    for x in cuts(P,[F(3,4),F(1)]):
        if p(x)<=F(3,4):
            assert norm(x)<=F(7,8)
        if F(3,4)<=p(x)<=1:
            assert norm(x)+2*(1-p(x))**2<=1
    assert F(9,100)*(F(1,2)+F(159,100)**2)<F(3,8)
    assert F(1,4)*F(7,4)**2<1
    print('PASS:',len(Q),'sorted root-polytope vertices,',len(P),
          'half-polytope vertices,',crossed,'polar crossings and two clearance bounds')


if __name__=='__main__':
    main()

# The reduced-radius endpoint in RADIUS_EXTENSION.md.
from fractions import Fraction
assert Fraction(23,2)**2 * Fraction(29,50)**9 < 1
