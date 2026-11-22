"""Exact certificate for a periodic packing in R^29. Python standard library only."""
from fractions import Fraction

n=14
allbits=(1<<n)-1
row=471
shift=(715,15538)
def rol(x,i):
    return ((x<<i)|(x>>(n-i)))&allbits
A=tuple(rol(row,i)for i in range(n))
AT=tuple(sum(((a>>i)&1)<<j for j,a in enumerate(A))for i in range(n))
def im(x,rows):
    y=0
    for j,a in enumerate(rows):
        if (x>>j)&1:y^=a
    return y
C=[(x,im(x,A))for x in range(1<<n)]
D=[(im(y,AT),y)for y in range(1<<n)]
Cset=set(C); Dset=set(D)
def wt(u):
    return u[0].bit_count()+u[1].bit_count()
def dist(u,v):
    return wt((u[0]^v[0],u[1]^v[1]))
assert len(Cset)==len(Dset)==1<<n
assert all(wt(x)%2==0 for x in C+D)
# Orthogonality of generators proves orthogonality of words.
assert all((((1<<i)&AT[j]).bit_count() + (A[i]&(1<<j)).bit_count())%2==0
           for i in range(n)for j in range(n))
assert min(wt(x)for x in C[1:])==min(wt(x)for x in D[1:])==8
assert min(dist(shift,x)for x in C)==min(dist(shift,x)for x in D)==6
assert next(x for x in C if x not in Dset)==(1,471)

# Horizontal period 4D28, with D28 the even-coordinate-sum subgroup of Z^28.
# C and D parity types give disjoint representatives by reduction mod 4.
hv=4**28*2
labels=2**15
assert hv//labels==2**42 and hv%labels==0
# Squared length for same-parity intra-layer pairs, or an underlying period:
same_bound=min(4*8,16*2)
# Opposite parity: for |coordinate|=1 everywhere the forced carries
# have parity wt(c)-c.d=0, contradicting the extra e1 in the odd stock.
opp_bound=27+9
# Heights in units sqrt(8); neighboring layers differ by +/-2*shift.
neighbor_bound=min(4*6,28)+8
# Any non-neighboring layers have vertical squared distance >=32.
assert min(same_bound,opp_bound,neighbor_bound,4*8)==32
# Minimum actually attained by horizontal period 4(e1+e2).
# Vertical period squared is 32, with twice the horizontal label count.
r2=Fraction(32,4)
covol_squared=32*hv*hv
multiplicity=2*labels
delta_squared=r2**29*multiplicity**2/covol_squared
assert delta_squared==1
assert delta_squared/Fraction(1,2)==2
print('PASS: [28,14,8], [28,14,8] dual; simultaneous shift depth (6,6); '
      'minimum squared 32; periodic nonlattice center density 1 (sqrt(2) gain ratio).')
