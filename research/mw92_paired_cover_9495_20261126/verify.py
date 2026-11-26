#!/usr/bin/env python3
"""Exact check, using MW92 (Schuett--Shioda) and equal-ball Coxeter--Few--Rogers.
Dirichlet coefficients follow research/mw92_cfr_taylor13_20261123/moments.py.
"""
from fractions import Fraction as F
from math import factorial, prod, comb, isqrt
import json
assert __debug__
def moments_cover(n=92,k=13):
    B=[factorial(2*j)//factorial(j) for j in range(k+1)]
    C=[1]+[0]*k
    for _ in range(n):
        C=[sum(C[j-i]*B[i] for i in range(j+1)) for j in range(k+1)]
    raw=[F(factorial(j)*C[j],prod(range(n,n+2*j))) for j in range(k+1)]
    m=raw[1]
    central=[sum(comb(j,i)*(-m)**(j-i)*raw[i] for i in range(j+1)) for j in range(k+1)]
    a=n//2
    return (n+1)*F(n,n+1)**a*sum((-1)**j*comb(a+j-1,j)*central[j]/(1+m)**(a+j) for j in range(k+1))
def atan(z,m):
    return sum(F((-1)**k,(2*k+1)*z**(2*k+1)) for k in range(m))
pi_lo=16*atan(5,40)-4*atan(239,41)
pi_hi=16*atan(5,41)-4*atan(239,40)
parent=47**7
A=F('11.00098526'); D=F('10.78402070')
s=16-A; t=16-3*s
assert A-t/3==F(32,3) and 4*s>16 and t>0
assert F(5,2)*D**46<A**46
assert pi_hi**46*A**46<factorial(46)*parent*moments_cover()
p=F(14,3); W=D-2*p; q=16-D
# 2+sqrt(W)<sqrt(D) iff D-W-4>0 and (D-W-4)^2>16W
assert W>0 and D-W-4>0 and (D-W-4)**2>16*W
assert 4*(q-p/4)>16 and 4*p>=16 and p+F(34,3)==16
assert q+D==16 and q+2*p+W==16 and q+6*p>16 and 1%3==(-2)%3
generic=3*s*s/4; fixed=p*(q-p/4); old=F(859313833,45000000); ext=F(192)
assert generic==F(187426112780529507,10**16)
assert 0<generic<fixed<old<ext
center2=F(4**94,parent**2)/fixed
occ2=center2*pi_lo**94/factorial(47)**2
def floor_sqrt_scaled(z,scale):
    t=z*scale*scale
    return isqrt(t.numerator//t.denominator)
oi=floor_sqrt_scaled(occ2,10**36)
occ=f'{oi//10**15}.{oi%10**15:015d}e-21'
gain=floor_sqrt_scaled(old/fixed,10**8)-10**8
extgain=floor_sqrt_scaled(ext/fixed,10**8)-10**8
assert F(10052535,10000000)**2<old/fixed
# Independent exact extension to dimension 95. All tail phases with last
# coefficient 1 are quotiented by the first-coordinate period three.
from itertools import product, permutations
def det(M):
    return sum((-1)**sum(o[i]>o[j] for i in range(len(o)) for j in range(i+1,len(o)))*
               prod(M[i][o[i]] for i in range(len(o))) for o in permutations(range(len(M))))
def short_positive(G):
    n=len(G); den=det(G)
    assert all(det([row[:j] for row in G[:j]])>0 for j in range(1,n+1))
    cof=[det([[e for l,e in enumerate(row) if l!=i] for j,row in enumerate(G) if j!=i]) for i in range(n)]
    bounds=[isqrt(int(16*c/den)) for c in cof]
    for zz in product(*(range(-h,h+1) for h in bounds)):
        if next((o for o in reversed(zz) if o),0)<=0: continue
        val=sum(z1*z2*G[i][j] for i,z1 in enumerate(zz) for j,z2 in enumerate(zz))
        if val<16: yield zz,val
d8=F('10.51475582'); d75=F('10.52951849'); z=16-d75
G=[[p,p/2,p/2],[p/2,q,q/2],[p/2,q/2,z]]
demands={}
for (i,j,k),v in short_positive(G):
    if not k:
        if not j: assert abs(i)==1 and v+F(34,3)>=16
        elif j==1:
            assert v+(W if i%3==1 else D)>=16
        else: assert False
    else:
        assert k==1
        phase=(i%3,j)
        demands[phase]=max(demands.get(phase,0),16-v)
tiny={h for h,f in demands.items() if f<=F(49,400)}  # radius <= 7/20
mid={h for h,f in demands.items() if F(49,400)<f<=F(121,100)}
large=set(demands)-tiny-mid
assert tiny=={(1,1),(0,-2),(1,-2),(0,1)}
assert mid=={(1,0)}
assert len(large)==6 and all(demands[h]<=d75 for h in large)
assert F(31,10)**2<d75  # dimer: 2+1.1<sqrt(d75)
assert F(17,6)**2>8 and (F(17,6)+F(7,20))**2<d75 # square circumradius 2*sqrt(2)
K=len(large)+F(len(mid),2)+F(len(tiny),4)
assert K==F(15,2) and K*d75**46<A**46
d_fixed=det(G)
# In the generic branch one six-centre covering plus six weak phases,
# each covered using one large centroid per three MW sites, replaces all
# former six-, seven-, nine-residue subcases.
Gg=[[F(5),F(5,2),-F(5,2)],[F(5,2),F(5),-F(5,2)],[-F(5,2),-F(5,2),16-d8]]
near={(0,0),(1,0),(0,1),(-1,1),(1,-1),(1,1)}
far={(0,-1),(-1,0),(-1,2),(0,2),(2,-1),(2,0)}
seen=set()
for (i,j,k),val in short_positive(Gg):
    if k==0:
        N=i*i+i*j+j*j
        assert N in {1,3} and val+(A if N==1 else t)>=16
    else:
        assert k==1 and (i,j) in near|far
        seen.add((i,j))
        assert val+(d8 if (i,j) in near else d8-10)>=16
assert seen==near|far and F(7,3)**2>F(16,3) and (F(7,3)+F(3,4))**2<d8
assert 0<d8-10<F(3,4)**2
assert (len(near)+F(len(far),3))*d8**46<A**46
d_generic=det(Gg); worst95=max(d_fixed,d_generic); old95=F(652349944595359843,9000000000000000)
assert 0<d_generic<d_fixed==worst95<old95
# A3 tail of minimum squared norm 16: scale root Gram by 8.
ext95=4*8**3 
def gain_pct(old,cur):
    i=floor_sqrt_scaled(old/cur,10**8)-10**8
    return f'{i//10**6}.{i%10**6:06d}'
def occ_floor(t):
    i=floor_sqrt_scaled(t,10**36)
    return f'{i//10**15}.{i%10**15:015d}e-21'
center95=F(4**95,parent**2)/worst95
ball95=2**48*pi_lo**47/prod(range(1,96,2))
print(json.dumps({'94':{'previous_det':str(old),'generic_det':str(generic),'fixed_det':str(fixed),
                         'coverage_centers_per_L_cell':'5/2','squared_full_radius':str(D),'squared_weak_radius':str(W),
                         'center_density_squared':str(center2),'occupied_density_floor':occ,
                         'gain_percent_floor':gain_pct(old,fixed),'derived_external_gain_percent_floor':gain_pct(ext,fixed)},
                  '95':{'previous_det':str(old95),'generic_det':str(d_generic),'fixed_det':str(d_fixed),
                         'worst_det':str(worst95),'generic_coverage_centers_per_L_cell':'8','fixed_coverage_centers_per_L_cell':str(K),
                         'fixed_large_dimer_square_counts':[len(large),len(mid),len(tiny)],
                         'center_density_squared':str(center95),'occupied_density_floor':occ_floor(center95*ball95**2),
                         'gain_percent_floor':gain_pct(old95,worst95),'derived_external_gain_percent_floor':gain_pct(ext95,worst95)}},
                 indent=2))

