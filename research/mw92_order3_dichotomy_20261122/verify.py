#!/usr/bin/env python3
"""Exact replay for the MW92 order-three 94--96D lattice existence bounds."""
import argparse
import json
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from itertools import permutations, product
from math import factorial, isqrt, prod
from pathlib import Path

PARENT = 47**7
A = Q("11.000980284")
DEPTH = {1: A, 3: Q("10.741"), 8: Q("10.51475106"),
         12: Q("10.422476"), 17: Q("10.3438567")}
A_DIAG = (16-A, 16-A, 16-DEPTH[12])
B_DIAG = (Q(14,3), 16-DEPTH[3], 16-DEPTH[8], 16-DEPTH[17])
HERE = Path(__file__).resolve().parent
BASE = "f35cf8724f816c5299b56911b8155ee2a7b89c94"

def det(G):
    n = len(G)
    return sum((-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
               * prod(G[i][p[i]] for i in range(n)) for p in permutations(range(n)))

def gram(diag):
    return [[diag[min(i,j)] / (1 if i==j else 2) for j in range(len(diag))]
            for i in range(len(diag))]

def short(G):
    n = len(G)
    assert all(det([r[:i] for r in G[:i]]) > 0 for i in range(1,n+1))
    D = det(G)
    C = [det([[x for j,x in enumerate(r) if i != j]
               for k,r in enumerate(G) if k != i]) for i in range(n)]
    B = [isqrt(int(16*c/D)) for c in C] # |z_i|^2 <= q(z)*(G^-1)_ii.
    for z in product(*(range(-b,b+1) for b in B)):
        j = next((i for i in reversed(range(n)) if z[i]),None)
        if j is None or z[j]<0: continue
        q = sum(G[i][k]*z[i]*z[k] for i in range(n) for k in range(n))
        if q<16: yield z,j,q

def atan(x,k):
    return sum(Q((-1)**i,(2*i+1)*x**(2*i+1)) for i in range(k))

def ball2(n,p):
    if n%2==0: return p**n/factorial(n//2)**2
    return 2**(n+1)*p**(n-1)/prod(range(1,n+1,2))**2

def floor_sqrt(sq,digits=15):
    with localcontext() as c:
        c.prec=60
        x=(Decimal(sq.numerator)/Decimal(sq.denominator)).sqrt()
        u=Decimal(1).scaleb(x.adjusted()+1-digits)
    t=sq/Q(u)**2
    k=isqrt(t.numerator//t.denominator)
    assert (Q(u)*k)**2<=sq<(Q(u)*(k+1))**2
    return f"{u*k:.{digits-1}E}".lower()

def percent(now,old):
    x=now/old*10**12
    k=isqrt(x.numerator//x.denominator)-10**6
    return ('-' if k<0 else '')+f"{abs(k)//10000}.{abs(k)%10000:04d}"

def covering():
    n=92
    mu=Q(2,n+1)
    second=Q(24*n+4*n*(n-1),prod(range(n,n+4)))
    third=Q(720*n+144*n*(n-1)+8*n*(n-1)*(n-2),prod(range(n,n+6)))
    var=second-mu**2
    skew=third-3*mu*second+2*mu**3
    h=n//2
    return (n+1)*Q(n,n+3)**h*(1+h*(h+1)*var/(2*(1+mu)**2)
                             -h*(h+1)*(h+2)*skew/(6*(1+mu)**3))

def eisenstein_gram():
    # O=Z[J] with coefficient norm e0^2-e0 e1+e1^2.
    p,b,c=Q(5),Q("5.95"),Q("4.999")
    G=[[p,-p/2,c/2,0],[-p/2,p,-c/2,c/2],
       [c/2,-c/2,b,-b/2],[0,c/2,-b/2,b]]
    J=[[0,-1,0,0],[1,-1,0,0],[0,0,0,-1],[0,0,1,-1]]
    assert all(sum(J[k][i]*G[k][l]*J[l][j] for k in range(4) for l in range(4))==G[i][j]
               for i in range(4) for j in range(4))
    assert det(G)==Q(9,16)*(p*b-c*c/3)**2
    return G

def verify_eisenstein(G,t):
    def j(z):
        return -z[1],z[0]-z[1],-z[3],z[2]-z[3]
    bank={}
    for z,_,q in short(G):
        d=16-q
        if z[2:]==(0,0):
            assert (q==5 and d<=A) or (q==15 and d<=t)
        else:
            o=(z,j(z),j(j(z)))
            rep=min(o+tuple(tuple(-i for i in v) for v in o))
            assert rep not in bank or bank[rep]==d
            bank[rep]=d
    weight=sum(x**46 for x in bank.values())
    assert PI_HIGH**46*weight < factorial(46)*PARENT
    return {"unit_orbits_with_nonzero_second_coefficient":len(bank),
            "union_measure_upper_ceiling_to_12_places":
                str(Decimal(1+10**12 * PI_HIGH**46*weight // (factorial(46)*PARENT))/10**12)}

def row(n,prior):
    r=n-92
    Gs=[gram(A_DIAG[:r]) if n<96 else eisenstein_gram(),gram(B_DIAG[:r])]
    qa=A_DIAG[0]
    # If d(Ta)^2 reaches this threshold, the A2 shell 3*qa works.
    t=16-3*qa
    assert t>0 and A-t/3==Q(32,3)
    assert Q(34,3)+B_DIAG[0]==16
    sa={j:set() for j in range(2,r)}
    sb={j:set() for j in range(1,r)}
    caseA={}
    if n==96:
        caseA=verify_eisenstein(Gs[0],t)
    else:
        for z,j,q in short(Gs[0]):
            if j<=1:
                assert (z[:2] in ((1,0),(-1,1),(0,1)) and q+A==16 or
                        z[:2] in ((-2,1),(1,1),(-1,2)) and q+t==16)
            else:
                assert z[j]==1 and q+DEPTH[12]>=16
                sa[j].add(z[:j])
    for z,j,q in short(Gs[1]):
        if j==0:
            assert z[0]==1 and q+Q(34,3)==16
        else:
            d=DEPTH[(3,8,17)[j-1]]
            assert z[j]==1 and q+d>=16
            sb[j].add((z[0]%3,)+z[1:j] if j in (1,3) else z[:j])
    if n<96: assert [len(s) for s in sa.values()]==[12][:r-2]
    assert [len(s) for s in sb.values()]==[3,8,17][:r-1]

    D=max(map(det,Gs))
    csq=Q(4**n)/(PARENT**2*D)
    old=Q(prior[n]["center_density_squared"])
    unres=Q(prior[n]["repository_unrestricted_comparison_center_density_squared"])
    ext=Q(prior[n]["external_lattice_comparison_center_density_squared"])
    assert ext== (Q(2**184,PARENT**2*2**r*(r+1)) if n<96 else Q(3**96,2**48))
    assert csq>old and ((csq>unres)==(n<96))
    return {
        "dimension":n, "kind":"lattice existence",
        "minimum_squared":16,
        "branch_grams":[[[str(x) for x in row] for row in G] for G in Gs],
        "branch_vertical_determinants":list(map(lambda G:str(det(G)),Gs)),
        "worst_vertical_determinant":str(D),
        "max_gram_determinant":str(PARENT**2*D),
        "center_density_squared":str(csq),
        "occupied_density_floor":floor_sqrt(csq*ball2(n,PI_LOW)),
        "gain_over_repository_lattice_percent_floor":percent(csq,old),
        "gain_over_repository_unrestricted_percent_floor":percent(csq,unres),
        "gain_over_named_external_lattice_percent_floor":percent(csq,ext),
        "unrestricted_improvement": n<96,
        **caseA,
    }

PI_LOW=16*atan(5,60)-4*atan(239,61)
PI_HIGH=16*atan(5,61)-4*atan(239,60)
def certify(previous):
    assert Q(3141592653589793,10**15)<PI_LOW<PI_HIGH<Q(22,7)
    for k,d in DEPTH.items():
        assert k*PI_HIGH**46*d**46<factorial(46)*PARENT*covering()
    past={x["dimension"]:x for x in previous["rows"] if x["kind"]=="lattice existence"}
    periodic=next(x for x in previous["rows"] if x["dimension"]==96 and x["kind"]=="periodic existence")
    past[96]["repository_unrestricted_comparison_center_density_squared"]=periodic["center_density_squared"]
    return {"comparison_base_commit":BASE,
            "strict_covering_depths":{str(k):str(v) for k,v in DEPTH.items()},
            "rows":[row(n,past) for n in range(94,97)]}

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write",action="store_true")
    ap.add_argument("--previous",type=Path,default=HERE.parent/"mw92_extensions_20261122"/"results.json")
    args=ap.parse_args()
    out=certify(json.loads(args.previous.read_text()))
    path=HERE/"results.json"
    if args.write: path.write_text(json.dumps(out,indent=2)+"\n")
    else: assert out==json.loads(path.read_text())
    for row in out["rows"]:
        print(row["dimension"],row["occupied_density_floor"],
              row["gain_over_repository_lattice_percent_floor"],
              row["gain_over_repository_unrestricted_percent_floor"])
