"""Integer-interval replay of explicit 29D and 30D Bacher lattice lifts."""
from fractions import Fraction as F
from hashlib import sha256
from math import factorial,lcm,prod
from pathlib import Path
import json,os,subprocess,tempfile

H=Path(__file__).resolve().parent
G_SHA='5f330e140ab6e336628a363d1ee25df8b5cf3912d3ed142484d69db57e15d864'

def ldl(A):
    n=len(A); L=[[F(i==j) for j in range(n)] for i in range(n)]; D=[F(0)]*n
    assert all(len(r)==n and all(r[j]==A[j][i] for j in range(n)) for i,r in enumerate(A))
    for i in range(n):
        D[i]=F(A[i][i])-sum(L[i][j]**2*D[j] for j in range(i))
        assert D[i]>0
        for j in range(i+1,n):
            L[j][i]=(A[j][i]-sum(L[j][k]*L[i][k]*D[k] for k in range(i)))/D[i]
    assert all(sum(L[i][k]*L[j][k]*D[k] for k in range(n))==A[i][j] for i in range(n) for j in range(i+1))
    return L,D

def input_rows(G,L,D):
    A=[]; W=[]
    for i in range(28):
        k=lcm(*(L[j][i].denominator for j in range(28)))
        V=[int(L[j][i]*k) for j in range(28)]
        assert all(V[j]==L[j][i]*k for j in range(28))
        A.append(V)
        W.append(D[i]/k**2)
    assert all(0<w.numerator<1<<55 and 0<w.denominator<1<<55 for w in W)
    assert max(abs(x) for r in A for x in r)<1<<32
    assert all(sum(W[k]*A[k][i]*A[k][j] for k in range(28))==G[i][j] for i in range(28) for j in range(i+1))
    return "\n".join(" ".join(map(str,r)) for r in [[t for w in W for t in (w.numerator,w.denominator)],*A,*G])+"\n"

def arctan_bounds(q,n=45):
    s=sum((-1)**i*F(1,(2*i+1)*q**(2*i+1)) for i in range(n))
    t=s+(-1)**n*F(1,(2*n+1)*q**(2*n+1))
    return min(s,t),max(s,t)

def pi_interval():
    a,b=arctan_bounds(5); c,d=arctan_bounds(239)
    return 16*a-4*d,16*b-4*c

def ball_squared(n,pi):
    if n%2==0:return pi**n/factorial(n//2)**2
    return 2**(n+1)*pi**(n-1)/prod(range(1,n+1,2))**2

def main():
    assert __debug__
    b=(H/'gram.txt').read_bytes()
    assert sha256(b).hexdigest()==G_SHA
    G=[list(map(int,r.split())) for r in b.decode().splitlines()]
    L,D=ldl(G);det=prod(D)
    assert det==603979776 and all(G[i][i]==8 for i in range(28))
    mask=json.loads((H/'witness.json').read_text())
    x,y=mask['half_labels']; d=mask['direction']
    eps=F(mask['shift_parameter']); h=F(mask['height_diagonal']); off=F(mask['height_off_diagonal_30'])
    assert eps==F(1,12) and h==F(11,6) and off==F(-2,3)
    assert len(d)==28 and all(type(v)==int for v in d)
    gd=[sum(G[i][j]*d[j] for j in range(28)) for i in range(28)]
    assert gd==[0,0,0,0,-6]+[0]*23 and sum(v*w for v,w in zip(d,gd))==24
    assert 0<x<1<<28 and 0<y<1<<28 and x!=y
    X=[F(x>>i&1,2)+eps*d[i] for i in range(28)]
    Y=[F(y>>i&1,2)+eps*d[i] for i in range(28)]
    q=[(0,1,7,*([0]*28))]; info=[]
    for k in (29,30):
        mineig=h-(abs(off) if k==30 else 0)
        assert 9*mineig>8
        for a in range(3):
            for b in range(-2,3) if k==30 else (0,):
                if not a and b<=0:continue  # arithmetic sign symmetry
                v=h*(a*a+b*b)+2*off*a*b
                if v>=8:continue
                t=[a*u+b*w for u,w in zip(X,Y)]
                den=lcm(*(u.denominator for u in t))
                target=(8-v)*den**2
                cutoff=(target.numerator+target.denominator-1)//target.denominator-1
                id=len(q);q.append((id,den,cutoff,*(int(u*den)%den for u in t)))
                info.append((k,a,b,den,cutoff))
    assert len(info)==10  # pinned complete short-altitude roster, plus the body control
    payload="\n".join(" ".join(map(str,row)) for row in q)+"\n"
    with tempfile.TemporaryDirectory()as t:
        p=Path(t); ex=p/'enum'
        (p/'ldl_fast.txt').write_text(input_rows(G,L,D))
        subprocess.run([os.environ.get('CXX','g++'),'-std=c++17','-O3',str(H/'enumerate.cpp'),'-o',str(ex)],check=True)
        r=subprocess.run([str(ex)],input=payload,text=True,stdout=subprocess.PIPE,check=True,cwd=t)
        rows=[l.split()for l in r.stdout.splitlines()]
        assert [int(l[0])for l in rows]==list(range(len(q)))
        assert all(int(l[4])==0 for l in rows)
        print("short layer checks",info)
        print(r.stdout)
    determinants={29:det*h,30:det*(h*h-off*off)}
    pi_lo,pi_hi=pi_interval()
    results=json.loads((H/'results.json').read_text())['rows']
    assert [r['dimension'] for r in results]==[29,30]
    references={29:F(1,3),30:F(3**27,2**44)}
    unrestricted={29:F(1,2),30:F(1)}
    for r in results:
        n=r['dimension']; V=determinants[n]; sq=F(2)**n/V
        assert r['minimum_squared']==8 and F(r['gram_determinant'])==V and F(r['center_density_squared'])==sq
        gain=F(r['gain_percent_floor']); o=F(r['occupied_density_floor'])
        assert F(r['previous_lattice_center_density_squared'])==references[n]
        assert sq>references[n] and sq<unrestricted[n]
        assert (1+gain/100)**2<=sq/references[n]
        assert o**2<=sq*ball_squared(n,pi_lo)
        print(json.dumps(r))
if __name__=='__main__':main()
