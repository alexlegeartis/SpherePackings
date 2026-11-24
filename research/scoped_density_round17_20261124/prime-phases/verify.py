from fractions import Fraction as F
from math import comb

def mul(a,b):
    return [sum(a[i]*b[j-i] for i in range(j+1)) for j in range(len(a))]

def power(a,k):
    b=[1]+[0]*(len(a)-1)
    for _ in range(k):
        b=mul(b,a)
    return b

def theta_a2(k):
    return [sum(a*a+a*b+b*b == t for a in range(-k-1,k+2) for b in range(-k-1,k+2)) for t in range(k+1)]

def leading(k,r):
    s=3-r if r else 0
    e=4*k-1 if r else 4*k
    A=theta_a2(k)
    cs=[]
    rem=[1]+[0]*k
    for i in range(k+1):
        col=mul(eta(k,3*r+9*(e-i)-3*i,-r-3*(e-i)+9*i,3),power(A,s))
        cs.append(rem[i])
        for j in range(i,k+1):
            rem[j]-=cs[-1]*col[j-i]
    assert not any(rem)
    count=F(cs[-1],3**s)
    depth=F(2,3)*(e-k+F(r,3))
    return depth,count

def eta(k,a,b,p):
    f=[1]+[0]*k
    for i in range(1,k+1):
        t=a+(b if i%p==0 else 0)
        g=[0]*(k+1)
        for j in range(k//i+1):
            g[i*j]=((-1)**j*comb(t,j) if j<=t else 0) if t>=0 else comb(-t+j-1,j)
        f=mul(f,g)
    return f

def p_iso(k,p,K=None):
    h=24//(p-1)
    if K is None:
        K=k
    rem=[1]+[0]*k
    cs=[]
    for i in range(k+1):
        c=rem[i]
        cs.append(c)
        terms=eta(k,h*(p*K-i),h*(i-K),p)
        for j in range(i,k+1):
            rem[j]-=c*terms[j-i]
    assert not any(rem)
    return cs[-1]


def circle_count_48():
    # Type II rank 48 ordinary theta, q^j counts squared norm 2j.
    N=16
    E=[1]+[240*sum(d**3 for d in range(1,j+1) if j%d==0) for j in range(1,N+1)]
    delta=eta(N,24,0,3)
    coeff=[0]*(N+1)
    for i in range(3):
        col=mul(power(E,3*(2-i)),power(delta,i))
        c=int(i==0)-coeff[i]
        for j in range(i,N+1):
            coeff[j]+=c*col[j-i]
    assert coeff[:4]==[1,0,0,52416000]
    covered=sum(coeff[j] for j in range(1,16,3))
    anisotropic_classes=3**47-3**23
    n_deep=anisotropic_classes-covered
    assert F(n_deep,anisotropic_classes)>F(127,1000)
    return covered, n_deep, anisotropic_classes

if __name__=='__main__':
    expected=((36,2,9),(1836,105,450),(105084,6060,25632),(6334524,366750,1541727))
    for k in range(1,5):
        rows=[leading(k,r) for r in range(3)]
        assert tuple(count for _,count in rows)==expected[k-1]
        assert all(rows[r][0]==F(2*k)-F((6-2*r) if r else 0,9) for r in range(3))
        print('rank, squared depths, number at max',24*k,rows)
    iso={
        2:(48,2400,136704,8222400),
        3:(36,1836,105084,6334524),
        5:(30,1530,87600,5281650),
        7:(28,1428,81760,4929540),
        13:(26,1326,75920,4577430),
    }
    for p,counts in iso.items():
        assert 24%(p-1)==0 and (24//(p-1))%2==0
        assert tuple(p_iso(k,p) for k in range(1,5))==counts
        print('isotropic prime, mean contact count',p,counts)
    for p in (5,7,13):
        orbit={(pow(2,j,p)%p) for j in range(p-1)}
        assert orbit|{(-x)%p for x in orbit}==set(range(1,p))
        b=(p-1)//2
        assert sum(i*i for i in range(1,b+1))%p==0
        for k in range(1,5):
            K=k*b
            c=p_iso(k,p,K)
            assert c==F(p*K,p*K-k)*eta(k,-(24//(p-1))*(p*K-k),0,p)[k]>0
            if k>1:
                assert F(2*((p+1)*K-k),p*b)<F(4*(4*k-1),6)
    print('48D type-1: shorter vectors, guaranteed maximal classes, all classes',circle_count_48())
