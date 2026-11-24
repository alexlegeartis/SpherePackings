from math import factorial
from fractions import Fraction as F

def rank(mat,p):
    r=0
    for j in range(len(mat[0])):
        k=next((k for k in range(r,len(mat)) if mat[k][j]%p),None)
        if k is None:continue
        mat[r],mat[k]=mat[k],mat[r]
        u=pow(mat[r][j],-1,p)
        mat[r]=[(u*x)%p for x in mat[r]]
        for k in range(r+1,len(mat)):
            u=mat[k][j]
            mat[k]=[(x-u*y)%p for x,y in zip(mat[k],mat[r])]
        r+=1
    return r
def tabulate(p,k):
    q=p*p; ranks={}
    for j in range(1,k):
        jp=k-j
        A=(p*jp-j)//k
        assert A*k==p*jp-j and 0<=A<p
        h=(j*q-1)//k;g=(jp*q-1)//k
        vals=[]
        for b in (0,1):
            C=[[0]*h for _ in range(g)]
            for i in range(1,h+1):
                l=-i%p
                if l>A:continue
                for u in range(A-l+1):
                    v=A-l-u
                    n=(i+l)//p+p*u+v
                    w=factorial(A)//(factorial(u)*factorial(v)*factorial(l))
                    assert 1<=n<=g
                    C[n-1][i-1]=w*b**v%p
            vals.append(rank(C,p))
        ranks[j]=tuple(vals)
    return ranks
want={(5,2):(6,), (5,3):(6,6), (5,6):(4,6,6,6,4), (7,2):(12,), (7,4):(10,12,10)}
lower={}
for (p,k),R in want.items():
    D=tabulate(p,k)
    assert tuple(D.values())==tuple((r,r) for r in R)
    # F mod p acts only on the elliptic anti-holomorphic line.
    rF=p*p-1-((k-1)*p*p-1)//k
    # F² adds a constraint on the other elliptic line.
    rFF=D[k-1][0]
    lower[p,k]=rF+rFF
    print(p,k,'Cartier ranks',R,'lower sigma',lower[p,k])
assert lower=={(5,2):18,(5,3):14,(5,6):8,(7,2):36,(7,4):22}
cap={(5,2):F(10,3)**24/5**18,(5,6):F(5,3)**24/5**8,
     (7,2):F(98,9)**24/7**36,(7,4):F(7,3)**48/7**22}
assert all(v<1 for v in cap.values())
print('PASS all four exact density comparisons')
