from fractions import Fraction as F
from math import factorial,prod,comb
from decimal import Decimal, localcontext
def raw(n,k):
    B=[factorial(2*j)//factorial(j) for j in range(k+1)]
    C=[1]+[0]*k
    for i in range(n):
        C=[sum(C[j-t]*B[t] for t in range(j+1)) for j in range(k+1)]
    return [F(C[i]*factorial(i),prod(range(n,n+2*i))) for i in range(k+1)]
def cover(n,k,c=None):
    M=raw(n,k)
    if c is None:c=M[1]
    C=[sum(comb(i,j)*(-c)**(i-j)*M[j] for j in range(i+1)) for i in range(k+1)]
    a=n//2
    return F(n+1)*F(n,n+1)**a*sum(F((-1)**j*comb(a+j-1,j),1)*C[j]/(1+c)**(a+j) for j in range(k+1))
if __name__=='__main__':
 for n in (48,56,80,92):
    h=cover(n,3)
    for k in (1,3,5,7,9,11,13):
        q=cover(n,k)
        print(n,k,float(q),100*(float(q/h)-1))
