from fractions import Fraction as F
L=F(1,3500); U=F(128,1875)
def fg(l):
   # integers 1/1000 increments, comparison is exact
   a,b=0,4000
   while b-a>1:
    n=(a+b)//2
    if 2**40*F(n,8000)**24<=16800*l**4:a=n
    else:b=n
   return F(a,1000)
def limits(g,h,d0,c0=F(2),b0=F(2)):
 assert 0<g<12
 d0=max(d0,c0,b0)
 T=(3*g-36*h/(b0*c0*d0))/(12-g)-1/d0
 if T<=0:return None
 B=min(2/T,1/(T-1/c0)) if T>1/c0 else 2/T
 if B<=b0:return F(0)
 worst=F(0);u=b0
 while u<B:
  v=min(B,u+F(1,50))
  if g*g/64 <= v/(3*d0):return None
  C=(h/d0+v/3)/(g*g/64-v/(3*d0))
  if T>1/u:C=min(C,1/(T-1/u))
  if C>max(c0,u):
   p=384*(h+v*C/3)/d0
   worst=max(worst,p)
  u=v
 return worst
def good(g,h,d0,c0=F(2),b0=F(2)):
 q=limits(g,h,F(d0),F(c0),F(b0));return q is not None and q<=g**3
table=[(12,181,50722),(11,128,22510),(9,90,10053),(8,64,4516),
       (7,46,2033),(6,33,928),(6,24,429),(5,13,161)]
lo=L
for B,C,W in table:
 hi=min(2*lo,U)
 g=fg(lo)
 assert 2**40*(g/8)**24/lo**4<=16800
 assert good(g,hi,B,2,B)
 assert good(g,hi,C,C,2)
 assert good(g,hi,W,2,2)
 lo=hi
assert lo==U
print("PASS all eight determinant slabs and their three unbounded tails")
