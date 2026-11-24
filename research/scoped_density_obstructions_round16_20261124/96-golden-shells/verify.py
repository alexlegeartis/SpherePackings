import itertools, json
from pathlib import Path
from fractions import Fraction as F
D=48
class Q5:
    def __init__(self,a=0,b=0): self.a,self.b=F(a),F(b)
    @staticmethod
    def cast(x): return x if isinstance(x,Q5) else Q5(x)
    def __add__(self,y):
        y=self.cast(y); return Q5(self.a+y.a,self.b+y.b)
    __radd__=__add__
    def __neg__(self):return Q5(-self.a,-self.b)
    def __sub__(self,y): return self+-self.cast(y)
    def __mul__(self,y):
        y=self.cast(y);return Q5(self.a*y.a+5*self.b*y.b,self.a*y.b+self.b*y.a)
    __rmul__=__mul__
    def __truediv__(self,y):
        y=self.cast(y)
        det=y.a*y.a-5*y.b*y.b
        return self*Q5(y.a/det,-y.b/det)
    def sign(self):
        a,b=self.a,self.b
        if a>=0 and b>=0:return int(a>0 or b>0)
        if a<=0 and b<=0:return -1
        if a<0:return -(-self).sign()
        return 1 if a*a>5*b*b else -1
    def __abs__(self):return self if self.sign()>=0 else -self

def gegenbauer(t,K,D=48):
 r=[Q5(1),t]
 for k in range(2,1+K):
  r.append(((D+2*k-4)*t*r[-1]-(k-1)*r[-2])/(D+k-3))
 return r
def localdata(num=2):
  entries=[]
  for p in range(-4,5):
   for b in range(p%2,21,2):
    a=Q5(F(3*p-4,20),F(3*b,100))
    u=Q5(F(3*p-4,20),F(-3*b,100))
    if (Q5(1)-abs(a)).sign()<0 or (Q5(1)-abs(u)).sign()<0:continue
    G=((12,6,4,2,4,2),(6,18,2,6,2,6),(4,2,8,4,p,(p+b)//2),
       (2,6,4,12,(p+b)//2,(3*p+b)//2),
       (4,2,p,(p+b)//2,8,4),(2,6,(p+b)//2,(3*p+b)//2,4,12))
    def allowed(t):
     s=tuple(z for i in range(0,6,2) for z in (2*t[i+1]-t[i],2*t[i]+t[i+1]))
     Gt=[sum(a*b for a,b in zip(r,t)) for r in G]
     q=sum(a*b for a,b in zip(Gt,t))
     return q==0 or q>=8 and (q>10 or abs(sum(a*b for a,b in zip(Gt,s))) <= 5*(q-8))
    if all(allowed(t) for t in itertools.product(range(-num,num+1),repeat=6) if any(t)):
        entries.append((p,b,a,u))
  return entries
C=json.loads(Path(__file__).with_name('delsarte.json').read_text())
C=[(i,j,F(c)) for i,j,c in C]
assert all(i>=j>=0 and i+j>0 and c>0 for i,j,c in C)
D=47
for p,b,t,u in localdata():
    x=gegenbauer(t,max(i for i,j,c in C),D)
    y=gegenbauer(u,max(i for i,j,c in C),D)
    f=Q5(1)+sum(c*(x[i]*y[j]+x[j]*y[i])/2 for i,j,c in C)
    assert f.sign()<0,(p,b,float(f.a)+5**.5*float(f.b))
B=1+sum(c for i,j,c in C)
assert B < 11040000
print('PASS exact local upper bound',B,float(B))

assert len(localdata())==25
