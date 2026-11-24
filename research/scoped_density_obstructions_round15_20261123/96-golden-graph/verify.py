from fractions import Fraction as F
MAX=10

def add(x,y):
 return [a+b for a,b in zip(x,y)]
def sc(x,a):
 return [a*b for b in x]
def mul(x,y):
 return [sum(x[j]*y[i-j] for j in range(i+1)) for i in range(MAX+1)]
def pw(x,k):
 return [1]+[0]*MAX if k==0 else mul(x,pw(x,k-1))
def ser(d):
 return [d.get(i,0) for i in range(MAX+1)]
def row(fs,i):
 return [f[i] for f in fs]

def alph(a,b):
 n=len(a)
 x=[[F(z) for z in u+[v]] for u,v in zip(a,b)]
 for k in range(n):
  pivot=next(j for j in range(k,n) if x[j][k])
  x[k],x[pivot]=x[pivot],x[k]
  x[k]=[a/x[k][k] for a in x[k]]
  for j in range(n):
   if j!=k:
    v=x[j][k]
    x[j]=[a-v*b for a,b in zip(x[j],x[k])]
 assert [r[:n] for r in x]==[[int(i==j) for j in range(n)] for i in range(n)]
 return [r[-1] for r in x]

def forms(deg):
 pairs=[]
 X=sc(mul(ser({1:1}),pw(ser({0:1,2:1,6:1}),4)),16)
 Y=pw(ser({0:1,1:-2,4:2,9:-2}),4)
 Z=add(X,sc(Y,2))
 XX=sc(Y,-1); YY=sc(X,-1); ZZ=add(XX,sc(YY,2))
 for j in range(deg//2+1):
  g=mul(pw(Z,deg-2*j),pw(X,2*j))
  h=mul(pw(ZZ,deg-2*j),pw(XX,2*j))
  if deg==11:
   T3=pw(ser({0:1,1:2,4:2}),2)
   T4=pw(ser({0:1,1:-2,4:2}),2)
   T2scaled=sc(pw(ser({0:1,2:1,6:1}),2),4)
   g=mul(g,mul(T3,T4))
   h=sc(mul(h,mul(T3,T2scaled)),-1) # strips t^1/2
  pairs.append((g,h))
 return zip(*pairs)

def signed(deg,nlow,nhigh,exp=6):
 fs,gs=forms(deg)
 M=[row(fs,i) for i in (0,2,4)]+[row(gs,i) for i in (range(4) if deg==12 else range(3))]
 vals=[1,0,0]+([0,0,nlow,nhigh] if deg==12 else [0,nlow,nhigh])
 co=alph(M,vals)
 return sum(c*f[exp] for c,f in zip(co,fs))


if __name__=='__main__':
 for degree,expected in ((11,(-1536,917504,32768)),(12,(-12800,196608,4096))):
  a,b,c=[signed(degree,*w) for w in [(0,0),(1,0),(0,1)]]
  assert (a,b-a,c-a)==expected
  print("weight-24 twist", degree, a,b-a,c-a)
 # Coefficients of the untwisted weight-24 extremal theta series.
 E4=ser({0:1,2:240,4:2160,6:6720,8:17520})
 E6=ser({0:1,2:-504,4:-16632,6:-122976,8:-532728})
 Delta=sc(add(pw(E4,3),sc(pw(E6,2),-1)),F(1,1728))
 basis=[pw(E4,6),mul(Delta,pw(E4,3)),pw(Delta,2)]
 theta=add(basis[0],sc(basis[1],-basis[0][2]))
 theta=add(theta,sc(basis[2],-theta[4]))
 assert theta[6]==52416000==4095*12800
 f8=signed(12,0,0,exp=8)
 assert f8==304800 and (theta[8]+4095*f8)/4096==9828000
 assert (8+4+4) > (24+6)//2
 print("PASS: no norm-six-free and balanced norm-eight graph glue within F4 rank six.")
