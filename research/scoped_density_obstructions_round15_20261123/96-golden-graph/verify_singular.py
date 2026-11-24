from verify import *
a0,a1,a2=[signed(12,*u,exp=8) for u in [(0,0),(1,0),(0,1)]]
assert (a0,a1-a0,a2-a0)==(304800,12058624,-98304)
X=sc(mul(ser({1:1}),pw(ser({0:1,2:1,6:1}),4)),16)
Y=pw(ser({0:1,1:-2,4:2,9:-2}),4)
def basis(x,y):
 p=mul(mul(pw(x,8),pw(y,2)),pw(add(x,y),2))
 return mul(p,pw(add(x,sc(y,2)),2)),mul(p,pw(x,2))
u,v=basis(X,Y)
assert (u[8],v[8])==(2**34,0)
u,v=basis(sc(Y,-1),sc(X,-1))
assert (u[2],v[2],u[3],v[3])==(256,256,0,-16384)
for w in (4**6,4**7,4**8):
 def expression(a,x,y):
  b=F(25,8)*(w-4096)-48*a
  K8=F(39007332000+(w-1)*304800+12058624*a-98304*b,w)
  T2=4*x-F(4,5)*a
  T3=48*(a-x)+8*y-F(4,5)*b
  return F(2**34,w)*(T2/256+T3/16384)+F(64,5)*K8
 vals=[expression(*c) for c in ((0,0,0),(1,0,0),(0,1,0),(0,0,1))]
 base=vals[0]
 assert base==F(2048,w)*(256901120-1295*w)>0
 assert tuple(v-base for v in vals[1:])==tuple(F(2048,w)*v for v in (122880,106496,4096))
print("PASS: singular rank-seven and rank-eight obstructions")
