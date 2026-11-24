from verify import Q5,F,itertools,Path,json,gegenbauer
def gram(pairs):
 K=len(pairs)
 G=[[0]*(2*K) for _ in range(2*K)]
 for i,j in itertools.product(range(K),repeat=2):
  p,b=pairs[i][j]
  G[2*i][2*j]=p
  G[2*i+1][2*j]=G[2*i][2*j+1]=(p+b)//2
  G[2*i+1][2*j+1]=(3*p+b)//2
 return G
def valid(pairs,R):
 G=gram(pairs); K=len(G)
 for t in itertools.product(range(-R,R+1),repeat=K):
  if not any(t):continue
  GT=[sum(a*b for a,b in zip(row,t)) for row in G]
  q=sum(a*b for a,b in zip(t,GT))
  if q==0:continue
  if q<8:return False
  if q>10:continue
  s=[z for i in range(0,K,2) for z in (2*t[i+1]-t[i],2*t[i]+t[i+1])]
  if abs(sum(a*b for a,b in zip(GT,s)))>5*(q-8):return False
 return True
C=[(i,j,F(c)) for i,j,c in json.loads(Path(__file__).with_name('nested.json').read_text())]
assert all(i>=j>=0 and i+j and c>0 for i,j,c in C)
n=0
for p in range(-4,5):
 for b in range(p%2,21,2):
  a=Q5(F(5*p-12,28),F(b,28))
  h=Q5(F(5*p-12,28),F(-b,28))
  if (Q5(1)-abs(a)).sign()<0 or (Q5(1)-abs(h)).sign()<0:continue
  if not valid([[(8,0),(p,b)],[(p,b),(8,0)]],2):continue
  z=[[(12,0),(4,0),(0,0),(0,0)],
     [(4,0),(8,0),(4,0),(4,0)],
     [(0,0),(4,0),(8,0),(p,b)],
     [(0,0),(4,0),(p,b),(8,0)]]
  if not valid(z,1):continue
  n+=1
  x=gegenbauer(a,max(i for i,j,c in C),46)
  y=gegenbauer(h,max(i for i,j,c in C),46)
  assert (Q5(1)+sum(c*(x[i]*y[j]+x[j]*y[i])/2 for i,j,c in C)).sign()<0,(p,b)
B=1+sum(c for i,j,c in C)
assert n==14 and B<29700
print('PASS nested local bound',B,n)
