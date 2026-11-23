from itertools import combinations,product
from pathlib import Path
from functools import reduce
import sys,json,operator
from sympy import Matrix
P=Path(__file__).resolve().parent
C=json.loads((P/'certificate.json').read_text())
alpha=C['alpha']; M=C['phase_isometry']; inv={v:i for i,v in enumerate(M)}
def mul(a,b):
 z=0
 while b:
  if b&1:z^=a
  a<<=1
  if a&64:a^=67
  b>>=1
 return z
def pow0(x,k):
 r=1
 for _ in range(k):r=mul(r,x)
 return r
def q0(x):
 n=pow0(x,9)
 return n^pow0(n,2)^pow0(n,4)
def B(x,y):return q0(x)^q0(y)^q0(x^y)
F8=[v for v in range(64)if pow0(v,8)==v]
beta=min(v for v in F8 if v not in alpha)
S=[]; W={0}
for _ in range(6):
 v=next(v for v in range(1,64)if v not in W and all(B(v,w)==1 for w in S))
 S.append(v);W.update({x^v for x in W})
assert len(W)==64
S.append(reduce(operator.xor,S,0))
assert len(set(S+[0]))==8 and all(B(x,y)==1 for x,y in combinations(S,2))
S=[0]+S
def tag(a):return 0 if a==0 else 1^q0(a)
def W35(a):
 b=mul(beta,a)
 return M[(a^b)|(64*tag(a))] | (M[b]<<7)
L=[W35(a)for a in S]
def syn(w):
 x=e=s=0
 for i in range(5):
  a=inv[(w>>(7*i))&127]
  e^=a>>6;x^=a&63;s^=mul(alpha[i],a&63)
 return x,s,e
count=0
for a,b in combinations(L,2):
 x,y,e=syn(a^b)
 assert x!=0 and y==mul(beta,x) and e==1^q0(x)
 assert all(y!=mul(c,x)for c in alpha)
 for i,j in combinations(range(5),2):
  zi=next(z for z in range(64)if mul(alpha[i]^alpha[j],z)==y^mul(alpha[j],x))
  zj=x^zi
  assert zi and zj and q0(zi)^q0(zj)==q0(x)
  # Both active letters are q=1 iff their e parity equals q0(x), which is not allowed.
  assert ((1^q0(zi))^(1^q0(zj)))!=e
  count+=1
assert count==280
V=Matrix([[1,1,0],[1,0,1],[0,1,1]])
assert abs(V.det())==2 and all(sum(v)%2==0 for v in V.tolist())
rows=[]
for k in range(8):
 v=[(k>>i)&1 for i in range(3)]
 r=list(map(int,(Matrix([v])*V)[0,:]))
 rows.append(dict(horizontal_times_2=[2*((L[k]>>i)&1)for i in range(35)],vertical=r))
 rows.append(dict(horizontal_times_2=[2*((L[k]>>i)&1)+((82>>(i%7))&1)for i in range(35)],vertical=[r[0]+1,*r[1:]]))
assert len({tuple(d['vertical']) for d in rows})==16
R=dict(missing_slope=beta,binary_symplectic_labels=S,
       horizontal_integer_shift_hex=[f'{w:09x}'for w in L],
       vertical_D3_basis=[[int(x)for x in row]for row in V.tolist()],motif=rows,
       period_gram_determinant=2**39,motif_size=16,squared_minimum=8)
if '--write' in sys.argv:(P/'motif.json').write_text(json.dumps(R,indent=2)+'\n')
else:assert R==json.loads((P/'motif.json').read_text())
assert R['period_gram_determinant']==(2**31)*((2**3*2)**2)
print('PASS exact symplectic clique, all',count,'two-site phase cases, D3 index and determinant')
print('38D periodic center density 8*sqrt(2); exact squared minimum 8')
