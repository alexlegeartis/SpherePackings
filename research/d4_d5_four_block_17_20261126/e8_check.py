"""All E8/2E8 symbol distances from the even-sum model, without a table."""
from math import isqrt
R=[[1,-1,-1,-1,-1,-1,-1,1],[2,2,0,0,0,0,0,0]]
for j in range(6):
 v=[0]*8;v[j]=-2;v[j+1]=2;R.append(v)
# Check the 8 rows span the doubled coordinate realization of E8 (covolume 256).
def det(m):
 if len(m)==1:return m[0][0]
 return sum((-1)**i*x*det([r[:i]+r[i+1:] for r in m[1:]])for i,x in enumerate(m[0])if x)
assert abs(det(R))==256
d={}
def visit(x,r,odd):
 if len(x)==8:
  if sum(x)%4:return
  a=[sum(i*j for i,j in zip(t,x))for t in R];assert not any(y%4 for y in a)
  k=sum(((y//4)%2)<<i for i,y in enumerate(a));z=(16-r)//4;assert z in (0,2,4)
  d[k]=min(d.get(k,z),z);return
 for i in range(-isqrt(r),isqrt(r)+1):
  if i%2==odd:visit(x+[i],r-i*i,odd)
for e in range(2):visit([],16,e)
assert len(d)==256 and [list(d.values()).count(i)for i in (0,2,4)]==[1,120,135]
print('E8/2E8 exact squared depths before scaling: 0:1, 2:120, 4:135')
