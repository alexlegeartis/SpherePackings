# Place next to the release's replay.py.
from fractions import Fraction as F
from itertools import product
import replay as r
d=[r.d[O[0]] for O in r.orbits]
H=[[a*b for a,b in product(h,k)] for h,k in product(r.H,repeat=2)]
sz=[a*b for a,b in product(r.sizes,repeat=2)]
z={i:F(x) for i,x in [
    (1,'160/57'),(2,'11075/1368'),(11,'1221/304'),(12,'4075/456'),
    (13,'15'),(17,'325/912'),(21,'223/19'),(22,'10713/304'),
    (23,'1659/608'),(27,'513/16'),(28,'2349/38'),(29,'1371/304'),
    (31,'2741/912'),(32,'1205/2736'),(33,'1705/608'),(34,'1509/304')]}
assert all(x>0 for x in z.values())
v=[s+sum(x*H[i][j] for i,x in z.items()) for j,s in enumerate(sz)]
B=F(30213,38); C=F(2268,19)
assert v[0]==B
for j,(a,b) in enumerate(product(d,repeat=2)):
    if j and a+b>=8:
        k=(j%8)*8+j//8
        ax=min(a,b)==0 and max(a,b)==12
        assert (v[j]+v[k])/2<=C*ax
threshold=F(923)*(923-B)/C
print('PASS: M <=',B,'+',C,'* P_axis/M')
print('PASS: no axis repeats implies M<=795; M>=923 implies P_axis >=',threshold)
