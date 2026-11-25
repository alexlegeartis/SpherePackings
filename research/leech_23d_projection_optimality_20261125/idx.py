from sympy.polys.matrices import DomainMatrix
from sympy import ZZ
import random,math
A=[[int(z) for z in l.split()] for l in open('triangles')];P=[[a[j]-a[j+24] for j in range(23)] for a in A];r=random.Random(16);g=0;kee=[]
for i in range(50):
 ids=r.sample(range(len(P)),23)
 B=DomainMatrix.from_list([P[z] for z in ids],ZZ)
 d=int(B.det());assert d%(8**11)==0
 ng=math.gcd(g,abs(d)//(8**11))
 if ng!=g:print('minor',ids,'det/P',d//(8**11),'newgcd',ng);kee+=ids
 g=ng
 if g==1:break
print('index divides',g,'witness triangles',len(set(kee)));assert g==1

