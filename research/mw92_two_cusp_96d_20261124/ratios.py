from verify import *
SQ=G[::96]
C=[set()for _ in range(8)]
for i in ('0','1','2','3','16'):
 A,B=get(i)
 w=[m(evaluate(A,(b,0)),inv(evaluate(B,(b,0))))for b in range(47)]
 for b,v in enumerate(w):
  ar=C[log[v]%8];f=inv(v)
  for j,u in enumerate(w):
   if j==b:continue
   r=m(u,f);ar.update(m(r,z)for z in SQ)
print('ratios by projected second-residue class modulo 8:',[len(r)for r in C])
assert all(len(r)==2208 and O not in r for r in C)
print('PASS: both-cusp additive residues of norm 18 meet every projective line.')
