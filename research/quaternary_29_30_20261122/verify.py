#!/usr/bin/env python3
"""Integer verifier for this nonintegral Construction A4 lattice."""
from collections import Counter
from fractions import Fraction
from itertools import combinations,product
import json
from pathlib import Path
P=Path(__file__).resolve().parent
obj=json.loads((P/'witness.json').read_text())
n=obj['dimension']
positions=[j for j in range(32) if j not in obj['omitted_coordinates']]
assert n==len(positions)
retained=sum(1<<j for j in positions)
X=obj['residue_masks']; Z=obj['lift_masks']; k=len(X)
assert len(Z)==k and all(0<=u<=retained and u&~retained==0 for u in X+Z)
def col(a):return 1|a<<1
def syn(bits):
 s=0
 for i in positions:
  if bits>>i&1:s^=col(i)
 return s
# Nonzero even codewords correspond to binary torsion words of minimum
# weight at least 4. A zero sum of <=3 columns is excluded by the common
# leading 1 (odd case) and by distinctness (pair case).
a=[col(j) for j in positions]
assert len(set(a))==n and all(i&1 for i in a)
s={0}
for u in a:s|={v^u for v in s}
assert len(s)==64
assert all(syn(x)==0 for x in X) # doubles of all new generators lie in 2T
C=set(); minima=Counter(); residue=Counter()
def affine_weight(bits,h):
 # Minimize all twos of a word whose odd positions are fixed. On an odd
 # position, toggling the upper bit changes 1 to 3 at no cost.
 d=[n+1]*64; d[0]=0
 for j in positions:
  u=col(j); c=int(not (bits>>j&1))
  d=[min(d[t],d[t^u]+c) for t in range(64)]
 return bits.bit_count()+4*d[h]
for bits in product((0,1),repeat=k):
 r=0; h=0
 for i,b in enumerate(bits):
  if b:
   h^=Z[i]^(X[i]&r)
   r^=X[i]
 assert r not in C
 C.add(r)
 if r:
  m=affine_weight(r,syn(h))
  assert m>=16
  minima[m]+=1
  residue[r.bit_count()]+=1
# Every word is uniquely an odd-generator selection plus a word of 2T.
# The literal 4Z^n frame attains Euclidean squared length 16, with no shorter
# vectors over any verified codeword. For L=(preimage B)/2, det L is its
# Gram determinant, and min L = 4; hence its squared center density is 1/det L.
code_power=n-6+k
det=Fraction(2**(2*n),2**(2*code_power))
nn=Fraction(1,det)
# A codeword pair with non-integer scaled dot demonstrates nonintegrality.
def row(x,z):return [((x>>i&1)+2*(z>>i&1)) for i in positions]
assert any(sum(a*b for a,b in zip(row(x,z),row(y,w)))%4
           for (x,z),(y,w) in combinations(zip(X,Z),2))
# Restrict to the hyperplane y_2 = y_3 in the normalized lattice.
# The root e_2-e_3 lies in the dual because the two residue columns agree.
section_normal=(2,3)
i,j=section_normal
assert all((b>>i&1)==(b>>j&1) for b in X)
# Verify primitivity: a vector from the binary torsion code pairs to 1.
u=sum(1<<j for j in (2,4,8,14))
assert u&~retained==0 and syn(u)==0 and (u>>i&1)-(u>>j&1)==1
# Pairing map L -> Z is surjective; section covolume = sqrt(2)*covolume(L).
# The section inherits min >=4 and contains the norm-4 frame vector 2e_4.
assert 4 in positions and 4 not in section_normal
out=dict(rows=[
    dict(dimension=29,minimum_squared=4,gram_determinant=str(2*det),
         center_density_squared=str(nn/2),
         preceding_repository_center_density_squared="16/33",
         cohn_lattice_center_density_squared="1/3",
         cohn_unrestricted_center_density_squared="1/2"),
    dict(dimension=30,minimum_squared=4,gram_determinant=str(det),
         center_density_squared=str(nn),
         preceding_repository_center_density_squared="64/105",
         cohn_lattice_center_density_squared=str(Fraction(3**27,2**44)),
         cohn_unrestricted_center_density_squared="1")
], residue_weight_histogram=dict(sorted(residue.items())),
   residue_euclidean_minima=dict(sorted(minima.items())))
for r in out['rows']:
    q=Fraction(r['center_density_squared'])
    assert q==Fraction(r['minimum_squared'],4)**r['dimension']/Fraction(r['gram_determinant'])
    assert Fraction(r['cohn_lattice_center_density_squared'])<Fraction(r['preceding_repository_center_density_squared'])<q
    assert q==Fraction(r['cohn_unrestricted_center_density_squared'])
print(json.dumps(out,indent=2))
