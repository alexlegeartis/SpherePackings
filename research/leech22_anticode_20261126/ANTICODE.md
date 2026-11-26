# A twelve-point anticode proves optimality in the canonical projected 22-dimensional Leech carrier

Let Λ be the Leech lattice in its usual Golay coordinates, scaled by 1/√8. Set W={x₀=x₁=x₂}, L=Λ∩W and P=π_W Λ=L*. The normal section is √2 A₂, of Gram determinant 12, so covol(L)=√12 and covol(P)=1/√12 by unimodularity. These are the canonical carriers for the [Conway–Sloane 22D packing](https://arxiv.org/pdf/math/0207182) (original article, p. 3).

**Theorem.** For ANY subset Y of this fixed grid P, periodic or not, with squared separation μ≤4 in these fixed coordinates, its upper center density is at most δ₀=3/√12·(11/12)^11. Y can comprise arbitrary parts of cosets with unrelated selections, not just entire labeled layers. The three-layer Conway–Sloane example in this grid attains δ₀ at μ=11/3. Thus the inequality is sharp. The hypothesis μ≤4 cannot be normalized away inside the fixed carrier: no claim covers P at μ>4, other rank-two normal sections, continuous displacements, or unconstrained 22D packings. This is a restricted result, not a new record.

Here are twelve raw integer Leech vectors, before projection and division by √8. Put B=their projected set, and T=the first four. The six squared distances of T are all 3; the 66 squared distances of B are all ≤11/3. The twelve projections are distinct. Coordinates and all pair distances can be checked in the tiny stand-alone script below. Source for Λ: all raw coordinates have the same parity m∈{0,1}, the vector of ((x_i−m)/2 mod 2) is in the extended binary Golay code, and Σx_i ≡4m mod 8. The chosen cyclic length-23 generator is 0xae3, followed by an overall parity coordinate.

If a finite S⊂P with diameter²<μ is translated by any t∈P, its translate holds at most one point of Y. For any radius R, if diameters of vectors in S about zero are bounded by H, count pairs (y,s) with y∈Y∩(x+B_R), s∈S by t=y−s: t∈P∩(x+B_{R+H}) and can account for at most one pair. Thus |S| times the upper point intensity of Y is at most 1/covol(P)=√12. Lattice point counts in balls have the required leading term uniformly in the center (use any bounded fundamental parallelotope). No measurability or periodicity assumption is needed. Center density equals intensity times (μ/4)^11.

1. If μ≤3, even the full grid gives at most √12(3/4)^11 < δ₀.
2. If 3<μ≤11/3, the four-point T gives intensity at most √12/4. At μ=11/3 the resulting bound is exactly δ₀, and for smaller μ it is smaller.
3. If 11/3<μ≤4, B gives intensity at most √12/12, so center density is at most 1/√12 <δ₀.

This excludes an improvement by arbitrary puncturing, switching or supercell selection on the old projection under the stated normalization. A genuine candidate must leave the grid, change the normal, or operate above the inherited minimum 4.

```python
from collections import Counter
raw=[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[2,-2,0,0,0,-2,2,2,0,2,0,2,0,0,0,0,0,0,0,0,0,0,0,2],[0,2,-2,0,0,0,2,2,-2,0,2,0,2,0,0,0,0,0,0,0,0,0,0,2],[2,0,-2,0,0,0,2,0,-2,2,0,0,0,0,2,0,2,0,0,0,0,0,0,2],[-4,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[4,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-1,-1,1,-1,1,-1,3,1,-1,1,1,1,-1,1,-1,-1,1,-1,-1,1,-1,1,1,1],[1,1,-1,1,-1,-1,3,1,-1,1,1,1,-1,1,-1,1,1,1,-1,1,1,1,1,1],[1,-1,-1,-1,-1,-1,3,1,-1,-1,1,1,1,1,1,-1,1,1,-1,-1,-1,1,1,1],[-1,1,1,-1,-1,1,3,1,-1,1,1,1,1,1,-1,-1,1,1,1,-1,1,1,1,1],[1,-1,1,-1,-1,1,3,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,1,1,1],[-1,1,-1,1,1,1,3,1,-1,1,-1,1,1,1,1,-1,1,1,-1,1,-1,1,1,1]]
G=set()
for m in range(4096):
    c=0
    for j in range(12):
        if (m>>j)&1: c ^= 0xae3<<j
    G.add(c|((c.bit_count()%2)<<23))
assert Counter(c.bit_count() for c in G)=={0:1,8:759,12:2576,16:759,24:1}
w=[]
for x in raw:
    m=x[0]%2
    assert len(x)==24 and all((u-m)%2==0 for u in x)
    assert sum((((u-m)//2)%2)<<i for i,u in enumerate(x)) in G
    assert (sum(x)-4*m)%8==0
    w.append((sum(x[:3]),*x[3:]))
assert len(set(w))==12
# In these coordinates the physical squared distance is d24/24.
def d24(x,y):
    return (x[0]-y[0])**2 + 3*sum((a-b)**2 for a,b in zip(x[1:],y[1:]))
def profile(n):
    return Counter(d24(w[i],w[j]) for i in range(n) for j in range(i))
assert profile(4)=={72:6}
assert max(profile(12))==88 and sum(profile(12).values())==66
assert 4*9**11 < 11**11 and 12**11 < 3*11**11
print(profile(4),profile(12))
```
