# Exact 48D and 72D optima for Leech product periods

Let Λ be the Leech lattice of covolume 1, squared minimum 4, and squared covering radius 2. A center density does not include the unit-ball-volume factor. In the following two statements the motif may be arbitrary: it does not have to be a subgroup. All equalities are within the described families, not statements about global sphere packing.

## 48D

Consider a finite union of distinct translates of an orthogonal period
√2Λ ⊕ √(2r)Λ, with r≥1. Any pair of similar Leech periods can be normalized to this form. Let μ be the squared minimum between centers; the period itself implies μ≤8. The P48 density is an upper bound in either of these settings:

1. The shifts belong to half the period, with no restriction on μ.
2. The shifts are arbitrary real vectors and μ≥6.

Writing N for the number of translates, the ratio to P48 is N(μ/6)^24/(4r)^12. If μ≥6, the second projections are distinct modulo √(2r)Λ. Indeed the first squared distance on its torus is at most 4, so equality of the second projections would contradict μ. The periodic union of second projections has squared minimum at least μ−4: the first coordinates of any difference can always be put at torus distance at most 4. Differences within the same second-period translate have squared minimum 8r≥μ−4. The optimal Leech packing in 24D therefore bounds N by (8r/(μ−4))^12. It follows that the ratio is at most [μ²/(18(μ−4))]^12 ≤ 1, because μ²−18μ+72=(μ−6)(μ−12)≤0.

If 4<μ<6 and the shifts are half-period, injectivity of the second projection still holds, so N≤2^24 and the ratio is <1. If μ≤4, take fibers with identical first projection. There are at most 2^24 fibers. In each, the Leech 24D packing theorem gives no more than (8r/μ)^12 elements. Substitution gives ratio at most (2μ/9)^12<1.

The half-period class therefore attains equality only at r=1, μ=6 and N=2^24; equality forces both projection maps to be bijections of the respective half-period labels. In the arbitrary-shift case equality requires μ=6 and the projected 24D packing to attain the 24D optimum. An attaining packing exists: Nebe's Theorem 4.10 explicitly guarantees an extremal even unimodular rank-48 lattice with an involution whose two fixed eigensublattices are √2Λ and √2Λ. Its shifts are half-period since twice any lattice point is the sum of its two integral projected doubles. Consequently the theorem also treats any real lattice with an involution and two Leech-similar eigensublattices; no integrality or modularity hypothesis is imposed on that lattice.

## 72D

Consider a finite union of half-period translates of √2Λ ⊕ √(2r)Λ ⊕ √(2s)Λ, ordered so 1≤r≤s. Its center density is at most 2^36, with equality attained by Nebe's Γ72. Write μ≤8 for the squared minimum and N for the motif size. The ratio to Γ72 is N(μ/8)^36/[2^36(rs)^12].

* For μ≤4, partition by the first two labels (at most 2^48 fibers). Each fiber in the third factor has at most (8s/μ)^12 points by Leech optimality. The ratio is ≤2^12(μ/8)^24/r^12≤2^−12.
* For 4<μ≤min(8,4+2r), the last two labels jointly determine the point: otherwise only the first factor contributes, at squared torus distance ≤4. Thus N≤2^48 and the ratio is ≤2^12[min(8,4+2r)/(8r^(2/3))]^36. The expression in square brackets is nonincreasing for r≥1: on [1,2] differentiate (4+2r)r^(−2/3), and above 2 the numerator is constant. The resulting bound 2^12(6/8)^36 is strictly below 1.
* For μ>4+2r, partition by the third label (≤2^24 fibers). Inside a fiber the second projection is distinct modulo its period. Moreover no two second labels can differ by Leech type 2: that would contribute 2r in the second factor and at most 4 in the first. The following exact spectral argument gives at most 4096 elements in each fiber. Thus N≤2^36 and the ratio is ≤1, strictly unless r=s=1 and μ=8.

The attaining Γ72 contains three orthogonal √2Λ periods, explicitly from its polarized construction: in the integer three-Leech model with metric divided by two, the common period is (2Λ)^3 and its index is 2^36.

### Exact 4096 bound

In Λ/2Λ every nonzero element has a shortest representative of squared norm 4, 6, or 8. The 98,280 norm-four classes define an undirected Cayley graph. Its character indexed by v is x↦(−1)^(v·x) because Λ is unimodular. For nonzero v choose its least representative. If its norm is 4, a shortest Leech vector has absolute dot with v in {0,1,2,4}, with exactly two vectors of absolute dot 4. If the norm of v is 6, possible dots are ≤3; if 8, ≤4. The exclusions follow by applying the minimum 4 to the sums and differences, except for the two coincident rays at norm 4. The known 11-design moments of the shortest shell give the full antipodal dot counts. For absolute dot values in increasing order starting at 0, the full signed vectors have these counts:
* Norm 4: 93150,94208,9200,0,2.
* Norm 6: 75900,97200,22356,1104.
* Norm 8: 65780,94208,32384,4096,92.

The Cayley eigenvalues for these v are respectively 4072, −24, −24, and the trivial eigenvalue is 98280. The ratio bound for independent sets is 2^24·24/(98280+24)=4096. This uses the exact spectrum, not a numerical heuristic. The companion script checks the dot moment systems and finite comparisons using rational arithmetic.

## Inputs

* [Cohn, Kumar, Miller, Radchenko and Viazovska, 24D Leech optimality](https://arxiv.org/abs/1603.06518).
* [Borcherds, The Leech lattice, covering radius √2 and deep holes](https://math.berkeley.edu/~reb/papers/leech/leech.pdf).
* [Nebe, On automorphisms of extremal even unimodular lattices, Theorem 4.10](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/DIM48.pdf).
* [Nebe, An even unimodular 72-dimensional lattice of minimum 8, §2.2 and Lemma 3.1](https://arxiv.org/abs/1008.2862).

No catalogue entries or unrestricted bounds are changed by these theorems.
