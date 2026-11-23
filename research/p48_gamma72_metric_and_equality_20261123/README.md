# Quantitative shape rigidity and a half-period equality classification

These are exact, scoped restrictions, not new lower bounds, global upper bounds, or priority claims.

## Arbitrary metrics on a fixed extremal even unimodular lattice

Let L be a 48-dimensional even unimodular lattice of squared minimum 6, or a 72-dimensional one of squared minimum 8 (e.g. P48n or Γ72). Measure its vectors using an arbitrary positive definite Euclidean quadratic metric G. If the ratio of the largest to smallest eigenvalue of G is at most **2 in 48D**, or **19/10 in 72D**, then the center density is no larger than that of L in its original metric. Apart from uniform scaling, the decrease is strict. All diagonal and off-diagonal entries may vary; no coordinate factorization or reflection symmetry is required. This does not compare shapes outside the indicated condition-number balls.

### Design and normalized moments

For an even unimodular lattice, harmonic theta series of harmonic degree k have elliptic modular weight n/2+k. The minimum implies at least order 3 at infinity in 48D, or 4 in 72D. Dividing a positive-degree harmonic theta series by Δ³ or Δ⁴ yields a holomorphic form of weight k−12, hence it vanishes for 0<k≤10. The lattice is centrally symmetric, and its shortest shell, rescaled to radius one, is a spherical 11-design.

Normalize G to mean eigenvalue 1. Write its least and greatest eigenvalues 1−a and 1+b, and put s=a+b, p=a/s. The scalar case needs no proof, so 0<p<1. The transformed shortest-shell squared norm divided by the original minimum is 1+sY, where Y is the value on a unit shell vector of the quadratic form with eigenvalues x_i=(G_i−1)/s. These eigenvalues have mean zero and lie in [−p,1−p]. Set
```
w = mean x_i² = u p(1−p),  0<u≤1,
z_j = mean x_i^j / w.
```
The upper bound on w follows by averaging (x+p)(1−p−x)≥0. Using the 11-design and standard Gaussian radial moments, the first five moments m_j of Y obey m_1=0 and
```
m_2=2w/(n+2)
m_3/m_2=4 z_3/(n+4)
m_4/m_2=(24 z_4+6nw)/((n+4)(n+6))
m_5/m_2=(192 z_5+80nw z_3)/((n+4)(n+6)(n+8)).
```
For clarity, the Gaussian quadratic-form cumulants are 2^(r−1)(r−1)! sum_i x_i^r; with zero trace the ordinary fourth moment is the fourth cumulant plus three times the square of the second, and the fifth adds ten times the product of the second and third cumulants. Divide each by the exact Gaussian radial factor n(n+2)...(n+2j−2). This uses only shell polynomial degrees ≤10.

Cauchy--Schwarz, or the nonnegative average of (x²−z_3 x−w)², gives z_4≥w+z_3²; pointwise x⁵≤(1−p)x⁴ gives z_5≤(1−p)z_4. Also z_4≤max(p,1−p)². The sharper third-moment interval
```
−p+u(1−p) ≤ z_3 ≤ 1−p−up
```
follows by averaging (1−p−x)(x+t)²≥0 and optimizing t=w/(1−p), and applying the same argument to −x.

### Determinant comparison

Let g=(det G)^(1/n), and let Q be the allowed condition number (2 or 19/10). Its definition gives s≤L=(Q−1)/(1+(Q−1)p). Put A=pL. Summing the logarithm series, or applying its monotone remainder bound for any eigenvalue increment y≥−a, shows
```
log(1+y) ≥ y − (1/2+a/3+a²/4+a³/(5(1−a))) y² .
```
For a negative y, apply the geometric bound to the remaining sum of terms after degree four. For a positive y, the elementary inequality log(1+y)≥y−y²/2 is stronger. Averaging and using exp(−v)≥1−v gives
```
1−g ≤ s² w (1/2+a/3+a²/4+a³/(5(1−a))) ≤ s R w
R=L(1/2+A/3+A²/4+A³/(5(1−A))).
```

### Exact interval certificate

The verifier establishes on the full unit square (p,u) that some shell vector satisfies Y<−r, where r=Rw. Only w>0 is needed in the theorem. If all shell values were at least −r, every polynomial weight P would satisfy mean ((Y+r)P(Y)²)≥0. Testing linear P gives the necessary determinant inequality r(r+m_3/m_2)≥m_2. Its failure follows if
```
R (r + 4(1−p−up)/(n+4)) < 2/(n+2).
```
For boxes where this first test is insufficient, choose P(y)=y²−c y−d, with c an affine function of p,w and d being w times an affine function of p,w, as specified in the verifier. The corresponding mean, divided by m_2, is
```
m_5/m_2 + (r−2c)m_4/m_2 + (c²−2d−2cr)m_3/m_2
 + 2cd + r(c²−2d) + R(n+2)d²/2.
```
Replace z_5 by (1−p)z_4. The expression becomes K+Z_4 z_4+Z_3 z_3 with coefficients written explicitly in the verifier. If a box has Z_4<0, replace z_4 by its lower bound w+z_3². Maximizing the resulting concave quadratic on the allowed z_3 interval uses only rational arithmetic, including its exact vertex. Otherwise use the coarse upper bound on z_4. Either method proves the weighted mean is strictly negative on every accepted box. All subdivisions, inequalities, and rounding are checked using integers. This contradiction gives a shorter shell vector with squared norm below the original minimum times g. Since covolume is multiplied by g^(n/2), the new density is strictly lower for non-scalar G.

## Equality with a single scaled Leech period in 48D

Consider an arbitrary periodic packing with one pure horizontal period √2Λ, where Λ is the Leech lattice of covolume 1 and squared minimum 4. In the regime μ≥6 of the exact one-period bound, attaining the optimal density (3/2)^24 forces μ=6 and the projection to be the unique densest periodic 24D packing, a Leech lattice of squared minimum 2. Each projected position has precisely one horizontal coset. After scaling both coordinates by √2, this is a phase assignment on Λ with values in Λ_R/2Λ. It is convenient to divide the horizontal phases by two to view them in Λ_R/Λ. Each shortest vertical step must change its horizontal phase by a Leech deep hole (squared depth 2 in this last normalization).

Every such deepest hole has 24 independent differences between closest lattice points. Otherwise a small displacement perpendicular to their span strictly increases all their squared distances, and all other lattice points remain strictly farther by local discreteness. After selecting one closest point as origin, the dot products of the hole with these 24 differences equal the corresponding squared lengths divided by two. Because the Leech lattice is integral and unimodular, Cramer's rule puts the hole in (1/D)Λ for a nonzero integer determinant D of the 24 difference rows. Each difference has squared length at most 8 by the triangle inequality, so |D|≤8^12 by Hadamard. The 80-root basis-generation check in leech_quad.py proves the shortest shell generates Λ. Thus after one common horizontal translation, **every phase in every attaining periodic motif belongs to a finite universal torsion alphabet**. For a lattice the phase assignment is a homomorphism, so this gives a finite search without assuming beforehand that the phases are half-periods. No enumeration of all torsion alphabet assignments is claimed.

### Full classification conditional on additive half-periods

In particular, suppose the attaining packing is a lattice and every horizontal label is a half-period. Its physical coordinate form is
```
   {(x/√2, y/√2) : x,y∈Λ,  x mod 2Λ = F(y mod 2Λ)}
```
for a binary linear map F on V=Λ/2Λ. The vertical projection is Λ/√2 as above. The opposite pure vertical intersection contains √2Λ because F has exponent two; it cannot have any other cosets. Indeed the same covering argument in the reversed direction with period √2Λ and μ=6 says that every fiber is just one coset. Thus F is injective and therefore invertible. Equivalently this follows from the two half-period projections being bijective at equality.

Let q(v)=||v||²/2 mod 2 be the nondegenerate quadratic form on V. A vertical minimum-vector class (least raw norm 4) must map to a horizontal class of least raw norm 8, by μ=6 and the Leech covering radius; both kinds have q=0. The script leech_quad.py constructs a basis from 80 minimal Golay vectors using column Hermite normal form; its determinant proves it generates the standard Leech lattice. It then uses 349 exact shortest vectors to certify that the binary quadratic evaluations have rank 299 in the 300-dimensional space spanned by {x_i x_j : 1≤i≤j≤24}, where diagonal monomials are read as x_i. The only nonzero quadratic vanishing on all shortest-vector classes is q. Since F is invertible, q composed with F is nonzero. It must therefore equal q. Consequently every such attaining lattice is an **extremal even unimodular lattice**, not merely an arbitrary real lattice.

Conversely, a q-isometry F gives equality in this form if and only if it sends the set of least-norm-4 classes disjointly from itself. Indeed nonzero Leech half-period classes have least raw norms 4, 6, or 8 (evenness, minimum 4, covering bound 8); the 6 classes are exactly the anisotropic ones and thus map to themselves as a set. Exclusion of a pair of norm-4 components makes all nonzero coset minima at least 12 in raw sum, so the physical minimum is at least 6. Evenness follows from q preservation and covolume one follows from the graph index. This gives a finite orthogonal-group characterization, not a count of all isometry types.

## Replay and sources

Run `g++ -std=c++17 -O3 -o /tmp/metric-checker verify.cpp && /tmp/metric-checker`, or `python verify_fast.py` for the same directed-integer interval bounds. Do not compile with NDEBUG. The exact 48 equality check is `python leech_quad.py` (SymPy).

- [Cohn–Kumar–Miller–Radchenko–Viazovska, unique optimal periodic 24D packing](https://arxiv.org/abs/1603.06518)
- [Nebe, extremal even unimodular Γ72](https://arxiv.org/abs/1008.2862)
- [Nebe, extremal 48D lattice and Leech eigensublattices](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/DIM48.pdf)
- [Borcherds, Leech lattice and covering radius](https://math.berkeley.edu/~reb/papers/leech/leech.pdf)

## Larger cross-shear chambers when the cross operator has low rank

There is an additional shape result for a fixed P48 model having the original Leech decomposition and a reflection in one factor. Change its metric to
```
||(x,y)||_T² = ||x||²+||y||²+2 xᵀT y.
```
For ranks at most 1, 2, 4, 8, 12, and 24, the original density remains maximal whenever the operator norm of T is at most 0.998, 0.981, 0.923, 0.819, 0.742, and 0.6, respectively. The strong low-rank limits are not covered by the preceding condition-number theorem; the latter also allows entirely different metrics.

For a shell vector write D=xᵀTy, S=TᵀT. The original 11-design gives
```
E[D⁴]/E[D²] = (tr S + 2 tr S² / tr S)/26.
```
For example use the symmetric quadratic-form matrix with off-diagonal blocks T/2 and Tᵀ/2 in the Gaussian moment formulas already given, and normalize the direction to radius √6. The reflection maps D to −D, so some shortest vector has D no greater than the negative square root of the displayed ratio.

Put k as an upper bound on rank, t as an upper bound on operator norm, and h=tr S/t² in [0,k]. The rank inequality tr S²≥(tr S)²/k gives a new squared minimum at most 6(1−a√(h/k)), where a=t√((k+2)/26)/3. Concavity of log(1−s) on [0,t²] gives det(I−S)≥(1−t²)^h. The density ratio is then at most
```
(1−a√(h/k))²⁴ / (1−t²)^(h/2).
```
The function −log(1−x)/x² decreases on [0,1/3]: the numerator controlling its derivative is x/(1−x)+2log(1−x)<0, using x/(1−x)≤(3/2)x and log(1−x)≤−x. Thus the bound is ≤1 whenever (1−a)^48≤(1−t²)^k. The short rational replay `python cross_ranks.py` checks precisely this endpoint condition using a lower rational enclosure on the square root and integer powers.
