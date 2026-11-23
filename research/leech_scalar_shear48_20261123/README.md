# Exact optimum for scalar-sheared Leech product periods in 48 dimensions

Let Λ be the Leech lattice with covolume 1, squared minimum 4 and squared covering radius 2. On the product period (√2Λ)² allow an arbitrary positive definite real 2-by-2 matrix on the factor index, i.e. measure squared length of (x,y) as a‖x‖²+2c⟨x,y⟩+b‖y‖². The maximum center density of a finite union of half-period translates is exactly (3/2)^24, the P48 center density. The same upper bound holds for arbitrary real translates if their squared minimum, after the normalization below, is at least 6. Neither assertion requires an additive motif. In the half-period problem equality forces the period metric to be scalar identity after integral factor-coordinate change, the minimum to be 6, and each projection to enumerate the entire half-period alphabet.

## Normalization and projection

An integral change of basis in the two factor coordinates preserves the period and its half-period subgroup. Gauss reduction and uniform rescaling give a=1, b≥1 and 0≤c≤1/2. Write t=b−c²≥3/4. In these coordinates the squared length is ‖x+cy‖²+t‖y‖², the determinant (covolume) of the period is 2^24 t^12, and the squared packing minimum μ is at most 8 by taking a period vector in the first factor. For a motif of size N, the ratio to P48 is
$$
\mathcal R=\frac{N}{2^{24}}\left[\frac{\mu^2}{36t}\right]^{12}.
$$
Here and below squared distances between phase labels mean the smallest over period translations. Shortest half-period representatives have squared lengths 0,2,3,4, respectively called zero and types 2,3,4.

For any prescribed representative y of a second difference, the squared distance can be made ≤4+t‖y‖² by choosing a first representative x nearest to −cy. Thus for μ>4 the second labels are distinct, and their packing in the projected lattice √(2t)Λ has minimum squared at least μ−4; period differences satisfy this as well, since 8t≥6 and μ≤8. Leech optimality gives N≤[8t/(μ−4)]^12. This alone closes the arbitrary-real case μ≥6:
$$
\mathcal R\le\left[\frac{\mu^2}{18(\mu-4)}\right]^{12}\le1\quad(6\le\mu\le8).
$$

## Half-period metrics with μ>4

Suppose first that some second-label difference is of type 2. If t≥1, the projection estimate yields μ≤4+2t. For 1≤t≤4 this is ≤6√t, by factoring the squared comparison or writing 2(√t−1)(√t−2)≤0. For t≥4 use μ≤8 instead. In all cases R≤1 because N≤2^24. The t≥1 equality possibility requires t=1, μ=6.

If 3/4≤t≤1, take that pair of second labels, and use separately shortest first and second representatives. If the first difference is not type 4 its contribution is at most 3, giving μ≤3+2b (signs of the representatives can be chosen independently). Otherwise the first coset contains 24 orthogonal antipodal pairs of length 2. Taking the inner products with a fixed second representative y of squared norm 2 gives 24 half-integral numbers whose squares sum to 8, so one absolute inner product is ≥1. Choosing the favorable sign gives μ≤4+2b−2c; this also subsumes the non-type-4 first difference. Throughout 1≤b≤1+c²,
$$
(4+2b-2c)^2\le36(b-c^2).
$$
The difference of the two sides is convex in b, so it suffices to check its endpoints. At b=1 the difference is −24c+40c², and at b=1+c² the unsquared left side is at most 6. Equality in this region forces c=0 and b=1. The same argument at t=1 eliminates c>0 from the earlier t≥1 case.

If no second difference is of type 2, the Hoffman bound in the appendix gives N≤4096. When t≥1, μ≤8 then implies R≤(8/(9t))^12<1. When t≤1, choose any pair if N≥2. Its nonzero second label has squared norm 3 or 4. The same first-frame argument makes one absolute dot ≥1 if the first label is type 4 (sum of squares is 12 or 16); if not, the first squared norm is at most 3. In either case μ≤4+4b−2c. The convex-in-b comparison
$$
(4+4b-2c)^2<72(b-c^2)
$$
follows at b=1 from −8−32c+76c²<0 (convex in c, negative at both endpoints), and at b=1+c² from 4+4b−2c≤8. Thus R<1. For N=1, simply use μ≤8 and t≥3/4.

## Half-period metrics with μ≤4

For fixed second label, the first labels form a 24D union of translates of √2Λ of squared minimum at least μ: realize the second difference by zero. Leech optimality bounds each fiber by (8/μ)^12, hence N≤2^24(8/μ)^12. When μ<(9/2)t this gives R<1.

It remains to handle (9/2)t≤μ≤4. Write H=Λ/2Λ and define f:H→ℝ by f(0)=24, f(v)=1 on the 98,280 type 2 classes, and 0 otherwise. The exact Fourier spectrum below makes f positive definite. Consequently
$$
F(x,y)=f(x)f(y)f(x+y)
$$
is positive definite on H²: compose with linear maps and multiply, using convolution of nonnegative Fourier transforms. The only nonzero pairs in its support have two type 2 entries and one zero among x,y,x+y (three pure patterns), or all three type 2. In the latter case two shortest root representatives u,v of the first two classes in Λ have |⟨u,v⟩|=2. Indeed the shorter of u±v has squared norm ≤8 and lies in a type 2 class. Such a class has no nonminimal representative of squared norm ≤8: together with an antipodal norm-4 representative, its plus and minus would be nonzero in 2Λ, with mean squared norm at most 12, below 16. Thus the shorter u±v has squared norm 4. Normalized in our torus by division by √2, the first two shortest representatives give a distance squared at most 2+2b−2c. The pure patterns give at most 2, 2b or 2+2b−4c, all no greater than the same bound. Also
$$
(9/2)t-(2+2b-2c)\ge \tfrac12+2c-(9/2)c^2>0,
$$
using concavity and the endpoint values 1/2 and 3/8. Hence all nonzero support points of F are forbidden differences.

The sum of F over H² is 24³+3·24·98280+98280·4600. The three pure patterns give the middle term; for each fixed type 2 x, exactly 4600 labels y have both y and x+y type 2 (see below). By averaging F over ordered motif pairs and using nonnegative Fourier coefficients,
$$
N\le \frac{24^3\,2^{48}}{24^3+3\cdot24\cdot98280+98280\cdot4600}.
$$
This rational cap is ≈8.4741·10^9, strictly below the ≈8.9463·10^9 labels needed even to tie P48 when μ≤4 and t≥3/4, namely 2^24(27/16)^12. The respective exact constants are 4398046511104/519 and 150094635296999121/16777216. This closes the final metric range.

## Leech alphabet facts and attainment

Two unequal, nonantipodal vectors of squared norm 8 in the same 2Λ class are orthogonal: both their sum and difference lie in 2Λ and have squared length at least 16. There are at most 48 norm-8 representatives per such class. Likewise the norm-4 and norm-6 shells each give exactly two representatives per class they meet and do not overlap each other or norm-8 classes, by the same sum/difference averaging. The three shell sizes are 196560, 16773120, 398034000 (the coefficients of E4³−720Δ). The covering radius gives a representative of norm ≤8 in every class of 2Λ. Since 398034000 = 48(2^24−1−196560/2−16773120/2), each norm-8 class has a complete frame of 24 orthogonal antipodal pairs, as used above.

The characters of H are indexed by v∈H, x↦(−1)^(v·x), since Λ is integral unimodular. Represent v by a least-length vector. Applying the minimum to its sums and differences with norm-4 vectors bounds the absolute dots at 2,3,4 for type 2,3,4 respectively, except for the two coincident vectors giving 4 in type 2. Even moments of the shortest shell, a spherical 11-design, now give the absolute dot counts indexed from 0:
* type 2: 93150, 94208, 9200, 0, 2;
* type 3: 75900, 97200, 22356, 1104;
* type 4: 65780, 94208, 32384, 4096, 92.

On the type-2 Cayley graph this gives respectively 4072, −24 and −24, and the trivial eigenvalue is 98280. Thus f has nonnegative Fourier coefficients. Hoffman's bound is 2^24·24/(98280+24)=4096. For the adjacency matrix A of this graph, the spectral identity (A−4072I)(A+24I)=552J follows from the four displayed eigenvalues. At a type 2 entry it gives (A²)_{0,x}=4072−24+552=4600. All shell and inequality arithmetic is replayed by the companion script.

P48 is attained in the half-period case by the extremal even unimodular lattice admitting an involution with both eigensublattices √2Λ (Nebe Theorem 4.10). Its phases lie in half the period because twice a lattice point is the sum of its integral projected doubles. In the equality case analyzed above c=0 and b=t=1, μ=6 and N=2^24. Both projections are injective since each factor's covering radius squared is 4. This proof does not identify all maximizing motifs or impose restrictions on general 48D lattices beyond the specified period representation.

## References
- Cohn, Kumar, Miller, Radchenko and Viazovska, [The sphere packing problem in dimension 24](https://arxiv.org/abs/1603.06518).
- Borcherds, [The Leech lattice and other lattices](https://math.berkeley.edu/~reb/papers/leech/leech.pdf), covering radius and deep holes.
- Nebe, [On automorphisms of extremal even unimodular lattices](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/DIM48.pdf), Theorem 4.10.
