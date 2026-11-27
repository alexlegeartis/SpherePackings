# Sphere 158: Rigidity of large ternary Coxeter–Todd layers

2026-11-27. No density record. This proof closes the nonlinear palettes in the fixed `K12² × A_k` ternary-palette construction, for k=1,2,3, that remained open in the exact [Sphere 65 scale-vertex reduction](https://github.com/alexlegeartis/SpherePackings/commit/415eaeb5c0b3cca957148471c82d6a18ad0228b8). Neither arbitrary horizontal shifts, different lattice carriers, nor a non-diagonal metric are covered. The new statement treats arbitrary membership of every palette, not merely additive or parallel layers.

## 1. The finite two-block model

Let V=F₃⁶, Q(x)=Σxᵢ², and σ(0)=0, σ(1)=1, σ(2)=−1. The Coxeter–Todd lattice at minimum squared length 4 has covolume 27. Its ternary phase quotient has minimum squared coset lengths 0 at zero and 4/3, 2, 8/3 for Q=+1, 0 (nonzero), −1, respectively. See the Coxeter–Todd construction [Conway–Sloane, §3](https://neilsloane.com/doc/Me107.pdf) and the [exact quotient verification and its literal model](https://github.com/alexlegeartis/KissingNumbers/tree/48af83ef2b009243848b5b8a9655606225fec72b/verifications/improved/dim36-ct-nonlocal). We choose the sign of Q so +1 is light. The standard elliptic form is Σxᵢ². In the left block make an anti-isometric change of labels: three diagonal copies of the 2×2 matrix with rows (1,1),(1,−1) multiply Q by −1. Thus the two-block site (x,y) carries costs with Q replaced by −Q on the left, not on the right.

A colour is a set of sites with squared two-block distances at least 4. No row x or column y repeats in one colour: a single block costs at most 8/3. Writing a colour as (x,f(x)), it is therefore an injection satisfying
`σ(Q(f(x)−f(x′))) ≤ σ(Q(x−x′))`
for all x,x′ in its domain. This equivalence follows from the fact that for a nonzero difference u the one-block cost is 2−(2/3)σ(Q(u)), with the opposite sign for the left block. For distinct colours at squared two-block separation ≥8/3, the only requirements are: right Q=−1 for a shared row, left normalized Q=+1 for a shared column, no repeated site. Both nonzero blocks already cost at least 8/3.

## 2. Rigidity with at least 440 sites

**Theorem.** In this normalization, every colour of size n≥440 is contained in a single affine isometry graph except possibly for one site. An outsider is possible only for n≤486. The n=486 upper for a colour with an outsider is attained.

Let N=729, m=N−n, S(u,v)=σ(Q(u−v)), and B(u,v)=Q(u+v)−Q(u)−Q(v). The matrix S is real symmetric, has zero row sums, and has operator norm 27. In fact, for ζ a primitive cube root of unity and character v↦ζ^B(t,v), its eigenvalue is −27(ζ^(−Q(t))−ζ^Q(t))/(i√3), by completing the square; the Gauss sum of Q is (1+2ζ)^6=−27, same for −Q. Consequently the spectrum is contained in {0,−27,27}.

In the common order of inputs x and images f(x), the difference P between the corresponding principal blocks of S is elementwise nonnegative. The entry sum in either principal block has absolute value at most 27mn/N: the centered indicator of n points, or its centered complement, has squared norm mn/N, and S kills the constant vector. Hence the sum of the row degrees d(x)=Σ_y P(x,y) is at most 54mn/N. Since defects are nonnegative integers, d(x) bounds the number of disagreeing distances at x.

There exists an affine basis of seven domain points whose total d is at most (7/n)Σd. For precision, order the points by nondecreasing d and use the usual greedy basis of the vectors (x,1) in F₃⁷. A prefix of rank r<7 is contained in an affine flat with at most 3^(r−1)≤rn/7 points; for rank 7 this prefix inequality is automatic. Thus the greedy basis picks at least 7j/n points in each initial segment of j points. Summation by parts proves the claimed weight bound. In particular at least
`n − floor((14/27)(729−n)) ≥ 291`
domain points disagree with none of the seven anchors.

For these common-agreement points, subtract the Q-distance equation for the first anchor from each of the other six. This gives six equations linear in x and f(x), with constant terms. The difference vectors of the seven image anchors must be independent: a nontrivial dependence would give a nontrivial linear equation on x, cutting the agreement set down to at most 243. The linear system therefore expresses f(x) as a single affine map g(x), and its linear part U is invertible by the same independence for the original anchors. For each direction v≠0, the >243 common-agreement points occupy some one of the 243 parallel lines in direction v twice. The original monotonicity thus implies σ(Q(Uv))≤σ(Q(v)) for every v; since U permutes V, the two sums over V are equal, so every inequality is equality. It follows that g is an affine isometry.

Keep all a sites agreeing with this g, and let k=n−a. Normalize g to the identity by applying its inverse to the output. An outsider at (0,s), s≠0, excludes precisely the following diagonal positions:
`
 E_s = {x : x=0 or x=s or σ(Q(s−x))>σ(Q(x))}.
`
The same formula translated by a₀ applies to an outsider (a₀,a₀+t), t≠0. A pair of outsiders is compatible only if their input and output are both distinct and their signed distances obey the colour condition in §1. The accompanying finite checker proves that the 728 E_s have size 244 (504 cases) or 245 (224 cases), and that any two compatible outsiders exclude diagonal subsets whose intersection has size at most 135. Translation fixes the first input at 0; orthogonal symmetry, by Witt's extension theorem, reduces its displacement to e₁, e₁+e₂, or e₁+e₂+e₃. The checker uses no quotient files, search heuristic, or floating-point arithmetic.

If k≥2, the excluded diagonal union has size ≥353, so a≤376. Over the at most 729−a missing diagonal sites sum their multiplicities in the k exclusions. This sum is ≥244k; the sum of squares is ≤245k+135k(k−1). Cauchy gives
`
  (135a−38879) k ≤ 110(729−a).
`
But k≥440−a, and on 291≤a≤376 the concave quadratic (440−a)(135a−38879) has both endpoint values ≥60494, while the right side is ≤48180. Hence k≤1. In case k=1 its exclusion has size at least 244, giving n≤486. Choosing a displacement with 244 exclusions and taking every unexcluded diagonal site alongside this outsider attains size 486.

## 3. Interaction between distinct colours

Let two distinct colours contain subsets of affine graphs of sizes a,b. Normalize the first reference graph to the identity as above, and call the second affine isometry h. The displacement norm d(x)=Q(h(x)−x) is h-invariant: h(hx)−h(x)=U(h(x)−x). Common rows R have |R|≥a+b−729 and d=−1. Moreover h(R) is disjoint from the domain of the first graph, because otherwise the sites share a column, prescribing displacement norm +1, incompatible with −1. Thus R and h(R) are disjoint inside d=−1. Common columns C similarly have |C|≥a+b−729 and d=+1, and C and h⁻¹(C) are disjoint inside d=+1. Therefore 4(a+b−729)≤729 and a+b≤911.

In addition, a single affine subset of size a forces every site in another colour to use a missing row or missing column of that graph. Indeed, otherwise sharing its row and its column prescribes both signs to the same displacement after normalization. Since the other colour has no repeated row or column, its size is at most 2(729−a).

Call a colour large if its size is at least 440. For ℓ≥2 large colours their affine subsets omit at most one site each, and their pairwise reference sizes sum to at most 911. Thus their total size is at most floor(911ℓ/2)+ℓ. Other colours contribute at most 439 each. If there is just one large colour of size ≤486, use the 439 bound for other colours. If it has size n≥487, it is all affine and the bound for a total of j colours is
`
  n + (j−1) min(439,2(729−n)).
`
For j=3,4 this is largest at n=509, by checking the two affine pieces on either side of the breakpoint. Combining those cases gives H₃≤1387 and H₄≤1826 for three or four colours at mutual squared two-block separation ≥8/3. Neither bound is claimed sharp.

## 4. Real packing and exact comparison

Use two K12 period blocks and normal A_k, independently scaled by squared factors (x,y,z)≥1 starting with root A_k additionally scaled by √2, so each period block starts at minimum squared length 4. Pick any subset of the `V² × A_k/3A_k` phases; a normal class v has displacement v/3. The normal cost is (2z/9) times the minimum root-normalized squared norm in its class. Radius 1 gives centre density equal to motif count divided by period covolume
`
  729 (xy)^6 2^(k/2) √(k+1) z^(k/2).
`
The formula for each relative distance is independent in the period factors, so this tests all translates, not just representatives. The [earlier scale-vertex proof](https://github.com/alexlegeartis/SpherePackings/commit/415eaeb5c0b3cca957148471c82d6a18ad0228b8) leaves (1,1,3) for k=1 and 3 and (1,1,1), (1,1,3) for k=2 as the only unresolved candidate vertices for a strict improvement against the respective pinned 25D, 26D and 27D records. The following resolves precisely those vertices:

- **25D**, z=3: the three A₁/3A₁ layers have pair cost 4/3, so H₃ bounds the motif. Centre density ≤1387/(1458√3)<1/√2.
- **26D**, z=1: in A₂/3A₂, the line J=span(1,−1) has two nonzero classes of minimum norm 6 in simple-root coordinates; the six other nonzero classes have a norm-2 representative. Thus the costs are 4/3 and 4/9 respectively. Recolour all nine labels using a nonzero linear map F₃²→F₃ injective on J. Distinct labels in one new colour have normal cost 4/9 and therefore squared two-block separation ≥32/9, which rounds up to 4 because costs are multiples of 2/3. Across colours the normal cost is at most 4/3. The whole motif is preserved: repeated sites are impossible since no distinct labels have cost 4. Hence M≤H₃. For z=3, select one label per each of the three cosets of J, uniformly. Differences outside J cost 4/3 so these sampled layers are a valid three-colour instance, and their mean motif is M/3. Hence M≤3H₃. Both vertices have centre density ≤1387/(1458√3), strictly below Bacher's 1/√3 by the factor 1387/1458<0.952.
- **27D**, z=3: take the four A₃/3A₃ offsets 0,e₁−e₂,e₁−e₃,e₁−e₄, which have all six pair differences of root norm 2 and normal cost 4/3. Their 27 translates average 4M/27 sites and each is a valid four-colour instance. Thus M≤floor(27·1826/4)=12325, and centre density ≤12325/(8748√6)<1/√3<1/√2. The last comparison is only for this vertex. The excluded vertices of the original theorem were tested against the unrestricted 27D target, not against the best lattice-only comparator.

The pinned comparisons are the [live table's](https://github.com/alexlegeartis/SpherePackings/blob/main/LOWER_BOUNDS_1_96.md) 25D laminated 1/√2, 26D [Bacher](https://doi.org/10.1007/s002220050180) 1/√3, 27D Bacher lattice 1/√3 and [Vardy's nonlattice packing](https://doi.org/10.1090/S0002-9947-99-02169-8) 1/√2, all **centre** densities. To compare with the occupied-fraction table, multiply by the unit ball volume in that dimension.

## Constructive consequence

At the critical 26D vertices this is a substantial fixed-architecture gap, not a failed choice of shifts. The exact structural leverage is that one-sided signed-distance comparison on a large partial matching forces almost all of it into an affine reference, and row/column conflicts then make different large layers incompatible. To improve by this mechanism requires changing the continuous metric, horizontal alphabet, or carrier, not optimizing more vertices of the closed family.
