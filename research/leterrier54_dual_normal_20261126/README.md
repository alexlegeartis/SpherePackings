# The scaled dual of Leterrier’s 54-dimensional lattice and consequences for 53-dimensional sections

No density record is claimed. This note determines the minimum and first-shell multiplicity of a specific scaled dual from the known first shell on the other side, and gives a constructive obstruction to a proposed new 53D packing. It does not rule out sections at higher minimum, other multilayers, translations, or other parents.

## Inputs

Let L be Leterrier’s 54D lattice with integral even Gram G, squared minimum 10, determinant 3^25, and exactly 15,309,000 vectors of norm 10. Here and below “norm” for a lattice vector denotes squared length. These are inputs from G. Leterrier, [Some Mordell–Weil lattices and applications to sphere packings](https://infoscience.epfl.ch/record/307330/files/EPFL_TH10417.pdf): Computational Theorem D (3.3.1) supplies the contact count and Computational Proposition 3.3.10 the full basis Gram. The [matrix transcription](../adaptive_layers_20261115/sources/leterrier54_gram.txt) has SHA-256 74a193d77f6faa0e94c79fa2795e57bcd7778d31b5a3840e276fc27acacde0f3. Direct integer multiplication verifies that H=3G^{-1} is integral with even diagonal. It is the Gram of J=√3 L*, not of the unscaled dual. Thus each of L and J has level dividing 3. Since det G=3^25, det H=3^29.

## A complementary coefficient identity fixes the entire short dual spectrum

Use q=e^{2πiz} and the theta convention θ_L=Σ q^{||x||²/2}. θ_L and θ_J have weight 27 on Γ₀(3) with quadratic character χ_{−3}. Let C=η(z)^9/η(3z)^3 and D=27η(3z)^9/η(z)^3. The eta-quotient criterion gives holomorphic weight-three forms of character χ_{−3}, their two cusp orders being (0,1) and (1,0). The monomials C^{9−j}D^j (0≤j≤9) span the weight-27 space: they have distinct leading powers 0 through 9, while a nonzero form cannot vanish to order 10 at infinity by the valence formula (index 4, maximal order 27·4/12=9). Write c_j for the resulting coefficients. The eta transformation interchanges C and D after removing the factor 3^{3/2}(z/i)^3; Poisson gives θ_L(−1/(3z))=3^{27−25/2}(z/i)^27 θ_J(z). Consequently c_j(J)=c_{9−j}(L)/3.

Write A_i for the number of vectors of norm 2i in J. Direct triangular power-series expansion to degree five gives
```
27^5 c5(L) = 398034 + 15309000 = 15707034,
27^4 c4(J) = 193914 + 16016 A1 + 1026 A2 + 45 A3 + A4.
```
By the reverse-coefficient identity, the first number is 81 times the second. Since 15707034=81·193914 and every A_i is nonnegative, all A1 through A4 vanish, not merely their signed sum. In the opposite direction,
```
27^5 c5(J) = 9·(27^4 c4(L)) = 9·193914 = 1745226,
A5 = 1745226 − 398034 = 1347192.
```
Therefore J has minimum exactly 10 with 1,347,192 minimal vectors; L* has minimum exactly 10/3. The supplied 54 integers separately exhibit a primitive shortest dual normal.

## Which 53D constructions this actually controls

For any primitive (unscaled) dual normal of squared length q, the perpendicular section has squared covolume 3^25 q. With the **inherited** packing radius √10/2 its center density is at most 20,890.80489885…, with equality exactly for q=10/3. The current independent 53D packing is the [P48n quaternionic lattice](../p48n_hurwitz_20261122/README.md), center density (3/2)^{53/2}/2=23,194.67914858…. A normal of norm 8/3 would have beaten it by 0.6982%, which is why we tried this construction; the theorem excludes it and every shorter candidate. Longer normals need not retain the original minimum.

For any harmonic quadratic polynomial, its shell-weighted theta series on L has weight 29 and order at least five at each cusp (the other cusp is the weighted theta of J, up to nonzero scale). The valence total is 29/3<10, so this series vanishes identically. In particular L’s shortest shell is a spherical 2-design. Its mean squared *integer* height with respect to a primitive dual normal is 10q/54. Every q<54/10 section therefore actually contains norm-10 vectors. In particular **every** shortest normal has at least 15,309,000·(1−50/81)=5,859,000 shortest vectors in its 53D section. This proves the exact inherited-minimum maximum is genuinely attained but does not exclude high-minimum sections at arbitrary large q.

## Several full sheets on an actual shortest normal

Let a be the supplied primitive integer label covector with aᵀHa=10, w=Ha primitive in L, and aᵀw=10. Orthogonal projections of L onto w-perp, modulo the bare integral section, have labels t=aᵀx mod 10. Every label occurs. The squared projected length is xᵀGx−3(aᵀx)^2/10. The supplied integer vectors for t=0,…,5 and their negatives bound the shortest nonzero projected differences, in tenths, by
```
100, 97, 88, 73, 72, 65, 72, 73, 88, 97.
```
These are **upper witnesses** for distances, not a claim of exhaustive closest-vector enumeration. A motif of M labels admitting squared minimum μ has relative center density M(μ/10)^{53/2} against the bare section. A direct exact check of all 512 subsets containing label zero shows that every nontrivial motif gives a ratio strictly less than one; the best multiple-sheet ceiling is 2·(97/100)^{53/2}=0.892238489…. It is attained by two adjacent labels: in the nearest lift of label ±1 one can take exactly that integer height by adding w, so inherited parent minimum gives μ≥9.7. Thus the one-sheet packing maximizes this specific ten-fiber model at every radius. This second assertion concerns **one** normal only.

## Construction screen outside this section (separate)

Starting at [Dorofeev–Sun–Wang’s six-sheet 43D P48m cross-section](https://github.com/wcgbg/sphere-packing/blob/main/packings_out/records/dim43-delta226/P48m-1f02996f.cert.txt), we also built small ordinary dual clusters after releasing one or two normal directions. With the published determinant-3072 normal K and B=48K^{-1} as in the verifier, release u=(0,0,0,0,1) in its dual: uᵀBu=12, the remaining rank-four primitive normal has determinant 768. Five projected classes (including zero) have squared diameter 7/24, hence yield 44D center density 5(137/96)^22/√768=450.99813389…. Release also v=(0,0,1,−1,0), with vᵀBv=12 and uᵀBv=0. The covector pair is primitive (unit minor); the residual rank-three normal has determinant 192. Four projected classes have squared diameter 1/4, giving 45D center density 4(23/16)^{45/2}/√192=1,015.28057628…. Orthogonal projection of K* gives the dual of the residual primitive normal: the basis covectors used need no extension assumed by dimension counting. Distinct chosen classes are distinct modulo a minimum-six lattice since the ordinary diameter is <6. Sources: [Chen–Hu–Li–Wang–Wu on the antipode construction](https://arxiv.org/abs/2505.02394), [existing density table](../../LOWER_BOUNDS_1_96.md). The actual comparators are 509.61898264… and 1,243.45967883… respectively. These are working constructions, not improvements or proven best projected clusters. Re-running these same packets does not bridge the gap. A subsequent [sharp class-wide analysis](../antipode44_five_sheet_20261126/README.md) rules out gains from at most five discriminant classes under the ordinary 44D normal-diameter certificate across all even min-6 unimodular 48D parents. The high-return next tests are genuinely six-sheet packets (for determinant 768 and diameter at most 7/24 this would give +6.2%) or a stronger certificate of the actual horizontal distances. The present packet does not establish either.