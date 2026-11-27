# Sphere 159: A rank-nine normal containing scaled D4 cannot improve the 39D incumbent

**Theorem.** Let L be an even integral positive-definite lattice of squared minimum at least 6. Every rank-nine lattice K contained in L that contains a copy of √3 D₄ has Gram determinant at least **21,870**. In case of equality K contains a copy of √3 E₈. This conclusion does not assume a particular root embedding, a parametrized K₉, an extremal ambient, or an exhaustive full-shell scan. It does not claim that the equality is attained only in one isometry class or one ambient embedding.

For the sphere-packing application take L to be any extremal even unimodular lattice of rank 48. If K is a primitive rational normal then its orthogonal section P has rank 39 and Gram determinant det(P)=det(K). This follows by taking a primitive lattice basis of the section and projecting a complementary basis: the projected lattice and the primitive normal are Euclidean dual, and volumes multiply. At inherited radius √6/2 the centre density is (3/2)^(39/2)/√det(K). The [Dorofeev–Sun–Wang construction, Theorem 14](https://arxiv.org/html/2607.20359v4#S3.SS3) is attained at det(K)=21870 in P48p, P48q and P48n; the [live project](https://github.com/alexlegeartis/SpherePackings/blob/main/LOWER_BOUNDS_1_96.md) has its occupied fraction 1.68028177365432e−7. In 39D a strict improvement at radius √6/2 requires det(K)<21870. The result therefore closes all bare sections at this radius whose normal contains √3D₄, not all 39D sections. In particular it does **not** control a section of actual minimum greater than 6 or a union of translated sections.

## 1. One gluing rule and four Hermite constants

For a seed S⊂K write R(S) for the **squared** glue radius, i.e. the largest squared norm of a shortest representative in S*/S. If R(S)<6, then S is primitive in K: any x∈K∩span(S) belongs to S*, and subtracting a shortest representative's period gives a contradiction to min(K)≥6 unless x∈S. Project K onto span(S)^⊥; it is a lattice Q of rank 9−rank(S), and det(K)=det(S)det(Q), by completing an S basis and taking its orthogonal Schur complement. For any nonzero projected point u, its lift's S-coordinate is in S*. If its class has least squared depth d, subtract a seed vector attaining d; this lift has squared length d+|u|²≥6. Since the ambient L is even, d+|u|² is an even integer. In each short-class split, any lift of length at least 8 is automatically outside the displayed cutoff.

For rank r in {2,3,4,5}, the classical Hermite values are det(Q)≥c_r(μ/2)^r for μ=min(Q), with c_r=3,4,4,4 respectively. They follow by scaling the exact minima-two constants of A₂, A₃, D₄, D₅; see the [Nebe–Sloane catalogue and references](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/density.html). Every comparison below is an exact rational inequality, repeated by the checker.

If a necessary short projection occurs, its shortest seed lift has length 6. Fixing the class and the representative fixes the cross inner products with S. Reversing the new vector and applying seed isometries are allowed. Thus each listed transition means the actual K contains the stated extension, not that every abstract Gram embeds in every possible ambient.

## 2. The root seeds

All lattices named in this section are scaled by √3. Exact upper glue levels and their denominators are:

| Seed | Gram determinant | largest squared depths (numerator / denominator) |
|---|---:|---|
| D₄ | 324 | 36, 28 / 12 |
| D₅ | 972 | 45, 37, 36, then ≤32 / 12 |
| D₆ | 2916 | 54, 46, 38 / 12 |
| E₆ | 2187 | 144, 112 / 36 |
| E₇ | 4374 | 54, 38 / 12 |
| E₈ | 6561 | 8/3 (already proved in [DSW, Table 3](https://arxiv.org/html/2607.20359v4#S3.T3)) |

For D_d the physical dual point is x/√3, where x∈D_d* has all coordinates integral or all half-integral. A representative with squared depth q/12 has integer coordinates a=2x and q=Σaᵢ². The checker verifies saturation from the exact maxima of all these tables. For E₆ and E₇ see §4.

Assume for contradiction det(K)<21870, or track a configuration achieving equality.

- **D₄ → D₅.** If the rank-five projection has squared minimum ≥11/3, the Hermite determinant bound already gives det(K)>21870. The only shorter value is 3, from a norm-six lift of a 36/12 class. There are three such classes. The one with x=3e₁ and the two with all |xᵢ|=3/2 are equivalent under signed permutations and D₄ triality: the normalized 4×4 Hadamard matrix with first column all 1/2 is an isometry of D₄ and maps the first point to one half-class; changing a sign gives the second. These are all three (checked). The representative 3e₁ allows the new root √3(e₁+e₅), giving D₅.
- **D₅ → E₆, D₆, or T₆.** With rank-four projection min≥10/3 the determinant is already ≥30,000. Shorter minima are 9/4, 35/12 or 3 from the three top depths. A representative of 45/12 has all |xᵢ|=3/2. The new generator can be matched to the half-integral E₆ root in the coordinate model below, giving E₆. The 36/12 case has x=3eᵢ and extends to D₆. The 37/12 case has one |xᵢ|=1/2 and four |xᵢ|=3/2, so signed permutations give x=s=(1,3,3,3,3)/2. This yields T₆.
- **D₆ → E₇ or T₇.** Squared minimum ≥17/6 in the rank-three projection would already give determinant >21870. The two shortest possible values 3/2 and 13/6 give, respectively, a √3E₇ half-root or T₇. The representative for 46/12 has one half-coordinate of magnitude 1/2 and five of magnitude 3/2; all such representatives are signed-permutation equivalents. The explicit identification with T₇ is in §3. 
- **E₆ → E₇.** The unique deepest pair of mutually negative E₆ classes has depth 4 and requires norm-six lift to get projected norm 2. Otherwise min≥26/9, and rank-three Hermite already gives det(K)≥26364. The norm-six lift coincides with the standard root extending E₆ to E₇, as checked below.
- **E₇ → E₈.** Exactly one depth-9/2 class exists; a norm-six lift gives projected norm 3/2, and in an even ambient every other value is ≥17/6. The rank-two bound then exceeds 21870. The norm-six lift gives E₈.
- **E₈.** Any nonzero orthogonal extension has min ≥6−8/3=10/3, so the rank-one Schur step gives det(K)≥6561·10/3=21870. This is the already published last step.

## 3. The non-root branch from D₅

In coordinates with √3D₅ in the first five dimensions let the first new vector of T₆ have horizontal component s/√3, s=(1,3,3,3,3)/2, and h₁²=35/12. On the D₅ simple basis e₁−e₂,…,e₄−e₅,e₄+e₅ its cross column is (−1,0,0,0,3) and its norm is 6. T₆ has determinant 2835.

**T₆ → E₆ or T₇.** Its largest squared depth is 1584/420=132/35; the next is at most 1376/420=344/105. If no norm-six lift of a deepest class occurs, the rank-three projected minimum is at least 286/105, which gives det(K)>21870. There are four deepest classes, with representatives described by the doubled first five coordinates a and the integer j₁ pairing with the new root:
`
 (a,j₁) = ((−6,0,0,0,0),−2), ((−6,0,0,0,0),1),
           ((−4,0,0,0,0),−3), ((2,0,0,0,0),−3).
`
A further vector v₂ of norm six has residual length squared 78/35 off T₆. Its projection off D₅ has inner product w/12 against the projected v₁, where w=12j₁−a·(1,3,3,3,3). In the last two cases w=−32,−38, and the squared length of the projection of v₂+v₁ off D₅ is
`
 78/35 + w²/420 + 35/12 + w/6 = 9/4.
`
A 9/4 projection off D₅ yields E₆ by §2. The first two cases are opposite modulo T₆: summing their dual representatives adds the period √3(−2e₁). Hence use the first representative, which appends the cross column (−3,0,0,0,0,−2) and diagonal 6. The extension T₇ has determinant 6318.

This T₇ also contains √3D₆: v₂ projects off D₅ to length squared 3 and has horizontal part −√3e₁, so D₅ plus v₂ can be represented using v₂=√3(−e₁+e₆). The other vector v₁ then has √3D₆ horizontal dual coordinate (1,3,3,3,3,−3)/2, as its pairing with v₂ is −2. This has exactly one half-coordinate of magnitude 1/2 and five of magnitude 3/2, and its residual squared norm off D₆ is 13/6. Thus the reverse T₇ identification closes the D₆ transition in §2.

**T₇ → T₈.** The largest squared depth is 44100/10920=105/26; all others except the two deepest are ≤41860/10920=23/6. Without a norm-six lift of a deepest class the rank-two projection has min≥13/6, so det(K)≥6318·(3/4)·(13/6)^2>21870. The two deepest classes have representatives (a,j₁,j₂)=((−6,0,0,0,0),−2,3) and ((−6,0,0,0,0),1,3). They are opposite modulo T₇ because they sum to √3(−2e₁), whose j₂=6. Use the first to append cross column (−3,0,0,0,0,−2,3) and diagonal 6, obtaining T₈. Its Gram determinant is 12393.

**T₈ is too large.** Its largest squared depth is 633360/185640=58/17. Therefore the residual rank-one min is at least 44/17 and det(K)≥12393·44/17=32076. This ends the non-root branch.

## 4. Exact finite certificate and identification of E₆/E₇

Run `node verify_d4.cjs` using the accompanying standard-library file. It enumerates dual vectors with squared norm at most 6 and assigns them to discriminant classes by `adj(G) m mod det(G)`, where m is the column of pairings with a lattice basis. Thus each observed exact depth is the least norm in its class: all shorter vectors are included because the searched box contains the radius-√6 ball. The number of classes with a representative equals the Gram determinant, so there are no unobserved classes. For types D and T, first five (or the appropriate number of D coordinates) are a_i/(2√3), and each j_i pairs with a norm-six vector. Consequently |a_i|≤√72<9, |j_i|≤6. All squared-norm formulas are integers after multiplication by their denominator, and the checker separately matches them to mᵀadj(G)m/det(G) at each enumerated point.

For the explicit T-coordinate formulas, put w=12j₁−a·(1,3,3,3,3), z=70j₂+35a₁+3w and y=52(j₃−j₂)+z. The successive orthogonal coordinates are w/(12h₁), z/(70h₂) and y/(52h₃), with h₁²=35/12, h₂²=78/35 and h₃²=51/26. Adding the successive squares gives the three integer numerator formulas in the replay. The orthogonal coordinates follow by taking the cross columns of §3 in order, or by directly computing their scalar products with the basis.

For E models use E₈ in `D₈ ∪ (D₈+(1/2)^8)`, with simple roots a₁=(1,−1,−1,−1,−1,−1,−1,1)/2, a₂=e₁+e₂, a₃=−e₁+e₂ and aᵢ=−e_(i−2)+e_(i−1) for i=4,…,8. The leading six or seven roots span E₆ and E₇. Their primal lattices are primitive in E₈ (subsets of a full root basis), so the projections of E₈ give their entire Euclidean duals. Here write x for the unscaled *dual* component, so the true point of (√3E)^* is x/√3.

- For E₆, x=(a₁′,…,a₅′, w/3,w/3,−w/3)/2 with common parity p of all aᵢ′ and w, and Σ(aᵢ′−p)/2+(w−p)/2 even. This is the projection of the integral and half-integral E₈ coordinates onto the hyperplane factor spanned by (1,1,−1). Classes in E₆*/3E₆ are exactly encoded by the five aᵢ′ modulo 3 and w modulo 9: if all vanish, dividing by 3 preserves the parity constraint, and makes the final three projected coordinates share the E₈ type. The numerator is 3Σ(aᵢ′)²+w² in units 1/36 after scaling. Complete enumeration through 144 reaches all 2187 classes. The deepest classes have all five residues zero and w≡3,6 mod9. The horizontal component of √3a₇ has x=3 proj(a₇), namely −3e₅+(e₆+e₇−e₈), which gives one of these classes, and the vertical squared length is 2 when scaled. For D₅→E₆ the simple roots a₂,…,a₆ are the usual D₅ on the first five coordinates, and √3a₁ has that projection of depth 15/4. For D₆→E₇ the same applies using a₂,…,a₇ on the first six coordinates, with squared depth 9/2 from √3a₁.
- For E₇, x=(a₁′,…,a₆′,w,−w)/2 with first six coordinates' numerators of common parity p, and Σ(aᵢ′−p)/2+w even. Projecting the last two E₈ coordinates gives this model. The class is the seven numerators modulo 3 and the bit w−p mod2. Vanishing residues allow division by 3 within E₇*, and vanishing bit makes the result a primal E₇ vector. The depth numerator is Σ(aᵢ′)²+2w² in units 1/12 after scaling. Complete enumeration through 54 reaches all 4374 classes. The unique deepest class has zero ternary residues and bit 1, represented by the x=3 proj(a₈)=−3e₆+(3/2)(e₇−e₈) projection of the standard next simple root. Its orthogonal squared length when scaled is 3/2.

Thus, after either norm-six exceptional lift, the known root extension is forced by the Gram. All gates away from √3E₈ are **strictly** worse than 21,870; tying requires the E₈ step. 
