# Sphere 159 addendum: exact determinant floor for every rank-nine even integral normal

**Theorem.** Every positive-definite even integral lattice K of rank nine and squared minimum at least 6 has Gram determinant **at least 21870**. The bound is attained. Every equality case contains a copy of √3E₈. No hypothesis on the ambient lattice, placement of the normal, or existence of any root sublattice is made.

Consequently, a bare 39-dimensional section of **any** extremal even unimodular 48-dimensional lattice cannot beat the Dorofeev–Sun–Wang packing at the inherited radius √6/2. A primitive normal K and its complementary section B have equal Gram determinants by unimodularity and orthogonal projection of a basis. The centre density is (3/2)^(39/2)/√det(K), so this closes the full bare-section determinant problem at that radius. [Dorofeev–Sun–Wang, Theorem 14 and Appendix A.2](https://arxiv.org/html/2607.20359v4#S3.SS3) gives a literal K₉ of determinant 21870 in P48p, P48q and P48n; the same [live 39D benchmark](https://github.com/alexlegeartis/SpherePackings/blob/main/LOWER_BOUNDS_1_96.md) has occupied fraction 1.68028177365432e−7. A section with a larger **actual** minimum could be used at a larger radius; neither that possibility nor translated unions, gluing, or an arbitrary 39D packing is excluded. This result establishes a sharp upper for the stated 39D construction class, not a new density lower bound.

The proof has two independent short replay files: this addendum routes every potentially competitive K to √3D₄; the companion root-seed proof handles all K containing that sublattice and shows that equality needs √3E₈. The former proof is a finite rank-one-through-rank-six extension argument, not a 48D shell search or a search through perfect rank-nine forms.

## 1. Shortest normal vector, primitive lifts, and exact gate

Assume det(K)≤21870. The squared minimum must be 6. For context, [Dutour Sikirić and van Woerden, arXiv:2508.20719v2, Corollary 1](https://arxiv.org/html/2508.20719v2#S1) prove γ₉=2. The present proof **does not need that classification**: classical Mordell transference gives γ₉⁷≤γ₈⁸=2⁸ and this already implies λ₁(K)²<8. Indeed, for rank n choose a shortest nonzero, hence primitive, dual vector u. The hyperplane section of K orthogonal to u has Gram determinant det(K)|u|² and squared minimum at least λ₁(K)². Apply γ_(n−1) to that section and γ_n to K*, then take the supremum to get γ_n^{n−2}≤γ_(n−1)^{n−1}. The exact inequality in use is 8^63>2^72·21870^7. Even integrality and min≥6 give the claim.

For each encountered seed S with Gram G, let D=det(G) and compute the least squared length q/D in each discriminant class S*/S. The largest is <6 in every seed visited. Therefore S is automatically primitive in K: any extra vector in K∩span(S) is in S*, and subtracting its shortest representative's period would give a nonzero vector of norm <6. The orthogonal projection Q of K off S is a lattice of rank r=9−rank(S), with det(Q)=det(K)/D. If u is a shortest nonzero vector of Q, choose a lift and subtract a period of S so that its S-component is a shortest dual representative of its class. The resulting lift has even squared norm N≥6, and the projected squared length is |u|²=N−q/D.

For r=1,…,8, the exact Hermite determinant bound is det(Q)≥c_r (|u|²/2)^r with (c₁,…,c₈)=(2,3,4,4,4,3,2,1): these are the determinants at squared minimum 2 of A₁,A₂,A₃,D₄,D₅,E₆,E₇,E₈. See the [Nebe–Sloane catalogue and references](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/density.html). Thus a candidate class and lift must pass the exact integer gate
`
  c_r (N D−q)^r ≤ 21870 D^(r−1) 2^r.
`
For every seed, even with q set to its maximum, N=8 fails this gate. Consequently all actual short lifts have N=6. Their Gram is obtained by adding diagonal 6 and their integer pairing column m; its determinant is 6D−q. Separately the checker confirms min≥6 of every enumerated extension: for coefficients k≥2 of the new vector, the least norm of that layer is k²(6−q/D) plus the exact depth of class k[m]. There is nothing to check after the first term alone reaches 6. This consistency check does not introduce any dependency on a floating shortest-vector result.

## 2. Every potentially optimal K contains √3A₃

The appendix checker gives the following short decision tree. A₃ and D₄ below denote their √3-scalings.

1. **No A₂ is impossible.** Start with one norm-six root; the gate from its Gram [6] forces pairing magnitude 2 or 3 for the next norm-six lift. Pairing magnitude 3 would give √3A₂, so up to sign the two-root Gram is C₂ with diagonal 6 and off-diagonal 2. If A₂ is absent, every later norm-six lift must avoid pairings of absolute value 3 with existing roots. Exactly one eligible class over C₂ can meet this necessary condition. Its extension is C₃, with diagonal 6 and all off-diagonal entries −2. Its three basis vectors and their sum have norm six. For each of the 18 eligible classes over C₃, either some pairing column entry has absolute value 3 or their sum has absolute value 3. The associated norm-six lift makes A₂ with one of those four roots. This is a contradiction.
2. **A₂ but no A₃ is impossible.** Start instead with the scaled A₂ with Gram diagonal 6 and off-diagonal 3. Its eight eligible classes give only A₃ or a seed T₃ with determinant 120, diagonal 6, and off-diagonals (−3,−3,2). If A₃ is absent, use T₃. Its eligible classes give either an A₃ inside the extension or two rank-four seeds of determinants 400 and 448 (their respective extension columns are (−1,−2,−2), (−1,−1,3)). Both in turn either give an A₃ or the same rank-five isometry type, determinant 1280, obtained from the 400 seed by column (−3,2,2,−2). Its only continuations avoiding A₃ are isometric to its extension by (−3,1,1,−1,1), determinant 3072. This last seed has exact squared glue radius 9728/3072=19/6. The remaining rank-three projection thus has min≥17/6, giving determinant floor 3072·(17/6)^3/2>21870. This contradicts the assumption.

At every branch “gives A₃” means three literal norm-six vectors with two absolute pairings 3 and third pairing 0, hence a rank-three A₃ after sign changes. It does not assume that all norm-six vectors form one root system.

## 3. The last branch: A₃, but no D₄

In the D₃ coordinates for A₃, choose the three simple roots e₁−e₂, e₂−e₃, e₂+e₃, scaled by √3. The 19 eligible classes give D₄ or exactly the following three rank-four seed choices, where the displayed integer column is appended with diagonal 6:
`
  (0,0,−3) : determinant 405;
  (0,−3,2) : determinant 405;
  (0,−3,1) : determinant 432.
`
Run the same Hermite gate at each seed, retaining only extensions not containing D₄ and identifying isometric lattices. “Contains D₄” is likewise literal: a norm-six central root has three norm-six neighbours with absolute pairings 3 to it and three mutual pairings zero. The checker reports:

| new seed rank | eligible class lifts | extensions already containing D₄ | remaining Gram determinants, up to isometry |
|---:|---:|---:|---|
| 5 | 60 | 22 | 1320, 1350, 1350, 1428 |
| 6 | 85 | 56 | 3375, 3648, 3795, 3920 |
| 7 | 0 | 0 | none |

The final four rank-six seeds have largest squared glue depths equal to 12150/3375, 13184/3648, 13398/3795 and 14112/3920 in that order. Subtracting these from 6 and using the rank-three Hermite constant strictly exceeds 21870 in every case. Thus every potentially optimal rank-nine K contains √3D₄. The companion root-seed proof proves det(K)≥21870 in that case, with equality requiring √3E₈. No uniqueness beyond that containment is asserted.

## 4. Why the finite checks cover all possibilities

Run `node verify_global.cjs` for this note and `node verify_d4.cjs` for the companion root-seed note. No third-party packages or data are needed. The global program prints all surviving Gram matrices through rank six. Its three finite subroutines are as follows:

- **Complete discriminant atlas.** For a dual vector u of squared norm ≤6, all integer pairings m_i=⟨u,b_i⟩ with the norm-six seed basis lie in [−6,6]. For a prefix of these pairings, the squared norm of the corresponding orthogonal projection is m_prefixᵀ G_prefix⁻¹m_prefix. This cannot exceed the norm of u, so the prefix can be omitted only when that exact value exceeds 6. At each leaf the class label is adj(G)m modulo D, since the coefficient vector is G⁻¹m and two coefficient vectors differ by an integral vector exactly when these labels coincide. Retain the least numerator in each class. The observed class count is checked equal to D, which proves that every class occurs in the enumerated ball and its least norm was found. Ties between shortest representatives are immaterial: they differ by a seed period and replacing the extension generator by it preserves the extension as a lattice.
- **Literal root shape check.** For an integral coefficient vector z with zᵀGz≤6, Cauchy gives z_i² D≤6 adj(G)_ii. Within these certified finite intervals enumerate coefficients in reverse order; the squared norm of the orthogonal projection off earlier basis vectors is computed using the exact Schur numerator and is a lower bound for full norm. At a leaf test the full Gram and retain one sign per root. Root motifs use absolute inner products, so the sign convention discards none.
- **Isometry deduplication.** The program uses a list of norm-six vectors and searches for a basis with the exact inner products of the other Gram. It reports success only if the coefficient change-of-basis matrix has determinant ±1. All constructed bases have norm-six rows, so this also suffices to discover any isometry; but completeness of the isometry finder is not needed to validate a deletion, because every successful deletion has a literal unimodular change of basis. A failure to find one simply keeps both cases. In the earlier decision tree its success on every asserted match is tested explicitly.

All enumerated Gram ranks are at most six. Matrices entering a tested lift can have rank seven, but the rank-six gate leaves zero such lifts. Determinants, pairings, labels, and Schur numerators have integer magnitude safely below 2^53 for the displayed small Grams, and Hermite gate powers use BigInt. All decisions are exact. Neither proof uses a numerical solver, approximate Gram matching, probabilistic test, or cached full-dimensional shell.


For [uniqueness of the equality lattice](EQUALITY.md) and [the actual minimum of its 39D sections](SECTION_MINIMUM.md), see the two separate addenda.
