# Sphere 159 follow-through: every equality 39D section really has minimum 6

Let L be any extremal even unimodular 48-dimensional lattice. Let K be a primitive rank-nine normal containing S≅√3E₈, and let B=L∩span(K)^⊥. Write h²=det(K)/6561 for the squared length of a primitive generator of K's projection off S. Then the number of norm-six vectors in B is at least
`
  1,092,000 · (1 − 3h²/20).
`
In particular, **if det(K)<43740**, then B has actual squared minimum 6. At det(K)=21870 every section has at least 546,000 norm-six vectors (both signs included). No assumptions on the ambient embedding are needed. This does not assert equality of kissing numbers and does not control higher-minimum normals with det(K)≥43740 or no √3E₈.

## A one-line conditional moment, made exact

Let X be the 52,416,000 norm-six vectors of L, a spherical 11-design by Venkov; see [Venkov, Even unimodular extremal lattices](https://www.mathnet.ru/php/archive.phtml?jrnid=tm&option_lang=eng&paperid=2270&wshow=paper) and the [DSW source's discussion of the shell](https://arxiv.org/html/2607.20359v4#S5). On the radius-√6 sphere set a(x)=|proj_S x|² and t(x)=⟨x,e⟩² for any unit e perpendicular to S. For x∈X not in S, the E₈ dual coordinate y=√3 proj_S x lies in E₈ and obeys |⟨y,r⟩|≤3 for every E₈ root r, by comparing x with both signs of √3r. This forces |y|²≤8: the integer-root pairs first give |y_i|+|y_j|≤3. For an integral y, a largest coordinate of size 3 is ruled out by parity, while if the largest is 2 there are at most four remaining unit entries, by the all-half root inequalities and parity. If all coordinates are at most 1, the assertion is immediate. For a half-integral y, a largest coordinate of size 5/2 forces the other seven to have size 1/2. Otherwise all magnitudes are 1/2 or 3/2; four or more coordinates of size 3/2 are impossible by the all-half root with matching signs adjusted, if necessary, at a 1/2 coordinate (or at any coordinate when all are 3/2). Thus a(x) belongs to {0,2/3,4/3,2,8/3} for x∉S. For the 240 roots x∈S instead a=6 and t=0.

Use the degree-four polynomial F(a)=(1−3a/2)(1−3a/4)(1−a/2)(1−3a/8), which has F(0)=1, F(6)=70 and vanishes at four nonzero small levels. On the sphere the conditional mean of t at fixed a is (6−a)/40. Since F(a)t and F(a)(6−a) both have coordinate degree 10, the design gives the exact identity
`
  Σ_(x∈X,a(x)=0) t(x)
      = Σ_(x∈X) F(a(x)) t(x)
      = (1/40) Σ_(x∈X) F(a(x)) (6−a(x))
      = (6/40) |{x∈X : a(x)=0}|.
`
For the count itself, the spherical moments are E[a^j]=6^j(4)_j/(24)_j; substitution gives E[F(a)]=11/520. Thus the a=0 fibre is 52,416,000·11/520 − 240·70 = 1,092,000, in agreement with DSW's explicit counts. This supplies a proof for every embedding, rather than a complete 48D shell enumeration.

Finally choose e to be the direction of the projection of a rank-one lift v∈K. Every x in the a=0 fibre pairs with v by an integer k, so t(x)=k²/h². Among that fibre, those not orthogonal to K contribute at least 1/h² each, proving the stated lower bound. For a 39D section with larger true radius, only those outside this determinant window (or outside the scaled E₈ family) remain possible.
