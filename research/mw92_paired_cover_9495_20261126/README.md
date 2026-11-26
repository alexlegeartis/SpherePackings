# Equal-ball covering repairs for MW92 lattice packings in 94D and 95D

A 94-dimensional **lattice** packing has occupied density at least **8.078244769088839e-21**. This raises the previous certified 94D lattice and unrestricted lower bound **8.036027386567117e-21** by **more than 0.52535%**. Both are existential lower bounds, not explicit 94D Gram matrices. The derived external comparison is Shioda's MW92 times a minimum-matched triangular lattice, with squared tail covolume 192; the improvement against that construction is more than 218.7546%. We do not attribute that derived 94D product to Shioda as a published numerical record or claim global optimality.

The same geometric replacement, using small triangular and square clusters, also raises the verified 95D lattice and unrestricted bound from **2.115959060080788e-21** to **2.128470378304569e-21** (**+0.591283%**, rounded down). This is a new improvement over the earlier result *of the same date*. It avoids the earlier 95D split by small torsion orders.

## Parent, literature, and competing constructions

Use [Shioda’s characteristic-47 construction (1991), Theorem 1.1](https://doi.org/10.2307/2374791), the even integral lattice L = MW92 with minimum squared norm 16 and covolume 47^7, as in [Schütt–Shioda, Theorem 13.35 and §13.4.2](https://link.springer.com/book/10.1007/978-981-32-9301-4) and [the project's parent audit](../mw92_extensions_20261122/). Its integral orthogonal automorphism J satisfies J²+J+I=0. Let d be the distance to L on its flat torus. The [degree-13 Dirichlet estimate](../mw92_cfr_taylor13_20261123/) for the **equal-ball** [Coxeter–Few–Rogers theorem](https://doi.org/10.1112/S0025579300002059) says that the volume of a 92-ball of squared radius A=11.00098526, divided by 47^7, is strictly below the necessary density of a Euclidean equal-ball covering. We reproduce the scalar estimates in the checker.

Three competing mechanisms make the distinction sharp. The minimum-matched product uses tail determinant 192. The published-on-main [three-full-coset covering](https://github.com/alexlegeartis/SpherePackings/blob/b90e0db1b446bb7a2344c20e8514e3ab2bbcbd7e/research/mw92_cfr_taylor13_20261123/result.json) has worst tail determinant 859313833/45000000. The construction below preserves the good generic Eisenstein shear, but replaces the limiting three-coset covering by a paired covering of only the two strong cosets and one weak coset. It uses precisely the same prior covering theorem, not an unproved bound for unequal balls.

## Exhaustive dichotomy

For completeness recall [the order-three argument](../mw92_order3_dichotomy_20261122/). Choose a with d(a)²>A, put s=16−A, t=16−3s, T=I−J. If d(Ta)²≥t, use shears a,−Ja with triangular vertical norm s(i²+ij+j²). The six norm-one Eisenstein labels have horizontal squared distance >A, and all norm-three labels have distance squared ≥t. Every other nonzero Eisenstein norm is at least 4, so this is a minimum-16 packing with tail determinant 3s²/4.

Otherwise let h be T⁻¹ applied to a shortest representative of Ta modulo L; let u=a−h. Then Tu∈L and ||h||²<t/3 because TᵀT=3I. The orbit of a closest lattice point to u under the affine lattice isometry v↦Jv+Tu gives three relative vectors from u of identical length summing to 0. Averaging squared lengths after translating by h yields d(a)²≤d(u)²+||h||². Thus d(u)²>A−t/3=32/3. Its values lie in (2/3)ℤ because 3||u−v||²=||T(u−v)||² is an even integer for every v∈L. Consequently d(u)²≥34/3. Since (I−J²)T=3I and u is nonzero modulo L, u has exact order three on the torus.

## A weak phase costs only half a full ball

Set p=14/3, D=10.78402070, W=D−2p, R=√D and r=√W. We need a b with d(b)²>D, d(b−u)²>D, and d(b+u)²>W. Its nonexistence would cover Euclidean space by radius-R balls centered on L and L+u, together with radius-r balls centered on L−u.

Let v be any minimum vector of L (length 4). It is primitive, so there is a homomorphism L→ℤ/2 taking v to 1. Its kernel H has index two. Partition the weak centers into pairs h−u, h+v−u, h∈H. One radius-R ball at h+v/2−u covers both radius-r balls by the triangle inequality because 2+r<R. This constructs a **periodic equal-ball covering** with at most five balls per H-cell: four from the original two strong cosets and one midpoint. Its density is at most (5/2)·vol(B₉₂(R))/47^7. But (5/2)D⁴⁶<A⁴⁶. This contradicts the strict equal-ball bound above. Thus the required b exists.

## Every tail vector

Use vertical Gram [[p,p/2],[p/2,q]] with q=16−D, with horizontal differences i u+j b for (i,j)∈ℤ². Completing the square gives
    Q(i,j)=p(i+j/2)²+(q−p/4)j².
For |j|≥2, Q≥4(q−p/4)>16. For j=0 the only nonzero possibly short differences have i=±1, and p+d(u)²≥16. Negation reduces |j|=1 to j=1. For i=0,−1 the vertical contribution is q, and the first two strong inequalities give a total >16. For i=1,−2 it is q+2p; both require the *same* third phase because 3u∈L, and the weak inequality gives a total >16. For every remaining i, Q≥q+6p>16. All horizontal fiber vectors are controlled by the 92D minimum 16, which is attained. Thus the combined lattice has squared minimum precisely 16 and squared covolume 47^14·p(q−p/4), independently of the (existential) b.

The repaired fixed-branch determinant exceeds the generic one but is strictly smaller than its predecessor. The density follows from vol(B₉₄(1))·2⁹⁴/(47⁷√det). The verifier uses exact rational arithmetic and rigorous Machin bounds for π; no finite floating-point search, large enumeration, or unquantified density comparison is part of the proof. Run `python3 -B verify.py` (assertions enabled).

## 95D: clusters of three and four parent points

Here are two exact periodic replacements of small balls used below.

* **Triangles.** If v is any minimum vector, 0,v,−Jv is equilateral of side 4. In L/3L the class of v is nonzero by primitivity. There is a linear functional χ taking v and Jv both to 1: if these classes are dependent, the scalar relating them cubes to 1 in the field ℤ/3 and hence is 1; if they are independent, extend the specified functional to a basis. Translates of the triangle by ker χ partition L, since its vertices have distinct χ-values. Every translated radius-r ball at one of its vertices lies in the ball about the centroid of radius 4/√3+r.
* **Squares.** There are orthogonal minimum vectors v,w in this particular parent. In the model y²=x³+1+t⁴⁸ take P=(−1,t²⁴) and Q=(−t¹⁶,1). They satisfy the equation and both have height 16: in the finite chart their coordinates are integral; at infinity, using s=1/t, x'=s¹⁶x, y'=s²⁴y, they extend with finite coordinates, so neither section meets O. The [height and irreducible-fiber proof for this exact model](../mw92_cusp_d4_20261122/README.md#why-specializing-to-zero-gains-two-in-the-horizontal-norm) gives ||P||²=16+2(P·O), and the same for Q. Alternatively, the base inversion t↦1/t interchanges P,Q after the displayed change of coordinates. The base isometry t↦ωt, with ω a primitive cube root in the constant field, fixes P and sends Q to JQ since 16≡1 mod 3. Therefore the three numbers ⟨P,Q⟩, ⟨P,JQ⟩, ⟨P,J²Q⟩ are equal, and they sum to zero. This proves orthogonality. Their classes in L/2L are independent: each is primitive, and (P+Q)/2 would have squared length 8, below the minimum 16. Choose a linear map to (ℤ/2)² taking them to the two basis vectors. Its kernel partitions L into translates of the four corners 0,P,Q,P+Q, with circumradius 2√2.
  
Both partition constructions work on any translate of L. Different cosets may use different types of clusters: intersect their period sublattices for a common finite period. If n_t ordinary phases use clusters of size t, the resulting equal-ball cover has at most Σ n_t/t centers per L-cell in this common period.

### Generic 95D branch (no torsion subcases)

Use the generic a from the 94D dichotomy. Let D8=10.51475582. The first two shears are a and −Ja; the third is c. Take

    [ 5,    5/2, −5/2  ]
    [ 5/2,  5,   −5/2  ]
    [ −5/2,−5/2, 16−D8 ].

The third-coordinate-zero section is protected by the norm-one and norm-three assumptions above (increasing the triangular scale to 5 causes no loss of minimum). For every nonzero third coordinate orient k>0. The verifier exhausts Q(i,j,k)<16 using the exact inverse-Gram bound z_l²≤16(G⁻¹)_ll. Only k=1 occurs. The six near phase labels (i,j) are
    (0,0),(1,0),(0,1),(−1,1),(1,−1),(1,1)
and each requires squared horizontal separation at most D8. The six outer phase labels are
    (0,−1),(−1,0),(−1,2),(0,2),(2,−1),(2,0);
each requires only D8−10<0.75². The center for each forbidden phase is −(i−jJ)a modulo L. One large ball of radius √D8 covers three outer-phase balls after triangle partition, since 4/√3+0.75<√D8. A failure to choose c would give an equal-ball covering with at most 6+6/3=8 centers per L-cell, impossible since 8·D8⁴⁶<A⁴⁶. This replaces the former [generic six-, seven-, nine-residue repairs](../mw92_eisenstein_repair95_20261126/) outright; we do not assume their conclusion in this verification.

### Fixed-order 95D branch

Retain the strengthened 94D fixed-order lattice with first two shears u,b, p and q as above. Set D*=10.52951849 and z=16−D*. Adjoin c with tail Gram

    [ p,   p/2,  p/2 ]
    [ p/2, q,    q/2 ]
    [ p/2, q/2,  z   ].

Again the third-coordinate-zero part is the 94D lattice, and every other oriented short vector (i,j,k) has k=1. The forbidden center for c is −iu−jb modulo L. Its phase (i mod 3,j) is well-defined since 3u∈L. Upon keeping the maximum demand per phase, exactly six phases need radii at most √D*, one phase (1,0) needs at most 1.1, and four phases (0,1),(1,1),(0,−2),(1,−2) need at most 0.35. The checker proves these inequalities, including that no further layer can be short. Handle the six large phases directly, the radius-1.1 phase with dimers, and the four radius-0.35 phases with squares: 2+1.1<√D* and 2√2+0.35<√D*. The resulting large-ball count is at most 6+1/2+4/4=15/2 per L-cell. Since (15/2)(D*)⁴⁶<A⁴⁶, a valid c exists.

The 95D squared tail covolumes are 572786627/8000000 (generic) and 644703354722796061/9000000000000000 (fixed). The second is larger, but it is strictly smaller than the previous [live 95D comparator 652349944595359843/9000000000000000](https://github.com/alexlegeartis/SpherePackings/blob/b90e0db1b446bb7a2344c20e8514e3ab2bbcbd7e/research/mw92_eisenstein_repair95_20261126/result.json). The new occupied density floor is **2.128470378304569e-21**. The gain over the derived minimum-matched MW92 × A3 product (tail determinant 2048) is **more than 434.6951%**. As in 94D the foliation and shear construction is a single lattice of squared minimum 16; the real horizontal offsets are certified to exist rather than computed explicitly.

