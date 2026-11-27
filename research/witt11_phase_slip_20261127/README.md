# Two exact tests of real, not just quantized, packing constructions

The two tests in this note are restricted obstructions. Their constructive escape is in [a separate note](PHASE_SLIP_11D.md): a rare phase slip at slightly reduced radius gives a strict 11D improvement. Here we study replacement at the exact old radius, or combining full Coxeter–Todd fibres with dense 9D/10D normal geometry. The attributed unrestricted occupied-volume benchmarks are 0.0662380270098011 (11D), 0.00292982910671570 (21D), and 0.00245103404412118 (22D); formulas in the proof are centre densities.

## 1. A movable 11D coset cannot be replaced by a denser packing over its fixed core

Use the [Leech–Sloane construction, §2.6 (1971)](https://neilsloane.com/doc/Me14.pdf): the 132 weight-six blocks of the Witt design are obtained from the xor of pairs of rows of the normalized order-12 Hadamard matrix and their complements. Add the six pairs of any perfect matching and their complements. For the explicit model here use the matching `(0,1), (2,3), …, (10,11)`, the zero row and 11 cyclic shifts of `11011100010` as the low 11 bits of the Hadamard rows, and shorten coordinate 11 to zero (indices start at 0). This produces 72 binary words in dimension 11, at minimum Hamming distance 4. Full translates of `2 Z^11` at these residues give the historical `P11a` centre-density benchmark `72/2^11 = 9/256`. We fix this explicit shortening; the proof does not assert a classification of all optimal binary codes.

Let `u=11111111111` and remove the word `r=u−e_10`, whose corresponding decimal bit mask is 1023. The remaining 71 complete cosets are called the core. **At squared minimum 4, any replacement packing retaining this fixed core has centre density at most `72/2^11`, even with nonperiodic replacement centres.** In particular the tempting 73-centre gain of 1/72 cannot arise this way. The result allows continuous motion: the point `u+s e_i` is valid for any coordinate `i` and `−1≤s≤1`, as every core word differs from `u` in at least five bits, so its squared torus distance from that point is at least `5+s²−2|s|≥4`. Taking one complete coset with any such shift ties `P11a`.

### Locating all new centres

Round an alleged new centre mod `2 Z^11` to a binary corner `b`. Let `x_i∈[0,1/2]` be its coordinate distances to `b`, and let `S_c` be the support of `b xor c` for core word `c`. Its squared distance to the nearest centre in the `c` coset is

```
Q + |S_c| − 2 sum_{i in S_c} x_i,      Q = sum_i x_i².
```

For a uniform distribution on any nonempty subset of core words, put `p_i=Pr(i∈S_c)`. The minimum is at most the average, which is at most `sum_i max(1/4,p_i)` by maximizing each convex quadratic on `[0,1/2]`. The checker uses the prefixes sorted by `|S_c|` (ties in generation order). Of the 2036 corners other than `u` and its 11 Hamming neighbours, 2032 have a strict certificate below 4. At the four remaining corners a prefix gives exactly 4, so equality requires *every* coordinate to be an endpoint of `[0,1/2]` (strict convexity). The literal half-coordinate evaluation gives squared distance at most 3 at every such vertex. Thus only `b=u` or `b=u−e_i` can round from a feasible new centre.

The 66 weight-five supports `u xor c` are the Witt `S(4,5,11)` design; the verifier checks every 4-subset appears once. Write `L=sum x_i`, and order `x_1≥…≥x_11` (temporary order statistics). When `b=u`, feasibility implies `2 M≤1+Q`, where `M` is the maximum sum on a weight-five support. The design gives `M≥(5/11)L` and `M≥S_4`, the sum of the first four. Also `Q≤L/2`, so `L≤(11/10)(1+Q)` and `Q≤11/9`. Taking `t=x_5` and applying the chord bound on `[t,1/2]` to the first four squares, and `x_i²≤t x_i` to all remaining squares, gives

```
Q ≤ S_4/2 + t(L−2).
```

If `L≤2`, then `M≥2Q`, giving `Q≤1/3`. If `L>2`, use `t≤S_4/4` to get `M≥4Q/L≥40Q/[11(1+Q)]`. This also forces `Q<1/3`, as `80Q>11(1+Q)^2` on `[1/3,11/9]` (check the increasing difference there). Thus in this rounding region all new centres have squared distance at most `1/3` from a lift of `u`.

When `b=u−e_i`, the distance from `u` in coordinate `i` is `s∈[1/2,1]`; all other 10 coordinate distances `x_j` are at most `1/2`. The 30 weight-five supports containing `i` become a `3-(10,4,1)` design after removing `i`. Put `a=(1−s)^2`, `Q=sum_{j≠i} x_j²`, `L=sum_{j≠i}x_j` and `M` equal to maximum sum on a design block. For those core words the formula for squared distance is `4+a+Q−2 sum_{j∈block}x_j`. Thus `2M≤a+Q` and `M≥(2/5)L`, with `Q≤L/2`. This gives `L≤(10/3)a≤5/6`. A block contains the top three coordinates, so `M≥S_3`; the same chord proof gives `Q≤S_3/2 + x_4(L−3/2)≤S_3/2`. Hence `Q≤a/3`. The total squared distance from `u` is at most `s²+(1−s)^2/3≤1`, with equality only at `s=1` and all transverse coordinates zero.

Every new centre therefore lies in one of the closed unit balls around `u+2 Z^11`. An interior centre uses that ball exclusively: its distance to any other point in the ball is less than 2. A boundary centre is an edge midpoint between adjacent ball centres, so assign it to the positive endpoint of the edge (according to the corresponding coordinate). Any two midpoints incident to a ball centre must be opposite along the same coordinate, since perpendicular midside differences have squared norm 2. Thus no ball is assigned two boundary centres; none can already have an interior centre. All new centres inject into `u+2 Z^11` with displacement at most 1. Counting in large balls (and changing their radius by at most 1) now shows that their upper centre density is at most `1/2^11`. This proves the claim without assumptions on replacement periodicity.

## 2. A coherent near-deep-hole obstruction to a 21D/22D fibre construction

Let `K=K12` have squared minimum 4 and covolume 27. Its dual is similar with scaling factor `1/sqrt(3)` (squared minimum `4/3`). In the quotient `K*/K ≅ F3^6`, `q(v)=(3/2)||v||² mod 3` is a well-defined quadratic form. The coset data of [Conway–Sloane 1983, §§3, 8–9](https://neilsloane.com/doc/Me107.pdf) give deep holes exactly in the `q=1` classes, with squared torus distance `8/3`. We also use the full Voronoi-vertex classification in [Pook-Kolb–Allen–Agrell, Table 2](https://arxiv.org/abs/2211.01987), in the same normalization: all vertices have squared norm `8/3` (deep) or at most `52/25`. In particular their exact nondeep gap is `g=44/75`.

**Lemma.** Among any nine real phases in `R^12/K` there are two at squared torus distance at most `64/27`.

To prove it, let `x` be in the closed Voronoi polytope of `K`, with `||x||²=8/3−epsilon` and `epsilon<8/27`. Write `x=sum lambda_v v` as a convex combination of vertices. The variance identity gives

```
epsilon = sum lambda_v[(8/3 − ||v||²) + ||v−x||²].
```

Because `g>epsilon`, some deep vertex `w` satisfies `d=||x−w||≤sqrt(epsilon)`: otherwise every summand with positive weight has bracket greater than `epsilon`. Set `A=sqrt(4/3)=min length(K*)`. Different deep vertices are at distance at least `A`; if `v` is nondeep, then `||v−w||²+g≥2 sqrt(g)||v−w||≥A||v−w||`. A second application of the variance identity, followed by the triangle inequality, gives

```
epsilon + d² = sum lambda_v[(8/3 − ||v||²) + ||v−w||²]
              ≥ A sum lambda_v||v−w|| ≥ A d.
```

First `d<A/2` since `epsilon<8/27<1/3=A²/4`. As `d(A−d)` is increasing there and its value at `A/3` is `8/27`, the inequality yields `d<A/3`. For any system of torus phases with all pairwise distances beyond `64/27`, pick real point representatives and round each difference to this nearby member of `K*`, choosing opposite signs for reverse edges. Every triangle of rounded differences sums to a member of `K*` of length less than `A`, hence is exactly zero. Thus these rounded differences really are differences of common quotient labels; there is no per-edge choice inconsistency. The rounded labels are distinct and all their nonzero differences have `q=1`. For nine labels, translate one to zero; their eight remaining labels have polar-form Gram matrix (using polarization divided by two) `2(I_8+J_8)` over `F3`, of rank seven. It cannot factor through a six-dimensional space. Contradiction.

Now take one complete real translate of `K` above each site of a fixed uniformly rescaled 9D `Lambda9` or 10D `P10c` normal packing. Horizontal phases need not be dual phases, additive or even periodic. Write `mu` for the total squared minimum, necessarily at most 4. The unscaled normal data `(minimum squared, covolume per centre)` are `(2,1)` for `Lambda9`, and `(4,2^10/40)` for `P10c`. For the first, `Lambda9=<E8,(e_1,1)>` with `E8` in the first eight coordinates; the nine sites `0, e_1+e_i (2≤i≤8), (1/2)^8` in its `E8` layer form a short-distance clique. For the second, the verifier supplies ten explicit sites in a binary translate of Best's code as obtained from [Best's original paper](https://ir.cwi.nl/pub/6831/6831D.pdf), all pairwise squared distance 4. If `s` is the squared transverse rescaling and `q0` the unscaled squared minimum, the lemma applied to one of these cliques implies

```
s q0 ≥ mu − 64/27.
```

For `m=9` or `10`, therefore, when `16/5≤mu≤4`, the overall centre density is at most

```
(mu/4)^((12+m)/2) (q0/(mu−64/27))^(m/2) / (27 · normal covolume per site).
```

This density expression applies also to a nonperiodic phase map: bounded representatives of phases can be taken in a fixed horizontal fundamental region, so the number of fibre points in large balls has the same leading term as the unshifted product. The bound is decreasing on this interval because the logarithmic derivative has the sign of `12mu−(12+m)64/27`, which is negative even at 4. At `mu=16/5`, it gives centre-density caps less than 0.18653 in 21D and 0.32380 in 22D. Both are strictly below the attributed benchmarks: `43^(21/2)/(2^41 3^(23/2))` from [Chen–Hu–Li–Wang–Wu 2025, §3](https://arxiv.org/html/2505.02394v1#S3) and `11^11/(2^23 3^(21/2))` from [Conway–Sloane 1996](https://neilsloane.com/doc/anti.pdf), respectively. Only fixed similar transverse stocks are excluded, and only for the indicated radii. No claim here for smaller radii, other transverse stocks or incomplete fibres.

Replay the finite checks and benchmark arithmetic with `node verify.js` (no dependencies). `PERIODIC_RIGIDITY.md` records a separate local theorem for the binary constructions.
