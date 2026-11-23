# All half-integral five-E7 stacking phases in 39–41D

Scoped upper bound for a construction family, **not** a general packing
or lattice upper bound. Let L=L35 be the five-E7 horizontal lattice
specified in certificate.json, of squared minimum 8 and covolume
2^(31/2). Place one complete translate of L above each point of an
arbitrary full four-dimensional lattice V. Each horizontal translate is
allowed any phase in (1/2)E7^5 / L. The phase assignment need not be a
homomorphism, need not be periodic, and is not otherwise restricted.
With squared packing minimum at least 8, this family cannot attain center
density 17.6, below the pinned 39D comparator 18.3592594841895.

## Exact horizontal facts

The fractional label F(v) of the horizontal phase above v is in
(1/2)E7^5/E7^5, encoded by five F64xF2 symbols. For a fractional
difference g, let d(g) be the maximum squared distance to L among all
integer translates of that fractional pattern. It is a half-integer
and is computed by the exact local mod-four table followed by the
syndrome min-plus convolution. The maximum over g is 7.5; d(0)=6.
Write H for the g with d(g)>6.5. Its members are the blanket radical
(all five symbols 64) and the one-defect labels (four symbols 64 and
one projected nonzero local type with minimum mod-two norm 4).
Two consequences of the completed finite census:
1. For any h in H and any g, d(g) and d(g+h) cannot both exceed 6.
2. For distinct h,h' in H, d(h+h')<=5.5.

The first is the prior complementary-radical maximum for the blanket;
for exceptional h, enumerate every deep (d>=6.5) g in symmetry
representatives and all 315 possible exceptional h relative to it. The
second follows because a sum has zero radical bits outside its
nonzero projected symbols, has projected support at most two, and is
nonzero. Different defects in the same slot have different projections
since local norm 4 fixes the radical bit for a nonzero projection.
Exact max depth is at most 5 for one occupied slot and at most 5.5 for
two when the outside radical bits are zero.

Common SL(2,8) and the simultaneous field-squaring slope permutation
preserve phase composition, depth, and the set H. Their verification
and full orbit-normal-form proof are inherited from the earlier
full-half census; deep.cpp augments its identical exact loop to print
all values >=26 in quarter units. Symmetry preserves fractional XOR;
integer shifts within a pattern are maximized over separately, always
an outer bound on the actual chosen translations. pair_exact.cpp works
on the printed deep representatives plus every support-zero pattern
and computes the other depth by independent fresh raw syndrome DP.
Support 3 can be skipped only because its already complete census
gives d<=6. All comparisons use integer costs.

## Geometry; no additivity assumption

Suppose covol(V)^2<=100/121. Let v be a shortest vector, a=||v||^2.
Global horizontal depth 7.5 implies a>=1/2, and the exact 4D Hermite
bound excludes a>=1.5. Hence each edge parallel to v of squared length
a has fractional phase difference in H, possibly varying by edge.

Project V along v to a ternary lattice T. For any projected nonzero
difference p, write q=||p||^2. Its two nearest longitudinal lifts have
offset magnitudes t and sqrt(a)-t, with min(t,sqrt(a)-t)<=sqrt(a)/2.
For any starting point of V, their two horizontal fractional
differences sum to a member h of H. Define
    R(a)=2-a for a<=1, otherwise 2-1/a,
    B0(a)=1.5-(a-.5)^2/(4a),    B(a)=min(B0(a),2.5-a).
If one fractional difference is zero, its vertical squared distance
is at least 2 (d(0)=6). The other vertical squared distance is at
least a (shortest-vector definition). Minimizing over the two offsets
gives q>=R(a). If neither is zero and at least one has depth>6.5,
it is in H and differs from h; their sum has depth<=5.5, so q>=2.5-a.
Otherwise the depths have upper bounds 6 and 6.5 in some order by
fact 1; minimizing max(2-t^2,1.5-(sqrt(a)-t)^2) gives q>=B0(a).
Therefore T has minimum squared at least min(R,B); and if q<B, for
**every** starting point one of the two lifts has zero fractional
difference. The smaller longitudinal lift cannot be that one since
q+a/4<B+a/4<2. Ties of the two offsets thus cannot occur for q<B.
The unique lift with larger longitudinal magnitude is a period of
the fractional labels independent of starting point. Its squared
length is less than B(a)+a<=2.5.

Take a 3D obtuse Selling superbase of T. If its cut vectors of squared
norm <B do not span the parity space, some nonzero binary functional
vanishes on them. Four parity labels then have squared norms >=B,
and the other three >=min(R,B). Exact enumeration of nonnegative
Selling conorm polytope vertices, with lower rounding and interval
bands for a, gives covol(V)^2>100/121 (geometry_final.py). Determinant
log-concavity and the nonnegative recession directions justify
minimizing over vertices. The same checker also applies the 3D or
4D Hermite bound where sufficient.

Otherwise three mod-two independent short cut vectors form a Z-basis
of T: their 0/1 coefficient determinant is odd, and a 3x3 0/1
determinant has magnitude at most 2. Lift all three by the unique
longer rule to W. Thus V=W+Zv, W has rank 3, all W shifts preserve
fractional labels, and the three basis lengths of W are <sqrt(2.5).
Every nonzero W vector has squared length >=2. By the 3D Hermite
bound p_W=covol(W)^2>=4. A KZ basis of W has squared projected (GS)
lengths x_i<=2.5 since there are three independent vectors with
these lengths. Hence its squared covering radius is at most
(sum x_i)/4 by nearest-plane rounding. Given 0<x_i<=2.5 and fixed
product p_W, sum x_i<=5+4p_W/25: write y=x_i/2.5 and use
1−y_1y_2y_3 <= sum(1−y_i) for 0<y_i<=1. There is therefore a
vector in 2v+W of squared length at most
    (5+4p_W/25)/4 + 4 covol(V)^2/p_W < 2.5.
The endpoint comparison follows from convexity in p_W on [4,125/8],
and is checked rationally. Its fractional phase difference from
any starting point is the sum of the two consecutive v increments
because W preserves fractional labels. If those increments differed,
fact 2 would require squared vertical separation >=2.5. Therefore
every two consecutive v increments agree: 2v also preserves labels.
The subgroup W+2Zv has index two in V and minimum squared at least
2. The 4D Hermite bound says its covolume >=2, so covol(V)>=1,
the final contradiction.

Covolume per occupied layer is covol(L)covol(V), so center density
is 16/covol(V). This proves the scoped claim (17.6), not a new record.
Different horizontal periods, sparse/disordered vertical site sets,
or shifts off the half-grid are not excluded.

## Exact replay

local4.cpp/local4.in regenerate the local table; joined field symbols,
local table, and slope data form table.in. The prior portable
e7_global_half_20261123.tar.gz gives the proof of the complete
symmetry reduction and the local table. A copy of the needed source
files is bundled here. Run verify35.py, verify_sl.py and replay.sh.
geometry_final.py uses NumPy for exact integer arrays, and Fraction
for all scaling/comparison arithmetic. deeplist.in and replay logs are
stored; verify_list.py checks strict agreement with all original
census counts. The exploratory LP shape search is not part of proof.

## The full 40D analogue

The same horizontal facts and deterministic longer-lift rule also
apply when V has rank five. In that case set U=61/100 and suppose
covol(V)<=U. The exact rank-five Hermite bound again excludes a>=1.5.
All four projected KZ basis vectors, after size reduction, have squared
length <B(a): geometry40.py checks the uniform bounds for
(x_i + (sum_{j<i} x_j)/4)/B(a) over the exact outward-rounded logarithmic
polytope of necessary KZ block Hermite inequalities. The largest bound
is 13529087/14375000 <1. The checker first excludes any interval
where the plain 4D projected Hermite bound already suffices. Thus
no Selling parity argument is needed here: the four KZ generators
all lift to periods of the fractional labels and generate W.

As before `p_W=covol(W)^2>=4`, since W has rank 4 and minimum
squared at least 2. Write its Gram as G + t t^T, where G is the Gram
of the basis in projection to v-perp. Then
    p_W = det(G)+t^T adj(G)t
        <= U^2/a + 4a (4 B(a)/3)^3 <= 2 U^2 + 32.
Here ||t||^2<=4a, trace(G)<=4B(a), and the operator norm of adj(G)
is the product of the largest three eigenvalues of G, at most
(trace(G)/3)^3. Pass from a B^3 to 27/8 using a+B<=2.5 and
B<=1.5, since (2.5-B)B^3 increases for B<=1.5.

The KZ projected squares of W are again <=2.5 since its four
constructed basis vectors have that bound. Apply the same product
inequality to four numbers to get sum x_i <= 7.5+p_W/(2.5)^3.
Consequently 2v+W contains a vector of squared norm at most
    (7.5+p_W/(2.5)^3)/4 + 4U^2/p_W < 2.5,
with the endpoint check in geometry40.py. Period two along v follows
as in 39D. Now the rank-five index-two integer-phase kernel has
minimum squared >=2 and covolume >=2, so covol(V)>=1, contradiction.
Hence center density <16 sqrt(2)/U < 3^16/2^20 in this 40D family,
also allowing arbitrary nonadditive full layer assignments.

The [P48 cut and antipode comparator](https://arxiv.org/abs/2607.20359)
and [Selling conorm geometry](https://neilsloane.com/doc/Me169.pdf)
supply comparison and reduction references respectively.

## Run receipt

The fresh support-5 extraction reached at least the printed index
647168 before its 520-second cap. The same normal-form loop resumed
at index 647168 and completed the remaining 9971 orbits, with maximum
depth 5.5 and no additional deep labels. The previous independently
completed support-5 histogram has exactly thirteen deep labels,
agreeing with the frozen list. A standalone uncapped replay may be
run later; replay.sh compares the deep records rather than requiring
the intentionally capped run's progress format. No dependence on
randomized or optimization scout output.

## 41D closure added after pruning and rank-four lift

This supersedes the unpruned 41D scout above. Assume rank-six vertical
V with covol(V)<=13/25. Direct rank-six and projected rank-five Hermite
bounds prune intervals of a. geometry41_pruned.py proves the first three
size-reduced projected KZ vectors have squared norm <B(a), maximum
ratio <.955. Lift them to fractional periods. The shortest vector in
the remaining rank-two projected quotient is primitive. Its squared
length is at most sqrt(4 det(T)^2/(3*x0*x1*x2)). Add the squared covering
radius of the initial rank-three projected span to lift that vector
to a fourth projected period candidate. The bound is (x0+x1+x2)/4,
improved by 7*min(R,B)/80 when x1<=.9*x0: the reduced angle's absolute
Gram-Schmidt coefficient is then >.3 by the KZ shortest condition;
splitting the first plane into two triangles improves its radius as
detailed in kz_reference.py and the prior KZ argument. Exact outward
log vertices in fourth41.py give the ratio <.986. Hence four projected
vectors lift to a rank-four period lattice W. Their projections
generate the primitive rank-four KZ prefix.

Use only the rank-five slice P=W+Zv. Its squared covolume is det(V)^2/x4
where x4 is the last GS squared length of T. The same exact KZ scan
gives covol(P)^2<331/500. By the rank-four adjugate estimate from 40D,
p_W=covol(W)^2<=2*(331/500)+32<34. Again W has minimum squared 2 or
more, and its constructed independent generators have squared lengths
<2.5. Its KZ squares lie in (0,2.5], x0>=2, consecutive ratios >=3/4,
and product in [4,34]. The exact vertices in escape41.py bound
sum(x_i)/4+4*(331/500)/product(x_i) by 3121103183/1250000000<2.5.
Thus 2v+W contains a vector forcing period two, even for nonlinear
assignments. The index-two subgroup W+2Zv of P has squared minimum
at least 2, hence P has covolume at least 1 by rank-five Hermite, a
contradiction. The resulting 41D scoped center cap is 32/(13/25) =
800/13, below 3^17/2^21. This gives all three 39–41 full-lattice
carriers, still not an upper bound for arbitrary packings.

Archived replay SHA-256: `afe77c62cc22a2d00e29882c6a4631d3b28e79e2cf8099d08b36eee3db33d685`.
