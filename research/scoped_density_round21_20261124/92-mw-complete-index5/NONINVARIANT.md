# Beyond invariant profiles: subspaces with a witnessed short row

Write K_j=D^j L. The complete fifth-jet proof uses vectors of
squared norm at most 24 **only in K_2**. Let A=L/K_5 and
B=K_2/K_5, of dimensions 10 and 6 over F_47. Their elementary
field structure follows from 47L=K_46. B is also a rank-two
module over R=F_47[D]/D^3.

Every primitive projective R-line P in B has a short row:
lift its three digits to an arbitrary primitive fifth-jet line
in A, and use the exact fifth closure. Its short witness in
K_2 restricts to P. Thus a sublattice M containing K_5, even
one not invariant under translation, inherits this witness
if W=M/K_5 and W intersect B contains such a primitive line.

## Exact count of that sufficient condition

The number of e-dimensional subspaces Q of B containing
some primitive R-line is C_e as follows (p=47):
- C_3=(p+1)p^2, the number of primitive projective lines.
- C_4=(p+1)p^2(p^2+p)+(p+1)p. For a fixed primitive line,
  there are p^2+p+1 one-dimensional extensions; just one
  is invariant because the quotient is cyclic of length three.
  A non-invariant 4-space cannot contain two distinct primitive
  lines: their sum would be invariant of dimension at least four.
  The invariant 4-spaces containing a primitive vector are of
  type (0,2) in B, counted by (p+1)p.
- C_5=(p^6−1)/(p−1). In fact any hyperplane given by a
  nonzero linear functional on the three jet pairs contains a
  primitive line. Choose the leading pair nonzero perpendicular
  to the last nonzero coefficient pair of the functional.
  The vanishing conditions on its generator and up to two
  D-shifts are then solvable successively for the two later pairs.

For fixed Q=W intersect B with dimension e, the quotient
projection of W in A/B has dimension d−e, where d=dim W.
For each such projection the graph of a linear map into B/Q
gives exactly one W. Thus the number of W is the sum over e
of C_e times the Gaussian binomial [4 choose d−e]_p times
p^((6−e)(d−e)). Here d=5,4,3 gives indices p^5,p^6,p^7.
The verifier checks the exact chart counts. This is an exact
count of the sufficient condition, not of all M that have a
short row.

For invariant profiles, if a≤2 and b≤5, extend a primitive
length-b line to length five and take its short row in K_2.
It lies in any type-(a,b) preimage of that line restricted
to K_a. This excludes all type (1,5) at index p^6, also (2,4)
there and (2,5) at index p^7. Types (3,3) and (3,4) contain
K_4 and its pinned norm-24 vector. Hence all 4,985,761 mixed
invariant profiles at index p^6 and 106,080 mixed profiles of
types (2,5),(3,4) at index p^7 fall short. At each such index,
24^46 < 16^46 p^index-exponent is the exact density comparison.

## Stronger count by ordinary field-lines

Projective field-lines spanned by known rows in B need not
generate the same cyclic R-module. For a row of depth d in
{2,3,4}, take its successive jet pairs v_0,v_1,v_2 as
available. Normalize the first unit leading coordinate to 1
and for d<4 use the unique translation T^j that makes the
corresponding next coordinate 0. The transformations of the
next two pairs are v_1−jv_0 and v_2−jv_1+binom(j,2)v_0.
This normalizes exactly the field-line and its T-orbit.
The orbit has 47 elements for d<4, since a fixed field-line
would have zero D-image, and one element at depth 4.
Close the normalized keys under the verified normalizer
U,C,F. The counts of keys at depths 2,3,4 are 68,724, 828, 12,
giving S=3,268,956 witnessed field-lines in B.

When e=1, use any witnessed line for Q itself. When e=2, let
c(Q) count witnessed lines in the plane Q. The incidence
sum I equals S [5 choose 1]_p, and the squared incidence
sum equals I+S(S−1), since any two distinct lines span
one plane. The two-term Bonferroni bound gives at least
I−binom(S,2)=10,955,198,305,026 such planes; the separate
Cauchy–Schwarz lower bound is 9,843,960,244,437.
These counts are disjoint from the earlier e≥3 count.
The same graph-lift enumeration thus excludes at least:

| Index in L | Sublattices containing K_5 obstructed by a short row |
| ---: | ---: |
| 47^5 | 9,179,147,596,354,897,233,629,066,132,555,953,562,124 |
| 47^6 | 4,184,518,740,099,439,915,762,768,085,636,498,388 |
| 47^7 | 857,658,426,082,247,424,333,241,185,192 |

The first bound is about 1.42% of the index-47^5 chart
containing K_5. These are scoped exclusions, not a bound
for the entire non-invariant chart nor a density record.

Replay: python -B verify_noninvariant_counts.py and
OPENBLAS_NUM_THREADS=1 python -B verify_short_lines.py
(the latter takes the same optional arguments as verify_fifth_complete.py).
