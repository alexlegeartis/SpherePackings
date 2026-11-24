# Fourth jet: partial short-vector obstruction at index 47⁴

For the same parent as in the main replay, at least **3,665,047 of
4,985,761** translation-invariant index-47⁴ sublattices have a nonzero
vector of squared norm at most 22. Since 22^46 < 16^46·47^4, none of
those improves the MW92 parent. The other 1,320,714 are not covered by
this check. No completeness of a short shell or positive construction
is claimed.

The local DVR argument in README still applies. The elementary divisor
types are (0,4), (1,3) and (2,2).

- Type (0,4): a primitive projective line modulo D⁴. There are
  (47+1)·47³ = 4,983,504 lines. Emit four successive residues to obtain
  the rank-eight jet matrix; D shifts as before. A depth-d witness hits
  all 47^d lifts of its modulo-D^(4−d) projective prefix.
- Type (1,3): D times one of the (47+1)·47 = 2,256 primitive
  index-47² lines. For a witness in DL, dividing by D simply
  shifts the jet by one, then keep the first two digits of the jet
  of the divided row. At depth one the original witness therefore
  supplies a two-digit ratio; at depth two it supplies one digit
  and covers all 47 possible next digits.
- Type (2,2): the unique D²L, containing the pinned norm-20 witness.

The bounded screen uses the 119 pinned rows plus the two final literal
rows in the main replay. Apply the same direct signed-pair and
cusp-inversion loops, but retain norms through 22. The 6,643 verified
short differences/sums and their images produce 317,147 distinct
nonzero initial fourth jets. Close their prefixes under the checked
normalizer action. Depth-2,1,0 give respectively 144, 3,318, and
3,528,666 prefixes, with no depth-three witness. Counting avoiding
overlap by taking each shortest witnessed prefix gives 318,096,
149,178, and 3,195,804 incremental exclusions, for 3,663,078
primitive lines. In addition, 1,968 mixed (1,3) lines and the
unique (2,2) are excluded.

Run **OPENBLAS_NUM_THREADS=1 python verify_fourth.py**, next to the
main verify.py; it first replays the published third-jet certificate.
The run checks actual integer norm witnesses and exact quotient
arithmetic. It does not treat the absence of a depth-three witness
as a lower bound on D³L.
