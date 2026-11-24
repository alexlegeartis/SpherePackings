# Equal-scale A2^5 third- and sixth-torsion

A literal nonlinear 18-site construction uses the ternary length-five
distance-three code recorded in `ternary-18.txt`. A symbol c labels the
third-period shift c(1,2)/3 in the simple-root basis of A2. Every nonzero
difference in this three-symbol line has squared depth 2/3. Thus three
differing positions give squared minimum 2, attained also by the period.
The center density is 18/(32 sqrt(243)) = 1/(16 sqrt(3)), tying the known
10D lattice density, not the unrestricted 5/128 record. There cannot be
a lattice among these 18 shifts because its subgroup index would be
a power of three. This is an explicit alternative carrier/packing.

The exact ordinary Fourier cap on arbitrary nonlinear third-torsion
motifs of the equal-scale A2^5 carrier is 1107/59, so the 18-site
construction is optimal in that fixed finite family. For the full
sixth-torsion carrier the corresponding bound is 144/7, allowing at
most 20 sites. Twenty would improve the record (about 2.64%); no 19-
or 20-site construction was found and the full sixth-torsion case
remains open. Exact replay with Python and NumPy:
`python -B verify.py`.

The stronger local constructive attempt used the 12-symbol subgroup
per block obtained from arbitrary 2-torsion plus the deep 3-torsion
line. Its (parity, ternary-symbol) pair depths in units 1/36 are 0 if
equal, 18 for parity change only, 24 for ternary change only, and 6 for
both changes. The normalized physical gate is total at least 72.
All 915 deletion-mask representatives with 6–14 deletions around the
18 construction were exhausted over the full 12^5 donor universe.
Their stabilizer group has 360 elements, obtained by enumerating
coordinate permutations and independent ternary symbol permutations;
these act isometrically on the 12-symbol metric and on the anchor set.
Every outside nondeep position blocks at least six anchors. For each
mask, retain the other anchors, include all donors whose exact blocker
set is contained in the mask (including removed anchors that can be
restored), and run the supplied branch-and-bound. The three pools larger
than 4096 used an 8192-bit universe; all other pools used 4096. All
node counts were below the stated cutoffs (max 2277). No gain and no
mixed-parity tie were found in that precise local window. This does
not bar much larger exchanges or positions outside the subgroup.

The full parity simplex16 palette with independently arbitrary deep
ternary shifts also had maximum 18 after complete exact branching.
Unrestricted 12-symbol greedy and restart searches reached no larger
motif. The 11D two-layer A2^5 third-torsion and A1 half-period branch
had a numerical Fourier aperture (37 cap versus 36 needed), but
constructive trials remained far below threshold.
