# 48D: pair of scalar E8 directions, with correlation

In the positive-definite 48D seven-slot zero-sum half-period E8
model with ordered weights (a,b,c,d,1,1,1), a<0 and 2≤b≤c≤d,
the center density is at most **16,800** for any motif if
every heavy weight is in [2,7/2], at **every** scalar determinant
D. This expands the previously closed [2,5/2] cube and remains
strictly below the [48D record](https://cohn.mit.edu/sphere-packing/).
This statement neither changes the record nor bounds unrestricted
48-dimensional packings.

## Two independent residue directions

Use the same normalization as the
[previous round](../../scoped_density_round21_20261124/48-e8-four-exceptions-scopes).
The carrier determinant parameter is
D=bcd+a(bc+bd+cd+3bcd); the squared minimum is g and the
center density is N(g/8)^24/D^4. From a scalar integer vector u
summing to zero, the doubled-period vector along a fixed E8
root has squared length p=4Q(u). It bounds g. If u is nonzero
modulo 2, it also embeds an E8/2E8 subgroup in the phase group.
For one direction, when g>p/4 the fibre has at most 16 labels
(differences are totally singular); when g>p/2 the fibre has one.
Consequently N g^24≤2^40 p^24.

Take two scalar directions u,v independent modulo 2, with
p=4Q(u), q=4Q(v), and K=|Q(u,v)|, where Q(u,v) is the
bilinear pairing of the scalar form. Their phase subgroup is
a direct sum of two full E8/2E8 factors with 2^16 labels. An
anisotropic phase has an E8 root lift of squared norm 2. Each
nonzero isotropic phase has 16 lifts of squared norm 4, forming
eight orthogonal axes. The standalone finite shell replay
verifies these facts using the doubled E8 coordinates and the
previous anisotropic label dictionary.

The phase difference lift in this subgroup has squared norm
Q(u)||x||²/2+Q(v)||y||²/2+Q(u,v)〈x,y〉. Either sign is
available for each representative. For a root and an isotropic
8-frame, the sum of squared inner products along that frame
is 8 and the pairings are integers, so one pairing has
magnitude at least 1. For two isotropic phases, take any
norm-four representative of the first; its squared pairing
sum in the second frame is 16, so some pairing has magnitude
at least 2. Thus choose the sign making the cross contribution
at most −K or −2K, respectively.

Set
    s = max(p/2, q/2, (p+q)/4, p/4+q/2−K, p/2+q/4−K),
    t = max(s, (p+q)/2−2K).
Above s, every nonzero allowed difference inside one paired
fibre must be nonzero isotropic in both factors; projection to
either is injective, and translating to contain zero makes
the projected span totally singular, of dimension at most four.
The fibre cap is therefore 16, giving N≤2^36. Above t all
nonzero differences are forbidden, giving N≤2^32. This is a
joint lift estimate, not a product of two marginal caps.
The single-direction bound handles the low range: since
s≥min(p,q)/2, if g≤s one has
N g^24≤2^40 min(min(p,q),s)^24 by considering
g≤min(p,q)/4, the middle single-direction range, and its
single-label range. The strict thresholds assign equality to
the lower range.

## Exact parameter coverage

At fixed D, a is separately convex and decreasing in each
heavy weight; for fixed weights a is increasing in D.
Therefore each period's maximum upper bound on a box is
evaluated at maximal D and the eight heavy corners.
For the cross pairing we use an outward rational interval
through monotonicity of a and the coefficient products; its
distance from zero is a safe lower bound for K.
Apply the three density inequalities above on every verified
bin, always using the box's lower D in the denominator.
Replay subdivides closed boxes and checks all endpoints as
rationals. No optimization solver or floating point is required.
Outside [1/3500,128/1875] the previous tail proof in D applies.

Run python3 -B verify.py. The receipt reports separately the
leaves where a single period at that leaf would not suffice,
and those where even the 16-label paired cap alone would not
complete the density comparison. In the latter leaves the
correlated lift forces the final one-label paired-fibre cap.

The only final low-determinant cluster missed by the shorter
period dictionary is covered by scalar coefficients
(-110,12,9,8,27,27,27). Its root period has squared length below
.439 on the entire enclosing box. This illustrates why failure to
find a short period with a bounded initial coefficient cutoff does
not certify a new packing candidate.
