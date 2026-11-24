# Every finite MW80 label cluster: exact optimum of the 79D loss guarantee

For the named MW80 parent of determinant 41^12 and inherited squared minimum
14, consider all primitive integer normals of rank one, and all finite
sets of at least two distinct integer labels. Consider the standard
projection construction when only the parent minimum and the largest
loss over label differences are retained. Its exact best **guaranteed**
density is the existing two-consecutive-label 79D construction with
dual-normal squared length 234/41. This is a global comparison of the
section/label formula, not of the actual minima when the worst loss fails
to coincide with a short vector.

Let b be the primitive normal, K=b G^(-1) b^T, M the number of labels, and
ell the difference between their largest and smallest labels. Then
ell is an integer at least M−1, and K=m/41 for a positive integer m:
the exact inverse identity G D=41I with D integral is also checked in
verify.py. Ignore clusters where ell²/K≥14, since these have no positive
guarantee. If x=ell²/K then the guaranteed square of the center density,
after dropping the common 4^79 det(G) denominator, is

    (M²/ell²) x (14−x)^79.

The unconstrained function f(x)=x(14−x)^79 is increasing before x=7/40
and decreasing after it. For M≥3, the multiplier is at most 9/4; for M=2
and ell≥2 it is at most one. Both are strictly inferior to four times
f(41/234), by exact rational comparison. For M=2 and ell=1, the allowed
x are 41/m. By monotonicity on each side of 7/40 only the integers
m=234 and 235 can maximize f(41/m). Exact arithmetic gives the former
strictly larger. This normal and the consecutive labels are attained
by the literal witness in sections.json and by the existing section
certificate, where the minimum-distance equality is proved.

This does not bound one-translate sections, replacements of the MW80
carrier, exact distance improvements from a stronger norm/label analysis,
or other relative shifts. For small arbitrary relative shifts of the
specific fixed 79D period see README.md. Comparator and exact parent:
[creative sections](https://github.com/alexlegeartis/SpherePackings/blob/6696a479a8a55e440ff8e65f10702702fe7d599c/research/creative_sections_20261116/PROOF.md)
and Schütt–Shioda, *Mordell–Weil Lattices*, Theorem 13.35.

Run python -B verify_projection.py for the exact comparisons.
