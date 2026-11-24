# Strict local optimum in the full cyclotomic trace family for two 88D parents

For the catalogued rank-88 parents X4049 and X4969, the unwarped
metric is a strict local sphere-density optimum among all positive definite
metrics invariant under its order-230 catalogue operator \(U\), modulo
scaling. This is a 43-dimensional shape family, not only the golden
one-parameter slice, but is not a global optimum and does not address
deformations that destroy the \(U\)-symmetry.

The catalogue operator is multiplication by \(-\zeta_{115}\) in the
cyclotomic ideal, so its 88 complex eigenvalues are distinct and pair in
44 conjugate rotation planes. A positive self-adjoint operator \(H\)
relative to the unwarped \(G\) is \(U\)-invariant exactly when it commutes
with \(U\); on each real rotation plane it is a scalar. Thus the real
space of such operators has dimension 44, with basis \(I\) and
\(U^k+U^{-k}\), \(1\le k\le43\). Independence follows since these are
polynomials of successive degrees in \(U+U^{-1}\), whose 44 distinct
eigenvalues are twice the rotation cosines.

The even unimodular minimum-eight theorem for these two parents is
taken from Gauvin's catalogue references. The norm-eight shell is a
spherical 3-design: for a harmonic quadratic \(P\), its harmonic theta
series has weight 46 and vanishes to order at least four, so division
by \(\Delta^4\) leaves negative weight and the series vanishes. The
shell is centrally symmetric, which also kills the odd terms. In
particular its covariance is scalar. The 44 covariants
\((\langle xU^k,x\rangle_G)_{k=0}^{43}\) from each choice of 44 norm-eight
witnesses span over \(\mathbb R\): the integer matrix has rank 44 modulo
1009 by exact elimination. This also spans the evaluations on the
displayed basis, since for nonzero positive \(k\) the basis's evaluations
are twice those covariants.

For any trace-zero nonzero \(U\)-invariant self-adjoint direction \(D\),
the norm-eight shell average of \(\langle xD,x\rangle_G\) is zero by
the 3-design; some value is nonzero by the spanning witness. Hence
some value is strictly negative. The minimum over the finite shell is
continuous in \(D\), so on the unit sphere of trace-zero directions it
is uniformly bounded above by a negative constant. For every sufficiently
small nonzero trace-free \(D\), the squared minimum of \(I+D\) is
therefore at most \(8-c\|D\|\) while the determinant correction is
\(1+O(\|D\|^2)\). This gives a strict local density loss. Every neighboring
metric can be scaled to this trace gauge without changing density.

Replay with the seven-parent mathematical input from the main packet:
python -B verify_commuting.py ../inputs.bin.gz, or pass its actual
path. The sparse release consists of this note, the verifier, and
vec.bin.gz with 44 literal norm-eight vectors in the original basis
for each lattice (X4049 first, then X4969, signed little-endian
int16). Both exact rank checks and all integer norms are verified.
