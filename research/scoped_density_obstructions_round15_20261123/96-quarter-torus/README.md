# Conditional rank-96 fourth-torsion bound

If an even unimodular lattice \(L\) of rank 96 with squared minimum 10 exists, then
\[
 \operatorname{dist}(v/4,L)^2\le 8\qquad(v\in L).
\]
This does not assert existence of such a lattice, or a bound on the ambient real covering radius. The hypothesis is extremality (minimum 10), not merely minimum 8.

The [half-torsion analysis](HALF_BOUND.md) provides \(H=\operatorname{dist}(v/2,L)^2\le8\), at most \(15/2\) for anisotropic classes, and the 2-design property at depth 8. The purely geometric argument of [PRIOR_GATES.md](PRIOR_GATES.md), with \(D=8\), then either immediately establishes the desired bound, or gives one of the 22 residue/depth rows listed in `certificates.txt`. The arguments below close all 22 by comparing near Fourier coefficients at two related cusps.

Use the notation of [PRIOR_GATES.md](PRIOR_GATES.md): \(r=v^2/2\bmod4\), and \(P=1,\alpha=0\) when \(r=0\), or \(P=\theta_4(2\tau)^r\theta_3(2\tau)^{4-r}, \alpha=(4-r)/4\) otherwise. A presumed counterexample with depth \(>8\) forces the scalar quarter-twist form to be of the shape
\[
 f=P Y^iW^j\sum_{u=0}^d c_uX^uY^{d-u},\qquad
 i=16+[r=0],\quad j=\lceil H/2-\alpha\rceil,\quad d=24-[r\ne0]-i-j.
\]
Its coefficients through exponent 4 at infinity are \(1,0,0,0,0\). At the secondary cusp \(M=ST^{-2}S^{-1}\), expand by substituting \((X,Y,W)\mapsto(W,-Y,X)\) and \(P\mapsto\theta_4(2\tau)^r\theta_2(2\tau)^{4-r}\), up to a fixed unit complex phase. Write \(a_j\) for the coefficients at exponents \(\alpha+j\) after removing this phase. They are real because the polynomial is fixed by conjugation at infinity before the formal substitution. Irrespective of that phase, Poisson gives \(\lvert a_j\rvert\le N_{2(\alpha+j)}\), where \(N_s\) is the unsigned number of norm-\(s\) points of \(v/2+L\). The unused potential leading cusp-zero coefficient also imposes \(c_d\ge0\); it may be zero if the vanishing order exceeds the one imposed.

In the half-torsion theta form, removing the multiplier \(\theta_3(\tau)^2\theta_4(\tau)^2\) in anisotropic cases leaves a polynomial divisible by \(Y^{\lfloor H\rfloor}W^{\lfloor H\rfloor}\) at argument \(\tau\). The quotient is a homogeneous polynomial in \(Z=X+2Y,X^2\) with weighted degrees 1,2 and total degree \(24-\mathbf1_{H\notin\mathbb Z}-2\lfloor H\rfloor\). Its first five original coefficients at integer exponents are again \(1,0,0,0,0\). The \(S\) cusp, now with \(X,Y\) interchanged and \(\theta_3^2\theta_2^2\) restored when required, gives affine formulas for \(N_s\) in the free coefficients. Ordinary unsigned constraints are \(N_s\ge0\) and \(N_s=0\) for \(s<H\). Also:
* If \(H<5\), then \(N_H=2\), and \(N_s=0\) for \(H<s<10-H\).
* If \(H=5\), then \(N_H\le192\).

Indeed non-antipodal vectors \(x,y\) in the half coset have nonzero \(x\pm y\in L\), so their mean squared norm \(x^2+y^2\) is at least 10. When both are of norm 5, this forces orthogonality and at most 96 antipodal lines.

Four rows give inconsistent initial coefficient equations alone. For each other row the verifier builds a short list of affine expressions that must be nonnegative by the above inequalities, and certifies a linear combination with explicitly supplied nonnegative rational weights equal to the negative constant \(-1\). Thus neither numerical feasibility decisions nor unknown shell classes are inputs to the proof. Run python verify96.py (Python 3 and SymPy).
