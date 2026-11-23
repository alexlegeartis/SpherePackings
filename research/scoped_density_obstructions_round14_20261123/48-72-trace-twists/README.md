# Exact density maxima for real-quadratic trace twists

This is an exact optimum within the specified twist families, not a global lattice-density bound or a new density lower bound.

## Scope and result

Let \(L\) be an even unimodular lattice of rank \(n\) with squared minimum \(m\). Let \(S\) be a self-adjoint integral endomorphism of \(L\) with \(S^2=D I\) for \(D=2,3,5\). Integral means \(S(L)\subseteq L\), not that \(S\) is an isometry. For \(D=5\), require \((I+S)/2\) to preserve \(L\), giving an action of the full ring of integers. Consider all positive forms \(aI+bS\) on \(L\) for real \(a,b\). This includes all positive trace twists of the action and all real interpolations of the two embeddings' weights.

- **\(D=2\):** for \((n,m)=(48,6),(72,8)\), the maximal Hermite invariant is \(m\). In rank 72 equality forces equivalence to the original via scale and a field unit. In rank 48 the proof allows additional equality only at two normalized parameters \(t=\pm4/9\). If \(\langle x,Sx\rangle\) is even on \(L\), these two possibilities are also strict.
- **\(D=3\):** for \((n,m)=(48,6)\), assume \(\langle x,Sx\rangle\in 2\mathbb Z\) on \(L\). The maximal Hermite invariant is \(m\). Only the original and endpoint forms can tie up to scale and unit equivalence; the endpoint ties exactly when it is itself extremal. The parity condition includes \(S=U+U^{-1}\) for any integral orthogonal automorphism \(U\) with \(S^2=3I\).
- **\(D=5\):** for \((n,m)=(24,4),(48,6),(72,8),(96,10)\), the maximal Hermite invariant is \(m\), attained only by scale and unit equivalents of the original. The final entry is conditional on existence of a rank-96 extremal parent and does not construct one.

The rank-48 center-density optimum in these families is \((3/2)^{24}\), and the rank-72 optimum is \(2^{36}\). The result works for every appropriate action on every extremal even unimodular parent, without a fixed coordinate representation. It does not apply to arbitrary SPD metrics or to lattices without such a parent structure.

## Normalization and shell moments

The eigenspaces of \(S\) each have dimension \(n/2\), since \(X^2-D\) is irreducible over the rationals and \(L\otimes_{\mathbb Z}\mathbb Q\) is a vector space over \(\mathbb Q(\sqrt D)\). After scaling, the form is \(I+tS\), \(|t|<1/\sqrt D\), with determinant relative to the base form \((1-Dt^2)^{n/2}\). For \(q=\langle x,x\rangle\) and \(s=\langle x,Sx\rangle\) the determinant-one squared length is
\[
\ell_t(x)=\frac{q+t s}{\sqrt{1-Dt^2}}.
\]
Equivalently the determinant-one form is \(B(r)=\cosh(r)I+\sinh(r)S/\sqrt D\), where \(t=\tanh(r)/\sqrt D\).

Choose the integral automorphism \(I+S,\ 2I+S,\ (I+S)/2\) for \(D=2,3,5\), respectively. Its eigenvalue absolute values are \(\varepsilon,\varepsilon^{-1}\) with \(\varepsilon=1+\sqrt2,\,2+\sqrt3,\,(1+\sqrt5)/2\), respectively. Pullback by its \(k\)-th power changes \(r\) by \(2k\log\varepsilon\). Hence it suffices to take \(|r|\le\log\varepsilon\). This is \(|t|\le1/2\) for \(D=2,3\) and \(|t|\le1/5\) for \(D=5\). Prove the inequalities for \(t\ge0\) and use \(-S\) for the opposite half.

The minimal shell is an 11-design in every parent above. For even harmonic degree \(j\) with \(0<j\le10\), the weighted theta series is a level-one modular form of weight \(n/2+j\) whose order at the cusp is at least \(m/2\). Division by \(\Delta^{m/2}\) produces negative weight \(j-12\), so it vanishes. Odd harmonics sum to zero by antipodality. Here \(n=24k,\ m=2k+2\).

Uniformly on the squared-radius-\(m\) sphere, \(s=m\sqrt D(2X-1)\), where \(X\) is \(\mathrm{Beta}(n/4,n/4)\) by projection onto either eigenspace. The 11-design therefore gives
\[
\mathbb E s=\mathbb E s^3=\mathbb E s^5=0,\qquad
\mathbb E s^2=\frac{m^2D}{n/2+1},\qquad
\mathbb E s^4=\frac{3m^4D^2}{(n/2+1)(n/2+3)}.
\]
On the even grid, if \(s\le-4\) never occurs, \(s^3-4s\) is nonnegative; its zero expectation restricts the support to \(\{-2,0,2\}\). This contradicts either \(\mathbb E s^2>4\) or \(\mathbb E s^4>16\). Apply also to \(-s\) for the positive sign.

The other witness thresholds use monic polynomials nonnegative on the remaining grid if \(s\le-h\) never occurs. Beyond the first root, all roots come in adjacent grid pairs. Their negative expectation is contradictory, and odd moments vanish, so the argument works also for \(-s\).

| \(D,n,m\) | Grid | \(h\) | Polynomial roots | Expectation |
|---|---|---:|---|---:|
| 2,48,6 | integer | 3 | -2,0,1,2,3 | -1152/25 |
| 2,72,8 | integer | 4 | -3,0,1,3,4 | -7040/481 |
| 5,72,8 | even | 6 | -4,0,2,4,6 | -286720/481 |
| 5,96,10 | even | 6 | -4,0,2,4,6 | -912000/833 |

## Proof on the whole interval

A norm-\(m\) witness with \(s\le-h\) gives \(\ell_t\le m\) whenever \(0\le t\le 2mh/(h^2+Dm^2)\), strictly inside this interval when \(t>0\). Compare \(m-ht\) with \(m\sqrt{1-Dt^2}\) after squaring; all the numerators used below are positive.

**Five.** The pairing \(s=2\langle x,((I+S)/2)x\rangle-q\) is even. The cubic yields \(h=4\) at ranks 24 and 48, the table yields \(h=6\) at 72 and 96. All four satisfy \(1/5<2mh/(h^2+5m^2)\), so any nonzero parameter within the fundamental interval is strictly worse.

**Two.** At \(t=1/2\) use the integral form \(A=2I+S\). It is 2-modular: relative to the original unimodular form the scaled dual has Gram operator \(2A^{-1}=2I-S=A(I+S)^{-2}\), an integral change of basis. Rains–Sloane give \(b=8\) or \(10\) as the squared-minimum bound in ranks 48 or 72, respectively, for both even and odd 2-modular lattices (Theorem 2 and the \(N=2\) assertion in Section 7). Pick an endpoint witness with \(2q+s\le b\). Its image under \(I+S\) has \(q'=3q+2s\ge m\). Hence \(q+q'\le2b\) and \(q\le2b-m\).

The shell witnesses cover \(0\le t\le4/9\) in both ranks. At \(4/9\), the endpoint affine witness satisfies \(q+4s/9\le(q+4b)/9\le(6b-m)/9\), while the target \(m\sqrt{1-2t^2}\) is \(7m/9\). At \(t=1/2\) it is strictly below the target since \(q+s/2\le b/2<m/\sqrt2\). Concavity of the target extends the bound through the full remaining interval. The join is strict in rank 72 because \(6b-m<7m\). In rank 48 the join can only tie if the everywhere-even parity hypothesis fails, since otherwise the cubic with \(\mathbb E s^4=576/25>16\) gives \(h=4\) and strictly improves the bound there.

**Three.** The cubic gives \(h=4\) since \(\mathbb E s^2=108/25>4\). Its strict interval extends beyond \(t=1/3\) because \(1/3<12/31\). In boost coordinates \(t=1/3\) is the midpoint \(r=(\log\varepsilon)/2\). The endpoint \(A=2I+S=B(\log\varepsilon)\) is even unimodular, hence has squared minimum at most 6. If it has minimum 6, repeat the same local argument at the endpoint, moving toward the middle: \(S\) is self-adjoint there and the new pairing \(2s+3q\) is even. The boost distance to that endpoint is at most \((\log\varepsilon)/2\) on the second half. Only the endpoint can then tie on the second half. The two endpoints are equivalent to one another but not asserted equivalent to the original.

If \(A\) instead has squared minimum \(b\le4\), an endpoint witness satisfies \(2q+s=b\). Its image under \(2I+S\) has \(q'=7q+4s\ge6\), whence \(q+q'=4b\) and \(q\le10\). At \(1/3\) this witness has unnormalized squared length \((q+b)/3\le14/3<6\sqrt{2/3}\). It is also strictly below the target at \(1/2\), since \(b/2\le2<3\). Concavity fills the second half strictly.


## Lower-minimum parents and a 96D dichotomy

The preceding scope in \(D=5\) extends to **all** even unimodular parents of ranks 24, 48 and 72, and in \(D=2\) to all even unimodular parents of rank 48. To see this, take a vector with squared norm \(q_0=\min(L)\) in the original metric and write \(z=s/\sqrt D\), \(R=\log\varepsilon\). For \(0<r<2R\), two available squared lengths are \(q_0\cosh r+z\sinh r\) and \(q_0\cosh(r-2R)+z\sinh(r-2R)\), the second coming from the inverse unit applied to the same vector. As functions of \(z\) these have opposite slopes. The maximum of their minimum is their intersection \(q_0\cosh R/\cosh(r-R)\), at most \(q_0\cosh R\); the endpoints also satisfy the bound. Now \(\cosh R=\sqrt5/2\) for \(D=5\), so any nonextremal parent at the stated ranks has \(q_0\le m-2\) and strictly lower density. The same is true in \(D=2\), rank 48, since \(4\cosh R=4\sqrt2<6\).

There is a sharper exact criterion in rank 96 over \(\mathbb Q(\sqrt5)\) when \(L\) has squared minimum 8. The norm-8 shell is a spherical 3-design: its degree-two harmonic theta series has weight 50, and division by \(\Delta^4\) leaves \(M_2(\mathrm{SL}_2(\mathbb Z))=0\). In particular the shell has mean \(s=0\).

1. If any norm-8 vector has \(s\ne0\), there are witnesses of both signs with absolute value at least 2. For \(0\le t\le1/5\) use the negative witness, obtaining \(\gamma\le (8-2t)/\sqrt{1-5t^2}\). Its derivative has the sign of \(40t-2\), so its maximum on that interval is the larger endpoint value \(19\sqrt5/5\). The positive witness supplies the opposite half. This value is below the repository lattice comparator (approximately \(8.6195\) in Hermite units).
2. If instead all norm-8 vectors have \(s=0\), the endpoint \(t=1/5\) has \(\gamma=4\sqrt5\) exactly. The unit \(U=(I+S)/2\) gives \(q'=(3q+s)/2\) and \(s'=(5q+3s)/2\). Thus \(q'\ge8\) and is even. If \(q=10\), then \(s\equiv2\bmod 4\) and \(s\ge-14\). The value \(-14\) would give \(q'=8,\ s'=4\), impossible by assumption, hence \(s\ge-10\) and \(q+s/5\ge8\). For \(q\ge12\), use \(s\ge16-3q\) to obtain \(q+s/5\ge8\). Norm-8 vectors attain 8 in the unnormalized form; applying \(-S\) proves the same at the other endpoint. Every intervening unnormalized minimum is likewise 8 by the affine inequality for each vector. The endpoint center density is therefore \(5^{24}\), and the packing fraction is about \(3.5040\times10^{-21}\), about \(5.9\) times the current lattice lower record quoted in the repository. This is conditional: no parent satisfying the norm-8 cancellation condition is established here.

## References

- E. M. Rains and N. J. A. Sloane, [The Shadow Theory of Modular and Unimodular Lattices](https://neilsloane.com/doc/Me231.pdf), Theorem 2 and Section 7.
- David Dursthoff, [Extremal Lattices and Hilbert Modular Forms](https://publications.rwth-aachen.de/record/685470/files/685470.pdf), chapters 6–7. Documents real multiplication examples on \(P_{48n}\) for both \(\mathbb Q(\sqrt2)\) and \(\mathbb Q(\sqrt5)\), and on \(\Gamma_{72}\) for \(\mathbb Q(\sqrt5)\).
- Nebe–Sloane Catalogue, [the nonuniform 11-modular twist of \(P_{48n}\) by \(4+\sqrt5\)](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/dim4811mod.html).
