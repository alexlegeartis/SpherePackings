# Distinct backup: rank-two Eisenstein tail

This appendix isolates one condition in a different 96D design. It is not a lower-bound claim. If verified it gives a **lattice**, and its vertical covolume \(573/40=14.325\) would improve the maintained unrestricted bound by about 0.785%.

Let \(R=\mathbb Z[\omega]\) with \(\omega^2+\omega+1=0\), \(N(i)=|i|^2\) and \(\omega\) acting as \(J\) on the torus of MW92. Put
\[
p=5,\qquad r=191/50,\qquad \theta=(2+\omega)/3.
\]
For \(i,j\in R\) use a vertical lattice with squared norm
\[
Q(i,j)=p\,N(i+\theta j)+r\,N(j),
\]
and horizontal offset \(i a+j b\). The initial depth \(d(a)^2>11\) is available by the one-ball covering result. To make \(j=0\) safe we also need \(d((1-J)a)^2\ge1\). The [94D dichotomy](https://github.com/alexlegeartis/SpherePackings/tree/main/research/mw92_paired_cover_9495_20261126) supplies it unless the order-three class needed in the main construction already exists.

The vertical covolume is \((\sqrt3/2)^2pr=573/40\): the shear by \(\theta\) is triangular and the two Eisenstein copies each have fundamental area \(\sqrt3/2\).

For \(j\) of norm 1, rotate by an Eisenstein unit to \(j=1\), which rotates the horizontal difference by \(\pm J^k\) as well. The only \(i\) requiring checks at \(Q<16\) have \(9N(i+\theta)=3,12,21\), in groups of 3, 3 and 6. This is a complete bounded integer check because \(N(x+y\omega)\ge (x^2+y^2)/2\). The strongest group requires horizontal squared depth at most
\[
D=16-r-5/3=1577/150,
\]
the three intermediate ones less than \(D\), and each of six weak ones at most \(77/150\). Use the [triangle partition from the existing 95D repair](https://github.com/alexlegeartis/SpherePackings/tree/main/research/mw92_paired_cover_9495_20261126) to cover the forbidden small balls for each weak offset: three centres form a triangle of circumradius \(4/\sqrt3\), and \(4/\sqrt3+\sqrt{77/150}<\sqrt D\). Thus failure to find a \(b\) satisfying all norm-one-layer conditions would imply an equal-ball cover with centre count at most \(3+3+6/3=8\) per MW92 cell. The existing exact eight-ball square cutoff \(10.51475582>D\) rules this out.

For \(j\) of norm 3 rotate to \(j=1-\omega\). Then \(\theta j=1\). All but \(i=-1\) satisfy \(Q\ge3r+5>16\). This leaves exactly
\[
d((1-J)b-a)^2\ \ge\ 16-3r=227/50.
\]
For \(N(j)=4\), rotate to \(j=2\); \(2\theta\equiv-\theta\bmod R\), so \(Q\ge4r+5/3>16\). For \(N(j)\ge7\) the \(rN(j)\) term alone suffices. The Eisenstein norm has no other positive values below 7.

The proof for the \(b\) in the unit layers does **not** prove the displayed norm-three condition for the same \(b\). Its preimage under \(1-J\) is not an ordinary extra torus ball of the same geometry: the multiplier has degree \(3^{46}\) on the torus. No equal-ball count has been claimed for these extra preimages. Finishing this simultaneous selection would give a genuinely new lattice witness, independent of the nine-site pattern.

