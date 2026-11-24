# Low-rank golden graph gluing cannot produce the high-density 96D twist

Let \(L\) be an even unimodular rank-48 lattice of squared minimum 6 with a self-adjoint \(S\) satisfying \(S^2=5I\) and \(U=(I+S)/2\) preserving \(L\). The inverse of \(U\) is \(U-I\). It makes \(V=L/2L\) a vector space over \(\mathbb F_4=\mathbb F_2[U]\). For an \(\mathbb F_4\)-linear map \(g\) on \(V\) preserving the binary quadratic form \(\overline Q(\bar x)=x^2/2\bmod 2\), consider the even unimodular lattice
\[
G_g=\{(x,y)/\sqrt2:x,y\in L,\ \bar y=g\bar x\}.
\]
Indeed \(g\) is invertible by nondegeneracy of the polar form, and the \(2^{48}\) graph classes extend \(\sqrt2 L\perp\sqrt2 L\), of covolume \(2^{48}\). The quadratic condition makes \(G_g\) even integral. The same \(S\) acts on both components.

If \(\operatorname{rank}_{\mathbb F_4}(g-I)\le6\), then either \(G_g\) contains a norm-six vector or it contains a norm-eight vector with nonzero golden pairing. Thus this family cannot realize the norm-eight cancellation hypothesis giving Hermite invariant \(4\sqrt5\). The same holds with \(I\) replaced by the reduction of any \(S\)-commuting isometry of \(L\), by applying that isometry inversely to the second factor. This does not exclude arbitrary gluing or arbitrary packings.

## Signed counting

Write \(F_j(v)=\sum_{x\in L:x^2=j}(-1)^{\langle v,x\rangle}\) and \(N_j(v)=\#\{x\in v+2L:x^2=j\}\). For every nonzero \(\bar v\in V\),
\[
 F_6(v)=
 \begin{cases}
 -12800+196608N_8(v)+4096N_{12}(v)& v^2=0\bmod4,\\
 -1536+917504N_6(v)+32768N_{10}(v)& v^2=2\bmod4 .
 \end{cases}
\]
In the first case, if \(N_8(v)=N_{12}(v)=0\), then \(F_8(v)=304800\). Run the accompanying verify.py with Python 3 for exact coefficient replay.

Put \(t=e^{\pi i\tau}\), \(X=\theta_2^4\), \(Y=\theta_4^4\), \(Z=X+2Y\). For the signed theta \(f=\sum_{x\in L}(-1)^{\langle v,x\rangle}t^{x^2}\), Poisson summation gives \(f|_{24}S_{\rm mod}=\sum_{w\in L+v/2}t^{w^2}\), where \(S_{\rm mod}\) is modular inversion, distinguished from lattice \(S\). The weight-two slashes act by \(T:(X,Y)\mapsto(-X,X+Y)\) and \(S_{\rm mod}:(X,Y)\mapsto(-Y,-X)\). For integral exponents in the shifted theta, \(f\) is a holomorphic weight-24 form on \(\Gamma(2)\): invariance under \(T^2,S_{\rm mod}T^2S_{\rm mod}^{-1}\) follows from the two expansions, and all three cusps are holomorphic by the expansions and \(T\)-invariance. A weight-\(2d\) form is a degree-\(d\) homogeneous polynomial in \(X,Y\), since these \(d+1\) monomials have distinct leading orders at infinity and span by the valence formula. The \(T\)-invariant weight-24 basis is \(Z^{12-2i}X^{2i}\) for \(0\le i\le6\). Match the first three even coefficients at infinity \((1,0,0)\), and the coefficients of \(t^0,\dots,t^3\) at the inverted cusp \((0,0,N_8(v),N_{12}(v))\). This determines \(f\) by a nonsingular exact linear system.

For half-integral exponents divide \(f\) by \(\theta_3^2\theta_4^2\), whose inversion is \(-\theta_3^2\theta_2^2\) of leading \(t\)-order \(1/2\), and which has no upper-half-plane zero. It has the same sign as \(f\) under \(S_{\rm mod}T^2S_{\rm mod}^{-1}\), and is \(T\)-invariant. Hence the quotient is a holomorphic weight-22 form with invariant basis \(Z^{11-2i}X^{2i}\), \(0\le i\le5\). Multiply the basis by the denominator and match the same first three even coefficients at infinity, together with \((0,N_6(v),N_{10}(v))\) at inverted orders \(1/2,3/2,5/2\). The exact system gives the second formula. The ordinary weight-24 theta with norm-two and norm-four coefficients zero gives \(F_6(0)=52,416,000=4095\cdot12800\) and \(F_8(0)=39,007,332,000\).

## Fixed-space consequence

Set \(C=\ker(g-I)\), let \(W=C^\perp\) in the binary scalar product, and let \(K\) be the full preimage of \(C\) in \(L\). A rank \(r\) over \(\mathbb F_4\) gives \(|W|=4^r\). Character orthogonality gives \(|K(j)|=|W|^{-1}\sum_{v\in W}F_j(v)\). Every nonzero summand at norm six is at least \(-12800\), so \(K(6)\) is nonempty when \(r<6\). At \(r=6\), if \(K(6)\) is empty then every summand reaches its lower bound. Thus each nonzero \(v\in W\) has \(v^2\equiv0\bmod4\) and \(N_8(v)=N_{12}(v)=0\), and all nonzero shifted classes \(L+v/2\) have minimum at least 4. The norm-eight count is
\[
 |K(8)|=(39,007,332,000+4095\cdot304800)/4096=9,828,000 .
\]

If \(P\) is homogeneous harmonic of positive even degree \(d\le6\), its theta on \(K\) is the character average of the weighted signed theta functions of \(L\). It is of weight \(24+d\), with \(t\)-order at least 8. By Poisson its inversion is, up to a common phase and scalar, a sum of weighted theta functions of \(L+v/2\) for \(v\in W\), with integral exponents and \(t\)-order at least 4 (the zero coset's constant vanishes). The same modular argument gives a homogeneous \(T\)-invariant polynomial in \(X,Y\) of degree \((24+d)/2\). The two vanishing orders and \(T\)-invariance imply divisibility by \(X^8Y^4(X+Y)^4\), of larger degree 16, so the form vanishes. Therefore \(K(8)\) is a 7-design (odd harmonics vanish by antipodality). Its mean \(\langle x,Sx\rangle^2\) is \(64/5\) by the spherical fourth moment, using trace zero and \(S^2=5I\). In particular there is a nonzero golden pairing. Every \(x\in K\) yields a diagonal \((x,x)/\sqrt2\in G_g\) with the same squared norm and pairing, which proves the claim.

References: Dursthoff, [*Extremal Lattices and Hilbert Modular Forms*](https://publications.rwth-aachen.de/record/685470/files/685470.pdf), Chapters 5–7 (theta transformation and concrete golden parents); Bannai et al., [*Configurations of Extremal Even Unimodular Lattices*](https://arxiv.org/abs/0706.3082) (weighted theta forms).
