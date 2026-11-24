# Rank seven and eight with a totally singular annihilator

In the notation of [README.md](README.md), the exclusion also holds for rank seven or eight if \(W=C^\perp\) is totally singular for \(\overline Q\). No assertion is made for the other rank-seven or rank-eight gluings.

Let \(w=|W|\) and \(M=L+W/2\). By self-adjointness of \(U\) the annihilator \(W\) is \(U\)-stable. The \(U\)-stable lattice \(M\) has integral norms and minimum at least two. For \(x\in M\), write \(q=x^2,\ s=\langle x,Sx\rangle\). Since both \(2x\) and \(2Ux\) have norms divisible by 4, the identity \((Ux)^2=(3q+s)/2\) gives \(s\in\mathbb Z\), \(s\equiv q\bmod2\). Applying both \(U\) and \(U^{-1}\) also gives \(|s|\le3q-4\) for \(x\ne0\).

Suppose \(K(6)\) is empty. Let \(a,b\) count \(M\)'s norm-two and norm-three vectors, and \(x,y\) count its \(q=2,|s|=2\) type and its \(q=3,|s|=3\) type. The norm-two types are \(s=0,\pm2\); norm three has \(|s|\in\{1,3,5\}\), with \(|s|=5\) count \(2(a-x)\) by the two unit bijections from \((q,s)=(2,0)\). The character sums at norms six and eight yield
\[
 b=25(w-4096)/8-48a,\qquad
 |K(8)|=\frac{39007332000+(w-1)304800+12058624a-98304b}{w}.
\]
The latter coefficients are from signed(12,nlow,nhigh,exp=8) in verify.py.

The quartic \(P=s^2-q^2/5\) is harmonic by the rank, trace and square of \(S\). Its sums over the first two shells of \(M\) are \(T_2=4x-4a/5\), \(T_3=48(a-x)+8y-4b/5\). On the other hand \(F=\Theta_{K,P}\) has \(t\)-order at least 8, and by character averaging and Poisson its inversion is \(w^{-1}\Theta_{M,P}\) (the degree-four phase is \(+1\)), hence has order at least 2 with integer exponents. As a \(T\)-invariant form on \(\Gamma(2)\) of degree 14, it is therefore spanned by
\[
 X^8Y^2(X+Y)^2 Z^2,\qquad X^8Y^2(X+Y)^2 X^2.
\]
The ordinary order-eight coefficients of the basis elements are \((2^{34},0)\). The inverted order-two and order-three coefficients are \((256,256)\) and \((0,-16384)\), respectively. Thus
\[
\begin{split}
 \sum_{v\in K(8)}\langle v,Sv\rangle^2
 &=\frac{2^{34}}{w}(T_2/256+T_3/16384)+(64/5)|K(8)|\\
 &=\frac{2048}{w}(256901120-1295w+122880a+106496x+4096y)>0
\end{split}
\]
for \(w=4^7,4^8\) since all counts are nonnegative. The diagonal then gives a noncanceling norm-eight vector as claimed.

Exact replay: run verify_singular.py with Python 3.
