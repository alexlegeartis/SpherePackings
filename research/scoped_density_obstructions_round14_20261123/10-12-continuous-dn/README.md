# A continuous D_n-period packing obstruction, and exact inherited-minimum motif sizes

**Scope.** For \(D_n=\{z\in\mathbb Z^n:\sum z_i\equiv0\pmod2\}\), consider any finite motif \(P\subset\mathbb R^n/D_n\), with no rationality, additivity, or grid restriction. Let \(q\) be the packing's squared minimum. In this normalization \(0<q\le2\), and its center density is \(M q^{n/2}/2^{n+1}\), \(M=|P|\). This is not an upper bound on other periods or on generic sphere packing.

- For \(n=10,11,12\), no such motif meets the respective unrestricted benchmark \(5/128,9/256,1/27\) in center density whenever \(q\ge5/4,137/100,59/40\), respectively. Benchmark source: [Cohn’s sphere-packing tables](https://cohn.mit.edu/sphere-packing/) and this repository’s [1–96 table](../../../LOWER_BOUNDS_1_96.md), which uses those literature benchmarks for 10–12.
- In fact the fixed-period upper bounds are strictly below (19/20,99/100,199/200) times those three respective benchmarks. Thus the same exclusion applies to a common linear distortion of the period and motif with Gram condition number at most (1011/1000,501/500,10009/10000), respectively, provided its squared minimum divided by the largest Gram eigenvalue meets the displayed cutoff.
- At inherited squared minimum \(q=2\), the exact motif capacities are \(1\) for \(2\le n\le7\), \(2\) for \(8\le n\le11\), \(4\) for \(12\le n\le13\), \(8\) in 14D, \(16\) in 15D, and \(32\) in 16D.
- For every \(2\le n\) and \(1<q\le2\), the Plotkin inequality below holds. In addition \(M\le2n\) for \(12\le n\le16\) and \(q\ge n/8\). This latter estimate is exact at \(n=12,q=3/2\) and \(n=16,q=2\), by the Hadamard codes of orders 12 and 16. For fixed \(0<q\le2\), the asymptotic statement from elementary torus volume bounds is \(\log M_{\max}=(n/2)\log n+O(n)\), not a useful high-dimensional density bound.

## Coordinate geometry and continuous Plotkin

For a pair of phases let \(r_i\in[0,1/2]\) be its unsigned coordinate difference modulo \(\mathbb Z\), let \(t=\max r_i\), and put \(S=\sum r_i^2\). Rounding to the closest integers gives the squared torus distance \(d=S\) if the integral rounding has even sum. Otherwise its parity can be changed at cost \(1-2t\) at a coordinate attaining \(t\); changing any odd number of coordinates costs at least this. Thus \(d=S+\delta\), where \(\delta=0\) or \(1-2t\). Ties have \(\delta=0\) if \(t=1/2\). If \(d>1\), then \(t>0\). In the odd case,
\[
 \sum_i r_i(1/2-r_i)\ge \frac{1-2t}{2t}\,S\ge 1-2t,
\]
because \(S=d-(1-2t)>2t\). The even case is immediate. Hence \(d\le\sum_i r_i/2\). For each coordinate, the probability that a uniformly rotated semicircle separates that pair is \(2r_i\). Each binary cut separates at most \(\lfloor M^2/4\rfloor\) pairs. Summation proves exactly
\[
 q\binom M2\le \frac n4\lfloor M^2/4\rfloor. \tag{1}
\]
This is the binary Plotkin proof applied to continuous random semicircle cuts, with a bound on the additional parity-sheet distance. In particular for \(q>n/8\), \(M\le q/(q-n/8)\). At \(q=2\), (1) proves the stated caps through dimension 15. Attainment uses half-grid codewords with Hamming distance 8: an 8-one difference, four repeats of the order-four simplex code (length 3), two repeats of the order-eight simplex code (length 7), and the order-sixteen simplex code (length 15). Phase coordinates are word coordinates divided by two; if words differ, parity has zero cost, so distance squared is one quarter the Hamming distance. Appending zero columns preserves the distance.

## A chordal embedding with a second harmonic remainder

Write \(C_0=369/256\). For \(0\le r\le1/2\), and \(113/128\le c<C_0\), the following pointwise bound holds:
\[
 \frac14\sin^2(\pi r)-r^2
 \ \ge\ c r^2(1-2r)+\alpha\sin^2(2\pi r),
 \qquad \alpha=\frac{C_0-c}{40}>0. \tag{2}
\]
For completeness set \(x=2r\). The polynomial \(p(x)=(25x-10x^3+x^5)/16\) is a lower bound on \(\sin(\pi x/2)\) on \([0,1]\). To see this, \(g(y)=\sin(\pi\sqrt y/2)/\sqrt y\) has \(g'''<0\) there: its differentiated sine series alternates, starting negative, with successive absolute term ratios at most \(\pi^2/72<1\). The order-two Taylor expansion of \(g\) at \(y=1\) is \(1+(1-y)/2+(3-\pi^2/4)(1-y)^2/8\); the remainder for \(y\le1\) is nonnegative. Since \(\pi^2<10\), the coefficient is at least \(1/16\), proving the claim. For \(0<x<1\), direct expansion yields
\[
 256\frac{(p(x)/x)^2-1}{1-x}
 =369+x(369-131x-131x^2+19x^3+19x^4-x^5-x^6)
 \ge369+143x.
\]
The polynomial in parentheses has value \(143\) at 1 and negative derivative, bounded above by \(-131-262x+133x^2<0\). Thus the left ratio without \(256\) is at least \(c+(C_0-c)(1+x)\). Euler's sine product gives \(\sin(\pi x)\le\pi x(1-x^2)\) on this interval, so its square is at most \(10x^2(1-x)(1+x)\). These bounds give (2), including the endpoints by continuity.

To absorb \(\delta\), in the odd case use \(d\le n t^2+1-2t\), and choose \(t_0\) such that both \(1\) and \(n t_0^2+1-2t_0\) are less than the required \(q_{\min}\). Convexity forces \(t>t_0\). If \(c(q_{\min}-1+2t_0)\ge1\), the first remainder in (2), summed over coordinates, is at least \(c\delta S\ge\delta\). In the even case no absorption is needed. Therefore
\[
 \frac14\sum_i\sin^2(\pi r_i)\ge d+\alpha h,\quad
 h=\sum_i\sin^2(2\pi r_i). \tag{3}
\]
This is strict relative to \(d\) if any \(r_i\) is strictly between 0 and \(1/2\). The parameters for 10–12 are in the exact verifier. For \(12\le n\le16\) at \(q_{\min}=n/8\), one may take \(t_0=3/10,c=15/16\): \(n t_0^2+1-2t_0<n/8\) and \(c(n/8-2/5)\ge1\).

Embed phases as \(u=(\cos(2\pi x_i),\sin(2\pi x_i))_{i=1}^n/\sqrt n\in S^{2n-1}\). From (3), writing \(b_{jk}=-\langle u_j,u_k\rangle\), define \(a=8q/n-1\), \(\beta=8\alpha/n\), \(e_{jk}=b_{jk}-a\), \(B=1-a\). All pairs satisfy \(0\le\beta h_{jk}\le e_{jk}\le B\).

## Frame coupling and exact rejection tests

Use \(L=\binom M2,\ R=M/2-aL,\ E=\sum_{j<k}e_{jk}\), \(H=\sum h_{jk}\). Positivity of the total Gram sum yields \(0\le\beta H\le E\le R\). Let \(F=\sum_{j,k}\langle u_j,u_k\rangle^2=M+2(a^2L+2aE+\sum e_{jk}^2)\). By discarding off-diagonal coordinate blocks in the frame operator,
\[
 F\ge\frac{nM^2+\sum_{\ell=1}^n |\sum_j e^{4\pi i x_{j\ell}}|^2}{2n^2}.
\]
The identity \(H=(nM^2-\sum_\ell|\sum_j e^{4\pi i x_{j\ell}}|^2)/4\) implies \(H\ge(nM^2-n^2F)/2\). Thus second-harmonic dispersion and frame potential cannot independently be small.

Put \(F_0=M+2[a^2L+(1+a)R]\), \(H_0=(nM^2-n^2F_0)/2\), and \(D=\sum e_{jk}(B-e_{jk})\), so \(F\le F_0-2D\). Since \(e\ge\beta h\), \(D\ge \beta B H-\beta\sum e h\).
1. Using \(h\le n\) gives \(D\ge\beta B H-\beta nR\). Thus \(H(1-\beta B n^2)\ge H_0-\beta n^3R\).
2. Jensen’s inequality gives \(h\le n(1-b^2)\le2n(B-e)\), using \(b\in[-1,1]\). It follows that \(D\ge\beta BH-2\beta nD\). Thus \(H(1-\beta B n^2/(1+2\beta n))\ge H_0\).

To reject either necessary inequality, maximize its left side over \(0\le H\le R/\beta\). This proves the two rational rejection tests in the verifier. At \(q=n/8,n=12,\ldots,16\), \(a=0\), \(H_0=nM(M-2n)/2\), and \(\beta n^2/(1+2\beta n)\ge1\) for \(c=15/16\), so test 2 rejects \(M>2n\). The 12D \(q=3/2\) witness is the rows and negatives of a Paley Hadamard matrix of order 12. The 16D \(q=2\) witness consists of all affine binary functions on \(\mathbb F_2^4\), equivalent to rows and negatives of Sylvester’s order-16 Hadamard matrix.

For bounded counts before running the rational tests: if all inner products are at most \(s=1-8q/n<1/(2n)\), the polynomial \((t-s)(t+1)\) yields \(M\le 2(1-s)/(1/(2n)-s)\). Its off-diagonal values are nonpositive; summing the constant, linear, and quadratic terms gives the bound from Gram positivity and the frame potential lower bound \(M^2/(2n)\). If \(q>n/8\), all off-diagonal products are strictly negative and \(M\le2n+1\). More generally a PSD matrix with nonpositive off-diagonal entries whose negative-edge graph is connected has corank at most one: subtract it from \(T I\) for \(T\) sufficiently large and use Perron–Frobenius simplicity of the maximal eigenvalue of the irreducible nonnegative matrix. At the 10D endpoint \(q=5/4\), any \(M>21\) therefore has a disconnected negative-edge graph; every cross-edge is orthogonal and equality (3) forces a half-integral phase difference. The cross-edges connect all points, giving a common half-grid. Different torus points with the same binary word would be separated by distance squared \(1\), so the motif is an injective binary length-10 code of min Hamming distance 5. Shortening on either value of any fixed coordinate gives a length-9 distance-5 code with at most 10 words by binary Plotkin. Consequently \(M\le20\), contradiction.

The verifier first caps \(M\), then lowers the minimum \(q\) to a rational value no greater than the larger of the cutoff and the minimum needed to match the benchmark. All tests thereafter are exact. Increasing \(q\) cannot make a rejected capacity feasible; no floating tolerances are used.

## Small nonisometric metric neighborhood

For a positive definite metric \(G\), write \(\lambda_{\max},\lambda_{\min}\) for its extreme eigenvalues, and \(q_G,q_I\) for the squared minima of the same period and motif in \(G\) and in the reference Euclidean metric. Then \(q_G\le\lambda_{\max}q_I\). If \(q_G/\lambda_{\max}\) meets the cutoff, the undistorted comparison therefore applies at \(q_I\). The ratio of center densities is
\[
 \frac{\delta_G}{\delta_I} = \frac{(q_G/q_I)^{n/2}}{\sqrt{\det G}}
 \le (\lambda_{\max}/\lambda_{\min})^{(n-1)/2},
\]
using \(\det G\ge\lambda_{\max}\lambda_{\min}^{n-1}\). Multiplying by the exact certified benchmark factors stays strictly below one at the condition numbers listed above; this check is part of the verifier.

## Asymptotic limitation

For \(0<q\le2\), disjoint torus balls of radius \(\sqrt q/2\) have ordinary ball volume because the injectivity radius is \(\sqrt2/2\). Thus \(M\le 2^{n+1}/(v_n q^{n/2})\). A maximal motif covers the torus with balls of radius \(\sqrt q\), each of torus volume no greater than its Euclidean volume, hence \(M_{\max}\ge2/(v_n q^{n/2})\). Stirling gives the stated order of \(\log M\) at fixed \(q\). Neither bound constrains different Euclidean periods.


## Replay

Run `python verify.py` with Python 3. The benchmark and integer capacity rejection tests use only rational arithmetic.
