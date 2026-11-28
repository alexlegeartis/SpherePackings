# A coupled Eisenstein tail improves the 96D lattice packing

A 96-dimensional **lattice** packing has occupied density at least **5.938937064720705e-22**, improving the verified lattice-only lower bound **5.9371093155920008e-22** by strictly more than **0.030785%**. This does **not** change the stronger [unrestricted periodic result](../mw92_96_refined_tail_20261125/), **6.367793872884866e-22**. The construction supplies a rational rank-four vertical Gram and proves the required horizontal shears exist. It does not give an explicit 96D Gram or an upper bound.

## Parent and known dichotomy

Let \(L=\mathrm{MW}_{92}\) have minimum squared norm \(16\) and covolume \(47^7\), from [Shioda's characteristic-47 construction](https://doi.org/10.2307/2374791) and [Schütt–Shioda, Theorem 13.35 and §13.4.2](https://link.springer.com/book/10.1007/978-981-32-9301-4). Its orthogonal integral automorphism \(J\) satisfies \(I+J+J^2=0\). Let \(T=I-J\) and \(d(x)=\operatorname{dist}(x,L)\). The [audited threshold](../mw92_cfr_taylor13_20261123/) \(A=11.00098526\) means a periodic equal-ball covering of \(\mathbb R^{92}\) is impossible at density \(\operatorname{vol}(B_{92}(\sqrt A))/47^7\), by the [Coxeter–Few–Rogers theorem](https://doi.org/10.1112/S0025579300002059). In particular there is \(a_0\) with \(d(a_0)^2>A\).

Write \(s=16-A\), \(t_0=16-3s\). The [order-three dichotomy](../mw92_paired_cover_9495_20261126/) gives one of two branches:

* **Generic:** \(a=a_0\) with \(d(Ta)^2\ge t_0\); use \(\sigma=s\) below.
* **Fixed:** If \(d(Ta_0)^2<t_0\), let \(h=T^{-1}e\), where \(e\) is a shortest representative of \(Ta_0\bmod L\), and put \(a=u=a_0-h\). Then \(Tu\in L\) and \(\|h\|^2<t_0/3\) since \(T^\top T=3I\). For a closest lattice point to \(u\), its orbit under \(v\mapsto Jv+Tu\) gives three equally near displacements whose sum is zero. Averaging their squared lengths after adding \(h\) yields \(d(u)^2>A-t_0/3=32/3\). Because \(T(u-v)\in L\) for \(v\in L\), \(3d(u)^2\) is an even integer, so \(d(u)^2\ge34/3\). Also \(3u\in L\) by \((I-J^2)T=3I\), and \(u\) has exact order three. Use \(\sigma=16/3\) in this branch.

## Rank-four tail and short-vector reduction

Take \(R=\mathbb Z[\omega]\), \(N(i+j\omega)=i^2-ij+j^2\), \(\theta=1-\omega\). In the orthogonal vertical space give a pair \(x,z\in R\) the squared length
\[
Q(x,z)=\sigma N(x+z/\theta)+\eta N(z),
\qquad
\eta=\begin{cases}514/125&\text{generic},\\387/100&\text{fixed}.\end{cases}
\]
The horizontal shear is \(F(x,z)=x(J)a+z(J)b\), with \(b\) chosen below and \(x(J)=iI+jJ\) for \(x=i+j\omega\). These define a single lattice in dimension 96: take the four vertical basis vectors with these real horizontal shears, and the rank-92 lattice \(L\) in the horizontal subspace. For explicit rational block matrices on the coordinates of \(x,z\), let \(H=\begin{pmatrix}1&-1/2\\-1/2&1\end{pmatrix}\) and \(S=\frac13\begin{pmatrix}2&-1\\1&1\end{pmatrix}\), since \(z/\theta=(2+\omega)z/3\). The vertical Gram is
\[
G=\begin{pmatrix}
\sigma H & \sigma HS\\
\sigma S^{\mathsf T}H & \sigma S^{\mathsf T}HS+\eta H
\end{pmatrix},
\quad
\det G=\frac9{16}\sigma^2\eta^2.
\]
This determinant follows from the block-triangular coordinate change \((x,z)\mapsto(x+Sz,z)\), of determinant one. The full lattice covolume is \(47^7\sqrt{\det G}\).

Every difference needs \(Q(x,z)+d(F(x,z))^2\ge16\); with \(x=z=0\), the horizontal lattice itself has minimum 16. Multiplication of both coordinates by the same unit of \(R\) preserves \(Q\) and applies one of the six orthogonal transformations \(\pm J^i\) to \(F\). Thus in analyzing \(z\neq0\) we may replace \(z\) by a chosen associate.

* If \(z=0\): norm-one \(x\) gives horizontal depth \(d(a)^2\). At norm three the horizontal depth is \(d(Ta)^2\), compensating \(3s\) in the generic case; it is unnecessary in the fixed case since \(3\sigma=16\). All larger norms give \(Q>16\), since both scales exceed 4.
* Norm \(N(z)=1\) reduces to \(z=1\). Write \(x=i+j\omega\). Then \(9N(x+1/\theta)=N(3i+2,3j+1)\), whose minimum is 3. The largest squared horizontal demand among these differences is \(D=16-\eta-\sigma/3\). The finite shell check is below.
* Norm \(N(z)=3\) reduces to \(z=\theta\). The then-minimal vertical direction is \(x=-1\) with \(Q=3\eta\), demanding \(d(Tb-a)^2\ge W=16-3\eta\). All other \(x\) add at least \(\sigma>W\) to the vertical norm.
* At \(N(z)=4\) we may use \(z=2\), since both integer coefficients must be even (mod-two norm of any nonzero pair is one) and \(z/2\) is a unit. The projected coset \(2/\theta+R=-1/\theta+R\) has squared norm at least \(1/3\), so \(Q\ge4\eta+\sigma/3>16\). Finally \(N(z)\ge5\) gives \(Q\ge5\eta>16\).

These cases are exhaustive: norm two is impossible since \(N(i+j\omega)\equiv(i+j)^2\bmod3\); for norm three this same congruence implies divisibility by \(\theta\), which has norm three and generates the kernel of \(i+j\omega\mapsto i+j\bmod3\).

## Choosing the second shear by moving several candidates together

Select a list of torus classes \(c_1,\ldots,c_m\) (explicitly below), choose any real representatives, and define \(\mu_i=T^{-1}c_i\) in the real parent space. We will choose \(b\) so **all** options \(b_i=b+\mu_i\) satisfy the norm-one conditions, then choose one option satisfying the norm-three condition as well.

Fix one \(b_i\). Its forbidden norm-one phases are centred at \(-x(J)a\), with squared radius \(16-Q(x,1)\) when positive. In choosing \(b\) rather than \(b_i\) these centres are simply translated by \(-\mu_i\). For a forbidden weak coset we use the following actual grouping, not an unequal-radius covering estimate: choose a minimum vector \(v_0\in L\) of length 4, take an index-two kernel sending the primitive vector \(v_0\) to 1 mod 2, and pair the centres along that translate of \(\{0,v_0\}\). Each pair of radius-\(r\) balls is contained in the ball of radius \(2+r\) at its midpoint. We intersect the period sublattices of all cluster types to count the final ordinary covering density.

* **Generic branch:** for \(z=1\), the signed shell labels \(N(3i+2,3j+1)=3,12,21\) have exact counts \(3,3,6\), respectively; no other \(x\) has \(Q<16\). The replay is complete because \(N(v,w)\ge(v^2+w^2)/2\) excludes any \(|i|\) or \(|j|\ge4\). Enlarge the first six labels to radius \(\sqrt D\). The six weakest have squared radius \(r^2=16-\eta-7\sigma/3\), and \(2+r<\sqrt D\) exactly. Use one midpoint for each two of their centres; their combined charge is 3. Hence nine equal balls per \(L\)-cell cover the entire forbidden union for one option. Choose \(c_1,c_2,c_3=0,a,-Ja\). With three options the covering count is at most 27, and \(27D^{46}<A^{46}\). If no \(b\) passes all norm-one conditions, these coverings would contradict the audited equal-ball threshold.
* **Fixed branch:** all \(x(J)u\equiv(i+j)u\bmod L\), so at a given \(b_i\) there are only three phase labels. Each has demand at most \(D=9317/900\). The checker finds only the 3 and 12 shells with counts 3 each. Choose a further class \(v\) with \(d(v+ju)^2>\beta=10.74\) for \(j=0,1,2\), possible by the three-full-ball inequality \(3\beta^{46}<A^{46}\). Take five target centres \(c_i=0,u,-u,v,v+u\). This gives at most fifteen balls per cell, and \(15D^{46}<A^{46}\), so some \(b\) passes all norm-one conditions.

Now suppose the \(m\) options simultaneously failed \(d(Tb_i-a)^2>W\). Choose a shortest real representative \(y_i\) of \(Tb+c_i-a\bmod L\), all with \(\|y_i\|\le\sqrt W\). Then
\[
  \sum_{i<j} d(c_i-c_j)^2
  \ \le\sum_{i<j}\|y_i-y_j\|^2
  =m\sum_i\|y_i\|^2-\Big\|\sum_i y_i\Big\|^2
  \le m^2 W.
\]
But in the generic branch the three pairwise torus distance squares all exceed \(A\) because \(I+J=-J^2\); exactly \(3A>9W\). In the fixed branch four within-group pairs contribute at least \(34/3\) each and six cross pairs exceed \(\beta\); exactly \(4(34/3)+6\beta>25W\) with \(W=439/100\). These contradictions give a single option passing both requirements. Apply the unit symmetry to it, giving full minimum squared norm 16.

## Density comparison and replay

The generic and fixed squared tail covolumes are
\[
D_{4,g}=\frac{37137921969123580223529}{156250000000000000000}
\quad\text{and}\quad
D_{4,f}=\frac{149769}{625}=239.6304.
\]
The larger is the fixed determinant. The [previous lattice certificate](../mw92_cusp_d4_20261122/) had minimum squared norm 18 and covolume \(2\cdot47^9\). To compare it at minimum squared norm 16, scale all dimensions by \(\sqrt{8/9}\), so its covolume divided by \(47^7\), squared, is \(4\cdot47^4(8/9)^{96}\), about \(239.77796396\). It is strictly above \(D_{4,f}\). The new radius-2 occupied density is at least
\(\pi^{48}2^{96}/(48!47^7\sqrt{D_{4,f}})\).
For the named historical external comparator, [Bos–Conway–Sloane (1982), Table 1 and Eq. (15)](https://neilsloane.com/doc/Me103.pdf), center density \(3^{48}/2^{24}\), the gain from the same formula is over **112.482560%**. This comparison does not establish literature priority beyond the stated benchmarks.

Run `node verify.js` for exact tail-shell bounds and counts, strict covering and variance inequalities, determinant ordering and normalized comparison, and the occupied-density floor via rational Machin bounds on \(\pi\). The [inherited analytic verifier](../mw92_paired_cover_9495_20261126/verify.py) checks the cited threshold \(A\) directly. Files for the prior unrestricted periodic packing stay unchanged.
