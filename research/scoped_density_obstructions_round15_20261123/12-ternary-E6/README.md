# Two sheared ternary E6 masks closed, and an orthogonal rigidity gate (12D)

This extends the [coupled E6 tensor-period result](https://github.com/alexlegeartis/SpherePackings/blob/d44afaa088b1c3b64f22dbb88e767b1d10cb925b/research/scoped_density_obstructions_round14_20261123/12-e6-coupled/README.md). Write the reduced two-slot Gram as
\(B=\left(\begin{smallmatrix} a&t/2\\t/2&c\end{smallmatrix}\right)\), normalize squared minimum to 2, and allow arbitrary motifs in \((E_6/3E_6)^2\) divided by three. The target is the 12D Coxeter–Todd center density (1/27), as in [Cohn's table](https://cohn.mit.edu/sphere-packing/) and [Conway–Sloane's construction](https://neilsloane.com/doc/Me107.pdf).

**Exact exclusions.** The sheared vertices \((3/2,3/2,3/2)\) and \((36/25,81/50,36/25)\) have **exact motif capacity 21** each, below the respective targets 35 and 43. The allowed difference sets at these vertices equal those of the previous cone relaxation (size 5,654 and 5,940). Hence the exclusion applies to all metrics reduced to either mask by that relaxation, not just the two precise matrices. This leaves the two orthogonal masks open for unrestricted motifs. At the anisotropic orthogonal mask \((9/5,18/5,0)\), the following rigidity theorem excludes a class of piecewise-quadratically described nonlinear motifs. No claimed global bound or record change.

## Exact sheared finite graph

Let (G) denote the six simple-root Gram, (u,v) integral lifts of a difference. Such a difference is forbidden precisely if a lift has scaled numerator
\(a u^TGu+c v^TGv+t u^TGv<18\).
Completing either square gives upper bounds on the squared norms of such lift components:
\[
 u^TGu<18/(a-t^2/(4c)),\qquad v^TGv<18/(c-t^2/(4a)).
\]
The inverse-Gram diagonal bounds each coordinate by Cauchy–Schwarz; the checker enumerates all lifts within these bounds for a representative of each orbit under the diagonal Weyl action and joint negation. Thus no lift omissions or numerical tolerances are involved. At either vertex the resulting exact allowed set is contained in the previous cone mask (which used only actual lifts). Counting gives identical set sizes, hence equality.

Every motif may be translated to contain zero. Weyl invariance and negation allow us to fix its second site (e) to one of the admitted orbit representatives. Every other site must be in the common neighbourhood of 0 and (e). Compute the graph there by the allowed differences. To decide whether the graph has a clique of size (k), greedily properly color the current graph, then branch in reverse color order. After eliminating later vertices, each prefix's number of colors is an upper bound for its clique size. Recursing on a vertex intersects the prefix with that vertex's neighbourhood. Exact bitsets implement the recursion. Checking each fixed second point excludes every 22-site motif. The script also verifies both 21-site matching witnesses, including an equality contact.

## Anisotropic orthogonal mask: piecewise-quadratic labels

Let (V=E_6/3E_6), (R=\operatorname{rad}(G\bmod3)), and (W=V/R) with nondegenerate form (Q(x)=x^TGx\bmod 3). The required gate is (d(x)+2d(y)\ge10) for every nonzero difference, where (d) is the minimum E6 squared length among integral lifts. Short-lift enumeration shows:

* (R) has three elements, and every residue outside (R) has a lift of norm at most 8. Thus each vertically fixed fibre has at most three sites, all with one common label in (W).
* For every horizontal (W)-label of nonzero (Q)-value, at least one of its three lifts to (V) has (d\le4).

Fix the (A_4/3A_4) vertical subspace (H) on the first four simple roots of the chain (0,1,2,3,4), with node 5 attached to 2. Suppose within each of the nine affine cosets of (H) the fibre labels where occupied are given by a degree-at-most-two polynomial map (F) into (W). Different cosets can have unrelated polynomials. The same result holds when this description applies to only a subset of sites of a motif: apply the theorem to that subset.

In a fixed coset, let (s(v)\in\{0,1,2,3\}) be the occupancies and (M=\sum_v s(v)). For any vertical root (r\in H), if
(P_r(v)=Q(F(v+r)-F(v))\ne0),
then (s(v)+s(v+r)\le3). Indeed if the sum is at least four, both sheet subsets in (R) are nonempty and their difference covers all three elements of (R). The corresponding pair with horizontal depth at most 4 has total depth (4+2\cdot2=8<10). The polynomial (P_r) has degree at most two. If it is nonzero, Schwartz–Zippel gives at least 27 bad inputs among the 81, and summing the occupancy bound over all translated pairs gives (2M\le27\cdot3+54\cdot6=405). Thus (M\ge203) forces all these polynomials to vanish.

Write (F(v)=A v+S(v,v)+c) with (S) symmetric bilinear, and (J_m(r)=A r+2S(m,r)). The midpoint identity (F(v+r)-F(v)=J_{v+r/2}(r)) shows that for every fixed midpoint (m), the quadratic form (Q(J_m(r))) vanishes on the 20 (A_4) roots. Their evaluation rows span all ten quadratic coefficients over (\mathbb F_3) (checked exactly). Therefore (Q(F(v')-F(v))=0) for every two vertical positions in this coset. By polarization the labels lie in an affine totally isotropic subspace of (W). It has at most nine points since (W) is nondegenerate of dimension five.

Fix one of those nine labels and one of its three horizontal residue representatives. Vertical differences between retained positions must have (d\ge6). Their (A_4) lifts suffice to bound the code by 7: in the zero-sum five-coordinate model, center each coordinate mod 3 and correct a possible sum (pm3) by replacing one corresponding (pm1) with (mp2). The resulting short-lift depth is number of nonzero coordinates plus absolute value of the centered sum. Only profiles ((2,0,3)) and ((0,1,4)) (zero count, then the sorted two nonzero counts) can reach depth 6. Their sizes are 20 and 10. Let (f_1,f_2\ge0) be the normalized autocorrelation averaged under permutations and sign. The characters represented by ((1,0,0,0,0)) and ((1,1,0,0,0)) have twice-real-part sums ((4,-10)) and ((-14,2)); both Fourier coefficients including the identity term (2) are nonnegative. The respective multipliers (15/11,20/11) give (1+20f_1+10f_2\le81/11<8). Consequently (M\le9\cdot3\cdot7=189), contradicting (M\ge203).

The nine cosets therefore contribute at most (9\cdot202=1818) sites. At least 1935 sites are needed to exceed (1/27) at this vertex, so any improving motif must disagree with every such piecewise-quadratic graph at **at least 117 sites**. This does not assume the actual packing extends to a quadratic graph. Without that hypothesis, the orthogonal mask remains open.

### Sharper affine subcase

If every piece is affine, the cap improves to 1089 sites. A bad root test is then constant on the whole coset and bounds its mass by \(\lfloor3\cdot81/2\rfloor=121\). If all roots pass, the linear part has image in an isotropic plane and its kernel has at least nine points. Since the vertical depth-6 code cap is 7, some nonzero kernel difference has depth below 6. On every three-point line in that direction the quotient label is constant and each horizontal sheet can appear at most once. Thus the coset has at most 81 sites in this case; summing the worse cap gives \(9\cdot121=1089\).

## Replay

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
g++ -O3 -std=c++17 -o cliques cliques.cpp
python verify_sheared.py
python verify_quadratic.py
```

Python needs NumPy only. The proof scripts do not use an LP solver.
