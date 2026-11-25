# One-third-grid product lattices in dimension 12

The comparison is center density \(1/27\) of the [Coxeter–Todd lattice](https://neilsloane.com/doc/Me107.pdf). These are exact exclusions of **specified construction families**, not a new packing, an upper bound for continuous shifts, or an upper bound for arbitrary packings in dimension 12.

- For three independently scaled \(D_4\) factors, no additive motif on the one-third grid improves the comparison. An additive motif is a translate of a subgroup of \((\tfrac13 L)/L\).
- For arbitrary motifs on that grid, if a motif can improve the comparison, then it can do so at the equal squared-scale shape \((9/5,9/5,9/5)\) after normalizing the minimum squared distance to 2 and minimizing covolume with the motif fixed. The residual problem is a nonlinear motif with depth at least 5 and between 645 and 729 points. No such motif is claimed here. A [separate exact test](NONLINEAR_ISOTROPIC_DIAGONAL.md) also closes an invariant two-block phase subfamily and a correlated three-block phase ansatz.
- For four equally scaled \(A_3\) factors, no one-third-grid motif improves the comparison, including nonlinear motifs.

## \(D_4\): distances, scale reduction and Fourier bounds

Coordinate residues identify \(D_4/3D_4\) with \(\mathbb F_3^4\). The centered residue has squared norm equal to its support size \(s\); if \(s\) is odd, adjusting one nonzero coordinate by 3 restores even parity at minimum cost 3. The local depth is half the minimum squared norm, giving the list \(0,2,1,3,2\) for \(s=0,\ldots,4\) and multiplicities \(1,24,24,32\) for depths \(0,1,2,3\). At squared scale \(w_i\), the physical squared quotient contribution is \(2w_i t_i/9\). The period translates have minimum squared length \(2\min_i w_i\) and the covolume is \(8\prod_i w_i^2\).

Normalize the minimum squared distance to 2. For a fixed motif we need \(w_i\ge1\) and \(t\cdot w\ge9\) for each occurring nonzero difference profile \(t\). Capping \(w_i\) by 9 preserves all constraints. Concavity of \(\sum_i\log w_i\) gives a minimizer at a vertex. The verifier enumerates intersections of three independent equalities of type \(w_i=1\) or \(t\cdot w=9\), including \(w_i=9\) via \(t=e_i\); it keeps solutions in \([1,9]^3\), up to permutation. This may overcount vertices, which is harmless for an exclusion. A motif of size \(M\) improves the comparison only if \(27M>512\prod_i w_i^2\).

For each local character class, the average Fourier coefficient on each displacement class, with denominator 32, gives the rows
\[
(32,32,32,32),\quad(32,8,-4,-4),\quad(32,-4,8,-4),\quad(32,-4,-4,5).
\]
The verifier checks this by enumerating all residues and characters, including the equivalent averages over characters at each fixed displacement. For a motif the normalized ordered difference distribution has origin mass one, total mass \(M\), and nonnegative Fourier transform, since each Fourier value is a squared character sum divided by \(M\). Average by the characterization above and by coordinate permutations preserving the allowed profile mask. This averaging needs no invariance of the motif itself. With common denominator \(S\), let \(A_{\chi,t}\) be the resulting integer numerators. Each stored nonnegative dual vector \(W\) satisfies
\[
q=\min_{t\ne0\text{ allowed}}-\sum_{\chi}W_\chi A_{\chi,t}>0.
\]
The nonnegative Fourier average gives \(M\le1+S\sum_\chi W_\chi/q\); integer floors and cross-multiplied density comparisons finish the test. [Standalone replay](verify_d4_scaled.js): 346 shapes, four closed by alphabet size, 341 by the exact duals, and one surviving shape. The replay verifies 10,216 pointwise inequalities over 303 masks and 342 density comparisons. Numerical optimization was used only to propose witnesses, not in verification.

## Exact additive closure at depth five

A subgroup in the residual interval must have \(3^6=729\) elements. View it as the image of a six-dimensional message space \(V\). The three projection kernels \(K_i\) have pairwise zero intersection; a nonzero word supported in just one block cannot reach depth five. The other two projections inject on \(K_i\), with no image of depth 1 (a remaining block contributes at most 3). The local depth comes from the split quadratic form \(-\sum x_j^2\) over \(\mathbb F_3\), with residue zero read as depth 3 except at the origin. A quadratic restriction avoiding residue 1 has rank at most one: diagonalize it, and either a coefficient is 1 or two coefficients 2 sum to 1. In a nondegenerate quaternary space a subspace of dimension at least three has restriction rank at least two; its radical is contained in its orthogonal complement. Thus each \(K_i\) is a plane and each projection has rank four.

Up to message-space basis change the three kernel planes are classified by \(\ell=\dim(K_3\cap(K_1+K_2))\in\{0,1,2\}\): choose bases of \(K_1,K_2\) so the intersection is a partial diagonal, then extend using the remaining basis vectors of \(K_3\). With zero-based basis indices, the third kernel representatives are \(\langle e_4,e_5\rangle,\langle e_0+e_2,e_4\rangle,\langle e_0+e_2,e_1+e_3\rangle\). In each quotient enumerate the 21,060 symmetric four-by-four matrices of determinant 1. These are exactly the nonsingular split forms: after diagonalization, two diagonal coefficients 2 can be changed to two coefficients 1 by sum and difference, and the square determinant means an even number of 2 coefficients. The no-shallow condition on the other kernels leaves 624, 810 or 3,888 forms per quotient.

Reduce the first quotient to orbits under the invertible maps \(I+u f^\top\) whose construction preserves each kernel: for each kernel either \(u\) lies in it or \(f\) annihilates it. Any such subgroup suffices, as the other two forms still range over all admissible choices. There are respectively four, four and five first-form orbit representatives. For each, test all compatible second and third forms, evaluating their total depth on all 364 projective message points. The replay examines 393,984, 554,040 and 13,716,864 form triples and finds none. [Standalone replay](verify_d4_additive.js).

## \(A_3^4\): arbitrary equal-scale motifs

For the integral sum-zero lattice in four coordinates, centered ternary residues sum to 0 or \(\pm3\). In the latter case the three nonzero trits are equal and a correction by \(\mp3\) costs 3 in squared norm. The half-norm depths \(0,1,2,3\) have multiplicities \(1,12,6,8\). Character representatives modulo common translation, coordinate permutation and sign are 0000, 0001, 0011 and 0012. Use their exactly enumerated Fourier averages in the same dual bound, symmetrized over the four factors. The unit-scale covolume is 16. For positive total depth \(T\) the comparison is \(27MT^6\le16\cdot18^6\), and the period minimum forces \(T\le9\). The verified size caps for \(T=1,\ldots,9\) are \(531441,106288,17971,2952,850,298,77,36,17\), respectively; each suffices. [Standalone replay](verify_a3_equal.js).

Each verifier runs with Node.js without packages, e.g. `node verify_d4_scaled.js`.
