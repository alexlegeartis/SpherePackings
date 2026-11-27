# A binary contact criterion for strict periodic jamming in every supercell

This is a local result, not a global sphere-packing bound. It proves that neither Best's `P10c` nor the fixed 71-coset core of the 11D packing in `README.md` can be densified by any sufficiently small perturbation of its cell and motif, for **each** fixed finite superperiod. The 72nd coset of `P11a` need not be locally jammed; the core still prevents a local change of cell volume. No uniform deformation radius over all superperiods is asserted.

## Criterion

Let `C` be a binary length-`n` code of minimum Hamming distance at least four, and write `P=2 Z^n+C`. For each coordinate `i` form the graph `G_i` on `C` whose edges are weight-four differences involving `i`. Assume:

1. Every pair of distinct coordinates occurs together in the support of a weight-four difference.
2. Each of the graphs `G_i` is connected.

Then in any fixed finite-index sublattice of the `2 Z^n` translation lattice used as a superperiod, `P` is a strict local optimum modulo Euclidean isometry for motif size fixed and squared minimum at least four. In particular, there is no small non-isometric flex at the old volume either.

## Proof

Write the superperiod as `2H`, with `H` a fixed finite-index subgroup of `Z^n`. Use strain Gram matrix `I+G` in the original coordinates, and perturbations `u_c(z)` for each pair `(c,z)` with `c∈C` and `z∈Z^n/H`, after removing a common translation. This parametrizes all nearby cells modulo isometry (polar decomposition), not just axis-aligned cells. An infinitesimal contact constraint for a difference `s` from `(c,z)` to `(d,z+t)` reads

```
s^T G s + 2 s dot [u_d(z+t)−u_c(z)] >= 0.
```

If `s=2e_i` and `d=c`, averaging over `z` shows `G_ii≥0`. The volume derivative is proportional to `tr G`. If it is nonpositive then all `G_ii=0`, and the raw axis-contact inequalities, whose average is zero, show that `u_{c,i}(z+e_i)=u_{c,i}(z)` everywhere.

For a code contact `c,d` the original separation can have any independent signs on its four-coordinate support `S`. For each `s` it is represented by exactly the shift `t_s=(s−d+c)/2` in `Z^n` (zero off the support). Average the contact derivative over all `z`; it becomes `s^T G s+2s dot (mean u_d−mean u_c)` and is nonnegative. The sum of these averages over all sign choices is zero because the diagonals of `G` are zero. Hence every average and therefore every individual original contact derivative is zero. Extracting the `s_i s_j` coefficients on the sign cube gives `G_ij=0` for all pairs co-occurring in `S`; the coverage assumption proves `G=0`.

Now take Fourier modes for displacement functions on the finite abelian group `Z^n/H`. For a character `chi`, let `z_{c,i}` be its Fourier coefficient. The axis contacts say it vanishes unless `chi(e_i)=1`. For each pair `c,d` at Hamming distance four, the code-contact equations at this mode are

```
sum_{i in S} s_i [chi(t_s) z_{d,i} − z_{c,i}] = 0     for all sign choices.
```

Choose an active component `i∈S` with `chi(e_i)=1`, and compare the two sign choices differing only in `s_i`. Since the shift changes by `e_i` up to sign, `chi(t_s)` is unchanged. Subtraction forces `z_{c,i}=chi(t_s)z_{d,i}` for every choice of the other signs. If `S` also contains an inactive `j` with `chi(e_j)≠1`, changing `s_j` in this identity forces `z_{c,i}=z_{d,i}=0`.

For a nontrivial character some generator `e_j` is inactive. Every active `i` shares a support with this `j` by assumption, so there is one zero edge in `G_i`; the equality along every edge of connected `G_i` propagates zero to all its vertices. All nontrivial modes vanish. The trivial character gives equality of the `i`-components along every edge of `G_i`, hence just a common translation, which was removed. This proves there is no nonzero infinitesimal contact-preserving direction with nonpositive volume derivative.

For true local strictness, suppose a sequence of distinct nearby cells/motifs had squared minimum at least four and volume at most the old volume. Remove common translations and rotations as above, divide their real perturbation parameters by their Euclidean norm, and pass to a convergent subsequence on the unit sphere of this finite-dimensional parameter space. The exact contact and volume inequalities give a nonzero limiting direction of the forbidden type. Contradiction. The argument works separately for every fixed finite superperiod; it does not give a size-independent quantitative local radius or rule out new contact topologies far from the old packing.

## Exact checks

In `verify.js`, Best's four seed words and all 10 cyclic shifts give the literal length-10 code, including all 45 coordinate pairs and connected `G_i` of order 40. This verifies the criterion for the known 10D packing with centre density `40/2^10=5/128` and ties the first result to [Best 1980](https://ir.cwi.nl/pub/6831/6831D.pdf). The verifier also checks all 55 coordinate pairs and connected `G_i` of order 71 for the shortened-Witt 11D core. For the full 72-site stock precisely one single vertex fails connectivity, on coordinate 10; an explicit free slide of that coset is proved in `README.md`. Its other contacts do not contradict strict jamming of the core.
