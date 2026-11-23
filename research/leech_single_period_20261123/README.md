# Exact projected packing bounds with one Leech period

Suppose a periodic packing in dimension 24+k is invariant under a horizontal copy H of √2Λ, where Λ is the Leech lattice with covolume 1 and squared minimum 4. There is no assumption on the transverse periods or on their horizontal components. Write μ for the squared minimum of the full packing in this normalization. The period H implies μ≤8. The following exact optima concern these families, not general packings.

* In **48D**, for every periodic packing of this kind with μ≥6, the maximum center density is (3/2)^24, attained by P48n.
* In **25D–36D**, for lattices of this kind with μ≥5, the exact center density is 2^(k/2) times the optimal k-dimensional lattice center density. In 25D–33D its values, in order, are 1/√2, 1/√3, 1/2, 1/2, 1/2, 1/√3, 1/√2, 1, 1. The remaining three are exact reductions to the Hermite constants in dimensions 10–12, not numerical determinations of those constants.
* The same optimal values apply to **arbitrary periodic** packings in 25D–27D and 32D. It suffices to assume μ≥4+1/64 for 25D–27D or μ≥35/8 in 32D. These improvements of scope do not change the density in the unrestricted catalogue.

## Projection inequality with skew transverse periods

Orthogonally project all centers to the complement of H. For any distinct projections the projected squared distance is at least μ−4: take the two original centers and change the horizontal difference by a vector of H so that its square is at most 4, the squared covering radius of √2Λ. If two original centers have the same projection, the same argument implies their horizontal difference belongs to H, since otherwise the residual is a nonzero difference of centers with square at most 4. Thus each projection fiber is a single H-coset, not a denser motif.

In these results μ>4, so the projected center set is discrete and separated. If K is a full period lattice of the original packing, we may take K+H as period group. It is discrete because it acts by translations of a uniformly discrete nonempty set. Its horizontal subgroup is exactly H: any additional horizontal period would give a different coset in the same fiber. Its projected group is a discrete full-rank lattice, because an orbit of one projected center lies in the discrete projected set, and K spans the ambient space. Choosing in K+H lifts of a basis of that projected period gives a block triangular period basis with covolume the product of covolumes. Motifs of the two packings are in bijection because fibers are precisely H-cosets. If the starting center set is a lattice, its projection itself is a lattice.

Let δ_k be the appropriate optimal transverse center density, lattice or unrestricted. If the per-center transverse covolume is V, then V≥((μ−4)/4)^(k/2)/δ_k. The full per-center covolume is 2^12V, so its density divided by 2^(k/2)δ_k is at most
```
R_k(μ) = (μ/8)^((24+k)/2) · (4/(μ−4))^(k/2).
```
The log derivative has the sign of μ−(4+k/6), so the maximum on a closed interval ending at 8 is at an endpoint. R_k(8)=1. Exact rational arithmetic checks the other endpoint for k≤12 at 5, for k≤3 at 4+1/64, and for k=8 at 35/8. 

To be explicit about that last comparison: in k=24 the projected bound δ_24=1 gives full center density at most μ^24/[2^36(μ−4)^12]. Division by (3/2)^24 gives [μ²/(18(μ−4))]^12, equal to 1 at 6 and decreasing up to 8.

## Attainment

For k≤12, take a densest k-dimensional lattice V rescaled to squared minimum 4. Such an optimum exists by compactness up to similarity. In the Leech half-period alphabet, use a 12-dimensional binary subgroup with every nonzero label of squared depth 4, provided by a polarization of Λ; see Nebe §2.2 and Lemma 3.1. Inject V/2V into this subgroup by a binary-linear map. Concatenate each vertical point of V with the corresponding horizontal coset of H. This is a lattice, and its per-center covolume is 2^12 covol(V). A vertical nonzero difference not in 2V contributes at least 4 in both factors; a nonzero difference in 2V contributes at least 16 by itself. Horizontal H differences attain 8. Thus the center density is exactly 2^(k/2)δ_k. In the unrestricted cases listed above the transverse optimum is attained by a lattice. The k=1–8 Hermite constants are classical; the 9D value is 2 by Dutour Sikirić and van Woerden.

In 48D, P48n has period H as an eigensublattice of an involution, determinant 1 and squared minimum 6. See Nebe Theorem 4.10 and Table 2. The 48D part extends the other [Leech scalar-shear result](../leech_scalar_shear48_20261123/) in the high-minimum regime, not its all-minimum half-period statement.

## Sources

- [Nebe: polarized Leech data and three-Leech Γ72](https://arxiv.org/abs/1008.2862)
- [Nebe: 48D extremal automorphism and eigensublattices](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/DIM48.pdf)
- [Borcherds: Leech covering radius](https://math.berkeley.edu/~reb/papers/leech/leech.pdf)
- [Cohn et al.: exact 24D packing](https://arxiv.org/abs/1603.06518) and [Viazovska: exact 8D packing](https://arxiv.org/abs/1603.04246)
- [Dutour Sikirić and van Woerden: exact 9D lattice packing](https://arxiv.org/abs/2508.20719)

The companion replay checks the normalization, endpoint inequalities and numerical lattice center-density values with rational arithmetic.
