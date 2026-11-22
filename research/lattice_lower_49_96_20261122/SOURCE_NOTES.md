# Source notes for the 49–96D lattice catalog

The bounds concern one lattice at a time. A minimum inherited from a parent
is marked as a guarantee unless an attaining vector is checked. The catalog
does not assert that every historical lattice has been located.

## 60D Coxeter–Todd section

Nebe's [Γ72 construction, §4.2](https://arxiv.org/abs/1008.2862) contains
a length-scaled Leech block √2Λ24. The rank-12 Coxeter–Todd lattice of
Gram determinant 3⁶ occurs as the fixed lattice of an order-three Leech
automorphism. See the [Conway–Sloane construction](https://neilsloane.com/doc/Me107.pdf)
and [Scheithauer's fixed-lattice table](https://www2.mathematik.tu-darmstadt.de/~scheithauer/papers/gkm1.pdf).
The Leech block is the intersection of Γ72 with its coordinate subspace: the two complementary binary polarizations force a lone occupied block to be even in integral Leech coordinates. Within this block the fixed sublattice is primitive. Thus Γ72 contains a primitive √2K12 of determinant 2¹²3⁶. In a unimodular lattice a primitive
sublattice and its orthogonal complement have equal Gram determinants.
The 60D complement has squared minimum 8, witnessed in a different Leech block, giving center
density at least 2²⁴/27, stronger than the historical Elkies E60 guarantee
quoted by [Vardy, §4](https://doi.org/10.1090/S0002-9947-99-02169-8). This
is a new derived bound from published inputs, not a density printed by Nebe.

## 59–63D bare Ne64 sections

Delete the 40th row and column of the [published Ne64 Gram](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/CQ64b.html).
The corresponding inverse diagonal entry is 4, which multiplies the
parent Gram determinant 3³² by 4 by the cofactor identity. The section
inherits squared minimum 12 and contains the first basis vector of
norm 12. Its density does not beat the higher reported Elkies E63 floor. The same principal-minor method supplies the alternative 59D, 61D and 62D lattices below the cited prior bounds. [Exact determinant replay](verify_ne64_sections.py).
There is no new 63D record claim.

## 77D Bachoc–Nebe section

In the [Bachoc–Nebe L80 catalogue Gram](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/L80.html), basis indices 1, 0 and 7 (zero-based) have Gram [[8,-4,0],[-4,8,-4],[0,-4,8]]. They span 2A3. The pairing map to Z³ is onto, as the exact HNF shows, so the orthogonal section has determinant 256 and inherits squared minimum 8 from the [published parent](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/extremm22.pdf). The ambient coefficient vector e3−e13 (zero-based) is in the section and has norm 8. This replaces the weaker determinant-300 prior comparator. [Literal normals and witness](published77.json) · [exact replay](published77.py).


## 79D bare MW80 section

The [certified 78D MW80 construction](../packing78_verified_20261117/README.md)
uses two primitive joint normal labels of squared norm 162/41. Each
row alone is primitive and defines a rank-79 kernel. For a primitive
integer label row B and parent Gram G, the section Gram determinant
is det(G) times B G⁻¹ Bᵗ. Therefore the determinant is 162·41¹¹. The
parent squared minimum 14 is inherited from [Schütt–Shioda](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf)
and attained by an original norm-14 basis vector in this kernel. The
center density is (7/2)^(79/2) / sqrt(162·41¹¹), giving an occupied
density of at least 1.10606606099445e−16, about 9.54% above the
previous half-root lattice. The unrestricted two-translate packing
at 1.11920881328318e−16 remains denser. This calculation does not add
its shifts to the lattice. [Literal witness](mw80_79.json) · [exact replay](verify_mw80_79.py).

## 86–91D MW92 bare kernels

For retained primitive dual normal rows B, the kernel of the label map
on MW92 is a single lattice, of Gram determinant 47¹⁴ det(B G⁻¹ Bᵗ)
and minimum 16, attained in each case by a saved parent atom. Use the bare kernel without any translates.
This may have a stronger density than the older lattice control even
though in most dimensions the corresponding periodic packing is denser still. The 91D normal has squared length 182/47 and the section exceeds the preceding unrestricted 91D bound by more than 0.18453%. [Literal labels](section_witnesses.json) · [exact replay](verify_mw92_bare.py).
[Primitive-kernel certification and coarea identity](../packing_release8590_20261117/PROOF.md).
The 89–91 comparison intervals are from the
[explicit Rogers-inequality derivation](../covering_and_deep_cosets_20261115/SECTIONS_PROOF.md);
their upper endpoints only control rounding of the quoted construction
guarantee, not the optimal lattice density.

## 93–96D lattice shears

These are [lattice existence proofs](../mw92_extensions_20261122/README.md#lattice-shears).
The horizontal generators and the individually chosen shear vectors
generate a single lattice; the 96D nine-offset motif is intentionally
excluded. In 96D the published [Bos–Conway–Sloane comparator](https://neilsloane.com/doc/Me103.pdf)
is Table 1 and equation (15), center density 3⁴⁸/2²⁴.
