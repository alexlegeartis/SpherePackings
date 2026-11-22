# A 96D lattice from MW92 additive reduction and D4

The occupied density is at least **5.9371093155920008e−22**, a **5.1657%**
improvement over the preceding 96D *lattice* bound and **112.4171%** over
the 96D lattice in [Bos–Conway–Sloane (1982), Table 1](https://neilsloane.com/doc/Me103.pdf).
The stronger unrestricted lower bound remains unchanged.
This is an explicitly specified lattice, not a union of unrelated shifts
and not an optimality claim.

Let k=F47[ω], ω²+ω+1=0, and L be the Mordell–Weil lattice of
E: y²=x³+1+t⁴⁸ over k(t), with Shioda's height pairing.
[Shioda (1991), Theorem 1.1](https://doi.org/10.2307/2374791)
gives rank 92, Gram determinant 47¹⁴, and squared minimum 16 for this
exact equation over any field containing F47². The theorem for this
family and its normalization are also written out as Theorem 5.6 in
[Muraz–Verger-Gaugry](https://www-fourier.univ-grenoble-alpes.fr/sites/default/files/REF_580.pdf).

## Why specializing to zero gains two in the horizontal norm

At all 48 zeros of 1+t⁴⁸ the vanishing order is one (derivative
48t⁴⁷), so each is an irreducible Kodaira type-II cusp. There are no other
bad finite fibres; at infinity (s=1/t), x'=s¹⁶x, y'=s²⁴y give
y'²=x'³+1+s⁴⁸, smooth at s=0. The minimal discriminant has degree 96,
so χ(O_S)=8. There are no reducible-fibre local corrections in the height
formula. In the convention of [Shioda's height pairing](https://doi.org/10.14992/00009994),
every nonzero section P obeys

    ||P||² = 2χ(O_S)+2(P·O) = 16+2(P·O).

This formula also rules out nonidentity torsion, since intersection
multiplicities of two distinct sections are nonnegative.

Take t₀=3+23ω, which satisfies t₀⁴⁸=−1, and use the specialization
homomorphism ρ:L→(k,+) into the smooth part of that type-II fibre.
The homomorphism is reduction in the smooth group scheme at the
minimal Weierstrass fibre (all sections meet its smooth locus). A
finite integral section cannot meet the cusp: both x,y vanishing would
make y²−x³ vanish to order at least two, contrary to the simple zero.
If coordinate functions have poles at t₀, their projective reduction
is O. The smooth cuspidal cubic y²=x³ is the additive group under
parameter x/y, with value 0 at O: its affine smooth points have
coordinates (u⁻²,u⁻³); three collinear parameters sum to zero by
the cubic equation of incidence, and negating y negates u.
Thus ρ(P)=0 implies P meets O at t₀, and every nonzero such P
has squared norm at least 18.

For an explicit section attaining height 16 take P₀=(-1,t²⁴).
It has no poles in any minimal chart, including infinity where
x'=−s¹⁶,y'=1. Hence P₀·O=0. Its specialization is
c=ρ(P₀)=t₀²⁴=4+8ω, since (t₀²⁴)²=−1.
The ordered pair (c,cω) is a basis of the additive group as an
F47-vector space.

## Glue to the rank-four D4 lattice

Use the orthogonal product metric on L⊗R and R⁴. Let
D₄={z∈Z⁴:∑z_i even} with the standard Euclidean metric. Specify
φ:D₄→(k,+) by the following columns in the basis (c,cω), using
the ordered D₄ basis (e₀−e₁,e₁−e₂,e₂−e₃,e₂+e₃):

    (1,0), (26,33), (36,35), (35,37)     mod 47.

The map is onto and exact inspection of all 624 nonzero signed vectors
of D₄ with squared norm at most 16 shows they do not lie in ker φ.
D₄ is even, hence the vertical zero-colour kernel has squared minimum
at least 18.

Set

    M = {(P,z)∈L⊕D₄ : ρ(P)=φ(z)}.

This is the kernel of the surjective difference homomorphism, because
φ alone is surjective. Its index in L⊕D₄ is exactly 47², and its
covolume is 47²·47⁷·2=2·47⁹.

For nonzero (P,z)∈M, there are three cases. If z=0 the horizontal cusp
kernel gives squared norm at least 18. If P=0 the vertical kernel
gives the same lower bound. If both are nonzero, the separate
minimum norms give 16+2=18. The pair (P₀,e₀−e₁) lies in M and attains
18, proving that the minimum is exact.

The center density is (9/2)⁴⁸/(2·47⁹), and occupied density is the
center density times π⁴⁸/48!. The verifier checks the finite-field
data, every vertical case, and the exact density comparison with the
preceding repository lattice and BCS96, with a rational Machin-formula
interval for π.

    python research/mw92_cusp_d4_20261122/verify.py

A separate integer and rational check can be run with
`python research/mw92_cusp_d4_20261122/audit_numbers.py`. Both
scripts check the finite arithmetic; the published MW92 theorem and
the height and reduction arguments above supply the geometric proof.
