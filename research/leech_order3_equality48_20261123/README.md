# Ternary phases cannot attain the 48D one-Leech-period optimum

Suppose a 48-dimensional **lattice** has a pure period √2Λ, where Λ is the Leech lattice, and squared minimum μ≥6 in this normalization. The maximum center density in this class is (3/2)^24. If this density is attained, its projected 24D lattice is Λ/√2 and μ=6, by the one-Leech-period projection argument and the uniqueness theorem for periodic 24D packings. This proof shows that the horizontal phase homomorphism cannot have exponent dividing three. It does not exclude larger torsion exponents, mixed primary components, or nonlattice motifs.

Scale the horizontal coordinates by 1/√2, and scale the vertical coordinates by √2. The lattice is a lift of y∈Λ with horizontal phase F(y)∈Λ_R/Λ. Its squared clearance over y is
```
||y||²/2 + 2 dist(F(y),Λ)².
```
Suppose 3F=0. Identify its image with a subspace W⊂Λ/3Λ by mapping a phase h to the class of 3h. The form q_3(v)=||v||²/2 mod 3 is a nondegenerate quadratic form because Λ is even unimodular. For any shortest y∈Λ, clearance μ=6 forces F(y) to have squared depth 2 (the squared covering radius). A nearest representative of h then has norm²2, so 3h has norm²18 and q_3(3h)=0. The homogeneous quadratic form y↦q_3(3F(y)) vanishes on all shortest vectors.

The verifier constructs a Leech basis from 80 exact Golay-model minimal vectors, with Gram determinant one and exact membership. It then checks that the evaluations on 600 further shortest vectors have full rank 300 among the 300 homogeneous quadratic monomials modulo three. Thus the quadratic form above is identically zero. Nondegeneracy gives dim W≤12, since W is totally isotropic and contained in its orthogonal complement.

The vertical kernel K of F must have squared minimum at least 12: its nonzero vectors can be lifted with horizontal coordinate zero. Its covolume is [Λ:K]=3^(dim W)≤3^12. Rescaling it by 1/√3 gives a 24D lattice of squared minimum at least 4 and covolume at most 1. Leech optimality forces equality, with K similar to √3Λ by uniqueness. In particular K has squared covering radius 6, so every nonzero coset in Λ/K admits a vector of squared norm 4 or 6.

Every phase h in the image has a representative x/3 with q_3(x)=0, so every possible squared length of that phase is a nonnegative multiple of 2/3. A shortest vertical representative of norm 4 or 6 forces the squared depth of every nonzero image phase to be at least 3/2. Since the squared covering radius is 2, the depth must equal 2. However this would make Λ plus all these phases a 24D lattice of squared minimum 2 and covolume 3^(-12). Scaling by √2 would then give a lattice of minimum 4 and covolume (2/3)^12<1, contradicting Leech optimality again.

Run `python verify.py` with NumPy and SymPy. Modulo 3 is the only rank test needed; the replay also checks the same rank modulo 5, 7, 11, 13 and 23.

Sources: [Cohn–Kumar–Miller–Radchenko–Viazovska, sphere packing and periodic uniqueness in 24D](https://arxiv.org/abs/1603.06518) · [Borcherds, Leech lattice and covering radius](https://math.berkeley.edu/~reb/papers/leech/leech.pdf).
