# 88D antipode certificates in 83–87D: uniform-loss obstruction

Let \(L\) be an integral self-dual rank-88 lattice of squared minimum 8,
such as either of the [88D lattices of Sakamoto](https://arxiv.org/abs/2609.15633v2).
The result does not apply to a merely covolume-one lattice of unknown
dual minimum. Let \(W\) be a rational normal space of dimension
\(k=88-n\), \(1\le k\le5\), put \(K=L\cap W\), \(T=L\cap W^\perp\), and
take any nonempty finite set \(S\) of distinct projected normal labels.
The usual certification using only the parent minimum and
\(b=\operatorname{diam} S\) guarantees squared separation
\(m=8-b^2\) when positive. The certified center density is
\[
 \delta_{\rm diam}=|S|\left((8-b^2)/4\right)^{n/2}/\operatorname{covol}(T).
\]

| Dimension | Upper cap on \(\delta_{\rm diam}/2^{n/2}\) | Pinned lower benchmark |
|---:|---:|---:|
| 83 | 0.158 | 0.5 |
| 84 | 0.169 | 0.5 |
| 85 | 0.202 | 0.5 |
| 86 | 0.277 | > 0.9 |
| 87 | 0.450 | > 1.9 |

Thus in 85D this uniform-loss certificate always falls below 41% of
the BN85 benchmark. This is a restriction on the **certificate**, not
on the actual packing minimum unless an active diameter difference
is attained by a norm-8 parent vector. Layer-specific bounds are not
limited by this argument.

## Argument

Primitivity of \(K\) implies that the natural map from \(L^*\) onto
\(K^*\) is surjective: an integral functional extends from a primitive
sublattice. Since \(L=L^*\), the projected label lattice is \(K^*\).
Computing covolume in a basis of \(T\) followed by lifts of a basis
of \(K^*\) gives \(\operatorname{covol}(T)=\operatorname{covol}(K)\).

Take an HKZ basis of \(K\) with Gram–Schmidt lengths \(a_i\). Then
\(a_1\ge\sqrt8\) and \(a_{i+1}\ge(\sqrt3/2)a_i\): the size-reduced
projection of basis vector \(i+1\) onto the complement of the first
\(i-1\) vectors has component at most \(a_i/2\) in the \(i\)-th axis
and component of length \(a_{i+1}\) in the next, but its total length
is at least \(a_i\) by minimality. In the reversed dual basis the
Gram–Schmidt lengths are the reversed reciprocals, by the triangular
inverse transpose. Nearest-plane rounding in that basis tiles
by translates of an orthogonal box \(B\) of side lengths
\(\ell_i=1/a_i\), after reversing axes. The interiors are disjoint
because a difference of lattice points has last nonzero triangular
coordinate at least the corresponding full width.

Put \(C=\operatorname{conv} S\). The boxes at the sites of \(S\)
lie in \(C+B\). Adding an axis segment of length \(\ell\) adds
\(\ell\) times the measure of the perpendicular projection, by
Fubini on line fibres. Iteration gives
\[
 |S|/\operatorname{covol}(K)\le
 \sum_{J\subseteq[k]}
 \operatorname{vol}_{|J|}(\operatorname{proj}_J C)
 \prod_{i\notin J}\ell_i .
\]
The zero-dimensional volume is 1 and the formula extends by
continuity to lower-dimensional \(C\). Each projection has diameter
at most \(b\), so its volume is at most \(V_{|J|}(b/2)^{|J|}\),
where \(V_j\) is the \(j\)-ball volume. This is the isodiametric
inequality, also immediate from Brunn–Minkowski on the difference
body. Substitution of the side-length bounds and multiplication by
\((1-b^2/8)^{n/2}\) bounds the normalized certificate density.

The exact replay uses rational side-length and \(\pi\) bounds. On
intervals of width \(1/1000\), evaluate the increasing polynomial at
the right end and decrease the loss factor from its value at the
left end; squaring removes all irrational exponents. Beyond 3 no
positive separation is certified. The 83–85 comparisons use the
[exact BN parity lifts](https://github.com/alexlegeartis/SpherePackings/tree/main/research/packing83_85_deep_verified_20261117);
86–87 use the [MW92 simplex certificates](https://github.com/alexlegeartis/SpherePackings/blob/main/research/mw92_small_simplexes_20261122/verified_results.json).
Even if one artificially replaces 8 by \(8+1/6\) in the scalar
separation formula, keeping the same normal restrictions, the 83–85
caps are 0.380, 0.417, 0.497, all below 0.5.

Replay: python -B verify.py. For independent shell information see
[SHELL.md](SHELL.md).
