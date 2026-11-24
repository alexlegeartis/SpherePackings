# A single skew similarity dualizes every translation-invariant MW92 profile

Let $L$ be [Shioda's MW92 lattice](https://doi.org/10.2307/2374791)
with the pinned Gram $G$, cusp isometry $T$ of order $p=47$, and
$D=1-T$. All lattices are row-coordinate lattices. For every
translation-invariant sublattice $M$ of index $47^{a+b}$,
with local elementary divisors $(a,b)$, **one fixed isometry** $J$
satisfies
$$
                  M J=D^{a+b+7}M^* .
$$
Here $D^kM$ denotes the image by the row operator. In particular
**every translation-invariant index-$47^{16}$ sublattice is
47-modular**, not merely the stem $D^8L$. The same explicit
similarity is available for all of them. These are structure
theorems, not improved density bounds for arbitrary 92D packings.

## Exact input identities

Use the basis of the [pinned Round 17 certificate](../../scoped_density_round17_20261124/92-MW-cusp-jets).
Reconstruct candidate integer matrices $Y,V$ modulo 16,777,213,
then check $GY=D^7$ and $GV=47I$ over the integers. The
determinants are $\det G=47^{14}$ and $\det D=47^2$; hence $Y$
is unimodular. Since the dual coordinate rows are $\mathbb Z^{92}G^{-1}$,
this proves $D^7L^*=L$. Also $V=47G^{-1}$ is an even integral
Gram, with maximal entry magnitude 456,438. Since $47G^{-1}$
is integral, the discriminant group is $(\mathbb Z/47)^{14}$.

The named isometry $U$ from the same input gives $J=U^{23}$.
The checker verifies $J^2=-I$, $J T J^{-1}=T^{-1}$,
and $J G J^\mathsf{T}=G$. Write
$$
     \pi=(1-T)(1+T)^{-1}.
$$
The inverse of $1+T$ is integral (a negative sum of the odd
powers of $T$, because the full cyclic trace vanishes).
Thus $\pi$ is integral with $\pi^*=-\pi$, and
$\pi J=-J\pi$. The factor $(1+T)^{-1}$ and its inverse
preserve every $T$-invariant lattice, so $D^dM=\pi^d M$
for all integral $d$, including negative $d$.

## The mechanism: an alternating, balanced pairing

The completions at 47 are over the Eisenstein ring
$\mathcal O=\mathbb Z_{47}[\zeta_{47}]$, with uniformizer $\pi$,
and $L$ is free of rank two. Suppose first that $M=N$ has
elementary type $(0,k)$. Then
$\pi^k L\subset N\subset L$, and $N/\pi^k L$ is cyclic
over $\mathcal O/\pi^k$. Consider the linking pairing
$$
  \beta_k(x,y)=\langle x J\pi^{-(k+7)},y\rangle\pmod{\mathbb Z}
      \quad (x,y\in L).
$$
Since $J^*=-J$, $\pi^*=-\pi$ and $J\pi=-\pi J$,
the operator $J\pi^{-(k+7)}$ is skew-adjoint.
Thus $\beta_k$ is alternating. It is also **balanced**:
moving $\pi$ from the first argument to the second introduces
one sign from anticommuting with $J$ and another from taking
the adjoint, which cancel.

The radical on $L$ is exactly $\pi^kL$, by $L^*=\pi^{-7}L$
and preservation of stems by $J$. The cyclic submodule
$N/\pi^k L$ is isotropic: for a generator $w$ the values
$\beta_k(\pi^i w,\pi^j w)$ equal
$\beta_k(w,\pi^{i+j}w)$, which is zero. In fact the same
equalities hold over $\mathbb Q$ before reduction for these
monomials, by skew-adjointness and balance. Polynomials in
$\pi$ generate the quotient ring: at 47,
$T=(1-\pi)/(1+\pi)$ and the inverse geometric series
truncates on the quotient. Lifting such a polynomial modulo
$\pi^k$ suffices. Consequently $NJ\subset\pi^{k+7}N^*$.
Both lattices have covolume $47^{k+7}$, so equality follows.

For general type $(a,b)$, use $M=\pi^a N$ of primitive
type $(0,b-a)$ after removing the common power. Here
$(\pi^aN)^*=\pi^{-a}N^*$ because $\pi^*=-\pi$.
The sign when commuting $J$ through $\pi^a$ does not change
a lattice. This proves the claimed identity, also at exponents
beyond the first 46 characteristic-47 digits. The same proof
applies to another prime-order rank-two module with such a
scalar dual stem and an order-four reversing isometry; these
hypotheses must be verified, not assumed.

## Modularity and jet consequences

The [Round 20 Gauss operator](../../scoped_density_round20_20261124/92-mw-gauss-period)
$B$ satisfies $B^2=-47 I$ and sends every invariant lattice
to its depth-23 sublattice with squared scale 47. For
$a+b=16$ the theorem reads $MJ=M^* B$.
The single skew-similarity $K=JB$ obeys $K^2=-47I$ and
$MK=47M^*$. It has squared scale 47, so equivalently
$\sqrt{47}M^*$ is isometric to $M$. For those profiles
$MK\subset M$ also follows directly from $b-a\le16<23$.
More generally if $a+b+7=23j$, $47^{j/2}M^*$ is isometric
to $M$. In particular the index-$47^{16}$ stem $D^8L$
is 47-modular and has exact squared minimum 32 by Round 20.
The depth-16 stem, separately, has the equivalent Gram $47G^{-1}$;
no exact minimum at depth 16 is claimed here.

On the first specialization jet $A$, $J$ sends $(x,y)$ to
$(-y,x)$ over $\mathbb F_{47}$. The exact identity
$D^{38}G\equiv4AA^\mathsf T\pmod {47}$ recovers that plane's
orthogonality. As $-1$ is nonsquare, every line differs from
its perpendicular; when $a<b$, the partner and original
lattice have sum $D^aL$ and intersection $D^bL$.
In particular both Norm-24 witness endpoints of the former
fifth-jet gap keep their own orbits under duality. The
three-digit ratios $(2,27,21)$ and $(6,24,0)$ map to
$(23,15,14)$ and $(39,18,29)$, respectively, via the same
operator $J$; direct pairing and the old $U$ action verify
these four labels. Duality did not itself close these orbits:
see the [independent literal witnesses](../92-mw-complete-index5).

Replay: `OPENBLAS_NUM_THREADS=1 python -B verify.py` with
NumPy. Supply the pinned Round 17 `input.json.gz` as an
optional explicit argument. Acceptance uses bounded integer
matrix arithmetic and exact determinant checks.
