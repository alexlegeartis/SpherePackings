# A mixed-discriminant alphabet in dimension 18

This is a **fixed-stock result, not a new density record or a global upper
bound**. It analyzes full translates of `K = (sqrt(2) E6)^3` at squared
minimum 4 whose translation classes belong to the product of three copies of
`E6*/2E6`. Unlike a half-period search, this alphabet retains both the binary
quadratic forms and the ternary discrimination of the three factors.

The center-density target is Bierbrauer–Edel's 18D value `(3/4)^9`, hence the
motif must have at least 200 elements: `(3/4)^9 * (8 sqrt(3))^3 =
199.75755495338794...`. Multiply by the unit 18-ball volume to compare
occupied fractions. ([Bierbrauer–Edel, 2000, Theorem 1 and §2](https://doi.org/10.1023/A:1008733715001);
[authors' paper](https://www.yvesedel.de/Papers/Sphere.pdf)).

For this alphabet:

- With no ternary part there can be at most **146** translates, by a
  four-variable character inequality.
- Allowing all ternary labels, but imposing squared minimum 4, there can
  be at most **206** translates, by an exact product character certificate.
- Repeating the same binary code with the three common ternary shifts
  `000,111,222` gives at most **192**; an explicit 192-phase lattice attains
  the bound and has the same density `1/(8 sqrt(3))` as classical `Λ18`.
  The upper proof uses the tripartite stabilizer normal form to handle every
  nonlinear binary code, not just linear codes or subcodes of the witness.
- Even if each binary point may independently carry several **complete**
  ternary lines with all three coordinates changing, the sharp maximum
  remains **192** phases. Thus within this alphabet any improvement requires
  partial lines, or an arrangement not made of those lines at all.

## Exact metric and character polynomials

In the normalization where the roots of `E6` have squared length 2, set
`V = E6/2E6`, `q(u) = (u,u)/2 mod 2`. The group `E6*/2E6` splits as
`V × F3`; a phase in physical coordinates is represented by
`sqrt(2) (u/2 + t w)` for one minimal dual vector `w` generating `E6*/E6`.
The binary labels 0,1,2 respectively mean zero, nonsingular (`q=1`), and
nonzero singular (`q=0`). Enumerating shortest vectors gives the following
table (all entries multiplied by 3):

| ternary difference | zero binary | nonsingular | nonzero singular |
|---|---:|---:|---:|
| zero | 0 | 3 | 6 |
| nonzero | 8 | 5 | 2 |

The unreduced phase distances add over the three orthogonal factors. The
verifier uses the integer `E6` Cartan matrix and exact radicand enumeration
up to `16/3` in both the lattice and nonintegral dual cosets. This finds
all 192 classes, with binary class sizes `1,36,27` in each ternary class.
The covolume of one stock is `8 sqrt(3)`, and in-stock squared minimum is 4.

For characters with binary type `h=0,1,2`, the average of that character on
binary relation `x=0,1,2` is the corresponding matrix entry divided by 27:

```
 27  27  27
 27  -3   3
 27   3  -5
```

Indeed the binary quadratic Gauss sum is −8 (minus type); for a nonzero
character `a`, the transform of `(-1)^q` equals `-8 (-1)^q(a)` by
translation. Subtract the identity, divide by `36` or `27`, and the shown
entries follow. The ternary factor contributes −1/2 if both the character
and relation are nonzero, and 1 otherwise. For three factors take products
and average over coordinate permutations of the relation. Each resulting
average `F_h` is positive definite and `F_h(000)=1`.

Label binary/ternary relations `0,1,2,3,4,5` in the order of the two rows
of the depth table. For a code of size M, let `D_s` denote the number of
ordered pairs with relation shape s, divided by M. Character positivity
gives `1 + sum_s F_h(s) D_s >= 0`. The verifier checks that with the
following nonnegative weights `y_h`, `-sum_h y_h F_h(s) >= 1` for every
allowed nonzero shape (depth at least 12):

| h | y |
|---|---:|
| 001 | 216/17 |
| 003 | 64/221 |
| 005 | 10368/221 |
| 025 | 54675/1547 |
| 055 | 2694141/24752 |
| 235 | 35235/24752 |

Thus `M <= 1+sum y = 24576/119 <207`. For the half-only alphabet use
character types `001,002,022` with weights `72/7,216/7,729/7`; this gives
`M <= 1024/7 <147`.

## Exact repetition maximum, using a symplectic normal form

In a full repetition motif the ternary difference is either zero or all
three nonzero. The depth table says the binary difference must therefore
have unordered type `022` or `112` (for distinct binary points). Each
allowed difference has total quadratic value zero; the span of the
translated binary set is totally singular for `Q = q1+q2+q3`. This is a
minus-type 18-dimensional binary quadratic space, so its singular span
extends to a maximal totally singular space `W` of dimension 8.

The underlying alternating form is symplectic. Extend W to a Lagrangian
`L` of dimension 9. The quadratic form restricts to a nonzero linear
functional on `L` (Witt index 8), so `W = ker(Q|L)`. Regard the three
6-dimensional factors as three groups of three binary Pauli coordinate
pairs. Every binary Lagrangian gives a stabilizer state by choosing
compatible Pauli signs. The [tripartite stabilizer normal form of
Bravyi–Fattal–Gottesman](https://arxiv.org/abs/quant-ph/0504208), stated
also as [Looi–Griffiths, Theorem 7](https://arxiv.org/abs/1107.1761), puts `L`
under local symplectic changes of basis into a direct sum of single-party
generators, Bell pairs between two parties, and GHZ triples. In particular,
the theorem proved there gives an LC normal form, so the maps on the
Pauli coordinates really are symplectic, not arbitrary local linear maps.

There are sixteen shapes after permuting the parties: choose g GHZ triples
and nondecreasing pair counts `e01 <= e02 <= e12`, with each party using
at most 3 coordinate pairs; fill any remaining slots with local generators.
The local coordinate changes need not preserve `q_i`: pull each `q_i`
back through its symplectic change. The possible pullbacks are exactly the
28 forms `q0(x+t)+1`, `q0(t)=1`, where `q0(x)=x0*x1+x2*x3+x4*x5`.
This enumeration permits all combinations of pullbacks, so it cannot miss
a shape needed by a packing.

An exact finite check of `16*28^3 = 351232` cases finds, inside each `W`,
a **2-dimensional subgroup all three of whose nonzero differences are
forbidden** (i.e. not 022 or 112). Therefore each of the 64 cosets of that
subgroup contains at most one permitted binary member, giving `M<=64`.
The checker tests this subgroup certificate, not maximal clique instances.

For attainment put `q(a,b,c,d,e,f)=ab+cd+e+f+ef` and take the 64 binary
triples

```
(x, (a,b+a+c,c,d+c,e,f), (0,e,0,f,a,c)),  x=(a,b,c,d,e,f).
```

The first two maps are injective. The sum of the three quadratic values
is zero, and all three values vanish exactly when `a=c=e=f=0`; then the
third component is zero. Thus the nonzero differences are 022 or 112.
Multiply by the common 3-line to give a 192-phase additive subgroup,
hence a lattice with center density `192/(8 sqrt(3))^3 = 1/(8 sqrt(3))`.
This is a construction at the classical lattice density, not a claim of
isometry to the catalogued lattice. It is below the unrestricted 18D
benchmark.

## Different complete ternary lines

For complete ternary lines (every coordinate changes), test the depth
table on all 3×3 pairs of their ternary labels. Lines must be parallel;
nonparallel full-support directions have no compatible binary relation.
Normalize this direction to `111` by changing signs of independent ternary
coordinates. Two parallel line syndromes can be compatible only if their
shift is zero modulo `111`, or has three distinct coordinates. Thus only
three syndromes can occur, with possible common plane `t0+t1+t2=0` up to
translation. Within one syndrome the binary relations are 022/112 as
above; across the three syndromes they are 000/001/111.

Count lines rather than individual phases. Put `C` for the number of
ordered cross-syndrome pairs divided by M, where M now counts the lines.
Use the same three binary character factors and a 3-letter syndrome factor
(average −1/2 if both difference and character are nonzero). The three
nonnegative test weights `(binary character, syndrome character, weight)`
are `(001,0,216/17)`, `(022,0,27/238)`, `(022,1,13095/238)`.
Checking on the five possible nonidentity relation types gives

```
M + (3/4) C <= 8192/119.
```

Each syndrome bank has size ≤64 from the preceding theorem. The inequality
already gives M≤68. For M≥66, convexity of squared bank sizes gives
`C >= 2*64*(M-64)/M`, and the left side is at least
`66+(3/4)*256/66 > 8192/119`. If M=65 and every bank has size ≤62,
the same argument gives `65+(3/4)*372/65 > 8192/119`. Thus the only
possibility above 64 lines has a bank of size 63 or 64 and an outside line.

**Anchor bound.** Translate the binary point of an outside line to 0.
The large bank's words have type `001` or `111` relative to it. (If a word
equals 0 then any second word of this bank would need a type both within
`022/112` and within `000/001/111`, which is impossible.) Their total
quadratic values are 1, while their pairwise differences have total
quadratic value 0. Thus their vectors, not just their differences, are
pairwise symplectically orthogonal. Extend their span to a Lagrangian L.
The words must lie in the odd coset of `W = ker(Q|L)`.

Run the same 351232 symplectic-normal-form cases. In each, either there
are at most 62 candidate words of the specified types in the odd coset,
or a disallowed 2-dimensional subgroup of W meets the odd candidate set in
at most 62 distinct **cosets of that subgroup**, giving the same capacity.
The checker makes exactly these local quotient statements except in 64
cases, all equivalent on the canonical L to one GHZ triple and two single
generators per party, with q anisotropic on each GHZ pair and trivial on
the single-generator bits. In these exceptions write `x` for the common
GHZ X bit, `zi` for the three GHZ Z bits with even parity, and `ui` for
the two local bits at party i. On this L, `Q=x`. Fixing the three `zi` in
the odd coset fixes `x=1`; differences of candidates at the same `z`
must therefore agree in exactly one of the three four-symbol `ui`
coordinates, by the `022/112` rule.

Such a three-coordinate code over four symbols has size ≤4. Indeed,
normalize a pair to 000,011. A word not starting 0 must have last symbols
01 or 10, with at most one word of each kind; if any such word occurs, no
further word beginning 0 can agree once with it and with 000 and 011. If
all words start 0 their middle symbols are distinct. There are four
possible even-parity `z`, giving at most 16 words in the exception.
This excludes the large bank and proves the claimed sharp bound of
64 lines or 192 phases.

No claim is made here about partially filled lines or continuous torus
phases outside this alphabet.

Replay all exact finite claims with `node verify.mjs` (no packages or
downloads). All density comparisons above concern unit-radius packings.
