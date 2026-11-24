# A radical-translation obstruction for the 12D anisotropic ternary mask

At the orthogonal metric \((a,c,t)=(9/5,18/5,0)\), an arbitrary ternary
\(E_6\oplus E_6\) motif invariant under any radical translation with
nonzero vertical component has at most **1,893 sites**. The center
density of such a motif is therefore at most
\(1893/[192(162/25)^3]\), strictly below the
[Coxeter–Todd benchmark](https://doi.org/10.4153/CJM-1953-043-4)
\(1/27\), which requires at least 1,935 sites at this metric.
The proof does not assume the motif is additive or polynomial.
Without the radical-translation symmetry the orthogonal mask remains
open.

## From phase fibres to affine polar labels

Write \(V=E_6/3E_6=\mathbb F_3^6\), \(G\) for the simple-root Gram
matrix, \(R=\operatorname{rad}(G\bmod 3)\), and \(W=V/R\).
Put \(Q(v)=v^TGv/2\bmod3\). The form on \(W\) is nondegenerate
of dimension five. The exact script gives literal lift lengths in
each quotient class: zero for \(0\), 2 for type \(Q=1\), 4 for
type \(Q=2\), and 6 for nonzero type \(Q=0\).
At this metric, with squared packing minimum at least 2, every
nonzero motif difference satisfies the necessary gate
\(d(u)+2d(v)\ge10\), where \(d\) is minimum squared lift length in
a residue of \(V\).

At a fixed actual vertical residue, horizontal sites must all have
the same class modulo \(R\). Otherwise a horizontal difference has a
lift of squared norm at most 8. To check the last statement it
suffices to enumerate the integer lifts of norm at most 8; these
cover all residues except the two nonzero elements of \(R\).
The corresponding fibre size is at most three.

Invariance under \(g=(r_1,r_2)\) with \(r_1,r_2\in R\) and \(r_2\ne0\)
makes both the fibre size \(s_v\in\{0,1,2,3\}\) and its horizontal
quotient label \(f(v)\in W\) well defined on each of the 243 vertical
quotient classes. Let \(S\subset W\) be the classes with \(s_v\ge2\).
Whenever \(v,w\in S\), a desired short representative of the vertical
difference can be selected by choosing vertical sheets. Then the
difference of the two selected horizontal subsets of size at least
two covers all three lifts of the horizontal quotient difference,
because subsets of size two in a group of order three have full
difference set. Therefore:
- \(Q(v-w)=1\) forces \(f(v)-f(w)\) to be nonzero isotropic, since
  vertical depth 2 requires every horizontal difference to have
  depth at least 6.
- \(Q(v-w)=2\) forces \(f(v)\ne f(w)\), using vertical depth 4 and
  horizontal depth zero if the labels agree.

A label has at most nine preimages in \(S\): if two labels agree
their vertical difference is isotropic. Translating such a set to
zero and polarizing shows its span is totally isotropic and has
dimension at most two. Consequently at most 81 vertices of \(S\)
can have labels in any single affine totally isotropic plane.

## A robust \(K_6\) seed and exact propagation

In the graph on \(W\) with edges \(Q(v-w)=1\), the verifier finds
432 six-cliques through zero. The six simple reflections and
negation act transitively on them; translations give transitivity
on all affine six-cliques. For any seed six-clique \(K\), let \(B\)
be \(K\) plus the 30 vertices having at least four root-graph
neighbors in \(K\). This count is verified in exact arithmetic.

Suppose \(|S|\ge195\), so at most 48 quotient vertices are missing.
The whole affine symmetry group is vertex-transitive. Averaging
the weighted missing count \(5|K\setminus S|+
|(B\setminus K)\setminus S|\) over all seeds is at most
\(48(30+30)/243<12\). Choose a seed of weight at most 11.
It must have 6, 5 or 4 existing vertices, with at most 11, 6 or
1 additional omissions in \(B\setminus K\), respectively.

The existing seed vertices have distinct labels and are pairwise
isotropically separated. At least four such points in \(W\)
must lie in, and affinely span, one totally isotropic plane \(P\).
For any other \(u\in S\), if it has four neighbors already known in
\(P\) with pairwise nonsingular differences in the *source* \(W\),
their labels are distinct and hence affinely span \(P\). Orthogonality
of the four isotropic differences \(f(u)-f(v_i)\), by subtraction of
their norm equations, puts \(f(u)-f(v_1)\) in the orthogonal
complement of the direction space of \(P\). That vector is isotropic,
so lies in the direction space: otherwise it would enlarge a
maximal totally isotropic plane. Thus \(f(u)\in P\).

The replay generates every valid four-neighbor witness in \(B\),
closes under this implication inside the occupied locations in \(B\),
then tests every outside candidate once. For each allowed deletion
pattern it finds at least 150 potential locations in \(P\).
Some outside candidates can themselves be missing, but deleting
all 48 gives at least 102 actual vertices in \(P\), more than 81.
Thus \(|S|\le194\), and the full motif count is at most
\(3(3\cdot194+49)=1893\).

The same proof bounds by 194 the number of vertical radical
classes in an arbitrary motif with the following two properties:
all three sheets have size at least two and their horizontal
quotient label is constant across the three sheets. For a putative
1,935-site motif and any radical \(g\) with nonzero vertical
component, deleting all points not in full \(g\)-orbits gives an
invariant motif of size at most 1,893. Thus at least 42 points
lie in partial \(g\)-orbits. In particular the overlap deficit
\(|M|-|M\cap(M+g)|\) is at least 21, since each partial orbit
contributes at most two points and contributes one to the deficit.

## Reproduction

From this directory run:

    g++ -O3 -std=c++17 closure.cpp -o closure
    python -B verify_vertical_radical.py

Standard library only. The C++ test uses 36-bit masks enumerating
the omitted first-shell locations. To cover fewer-than-maximum
omissions it suffices to test the maximum number: with fewer holes,
more first-shell and outside positions can be reached. The output
is an exact finite check, not a numerical optimizer or a catalogue
density change.
