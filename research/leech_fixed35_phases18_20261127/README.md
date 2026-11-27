# Exact 18D antipode stock from the rank-six Leech fixed lattice
Date: 2026-11-27. A sharp result about a specified period and phase
alphabet; not a new unrestricted density result.

Let \(\Lambda\) be the Leech lattice, minimum squared norm 4, and take
its rank-six fixed-point sublattice \(K\) in orbit No. 35 of
[Höhn–Mason, Table 1 and §4](https://arxiv.org/pdf/1505.06420).
Let \(L=\Lambda\cap K^\perp\) be the rank-eighteen complementary lattice.
The following theorem permits **any** subset of the torus phases
\(L^*/L\), not just a simplex chosen in the normal dual or a linear code.

**Theorem.** Every union of cosets of \(L\) taken from \(L^*/L\) has
center density at most
\[
\frac6{\sqrt{243}}\left(\frac56\right)^9
=\frac{5^9}{2^8 3^{21/2}}.
\]
This bound is attained by six phases. It equals the value in the
[18D remark of Chen–Hu–Li–Wang–Wu](https://arxiv.org/html/2505.02394v1#S3);
we make no claim that the remark uses this same embedding. It is
strictly smaller than the \((3/4)^9\) of
[Bierbrauer–Edel, Theorem 1 and its 18D application](https://www.yvesedel.de/Papers/Sphere.pdf).
This theorem does not cover arbitrary real off-alphabet shifts, changes
in the period or horizontal metric, or other Leech sections.

## 1. The normal discriminant metric
No. 35 has discriminant symbol \(3^{+5}\), rank six and one class in its
genus. It is \(K=\sqrt3 E_6^*\), where \(E_6\) has root-length-squared 2.
For verification, let \(A\) be the \(E_6\) Cartan matrix displayed in
the [Nebe–Sloane lattice entry](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/E6.b.html).
The even Gram \(B=3A^{-1}\) has determinant \(3^5\) and rank one modulo
3. Its unit block over \(\mathbb Z_3\) can be taken as \([4]\);
the remaining 3-modular block after dividing out 3 has determinant
\(1/4\), a square unit. This verifies the symbol, and the signature
agrees. The polar reverse index \(i^G=1\) and the normalizer index
\(\mathrm{ind}=1\) of No. 35 in the same table imply that every
isometry of its finite discriminant form is induced by a Leech
automorphism preserving \(K\) and \(L\).

The quotient \(K^*/K\) is \(E_6/(3E_6^*)\). In \(E_6\) basis
coordinates label the class of \(v\) by \(Av\bmod3\); the kernel
matches \(3E_6^*\). This is a 5-dimensional vector space over
\(\mathbb F_3\). Its nondegenerate quadratic form is
\(Q(v)=v^\top Av/2\bmod3\). The checker verifies that the minimal
squared **normal** norms in the 243 labels are \(0,2/3,4/3,2\)
with counts \(1,72,90,80\), respectively; the three nonzero classes
of types \(Q=1,2,0\) take norms \(2/3,4/3,2\) in that order. This
bounded search is exact: \(v^\top Av\le6\) implies
\(v_i^2\le6(A^{-1})_{ii}=2B_{ii}\) by Cauchy–Schwarz, and all 243
classes are reached. The nonzero \(Q\)-level sets are orbits since
any isometry between 1-dimensional subspaces mapping one chosen
nonzero vector to another of the same \(Q\) extends by Witt's lemma.

## 2. The Leech design determines every horizontal coset minimum
Self-duality and primitivity give \(\pi_K(\Lambda)=K^*\),
\(\pi_L(\Lambda)=L^*\), with compatible classes matched by an
anti-isometry \(K^*/K\cong L^*/L\). Since \(v\in K^*\) has squared norm
in \((2/3)\mathbb Z\), the projection of a norm-four Leech vector has
squared norm \(2i/3\) for \(0\le i\le6\). Let \(N_i\) count these
vectors. The 196560 minimal vectors of the Leech lattice form a
[spherical 11-design, e.g. Cohn–Kumar, Table 1](https://www.ams.org/journals/jams/2007-20-01/S0894-0347-06-00546-7/S0894-0347-06-00546-7.pdf).
The squared length of the rank-six projection of a uniform unit
vector in 24 dimensions is \(\mathrm{Beta}(3,9)\), so
\[
\sum_{i=0}^6 N_i i^j=196560\cdot6^j(3)_j/(12)_j
\quad (0\le j\le5).
\]
Here \((a)_j\) is a rising factorial. Also \(N_6=54\): such a vector is
in \(K\), and all vectors of \(B\)-norm at most 4 have basis coordinates
of absolute value at most 1 since \((B^{-1})_{ii}=2/3\). The unique
solution (Vandermonde on nodes 0 through 5 once \(N_6\) is fixed) is
\((6480,104976,65610,19440,0,0,54)\); the checker verifies all exact
integer equations. In particular the horizontal lattice has minimum
squared norm 4, since \(N_0>0\).

Fix a nonzero normal quotient class \(c\) and a shortest normal vector
\(v\) in it, of norm \(2i/3\) with \(i=1,2,3\). Every horizontal
vector \(z\) in the glued coset gives a Leech vector \((v,z)\), so
\(\|z\|^2\ge4-2i/3\). For the reverse bound, \(N_i>0\) gives a norm-four
Leech vector whose normal component has norm \(2i/3\). Its class is
nonzero (since \(K\) has minimum 4) and is of the right \(Q\)-type.
Use the orbit transitivity and the lifted automorphisms in §1 to map
that class to \(c\); its horizontal component has squared norm
\(4-2i/3\). Thus the true horizontal torus squared distances on the
three nonzero relations \(Q=1,2,0\) are exactly \(10/3,8/3,2\).

## 3. Sharp packing within the stock
If \(M\) phases are used at squared separation \(\mu\le4\), center
density is \(M(\mu/4)^9/\sqrt{243}\), because
\(\det L=\det K=243\). The exact distance alphabet gives these gates:

* When \(\mu\le2\), \(M\le243\).
* When \(2<\mu\le8/3\), distinct phases cannot differ by an isotropic
  vector. The nonzero-isotropic Cayley graph has degree 80, minimum
  eigenvalue \(-10\), and 243 vertices; Hoffman's ratio bound gives
  \(M\le27\). The checker computes every character sum by representing
  all quotient vectors with first coordinate zero (the remaining
  principal Cartan minor has determinant 4, nonzero modulo 3).
  On the 80 nonzero isotropic vectors a character of frequency \(h\)
  has eigenvalue \((3N_0(h)-80)/2\), where \(N_0(h)\) counts those
  with \(h\cdot v=0\); symmetry makes the other two counts equal.
* When \(8/3<\mu\le10/3\), a phase clique has \(Q(x-y)=1\) for each
  distinct pair. Translate one point to zero. For seven total points
  the other six would have Gram \(2(I_6+J_6)\) for the polarization
  \(B_Q(x,y)=(Q(x+y)-Q(x)-Q(y))/2\), which has rank 6 over
  \(\mathbb F_3\), impossible in rank 5. Hence \(M\le6\).
  A size-six clique is obtained from zero and the five successive
  sums of simple roots along the \(A_5\) chain \(0,2,3,4,5\) in the
  chosen \(E_6\) diagram: every difference is a root of norm 2.
* When \(10/3<\mu\le4\), \(M\le1\).

The three other gates are strictly below the attained six-phase
density, by comparing \(243\cdot6^9\), \(27\cdot8^9\) and \(12^9\) to
\(6\cdot10^9\). The six-phase value is strictly below the 18D
unrestricted benchmark by squaring both positive densities. Replay
the finite checks with Node on the companion verify.js; it uses exact
integer norm, character and moment arithmetic.
