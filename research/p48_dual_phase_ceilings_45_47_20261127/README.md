# Actual phase distances for three small normals of extremal 48D lattices

**Result.** The external 45D, 46D and 47D antipode densities are sharp in the following complete dual-phase families. In any even unimodular rank-48 lattice of minimum 6, take a primitive normal lattice of the stated Gram and let \(L\) be its perpendicular section. Permit an arbitrary subset of the **entire** discriminant group \(L^*/L\) as translates of \(L\), without requiring a simplex, a small-diameter choice of normal representatives, or a subgroup.

| Dimension | Normal Gram | Sharp center density in the family |
| --- | --- | --- |
| 45 | \(8I_3-2J_3\), determinant 128 | \(23^{45/2}/2^{183/2}\) |
| 46 | \(\sqrt3A_2\), determinant 27 | \(13^{23}/3^{93/2}\) |
| 47 | squared length \(6,8,\) or \(10\) | \(47^{47/2}/2^{118}\), using length squared 8 |

The attaining constructions are credited to [Chen, Hu, Li, Wang and Wu (2025), §3](https://arxiv.org/html/2505.02394v1) in 45D and 47D, and [Conway–Sloane (1996)](https://arxiv.org/abs/math/0207182) in 46D. This is not a global sphere-packing bound, and no new density record is claimed.

## Projection and the moment inputs

Let \(N\) be the normal lattice, \(G\) its Gram, and \(L\) its perpendicular section. Even unimodularity and primitivity identify their discriminant groups, with opposite quadratic forms, by projecting vectors of the ambient lattice \(\Lambda\). Both determinants are \(\det G\). For a nonzero class choose a shortest dual normal representative, with label \(a\) and squared norm \(\gamma=a^T G^{-1}a\). For any horizontal vector in the corresponding class we can choose a lift with this normal component by subtracting a vector of \(N\). Its horizontal squared norm is therefore \(m-\gamma\), where \(m\) is the lift's squared norm. In particular \(6-\gamma\) is a lower bound and is **the true distance** if a norm-6 lift exists. If no such lift exists at any shortest representative, but a norm-8 lift exists at one, the true distance is \(8-\gamma\). The zero class is special: two distinct copies of the same phase are redundant, and the within-period minimum is checked separately.

Every shell is an 11-design by Venkov's harmonic-theta argument: for a positive even harmonic degree at most 10, the weighted theta series has weight below 36 and vanishes through q², so division by Δ³ would give a holomorphic modular form of negative weight. Odd moments vanish by central symmetry. The weight-24 ordinary theta series is E₄⁶ − 1440 E₄³Δ + 125280 Δ², from its constant term and its vanishing q and q² coefficients. In particular its shell sizes are 52,416,000 and 39,007,332,000 for squared norms 6 and 8. For a shell of norm \(m\), its projected moments of total degree \(2j\le10\) are the Wick moments for Gaussian covariance \(G\), multiplied by

\[
\frac{\#\Lambda(m)\, m^{j}}{48\cdot50\cdots(48+2j-2)}.
\]

This follows from the rotationally invariant moments of a sphere; odd moments vanish. The companion verifier derives the finite censuses from integer linear algebra and exact rational arithmetic.

## Dimension 46: every class has a norm-6 lift

Choose basis \(u_1,u_2\) of norm 6 and dot product \(-3\), so \(u_1+u_2\) also has norm 6. The six pure normal roots have known labels. Any other norm-6 vector satisfies \(\lvert a\rvert,\lvert b\rvert,\lvert a+b\rvert\le3\), by comparing its distance to those six roots. This leaves 37 sites. For a fixed row \(a\ne0\), interpolate a site with a degree-six indicator for \(a\) and a degree-at-most-five indicator for \(b\) within that row. The 11-design suffices. The row \(a=0\) then follows by subtracting all other rows from the univariate moments in \(b\).

For shortest representatives the census is:

| \(9\gamma/2\) | Number of classes | Norm-6 lift count per representative |
| ---: | ---: | ---: |
| 0 | 1 | 12,309,600 |
| 1 | 6 | 5,365,440 |
| 3 | 6 | 922,185 |
| 4 | 6 | 361,584 |
| 7 | 6 | 16,767 |
| 9 | 2 | 1,656 |

At squared minimum above \(16/3\), only the shortest differences of the triangular dual lattice are available. They form a cycle in the neighbors of one fixed phase, so at most three phases can be used, with minimum at most \(52/9\) if there is more than one. For minimum in \((46/9,16/3]\), a size-three subgroup whose nonzero classes have \(9\gamma/2=9\) limits the motif to at most nine phases, by at most one per subgroup coset. For minimum at most \(46/9\), the trivial bound of 27 phases suffices. The verifier compares each density exactly. The one-phase period has minimum 6 by the zero row census and does not match the three-phase incumbent.

## Dimension 45: one free moment, fixed at the next shell

Write the four normal tetrahedral vectors as \(u_1,u_2,u_3,u_4=-u_1-u_2-u_3\). All have norm 6 and mutual dot products \(-2\); their pair sums have norm 8. For a norm-6 vector outside these eight pure normal roots, the four inner labels satisfy \(\sum a_i=0,\ |a_i|\le3,\ |a_i+a_j|\le4\). Opposite labels have equal counts. The resulting 101 antipodal classes (including zero once) are sufficient for the entire projected shell.

Using the 11-design, linear reconstruction leaves one parameter \(\tau\), defined as the **combined** count at the label pair \(\pm(3,1,-1,-3)\). Nonnegativity gives \(0\le\tau\le368\). The shortest normal representatives have \(8\gamma=q\) in the following table. “Positive” means a lift at every shortest representative; the last two rows may vanish at a parameter endpoint.

| \(q\) | Number of classes | Norm-6 lift count at a shortest representative |
| ---: | ---: | --- |
| 0 | 1 | positive, so the period minimum is 6 |
| 2,4,6,8,10 | 12,6,24,12,24 | positive |
| 12,14,16 | 4,24,3 | positive |
| 18 | 12 | \(368-\tau\) |
| 20 | 6 | \(\tau/2\) |

For norm 8, remove the six pure normal pair sums (and negatives as distinct vectors, six in total). The projection bound then has \(|a_i|\le4,\ |a_i+a_j|\le5\). If \(|a_i|=4\), subtract the appropriate normal root; if a pair sum has magnitude 5, subtract the matching normal pair sum. Each is an exact bijection to a norm-6 slice whose labels are already controlled. After subtracting these completely known boundary counts, the norm-8 problem has the same 101 interior sites and becomes affine in \(\tau\) and a new scalar \(\sigma\), the combined norm-8 count at \(\pm(3,1,-1,-3)\).

The extra equation is a **cross-shell** identity. For a harmonic polynomial \(H\) of degree 12, its weighted theta series has weight 36 and vanishes through \(q^2\), so it is a constant times \(\Delta^3=q^3-72q^4+\cdots\). For any degree-12 homogeneous polynomial, the sums of its top harmonic part are exactly the deviations from the spherical mean on these two shells, since all lower degrees are already a design. Apply this to \(a_1^2a_2^4a_3^6\). Subtract the two spherical means and use the weighted-theta coefficient ratio \(-72\); the verifier obtains

\[
\sigma=1,887,360-40\tau.
\]

The combined norm-8 count at a shortest \(q=18\) representative is \(8,662,080+4\tau-2\sigma\), and at \(q=20\) it is \(\sigma\). Both are strictly positive for every allowed \(\tau\). Thus the two potentially missing norm-6 slices can grant a phase separation of at most \(23/4\) and \(11/2\), respectively. In particular the attractive record mechanism “a far class has neither norm-6 nor norm-8 lift” is ruled out. This step cannot be replaced with a norm-6 check alone.

For a uniform bound it is safe to allow **both** special phase types at their larger distances at once. The resulting graph has clique sizes at most four above \(11/2\) and at most six above \(21/4\), found by a small exact clique recursion with one item fixed at zero. For \(\mu\in(9/2,21/4]\), the eight-element 2-torsion subgroup has all nonzero phase distances at most \(9/2\), so the motif has at most 16 items. Below that, use 128. Density computations at the interval endpoints are exact. Above \(23/4\) only one phase can occur, which also loses. The published four shortest phases attain the displayed 45D ceiling.

As an attribution check, the external [Sun–Wang exact kissing certificates](https://github.com/wcgbg/sphere-packing/tree/main/packings_out/records/dim45-kissing7379838) give 7,379,838 contacts in the particular \(P_{48m},P_{48q}\) embeddings at this density. The projected formula gives \(7,380,720-9\tau\), i.e. \(\tau=98\) for those embeddings. This observed interior parameter is not an assumption in the universal proof.

## Dimension 47: three small normal norms

For a primitive normal \(u\) of squared norm \(t=6,8,10\), the normal labels of norm-6 vectors have magnitude at most \(t/2\), apart from the two pure normals at \(t=6\). After subtracting those known points when necessary, at most six squared label values remain, and degree-10 moments fix their counts. They are all positive. Consequently the exact nonzero phase distances are \(6-k^2/t\), where \(k\) is the shorter cyclic distance modulo \(t\). Checking subsets containing zero in these groups of sizes 6, 8 and 10 gives the displayed sharp ceiling at \(t=8\). Larger normal norms are not covered.

## Remaining construction scope

The proofs concern the displayed fixed normals and dual phases. They do not exclude continuous phases outside the dual group, different normal Gram types or deformation of the period metric. The most promising transferable ingredient is the degree-12 **cross-shell** constraint: a suspected missing projected contact at one shell can force a positive contact count at the next before any full shell enumeration.

## Verified 44D addendum

For the particular published determinant-605 normal realized in the literal P48p stock, [the 44D proof](CHEN44.md) gives exact explicit lifts in all 605 dual classes and a sharp density ceiling for every dual-phase union. This addendum is specific to that embedding.
