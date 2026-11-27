# Effective-resistance formulation of six Ne64 sheets in 59D

**Result.** The current matched 59D Ne64 construction is the unique best
six-label **diameter certificate** among primitive rank-five normal stocks
with a six-vector obtuse superbase consisting entirely of minimum-norm
vectors. This does not bound the true packing density of an arbitrary
normal section: an exact horizontal minimum can be greater than the
diameter guarantee. Other phase clusters, non-obtuse superbases, varying
normal norms, and noncanonical shifts are open. No new record is claimed.

## The network reduction

Let \(L\) be the 3-modular Ne64 lattice of minimum squared length 12 and
determinant \(3^{32}\). In the normal space use a primitive rank-five
sublattice of \(L^*\) with superbase \(b_0,\ldots,b_5\), where
\(\sum_i b_i=0\), \(b_1,\ldots,b_5\) are a basis,
\(\|b_i\|^2=4\), and \(b_i\cdot b_j\le0\) for \(i\ne j\).
Set \(w_{ij}=-3b_i\cdot b_j\), integers by 3-modularity.
The resulting multigraph on six vertices is 12-regular and its reduced
Laplacian \(H\) (deleting vertex zero) is three times the normal Gram.

Every nontrivial cut has weight at least 12: a cut \(S\) has weight
\(3\|\sum_{i\in S}b_i\|^2\), and the indicated sum is nonzero by the
superbase basis condition. Equivalently, among 12-regular graphs, all
edge weights are at most 6 and all internal triangle weights at most 12;
these are the only two added cut inequalities up to taking complements.
Conversely these cut conditions give a positive definite \(H\), and each
nonzero integer vector of its lattice has norm at least 12: write its
coordinates on the graph with the anchored coordinate zero. Adjacent
integer-coordinate differences satisfy \((x_i-x_j)^2\ge |x_i-x_j|\);
the sum of weighted absolute differences decomposes into the weights
of the intervening level cuts, at least one of which is nontrivial.
Thus this graph space neither omits a permissible normal Gram nor assumes
that every admissible Gram actually embeds in Ne64.

Take the standard antipode labels \(e_0=0,e_1,\ldots,e_5\). Their
squared distances in the dual normal metric are
\(3(e_i-e_j)^{\mathsf T}H^{-1}(e_i-e_j)\), i.e. three times the
effective resistance between the corresponding graph vertices.
Primitivity ensures all of these labels lift to the ambient lattice.
The 59D period Gram determinant is
\(3^{32}\det(H/3)=3^{27}\det H\).
The antipode diameter therefore guarantees the squared minimum
\(\mu=12-\beta\), where \(\beta=3\max R_{ij}\). In the feasible graph
set \(\beta<12\): a positive-weight path between two vertices has at
most five unit-or-better resistors, and the
checker verifies \(3\max R_{ij}<12\) directly on all graphs.
Consequently all six projected translates are distinct, and the quantity
being optimized is precisely
\[
\delta^2_{\rm diam}
= 36\,( (12-\beta)/4)^{59}/(3^{27}\det H).
\]

## Exact bound and what would actually improve the packing

The attached [integer enumeration](verify.js) checks all 55,176
admissible graphs after sorting the incident weights of vertex zero.
At equality there are just three labelings under this convention,
all isomorphic to weights 2 everywhere plus 2 on a perfect matching.
This gives \(\det H=98304\) and \(\beta=7/16\), exactly the parameters
of the live optimum within the stated certificate class.

The superficially promising rounded graph with weights 3 on a six-cycle
and 2 otherwise is not new. Deleting a vertex gives the five-vertex path
Gram in the repository's [older 59D construction](https://github.com/alexlegeartis/SpherePackings/tree/eef6ba7/constructions/dim59).
It has \(\det H=101400\). Its three resistance classes, multiplied by
3, are \(209/520,28/65,45/104\) for cyclic distances one, two and three.
The diameter result is 0.3235% below the current matched construction.

There is a specific geometry-based escape. An embedding of this
cycle normal in which the three opposite-pair labels all lacked ambient
norm-12 vectors would have certified squared minimum
\(12-28/65\), not \(12-45/104\): opposite labels would then have even
ambient norm at least 14. At the same period this would give a
**strict 0.1665% gain** over the current match. The existing embedding
does not have this property, as its exact attained minimum in the
older certificate certifies a norm-12 contact at an opposite pair.
No alternative embedding with the needed exclusions is established.

### Comparators and replay

* [Nebe's Ne64 catalogue entry](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/CQ64b.html)
  gives the 3-modularity and minimum used here (Nebe, *Journal of Algebra*
  199 (1998), 472–498).
* [Current matched construction](https://github.com/alexlegeartis/SpherePackings/tree/eef6ba7/research/ne64_dual_matching59_20261123):
  centre density \(275413.6571763477\ldots\), occupied density
  \(2.646163932713027\ldots\,10^{-12}\).
* The historical proof audits a weaker **derived external** comparison
  \(2^{22}/(9\sqrt3)=269064.72705\ldots\) from
  [Nebe's \(\Gamma_{72}\)](https://arxiv.org/abs/1008.2862) and the
  [classical \(\kappa_{13}\) section](https://github.com/alexlegeartis/SpherePackings/blob/eef6ba7/constructions/dim59/README.md#literature-derived-comparison).
  This is not a printed 59D result of Nebe alone.
* Run node verify.js. It checks both explicit benchmark graphs and exact
  squared-density comparisons, including the conditional opposite-contact
  gain. No search over ambient minimal vectors or numerical eigenvalue test
  enters the bound.
