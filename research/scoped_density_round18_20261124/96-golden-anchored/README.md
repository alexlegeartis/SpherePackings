# Two anchored tests of the 96D golden extremal branch

Assume an even unimodular rank-96 lattice L of squared minimum 10,
with integral self-adjoint U satisfying U²=U+I. Let S=2U−I,
q(x)=x·x, and s(x)=x·Sx. This is the A=0 branch of the
[golden trace analysis](../../scoped_density_obstructions_round15_20261123/96-golden-trace/README.md).
Neither test establishes existence or nonexistence, or changes a density bound.

## Cross-shell constraint: the old endpoint witness fails

There are exactly **11403180502176153600000 ordered pairs of
norm-12 vectors with ordinary inner product +7**, if an extremal
even unimodular 96-dimensional lattice exists, with or without
the golden action. The A=0 witness from
[round 17](../../scoped_density_round17_20261124/96-product-two-point/README.md)
would put about **1.5384994 times that total** in the balanced
norm-12 subshell alone. The old relaxation did not include this
whole-shell budget.

For a direct coefficient calculation, the map (x,y) → (x−y,x)
is a bijection from the pairs above to pairs of norms 10,12
and dot 5. This is coefficient (6,5) of the index-5 Fourier–Jacobi
coefficient J5 of the ordinary degree-two theta series.
Write J5 in the degree-5 monomials of the standard index-one
weak Jacobi generators φ(−2,1), φ(0,1) with elliptic modular
coefficients adjusted to weight 48. Compute these generators by
product expansion and the normalized Jacobi heat derivative,
and use the basis Δ^j E4^u E6^v with v=0,1 in each elliptic
weight. At z=0 the Jacobi form equals the ordinary theta series
times the number of marked shortest vectors. Fourier coefficients
vanish for negative discriminant or when their positive binary
Gram has a nonzero integral vector of squared norm <10, tested
by exact Gauss reduction. These rows determine J5 uniquely in
a dimension-29 basis; symmetry with J5 also determines J6
in its dimension-35 basis and gives all the norm-12 dot counts.
This agrees with [Kominers' genus-two uniqueness theorem](https://arxiv.org/abs/2606.08733).

## A paired shortest-anchor histogram

There are 1,655,821,440 norm-10 vectors of each sign s=10,−10.
Fix x with s(x)=−10. Then Ux is orthogonal to x, also has norm
10, and has s(Ux)=10. For every norm-10 y, use labels
p=x·y, c=Ux·y, s=s(y). Apart from four known vectors ±x,±Ux,
both |p| and |c| are at most 5. The Gram on x,Ux,y,Uy is
the symmetric matrix with rows (10,0,p,c), (0,10,c,p+c),
(p,c,10,(10+s)/2), (c,p+c,(10+s)/2,(30+s)/2).
Require Schur positivity and no positive squared norm below 10
for coefficient vectors in [−2,2]^4. Negation of y, unit
multiplication from s=−10 to +10, and translations
y+u x+v Ux preserving norm 10 give exact equality of the
label counts within their corresponding orbits. Translation
is used for u,v in {−1,0,1} when 5(u²+v²)+pu+cv=0.

The exact remaining inputs are:
- Known signed shell sizes, and covariance within each orthogonal
  eigenspace, from the
  [split-shell designs](../../scoped_density_obstructions_round16_20261124/96-golden-shells/HILBERT.md).
  With b=2c−p, this gives conditional means
  E[5p²+b²]=(500−10s)/48 and E[pb]=(5s−50)/48.
- The ordinary norm-10 shell is an 11-design in R^96, and its
  degree-14 harmonic moments also vanish. The weighted theta
  series is divisible by Δ^5, leaving negative weight for
  positive even degrees up to 10 and weight 2 for degree 14.
  Antipodality handles odd degrees. The replay computes harmonic
  projections of the binary monomials p^i c^(14−i) by Laplacian.
- The doubly marked theta series has diagonal Jacobi index (5,5)
  and weight 48. Its odd-odd part vanishes through q^4.
  Divide it by Δ^5 and apply [Wang–Williams, Theorem 2.4](https://arxiv.org/abs/2102.11190).
  Odd-odd terms contain φ(−1,2) in each variable, then powers
  φ(−2,1)^a φ(0,1)^(3−a) and similarly exponent b in
  the other variable. The remaining elliptic weight is
  2(a+b)−10, so only (a,b)=(2,3),(3,2) occur. Their leading
  Laurent polynomials give two exact odd-odd directions.
  This constraint concerns both dots to the same y.

**Limit.** The 222-orbit anchored model admits a formal nonnegative
rational solution on 45 orbits satisfying all 93 rows. This is
not an integer histogram, not a lattice, and not a relaxation of
the whole genus-three series. We have not imposed the full
even-even sector or consistency across multiple marked pairs.

Run **python -B anchored_replay.py** (standard library) and
**python -B jacobi_replay.py** (SymPy for exact row reduction).
Neither verifier uses a numerical solver.


The replacement distribution in [refit-exact.json](refit-exact.json)
satisfies the round-17 constraints and all eight ordinary whole-shell
dot budgets at A=0. It uses less than the full dot+7 budget, with a
finite nonconstant harmonic margin exceeding 0.2 and tail below 0.0025.
Consequently these ordinary genus-two budget constraints alone do not
close the branch. Replay with python -B refit_replay.py in the repository
layout (it reuses the unchanged round-17 verifier).
