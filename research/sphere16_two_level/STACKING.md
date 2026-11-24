# The exact quarter-phase stacking optima survive unequal contractions

Fix any extremal even unimodular lattice L of rank n=48 or 72, min m=6
or 8, and write D=n/12=m-2. Use full horizontal L layers in rank n+d,
with each layer's phase belonging to (1/4)L/L. Only the horizontal
metric is changed, to G=I-A with A positive semidefinite of rank at
most p and condition ratio at most Q. The d-dimensional vertical
space stays Euclidean and orthogonal to the horizontal space.
The contracted eigenvalues and their directions need not agree.

Under the normalized squared-minimum hypothesis mu>=D+d/12, the
center density cannot exceed that of the best uncontracted matched
stack, and any nonzero A gives strict inequality. The same scopes
as in the uncontracted quarter-period theorem apply: arbitrary
periodic stacks for 1<=d<=3, and lattice-only for 4<=d<=5 in rank
48 or d=4 in rank 72. Sample certified Q (all d in the respective range):

| Contracted rank p | Rank 48 carrier | Rank 72 carrier |
|:--:|:--:|:--:|
| 1 | 1150 | 2920 |
| 2 | 45 | 82 |
| 3 | 14 | 22.7 |
| 4 | 7.37 | 11.24 |
| 5 | 4.74 | 6.88 |
| 6 | 3.38 | 4.76 |

For arbitrary extremal carriers these are common upper ceilings,
not attainment for every possible carrier. The uncontracted
examples from the earlier theorem attain them in the stated
dimensions. This does not give a general covering-radius or
unrestricted density bound.

## Vertical clearance under contraction

The previous exact quarter-period theorem gives squared torus
depth at most D for every quarter phase under the original metric.
Under contraction every vector has weakly smaller squared norm,
so all the same order-four phases still have depth at most D.
With mu>D, each vertical location can then carry only one layer
coset; projecting the centers vertically gives min at least mu-D.
If the original set is a lattice its projection is a discrete
lattice, by reducing in the horizontal period; for a periodic
stack its projection is periodic and the fiber count remains one.
By the argument of the earlier theorem, compare with the sharp
vertical packing constant in the required class. If delta_0 is
the matched uncontracted bound at min m, the squared density ratio
is at most

    (mu/m)^(n+d) / [ det G * ((mu-D)/2)^d ].

For fixed G the right side is nondecreasing for mu>=D+d/12 by its
logarithmic derivative. This is why we keep that lower cutoff.

## Endpoint reduction for unequal contraction rates

From the projection moment tests in README.md take the exact
increasing lower witnesses c_1,...,c_p for the maximum squared
projections of a unit direction of the shortest shell onto an
arbitrary j-dimensional space. If T=1-1/Q and alpha_j are the
decreasing eigenvalues of A, define h=max_j(alpha_j c_j)/T. The
horizontal period has squared minimum at most m(1-Th), and

    det G >= product_i (1-T min(1,h/c_i)).

Here 0<=h<=c_p, and the packing minimum hypothesis implies
h<=(24-d)/(12mT). Evaluating the density ratio at the upper bound
m(1-Th) therefore gives

    F(h) = (1-Th)^(n+d) /
           [ (1-mTh/2)^d * product_i (1-T min(1,h/c_i)) ].

On each interval cut out by the c_i and the upper cutoff, multiply
the log derivative of F by the positive factor 1-Th. The result is

    -(n+d)T + d(mT/2)(1-Th)/(1-mTh/2)
             + sum_(i:h<c_i) (T/c_i)(1-Th)/(1-(T/c_i)h).

Each variable summand is strictly increasing in h since m/2>1
and c_i<1. Thus F cannot have an interior maximum on such an
interval. Since F(0)=1, exact strict inequalities at every
positive endpoint show F(h)<1 for nonzero contraction. The
verifier applies this to every listed rank and every d in scope,
including the endpoint from the clearance cutoff when necessary.

Run python -B stacking.py, standard library only. Its exact
inequality checks are small and independent of the three-level
spectral verifier. The needed quarter-phase input and equality
benchmarks are the [prior exact fourth-torsion theorem](https://github.com/alexlegeartis/SpherePackings/blob/6696a479a8a55e440ff8e65f10702702fe7d599c/research/scoped_density_obstructions_round15_20261123/48-72-quarter-torus/README.md).
