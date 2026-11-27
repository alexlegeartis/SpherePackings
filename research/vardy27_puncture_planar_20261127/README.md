# A 27-dimensional obstruction to flattening the Vardy binary stocks

**Scope and outcome.** The explicit dual binary stocks used in the repository's
28D Vardy-type construction cannot improve the 27D Vardy centre density
\(1/\sqrt2\) by replacing any three coordinates with an arbitrary real planar
torus, with one planar site per three-bit label in each bank. The result holds
at every radius, with arbitrary codeword deletion and arbitrary horizontal
high-digit parity, and does not require the two banks to obey any particular
cross-bank rule. It is not an upper bound for all 27D packings; the incumbent
remains unchanged.

## Model

Let \(C=\{(x,Ax):x\in\mathbb F_2^{14}\}\) and
\(D=\{(A^{\mathsf T}y,y):y\in\mathbb F_2^{14}\}\), where the columns of \(A\)
are the 14-bit cyclic shifts of the integer 471. Let \(F\) be any three
coordinates and write a word as \((u,f)\) for its remaining 25 and selected
three bits. Bank \(b=0\) uses \(C\), and bank \(b=1\) uses \(D\).

Fix a full-rank planar lattice \(T\) with area \(V\). Independently in each
bank choose a planar site \(z_{b,f}\bmod T\) for each \(f\in\mathbb F_2^3\).
Delete any set of original codewords; for every retained word choose one
whole horizontal coset
\[
 b\mathbf1-2u+4\eta(u,f)+4D_{25},\qquad \eta(u,f)\in\mathbb Z^{25},
 \qquad D_{25}=\{x\in\mathbb Z^{25}:\textstyle\sum_i x_i\equiv0\pmod2\},
\]
at that word's bank/site. In particular a retained word retains an entire
high-digit parity half-fiber of \(4\mathbb Z^{25}\); this does not cover
thinning individual points inside such a coset or adding both parity cosets.
Take all product period translates. The period has covolume
\(2\cdot4^{25}V\), and its horizontal minimum squared length is 32. If
\(\mu>0\) is the squared minimum distance of a nonempty packing, then
\(\mu\le32\) and the minimum squared length of \(T\) is at least \(\mu\).
Planar lattice reduction gives \(V\ge\sqrt3\,\mu/2\): a reduced basis with
shortest first vector \(a\) and second vector \(b\) has
\(\|b\|\ge\|a\|\) and \(|a\cdot b|\le\|a\|^2/2\).

## What the binary stocks force

The attached integer checker verifies for each code that its 546 weight-eight
words realize **all eight masks on every one of the \(\binom{28}3=3276\)
coordinate triples**. Consequently puncture projection to \(F\) is surjective
on either linear code and each of its fibers has exactly \(N=2^{11}\) words.

Fix a bank and distinct labels \(f,g\). A weight-eight word \(w\) having
puncture projection \(f+g\) pairs the entire \(f\)-fiber bijectively with the
\(g\)-fiber by addition. If their retained fractions \(p_f,p_g\) sum to more
than one, some pair survives. Its horizontal support of differences has
size \(8-\mathrm{wt}(f+g)\ge5\). On unchanged coordinates choose the
high-digit lattice translate to give difference zero. On every changed
coordinate the differences \(+2\) and \(-2\) are both available and their
integer translate coordinates differ by one. Hence choosing one sign
appropriately realizes the required total translate parity and the
horizontal minimum squared distance is exactly
\(4(8-\mathrm{wt}(f+g))\), regardless of the two chosen high-digit parities.
The normal torus minimum can be realized independently because the period is
a product.

Write \(A_0=\mu-28\), for now assuming it is positive. If \(f,g\) have
different parity, a torus distance squared below \(A_0\) implies
\(p_f+p_g\le1\), since \(\mathrm{wt}(f+g)\) is 1 or 3.
If they have the same parity, their weight difference is 2 and the true
cutoff is \(\mu-24\ge2A_0\). Thus in this case too, a torus distance squared
below \(2A_0\) implies \(p_f+p_g\le1\).

## Two-colour area bound

If a planar periodic set has two colours, unlike distances at least 1, and
like distances at least \(\sqrt2\), the area of every Voronoi cell is at
least one. To see this, center a radius-\(1/\sqrt2\) disk at one site. Only
opposite-colour sites at distances \(d\in[1,\sqrt2)\) cut positive-area
caps. Set \(\theta=\arccos(d/\sqrt2)\in(0,\pi/4]\). The angular intervals
of half-width \(\theta\) around the cap directions have disjoint interiors:
if the smaller angle between two directions were less than
\(\theta+\phi\), the two opposite-colour sites would be closer than at
angle \(\theta+\phi\), where the squared distance is
\(2\sin^2(\theta+\phi)\le2\) by the cosine rule. This contradicts
their like-colour minimum. Thus \(\sum\theta\le\pi\).

Each cap has area \(g(\theta)=
(\theta-\sin\theta\cos\theta)/2\). The function \(g\) is convex on this
interval, so \(g(\theta)\le(\theta/(\pi/4))g(\pi/4)\). All caps together
have area at most \(\pi/2-1\), and the disk has area \(\pi/2\).
This proves the claim. Scaling and summing over a period, \(k\) sites on
a torus with unlike square minimum \(A_0\), like square minimum \(2A_0\),
and period minimum at least \(2A_0\), force \(V\ge kA_0\).

## Density contradiction

Let \(W_b=\sum_{f}p_{b,f}\), \(W=\max(W_0,W_1)\le8\) and \(t=\mu/32\).
At most \(N(W_0+W_1)\) motif sites are distinct; the radius is
\(\sqrt\mu/2\). Thus the ratio of centre densities to \(1/\sqrt2\) is at most
\[
\frac{N(W_0+W_1)(\sqrt\mu/2)^{27}\sqrt2}{2\cdot4^{25}V}
=\frac{2(W_0+W_1)}{V}t^{27/2}\le\frac{4W}{V}t^{27/2}.
\]
If the ratio exceeded one, \(W>V/(4t^{27/2})\ge4\sqrt3\,t^{-25/2}\).
Therefore \(W>13/2\) and \(t^{25}>3/4\), so \(t>49/50\), by
\(4\cdot49^{25}<3\cdot50^{25}\). In this range \(A_0>0\) and
\(A_0\ge4t^{27/2}\): the logarithmic derivative of
\((\mu-28)/(\mu/32)^{27/2}\) has the sign of
\(378-25\mu/2<0\), and the ratio is 4 at \(\mu=32\).

In the heavier bank draw a bad edge between any two labels violating
the weakened squared thresholds \(A_0\) and \(2A_0\) above (an empty
label can still be assigned its arbitrary site). Every bad edge requires
\((1-p_f)+(1-p_g)\ge1\). Since the total deficit \(8-W<3/2\), there
are neither two disjoint bad edges nor a bad triangle. A pairwise
intersecting family of edges with no triangle is a star: two distinct
edges sharing one endpoint force every other edge through that endpoint
unless it completes a triangle. Remove its center, or any one label if
there are no edges. The seven remaining labels now satisfy the
two-colour lemma, including repeated translates since the period
minimum is at least \(\mu\ge2A_0\). Hence \(V\ge7A_0\).
It follows that \(W>7\); the total deficit drops below one, so there
cannot be any bad edge at all. Applying the area bound to all eight
sites gives \(V\ge8A_0\ge32t^{27/2}\), forcing \(W>8\).

## References and open extensions

* [Vardy, Density doubling of periodic packings and new sphere packings in
  \(\mathbb R^{\!27}\) and \(\mathbb R^{\!28}\)](https://doi.org/10.1090/S0002-9947-99-02169-8),
  *Trans. Amer. Math. Soc.* 351 (1999).
* [Bierbrauer–Edel, A family of 2-groups and packing densities
  obtained from it](https://www.yvesedel.de/Papers/Sphere.pdf), §2 has the
  two-bank construction and both comparator centre densities.
* [Repository's B*-type dual 28D code and 29D stack](https://github.com/alexlegeartis/SpherePackings/tree/main/research/vardy28_stack29_20261122).
* [Pinned project table](https://github.com/alexlegeartis/SpherePackings/blob/a3d7357/LOWER_BOUNDS_1_96.md),
  27D packing fraction \(\approx0.000157594390727866\) (multiply the centre
  density by the volume of the 27-dimensional unit ball).

The code argument in fact holds for any two length-28 binary linear
dimension-14 stocks with the stated weight-eight projection property.
Ways outside this result include more than one normal site per puncture
label in a bank, a mixed horizontal/normal period, a different binary
stock without these shortest-word projections, and internal subcoset
thinning.
