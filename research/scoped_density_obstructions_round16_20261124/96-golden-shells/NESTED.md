# Degree of the tight-edge graph inside the norm-12 local layer

Use the [golden cancellation hypotheses and local notation](README.md).
For fixed balanced norm-12 x, the vectors r=x-y for balanced norm-12
y with x.y=8 form a set of size <11,040,000. In that set the graph
of edges r.t=4 has **degree less than 29,700** at every vertex.
This is a bound on a local subgraph, not a new density bound.

Indeed, fix r and map t to w=r-t. Then q(w)=8, s(w)=0, and
r.w=4, r.Sw=0, x.w=x.Sw=0. The two transformed off-diagonal
pairings vanish by cancellation applied to r-t and to r,t.
The span of x,r in either factor has dimension two, since their
factor squared norms and dot are 6,4,2 respectively. In each factor
subtract (3/5)(r-x/3) from w (taking the components in that factor).
The residual is orthogonal to x,r in that factor with squared norm
14/5, and the dimension there is 46. For two distinct w,w' write
p=w.w', b=w.Sw'. The normalized product-sphere dots are
(5p-12 ± sqrt(5)b)/28.

As in the previous proof, p is an integer in [-4,4], p=b mod 2,
and |b|<20. Opposites cannot occur since r.w=4. The symmetric
kernel permits requiring b>=0. The exact filter checks Cauchy
and the same minimum and cancellation constraints on the integral
span of w,Uw,w',Uw' with coefficients in [-2,2], and then on the
larger integral span with x,Ux,r,Ur as well with coefficients in
[-1,1]. These are only necessary conditions. The 14 survivors are
all strictly negative for the product Gegenbauer kernel in
nested.json, which has constant term 1 and other coefficients
positive. Its diagonal is 74172033582491/2500000000 <29,700.

Run python -B nested.py (stdlib only). The second local-degree
bound is potentially useful for refining the shell LP; no shell
exclusion is claimed.
