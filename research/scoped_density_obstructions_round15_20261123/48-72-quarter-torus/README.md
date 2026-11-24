# Sharp bounds on fourth-torsion depths of extremal 48D and 72D lattices

Let \(L\) be an extremal even unimodular lattice of rank \(n=48\) or \(72\), with squared minimum \(m=6\) or \(8\), respectively. For every \(v\in L\),
\[
 \operatorname{dist}(v/4,L)^2\le n/12.
\]
In particular the squared bound is **4 in 48D** and **6 in 72D**. These bounds are sharp across the stated families: the P48 and Gamma72 half-period examples from the corresponding stacking constructions lie in the fourth-torsion subgroup. This does not assert an attainment classification or a bound on general real phases. No new global density lower or upper bound is implied.

Put \(r=v^2/2\bmod4\) and \(H=\operatorname{dist}(v/2,L)^2\). The complete scalar reduction, including the half-period bound, is reproduced in [PRIOR_GATES.md](PRIOR_GATES.md) and [HALF_BOUND.md](HALF_BOUND.md); it confines a hypothetical counterexample to seven rows. Here we close all seven by using both the order and coefficient of a secondary cusp, compared to the actual **unsigned** nearby half-coset shell counts.

## Quarter and half coefficients

At the quarter-level cusp \(M=ST^{-2}S^{-1}\), the coefficient of \(t^{s/2}\) is, up to a fixed complex phase, a signed shell count on \(v/2+L\) at squared norm \(s\); in any case its modulus cannot exceed the count. In the notation of the scalar reduction, \(X,Y,W\) transform as \((W,-Y,X)\), and \(P\) transforms to \(\theta_4(2\tau)^r\theta_2(2\tau)^{4-r}\) up to a fixed unit. Expand exactly before taking any absolute values. In each row below \(A\) denotes the signed coefficient with the fixed common phase removed (that normalization need not be tracked).

| \(n\) | \(r\) | \(H\) | Required modulus, with squared norm \(s\) | Unsigned shell bound |
|---:|---:|---:|:---|---:|
| 48 | 1 | \(5/2\) | \(s=7/2:\ 516\) | 100 |
| 72 | 1 | \(5/2\) | \(s=11/2:\ 29648\) | 5680 |
| 72 | 1 | \(7/2\) | \(s=11/2:\ \lvert29648-66A\rvert\), \(A\) from norm \(7/2\), \(\lvert A\rvert\le2\) | 6272 |
| 72 | 1 | \(9/2\) | \(s=11/2:\ 29648\) | 6272 |
| 72 | 2 | \(3\) | \(s=5:\ \lvert17376-68A\rvert\), \(A\) from norm \(3\), \(\lvert A\rvert\le2\) | 5544 |
| 72 | 2 | \(4\) | \(s=5:\ 17376\) | 5696 |
| 72 | 2 | \(5\) | \(s=5:\ 17376\) | 5696 |

For the 72D \(r=1,H=5/2\) row, the same formula \(29648-66A\) applies, but the norm \(7/2\) shell is **empty**: if \(a\) realizes \(H\), a norm \(7/2\) vector \(x\) would give \(x\pm a\in L\setminus\{0\}\) with their mean squared norm \(6<8\). Likewise, a shortest half-shell of squared norm \(H<m/2\) is exactly one antipodal pair, which gives \(\lvert A\rvert\le2\) in the indicated rows.

For the four half-series gates, write \(N_s=\lvert\{x\in v/2+L:x^2=s\}\rvert\). In the half-torsion scalar series at \(z=\tau\), put \(X=\theta_2(\tau)^4\), \(Y=\theta_4(\tau)^4\), \(W=X+Y\), and \(Z=X+2Y\). Under \(T\), \(X\mapsto-X\) and \(Y\mapsto W\), so \(Z\) is invariant. Factoring \(Y^eW^e\) for the two equivalent shifted cusps leaves a polynomial in \(Z,X^2\), after dividing by the zero-free multiplier form \(\theta_3^2\theta_4^2\) in the half-integral case. Set \(e=\lfloor H\rfloor\) and normalize its first four integer \(q=\exp(2\pi i\tau)\) coefficients to \(1,0,0,0\). At the cusp \(S\), multiply back \(\theta_3^2\theta_2^2\) if applicable and exchange \(X,Y\); the common signs cancel in rank 72. Then the untwisted shifted theta identities are
\[
\begin{aligned}
 H=7/2:&\quad N_{11/2}+4N_{9/2}+78N_{7/2}=6272,\\
 H=9/2:&\quad N_{11/2}+4N_{9/2}=6272,\\
 H=4:&\quad N_6+24N_5+96N_4=136704,\\
 H=5:&\quad N_6+24N_5=136704.
\end{aligned}
\]
The displayed unsigned bounds follow by nonnegativity. The verifier uses \(t=\exp(\pi i\tau)\) in this half-level expansion.

## Unsigned counts from the parent shortest-shell design

For the other three unsigned bounds, take \(a\in v/2+L\) with \(a^2=H\), and set \(u=2a\in L\), \(s_0=u^2=4H\). For a norm-\(s\) point \(x\in v/2+L\) with \(s>H\), the nonzero \(x\pm a\) belong to \(L\), so \(H+s-2\lvert x\cdot a\rvert\ge m\). In these three cases \(H+s=m\); hence \(x\cdot a=0\) and \(x-a\) is a shortest vector of \(L\) with dot product \(-s_0/2\) against \(u\). This applies to the 48D row and the 72D rows \(r=1,H=5/2\) and \(r=2,H=3\).

For any shortest \(p\in L\), \(\lvert p\cdot u\rvert\le s_0/2\) follows from the nonzero \(u\pm p\) and their lower norm bound \(m\); the dots are integers. Let \(b_j\) be the count at absolute dot \(j\). The parent shortest shell is a spherical 11-design, hence for \(k=0,\ldots,5\),
\[
 \sum_{j=0}^{s_0/2} j^{2k}b_j
 = \lvert L_m\rvert(ms_0)^k\frac{(1/2)_k}{(n/2)_k},
\]
where the power for \(k=0\) is one. The design follows directly by the vanishing of positive-degree elliptic harmonic theta series for degrees \(2,4,6,8,10\): they vanish to order \(m/2\) and division by \(\Delta^{m/2}\) gives negative weight. Odd moments vanish by antipodality. The ordinary theta series similarly gives \(\lvert L_6\rvert=52416000\) in 48D and \(\lvert L_8\rvert=6218175600\) in 72D.

In these three cases \(N_s\le b_{s_0/2}/2\) by antipodality of the shortest shell. The bounds are recomputed as exact nonnegative moment polytope maxima by enumerating all six-site bases. The bounded polytope has six independent equality constraints, and every vertex has support of size at most six; zero coefficients are allowed in enumerated bases.

## Density consequence for a restricted family

Allow horizontal full layers of \(L\) stacked in \(\mathbb R^k\), with every layer shift in the fourth-torsion subgroup. For squared packing minimum \(\mu>D=n/12\), each vertical site supports exactly one phase, and vertical squared distances are at least \(\mu-D\). With \(c_k\) the optimal center density of the class of vertical packings, the center density is at most
\[
 \frac{c_k}{2^n}\mu^{(n+k)/2}(\mu-D)^{-k/2}.
\]
Indeed the vertical point density is at most \(c_k(2/\sqrt{\mu-D})^k\), and each horizontal layer has covolume one. For periodic stacks the vertical point set is periodic; for a full lattice it is a lattice (project modulo the fixed horizontal period, or extend an integer basis of that primitive period). Since \(L\) occurs in each layer, \(\mu\le D+2\). The bound is nondecreasing for \(\mu\ge D+k/12\), by logarithmic differentiation. It is therefore no larger than its value at \(D+2\); where the known half-torsion stacks attain that value, they are also optimal among fourth-torsion stacks in the same numerical interval. In particular this extends the previously matched P48 scope for periodic 49–51 and lattice 49–53, and the Gamma72 scope for periodic 73–75 and lattice 73–76. This is a partial packing-family result, not a new catalog entry.

## Replay and references

Run python half_bound.py and python verify.py (Python 3 + SymPy). Only exact integer and rational arithmetic is used. See [Nebe, harmonic-theta/design theory](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/venklat.pdf), [Nebe, rank-48 constructions](https://arxiv.org/abs/1312.3780), and [Nebe, extremal rank-72 lattice](https://arxiv.org/abs/1008.2862).
