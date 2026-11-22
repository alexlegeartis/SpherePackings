# Poisson bounds for short vectors in selected intervals

These density bounds apply to every full-rank lattice in the listed dimensions. The proof counts short vectors in the squared-length intervals that contribute most to the earlier parity count. It then recomputes that count with the new constraint.

| Dimension | Density upper bound | Reduction from Zubrilina | Reduction from preceding project certificate |
|---:|---:|---:|---:|
| 25 | 0.001286546083743192189889370674 | 0.00092984% | 0.00082474% |
| 31 | 0.000135410497501983337616757649 | 0.01411027% | 0.00257238% |
| 32 | 0.000095866360403730692795236277 | 0.01099298% | 0.00274740% |

## Proof

Normalize the squared lattice minimum to one. Write $V_n$ for the unit ball volume, $D$ for the covolume and $U_n$ for the preceding certified density upper bound. Set $T=0.99U_n$. Lattices of density at most $T$ already satisfy the new bound. On the remaining branch, the two density bounds give

$$
\frac{V_n}{2^n U_n}\le D\le \frac{V_n}{2^n T}.
$$

For a shortest reciprocal vector $y$, the section perpendicular to $y$ has covolume $D\lVert y\rVert$ and minimum at least one in length. The certified bound $U_{n-1}$ therefore gives $\lVert y\rVert\ge 2V_{n-1}T/(V_n U_{n-1})$.

The certificate gives a Gaussian times a radial polynomial $f$ which is nonnegative at every nonzero lattice length, and at least $\mu w_i$ on each chosen squared-length interval $I_i$. Its Fourier transform is nonpositive at every permitted nonzero reciprocal length. Let $M_i$ count opposite pairs assigned to $I_i$, counting each pair once. Poisson summation gives

$$
2\mu\sum_i w_i M_i + f(0)
\le \sum_{x\in L} f(x)
\le \frac{\widehat f(0)}{D}.
$$

The checker uses the covolume endpoint giving the largest right side. It bounds the weighted count by $C$. It then subtracts a rational price $\lambda w_i$ from the score for each pair in $I_i$ and adds the allowance $\lambda C/(2^n-1)$ in the original units. Replaying every parity class keeps all admissible longer companions of each short vector. Directed rounding gives the density bound for this branch. The published bound is the larger of the resulting upper bound and $T$.

## Replay

Use Python 3.11 or later, NumPy and python-flint 0.9.0. The verifier runs the earlier angular proof, checks polynomial signs in exact arithmetic, and encloses the Gaussian extrema and density factors with Arb. From this directory, supply the pinned ZIP or its extracted directory. [Pinned ZIP](https://github.com/AstroNinja389/sphere-packings/raw/7b63ee2094b7489aea3da83349c547c8da4cb388/proofs/n25-48-proofs-r9.zip) · [Earlier angular verifier](../angular_kernel_followup_20261122/)

~~~sh
python -B verify.py /path/to/n25-48-proofs-r9.zip
~~~

If the earlier verifier lives outside the default sibling directory, supply its path with the --angular-proof option.
