# Explicit lattices from Γ72 in dimensions 73–78

| Dimension | Occupied density at least | Gain over the preceding repository bound | Gram determinant |
|---:|---:|---:|---:|
| 73 | **3.81440098634152e-14** | **9.5445%** | 5/2 |
| 74 | **1.14405687613928e-14** | **20%** | 75/16 |
| 75 | **3.61534657312238e-15** | **31.4534%** | 125/16 |
| 76 | **1.31058292710898e-15** | **44%** | 625/64 |
| 77 | **3.98928481450320e-16** | **33.3180%** | 4375/256 |

The **78D lattice** gives occupied density at least **1.39319236921112e-16**, with Gram determinant `91875/4096`. It improves the named lattice comparison by **46.2857%**. That comparison is the perpendicular section of Bachoc–Nebe’s L80 by a primitive `2A2`, with determinant 48 and squared minimum 8. The [previously certified three-translate packing](../packing78_verified_20261117/) gives the stronger 78D entry in the unrestricted table. [Exact derivation of the L80 comparator](../packing78_verified_20261117/packing78_independent_20261117/external_result.json).

The 73–77D gains use the exact lower bounds on [commit `0974b5a`](https://github.com/alexlegeartis/SpherePackings/commit/0974b5afeac553a218a66a6f5340d926ecc8b01f). [Results with exact density squares and comparisons](results.json).

## Lattice construction

Let $G$ be the catalogue Gram of Nebe’s Γ72. The construction uses the published squared minimum 8. The verifier checks positive definiteness, evenness and determinant 1 against the archived catalogue. [Minimum theorem](https://arxiv.org/abs/1008.2862) · [Catalogue](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/neb72.html) · [Archived catalogue](../novel_routes_20261113/higher_parent/neb72.html).

For dimension $72+r$, use the first $r$ rows of $B$ and the leading rank-$r$ block of $V$. The input gives six integer numerator rows $B$ and a rational vertical Gram $V$. [Certificate](certificate.json).

~~~
[ G                G B^T/2          ]
[ B G/2            B G B^T/4 + V    ]
~~~

In an integer vertical layer $z$, the horizontal component lies in the half-coset $zB/2+\Gamma_{72}$. Its squared depth is one quarter of the least squared norm in $zB+2\Gamma_{72}$. The verifier checks all nonzero vertical vectors with squared norm below 8 through the coordinate bounds $z_i^2\leq 8(V^{-1})_{ii}$. Every nonzero vector of zero parity already has vertical norm at least 8. For every remaining short parity, the verifier calculates the largest Γ72 numerator norm that could violate the minimum and excludes it and all smaller norms. Evenness fixes the numerator norms modulo four. The embedded Γ72 attains squared norm 8. The determinant is $\det V$ by the Schur complement, and the squared center density is $2^{72+r}/\det V$.

The first four vertical ranks use the leading root forms of $D_4$, scaled by $5/4$. An order-three isometry of Γ72 cycles the nonzero colors in the first two layers. Adding two colors fixed by this isometry produces four orbits among the rank-four root colors. The decoder excludes numerator norms through 18 in one class per orbit. The corresponding depth is at least $11/2$. In rank six the vertical Gram uses an $E_6$ root form with an additional positive form on the last two coordinates. The same separation check handles every rank.

## Exact three-block decoder

From the catalogue Gram the verifier reconstructs an integral matrix $T$ and a reduced Leech Gram $F$ with the following identity. [Nebe, §§2.1–2.2](https://arxiv.org/pdf/1008.2862).

~~~
T diag(F,F,F) T^T = 2G.
~~~

The image of $T$ contains $2\mathbb Z^{72}$ and has index $2^{36}$ in the triple of Leech blocks. Two complementary rank-12 binary projections, $A$ and $I-A$, describe the image. The $A$ components of the three blocks agree, and the $I-A$ components sum to zero. These congruences have the same index as $T$. The verifier checks them on every basis row of $T$.

For a numerator class $w+2\Gamma_{72}$, write $wT=(r_1,r_2,r_3)$. A block candidate $y_i$ belongs to $r_i+2\mathbb Z^{24}$. Its glue labels are the projections of $(y_i-r_i)/2$ modulo two. The integer enumeration computes the least block norm for each pair of glue labels.

1. Let $J$ be the supplied integral frame with $JFJ^T=8I$. In frame coordinates $a=y_iFJ^T$, the squared block norm is $\sum_j a_j^2/8$.
2. Enumerate the 8,192 residues $\beta$ in the image of $FJ^T$ modulo four. For each $\beta$, enumerate every bounded integer vector $a\equiv r_i FJ^T+2\beta\pmod 8$.
3. Retain precisely those vectors for which $(aJ-8r_i)/16$ is integral. That quotient equals $(y_i-r_i)/2$ and supplies the two binary labels.

For a cutoff $c$ on the Γ72 numerator norm, the combined Leech norm budget is $2c$. Lower bounds on the other two norms limit each block enumeration. The join pairs the first two blocks with a common first label and uses the XOR of their second labels to look up the third block. This covers every triple within the cutoff since the norm is a sum of independent block norms. Frame enumeration and all glue and comparison arithmetic use integers or rationals.

## Replay

Run `python research/gamma72_leech_layers_20261122/verify_gamma_layers.py` with python-flint and a C++17 compiler. It rebuilds the catalogue binding, frame residues, projections and exact vertical cutoff list, compiles the decoder and reruns all exclusions. The controls test the same known norm-8 class at cutoffs 7 and 8. Use `--write-results` to regenerate the committed report.
