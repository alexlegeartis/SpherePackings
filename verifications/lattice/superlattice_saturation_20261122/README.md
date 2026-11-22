# Lattice packing bounds with forced parity classes

The superlattice argument gives stricter lattice-only upper bounds in dimensions 17, 18, 20, 21, and 22. A separate Fourier refit gives the strongest bound in dimension 18. The [13–48D comparison](COMPARISON.md) compares these bounds and the earlier follow-ups against Zubrilina’s `sphere-packings` benchmark. [Native proof sources](https://github.com/AstroNinja389/sphere-packings/blob/d7ffd0df47dfeb11814631858205f05ed73ec8ff/proofs/COMPARISON.md)

| Dimension | New upper bound | Reduction from native certificate |
|---:|---:|---:|
| 17 | 0.013626238385270738359807 | 5.7491e−8% |
| 18 | 0.01020979021930188741773 | 0.00018783% |
| 20 | 0.005812874640290415375654 | 1.0512e−7% |
| 21 | 0.004409193412310544839783 | 5.7810e−8% |
| 22 | 0.003330067292822832547769 | 3.8670e−8% |

All upper bounds are rounded upward. The verifier checks the exact values in [certificates.json](certificates.json) and [refit-18.json](refit-18.json).

## Forced classes in a densest lattice

Normalize the squared minimum to one. For a $k$-dimensional subspace $W\subseteq L/2L$, adjoining the half-vectors from the classes in $W$ gives an index-$2^k$ superlattice. If $s$ is the smallest squared length of a representative in a nonzero class in $W$, the new squared minimum is $\min(1,s/4)$. Densest-lattice extremality gives $2^k\min(1,s/4)^{n/2}\le1$, so $s\le4\,2^{-2k/n}$.

The classes at or below that radius meet every $k$-subspace and hence number at least $2^{n-k+1}-1$. To see this, choose a subspace maximal among those missing these classes. It has dimension $d<k$, and each nonzero coset of it meets the class set. This gives at least $2^{n-d}-1$ distinct classes.

The replay uses the original native proof's per-bin Gaussian bounds and signed prices. All residual scores are nonpositive. Each forced class has a shortest representative inside a certified bin. Every cutoff used in the deduction is below squared length four, so the representative is primitive. Count the $n(n+1)/2$ exact minimum classes already used by the parent first. The additional classes up to each successive cutoff give the stated numerator deduction. A class already counted at a smaller cutoff is paid once.

## Fourier refit in dimension 18

The new degree-81 Gaussian polynomial keeps the native population inequalities and changes their weights. The spatial sum over the first three multiples of each primitive representative is enclosed on every price bin. The positive residuals have separate short and long allowances. A parity class has at most one representative line below squared length two, so the short allowance enters once. The degree-four spherical bound from the parent controls the total class size.

For a contender of density above the recorded threshold, the rank-17 section bound gives a lower bound on the shortest dual squared length in Fourier coordinates. The verifier checks the rational dual polynomial from that lower endpoint, and it checks the physical tail from the last price bin. It covers the threshold case by taking the larger of the threshold and the conditional density bound.

## Offline replay

Download the [pinned native archive](https://github.com/AstroNinja389/sphere-packings/raw/d7ffd0df47dfeb11814631858205f05ed73ec8ff/proofs/native-proofs.zip). The verifier checks the archive and every manifest entry, replays the final parent proofs with their recorded census receipts, and checks the new deductions. The upstream full census replay is documented in its proof guide. [Upstream replay guide](https://github.com/AstroNinja389/sphere-packings/blob/d7ffd0df47dfeb11814631858205f05ed73ec8ff/proofs/README.md)

```sh
python -m pip install python-flint==0.9.0 numpy scipy sympy mpmath
python -B verify.py /path/to/native-proofs.zip
python -B verify.py /path/to/native-proofs.zip 18
```
