# Exact quarter-torsion gate after the continuous half/third-period obstructions

This companion check is an exact finite result, not a record. The nearby constructive escape was to go from a pure three-period normal architecture to a pure four-period one, with fourth-order Coxeter–Todd shifts. At that higher denominator the raw \(D_4\) geometry alone appears capable of a large improvement, so it is important not to inherit the ternary phase assumption.

**Exact phase radius:** with \(K=K_{12}\) of squared minimum 4, the largest squared torus distance among \((K/4)/K\) is **\(17/8\)**. A literal needle attaining it in the Eisenstein model below has numerator
\[
(0,0,1,1,5+\omega,6+3\omega).
\]
The distance is \(34/16\), substantially shallower than the true covering radius squared \(8/3\).

Use \(K=\{x\in R^6:x\bmod2\in H_6\}\), \(R=\mathbb Z[\omega]\). The hexacode formula \((a,b,c,a+b+c,a+\omega b+\bar\omega c,a+\bar\omega b+\omega c)\) gives 64 words over \(R/2R\). To compute a phase \(x/4\), use only \(x\bmod8R^6\) because \(8R^6\subset4K\). First divide the possible parity words \(x\bmod2\) into three weights \(0,4,6\). There is a monomial hexacode automorphism sending a chosen representative of each weight to any word of that weight: the checker exhibits one for every word by checking the images of just the three code generators. It lifts to an isometry of \(K\) by unit multiplication and permutation of the Eisenstein coordinates.

For each of these three parity words use a \(4K\) translation to make both coefficients of each of the first three coordinates lie in \(0,\ldots,3\bmod8\): the high bits can be set independently precisely because the first three hexacode coordinates are systematic. Independent sign changes identify opposite residues in the last three coordinates while preserving the first three. This gives 129,536 representatives in total. For each, the exact squared distance times 16 is
\[
 \min_{h\in H_6}\sum_{i=1}^6
 \min_{(a,b)\equiv x_i-4h_i\pmod 8}(a^2-ab+b^2).
\]
The inner minimum uses only the displayed residue within \([-8,8]^2\); a centered representative has norm at most 48, and the norm is at least \(3a^2/4\) and \(3b^2/4\). The three maxima are \(32,34,34\). Thus this is full quotient coverage via exact symmetry/translation, not sampling.

**16D consequence:** with one whole quarter-torsion Coxeter–Todd translate per occupied site of \(\sqrt{s}D_4\), pure normal period \(4\sqrt{s}D_4\), and squared minimum \(\mu>17/8\), **no packing reaches Barnes–Wall's \(1/16\)**. As in the [main proof](README.md), at most one translate occupies a site. If \(2rs\) is the least occupied normal squared separation, \(r\ge1\) and
\[
rs\ge\max(\mu/32,(\mu-17/8)/2)
\]
by same-site period repetition and the exact phase radius. The elementary 4D normal (not lattice-only) density bound gives center density at most \((4/3)F(\mu,rs)\), with \(F(\mu,u)=(\mu/4)^8/(54u^2)\) as above. The endpoint maxima with this lower bound on \(rs\) occur at \(34/15\) and \(4\); both are strictly below \(1/16\). This does not treat arbitrary phases in a four-period stack or multiple transforms at one site below the stated radius.

**22D consequence:** for quarter-torsion phases and squared minimum \(\mu\ge3\), **all** normal site geometries are excluded, even without a normal period and even with metrics outside the main stretching theorem. Distinct normal sites are separated by at least \(\sqrt{\mu-17/8}\). The same Blichfeldt gauge in 10D gives unnormalized normal center density at radius 1 at most \(45/(2\pi^5)<5/54\). Thus the 22D center density is at most
\[
\frac{5}{54\cdot27}\left(\frac{\mu}{4}\right)^{11}
                         \left(\frac{4}{\mu-17/8}\right)^{5}.
\]
It suffices to check endpoints 3 and 4, since the only possible interior critical point is a minimum. At both endpoints the bound is below the attributed \(R_{22}\). This claim is for complete \(K_{12}\) fibers, not arbitrary subsets.

Run `node verify_quarter.js`. All replay calculations are integers, with the squared comparator checked without decimal rounding.
