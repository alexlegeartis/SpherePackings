# Nonradical escape: a literal odd 15-cycle at depth 7

Keep the identical \(H\) carrier and notation from [the main note](README.md). There is an exact phase family outside the radical two-coloring. This is a geometric target, not yet a denser packing.

**Four radial blocks suffice at their absolute minimum.** In the table, a mask with four radical blocks and one zero has maximal four-times depth 24 and lower bound \(4\cdot6+0=24\); hence depth is exactly 6 at every syndrome. The minimum in the zero block can be zero only at its zero label. Therefore the four radical sites alone realize *every* syndrome at total four-times cost 24. Put an arbitrary coefficient mask \(\beta\) into the fifth site and fix it at any of its locally minimizing labels; the other four repair every desired syndrome at cost 24. The matching lower bound is immediate from the five local absolute minima. Thus, all integral tags have squared depth
\[
6+\frac14\min_{z\equiv\beta\bmod2}\|z\|_E^2.
\]
In particular:
- \(\beta=r\) gives depth \(15/2\);
- any nonzero E7 parity class with \(Q(\beta):=\|z\|_E^2/2\bmod2=0\) has local minimum 4 and depth 7;
- a nonradical class with \(Q=1\) has local minimum 2 and depth \(13/2\).

The local minima follow from the complete mod-4 table after minimizing over integral labels; the verifier checks all 128 masks. The radical vector is in the kernel of the polar form of \(Q\) and \(Q(r)=1\), which the verifier confirms from the Cartan matrix.

**Graph of explicit guaranteed depth-7 moves.** Define a graph on the fractional coefficient masks \((E/2E)^5\). Allow as difference either (a) \(r\) in all five positions, or (b) \(r\) in four positions and a nonzero \(Q=0\) mask in one exceptional position. Any allowed edge has exact depth at least 7, regardless of the integral tags chosen. This graph has odd girth exactly 15.

To prove the lower bound, consider any odd closed walk of length \(m\), so its fractional differences XOR to zero. In a fixed position let \(t\) count exceptions. If \(t=0\) the sum is \(r\); if \(t=1\) it is one nonzero exceptional mask. If \(t=2\) vanishing would require two \(Q=0\) masks differing by \(r\), impossible as \(Q(x+r)=Q(x)+1\). Thus each of five positions has at least three exceptions. Each edge contributes at most one exception, so \(m\ge15\).

For equality, use the three masks \(\{6,18,20\}\) in each position, each once as the exception. They are nonzero with \(Q=0\) and XOR to zero. Across the resulting 15 moves, each position has three exceptions summing to zero and twelve copies of \(r\), so the total sum is zero. The resulting odd walk has no repeated vertex except its endpoints, since that would give a shorter odd closed walk. This is therefore a genuine 15-cycle, not a fractional relaxation of the palette. All its vertices may be realized as actual torus phases by taking their literal half-integer representatives with any integral tags.

This odd-girth claim is about the graph of the indicated moves only. It does **not** classify all nonradical mask differences or claim that the full depth graph has girth 15. It gives an explicit way around the radical bipartition, while ruling out a triangular assembly based solely on these strongest verified moves. A denser 37D/38D construction still needs an actual normal arrangement and verification of every smaller-distance pair (not just the 15 marked edges).

