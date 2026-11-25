# Exact limit for translated whole \(P_{10c}\) layers in 11D

This result concerns stacks using identical, equally oriented full horizontal \(P_{10c}\) layers at minimum squared distance 4. It does not constrain different orientations, changed motifs or periods, other scales, or arbitrary 11-dimensional packings.

Use the Construction A representation in [Conway–Sloane, Quaternary constructions](https://neilsloane.com/doc/Me185.pdf): the binary Gray image \(C\) of their 40-word pentacode is a motif in \(2\mathbb Z^{10}\). Between identical layers all horizontal differences, modulo the period, are the binary set \(D=C\mathbin{\mathrm{xor}}C\) with 261 words. The maximum over all real shifts of the minimum squared distance to \(D+2\mathbb Z^{10}\) is \(5/2\).

For a shift, let \(a\) be its nearest binary corner modulo two and \(w_i\in[0,1/2]\) its coordinate distances to that corner. Its squared distance to a binary difference word \(d\) modulo the period is
\[
\sum_i w_i^2+\sum_{i\in a\mathbin{\mathrm{xor}}d}(1-2w_i).
\]
If \(a\in D\), choose \(d=a\). Otherwise the [exact corner replay](verify.js) finds four difference words whose mismatch supports are pairwise disjoint, for each of the other 763 corners. Averaging their costs bounds the minimum by \(\sum_i(w_i^2+(1-2w_i)/4)\le10/4\), since \(w_i^2-w_i/2\le0\). The all-half shift attains equality.

Thus neighboring layer heights differ by at least \(\sqrt{3/2}\). Alternating an unshifted layer and the all-half shift at this spacing attains it; nonadjacent layers already have vertical distance greater than 2. The maximum center density in this fixed-layer architecture is
\[
\frac{40}{2^{10}\sqrt{3/2}}=\frac{5}{128\sqrt{3/2}},
\]
whose squared ratio to the \(P_{11a}\) comparison \(9/256\) is \(200/243\). Improving that comparator at this scale with full 40-point layers would instead need squared horizontal cross separation above \(224/81\). Hence seek different layer motifs or a genuinely different architecture rather than tuning the same shift.

Run `node verify.js` for all 1,024 corners. No density table was changed.
