# All diagonal half-period E8 metrics: the high-minimum range is excluded

Take seven real squared metric weights on the zero-sum E8 carrier from the two-exceptional proof, now allowing **all seven weights to differ**. Normalize the **second-largest** weight to 1, let $a\le b$ be the two smallest, and suppose the squared minimum of an arbitrary nonempty half-period motif is $g>7$. This is impossible if $a\le0$: the period contains opposite root vectors in the two least-weighted slots, giving $g\le4(a+b)\le4$. More precisely $t=a+b>7/4$, and since $b\le1$, all the weights are positive.

Fix a label of largest weight in each phase. The remaining six weights are at most 1. Every uniform six-slot witness of length at most 7 remains at most 7, so the uniform Fourier dual in the verifier bounds the number of phases across all such fibres by $B=38{,}246{,}144$.

Write $D=\sum_i\prod_{j\ne i}w_j$, the determinant of the scalar zero-sum metric after eliminating one slot. By monotonicity in positive weights, first lowering the largest weight to 1, retaining one slot of weight 1 and lowering every weight other than $a$ and this 1 to $b$ gives
$$
 D \ge b^4(b+5a+ab)=b^4(5t-4b+bt-b^2).
$$
At fixed $t>7/4$, for $t/2\le b\le1$ the derivative of the right hand side in $b$ is
$$
 b^3[20(t-b)+b(5t-6b)]>0.
$$
Thus $D\ge x^5(6+x)$ where $x=t/2\in(7/8,1]$. The center density is at most
$$
 B\frac{(g/8)^{24}}{D^4}
 \le B\frac{x^{24}}{x^{20}(6+x)^4}
 =B\left(\frac{x}{6+x}\right)^4
 \le \frac{B}{7^4} < (3/2)^{24}.
$$
Here the period covolume is $2^{24}D^4$ as before. No motif symmetry or additivity is used. The conclusion is **only** for $g>7$ with the second-largest weight normalized to 1, but it places every remaining diagonal metric search in the range $g\le7$.
