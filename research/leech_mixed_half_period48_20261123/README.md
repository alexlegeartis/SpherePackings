# A density gap for mixed 24D half-period factors in 48D

Let Λ be the Leech lattice of covolume 1 and squared minimum 4. Let U be **any** covolume-one 24-dimensional lattice of squared minimum at most 3; it need not be integral. Consider an arbitrary motif of half-periods of the orthogonal product H=√2Λ ⊕ √(2r)U, for arbitrary real r>0. This need not define a lattice. Its center density is **strictly below 18/25 of the P48 center density** (3/2)^24. In particular, pairing an even Leech factor with an odd Leech factor, or with any rooted even unimodular 24D lattice, cannot attain or improve P48 in this family. This is not a statement about non-half-period motifs, other dimensions, or arbitrary 48D packings.

Let μ be the squared minimum and N the motif size. Then 0<μ≤min(8,6r) by the two period factors. The density ratio to P48 is
```
D = N/(2^24 r^12) · (μ/6)^24.
```

1. **μ≤4.** Partition the motif by its first half-period label (at most 2^24 choices). Within a fixed fiber the second coordinates form a periodic 24D packing with squared minimum at least μ and period √(2r)U. By 24D sphere-packing optimality, its number of labels is at most (8r/μ)^12, so D≤(2μ/9)^12≤(8/9)^12.

2. **4<μ≤11/2.** Two motifs with the same second half-period label would be at squared distance ≤4 (the squared covering radius of √2Λ) after a horizontal period translation; therefore the second labels are distinct and N≤2^24. Since r≥μ/6, D≤(μ/6)^12≤(11/12)^12.

3. **11/2<μ≤371/60.** Choose a shortest vector v of U. It is not divisible by 2 in U, since U has positive minimum, and its half-period in the second factor has squared depth ≤(3/2)r. If μ−4>(3/2)r, the occupied second labels cannot include both elements of any pair {h,h+v} modulo 2U: after appropriate period translations the total squared distance would be ≤4+(3/2)r. These disjoint pairs imply N≤2^23. Hence D≤(1/2)(μ/6)^12≤(1/2)(371/360)^12. Otherwise r≥(2/3)(μ−4); the trivial N≤2^24 gives D≤[μ²/(24(μ−4))]^12. This expression decreases for 4<μ≤8, so its value at 11/2 is an upper bound.

4. **371/60<μ≤8.** Project to the second 24 coordinates. The projected centers form a periodic packing with one coset per second label and squared minimum at least μ−4 by the horizontal covering radius. 24D optimality bounds their point density by [4/(μ−4)]^12. Dividing by the first-factor covolume and converting to center density gives D≤[μ²/(18(μ−4))]^12, again decreasing up to 8. Use its endpoint value at μ=371/60.

All five rational endpoint bounds above are strictly less than 18/25. Run `python verify.py` to check these comparisons. This argument needs no classification or computation of the odd Leech half-period alphabet.

Sources: [Cohn–Kumar–Miller–Radchenko–Viazovska, 24D sphere-packing optimum](https://arxiv.org/abs/1603.06518) and [Borcherds, Leech lattice covering radius](https://math.berkeley.edu/~reb/papers/leech/leech.pdf).
