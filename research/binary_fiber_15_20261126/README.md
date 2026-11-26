# Exact bounds for binary-fibre and two-weight constructions in 15D

Date: 26 November 2026. Comparator: both overall and lattice center-density incumbent 1/(16√2), occupied density 0.0168575706567626…, unchanged in alexlegeartis/SpherePackings at df55ff74ab95b8cea80baa359f8a1b3bc4c7713b. No improvement or new record is claimed.

## Exact optimum of the 14+1 dimensional construction

Write Dₙ={z∈ℤⁿ: Σzᵢ even}. Normalize the 15D minimum squared distance to 8, i.e. packing radius √2. In each horizontal plane at any real height, independently use any collection of integer translates of √x·2D₁₄, with x≥1. “Integer translate” means shifts √x·u, u∈ℤ¹⁴, in units of the unscaled core, not arbitrary real shifts. Fibers may be repeated, heights need not be periodic. Each individual core has center intensity 1/(2¹⁵x⁷). If the total number of fibers per unit axial length is λ, the center density is λ/(2^(15/2)x⁷): beating the incumbent needs λ>8x⁷.

Among these whole cores all heights, offsets and x are covered by a slab bound. Each core's integer offset has a 14-bit residue word, with two 2D₁₄ lifts. Different binary words at Hamming distance d have actual horizontal minimum² dx regardless of either parity lift: choose equal signs where the words agree and use a sign choice at one differing coordinate to meet the parity constraint. Opposite lifts of the same word have distance² 4x.

Set hₖ=√(8−kx) for positive radicand. Every same-width half-open height interval has at most Cₖ fibers. Its distinct words have Hamming distance at least k+1. For k≥4 even opposite lifts of the same word are incompatible in this interval, giving Cₖ=A₂(14,k+1); for k=2 or 3 one may allow both lifts but no repetition of an identical lift, giving Cₖ=2A₂(14,k+1). Finally C₀=32768. The cap implies λ≤Cₖ/hₖ, including upper asymptotic density. From the original authors’ current short binary-code table (https://aeb.win.tue.nl/codes/binary-1.html), including A(n,2e−1)=A(n+1,2e), use:

| k | Cₖ | interval for x |
|---|---:|---|
|7|8|[1,12/11]|
|6|16|[12/11,120/91]|
|5|64|[120/91,3/2]|
|4|128|[3/2,39/20]|
|3|1024|[39/20,5/2]|
|2|2048|[5/2,3]|
|0|32768|[3,∞)|

On each interval x⁷hₖ≥Cₖ/8 by log concavity and endpoint checks. For x=p/q reduced the required sign is that of (8q−kp)p¹⁴−(Cₖ/8)²q¹⁵. Endpoint signs in interval order: [0,958425688779805]; [3833702755119220,55102069536060921287360412596]; [881633112576974740597766601536,2685817]; [10743268,66940788776526970436164]; [272920103347664932188763,5566644713]; [34473610102,9500402]; last at 3 gives 21486536 and then increases. Thus no member of this family beats the comparator.

Literal equality: take both signs of all rows of any 16×16 real Hadamard matrix with normalized first row and column, expressed as binary words (32 words). Delete any two nonconstant columns giving four head groups of 8 (alternatively any two distinct columns, since the antipodal doubling balances them). Set x=1. Place the four groups at heights 0,1,2,3 modulo 4 in Gray order of their deleted two-bit heads; each tail uses its integer 2D₁₄ coset. Within a group distinct words differ in 8 horizontal entries. Adjacent heads have tail distance 7 and vertical difference at least 1. Opposite heads have tail distance 6 or 14 and height separation at least 2. Repeating a single word occurs at vertical distance 4 and the core itself has minimum² 8. Thus λ=32/4=8: this is an explicit equality packing, not a new best density.

## Allowing a second continuous coordinate weight for a fully specified binary motif class

Take the binary [15,4,8] simplex code S given by all 16 evaluations of linear functions on nonzero vectors of F₂⁴. Choose ANY union of full S-cosets inside F₂¹⁵; independently lift each binary word once to an integer 2D₁₅ coset. Apply any diagonal metric with two coordinate weights. Ignore global scale and write the squared weights as 1 and 1+z, with z≥0 on a chosen heavy-coordinate subset H of size r. When r=1 or 14, this is contained in the preceding 14+1 theorem, treating the distinct coordinate as the height. The all-equal case is covered there too. For 2≤r≤13 there are ≥2 light coordinates, so the horizontal base minimum² is 8. Nonidentical words whose XOR has total weight d and heavy weight j have separation² d+zj, regardless of parity-lift choices. Thus a union of K simplex cosets has density ratio K(m/8)^(15/2)/(1+z)^(r/2), where m is the minimum of 8 and all its differences' lines d+zj.

On a linear segment for m, this ratio’s derivative has the sign of 15j−rd+j(15−r)z. Its zero, if any, is a minimum. Suprema occur at z=0 or rational kinks z=p/q with reduced 1≤p,q≤7: below the ceiling 8 any contributing positive-slope line has 1≤d≤7, 0≤j≤d. For z≥7 only zero-slope differences or the ceiling remain. GL(4,2) acts on the 15 coordinate labels; its elementary transvections generate 46 orbits of subsets altogether. The quotient F₂¹⁵/S has 2048 cosets. Translate a candidate to contain zero. At each kink evaluate all attainable thresholds m on the quotient. Admissible differences are all nonzero cosets whose 16 vectors have metric weight ≥m. A union exists exactly when there is a clique of the required size containing 0 in this Cayley graph.

Two inexpensive upper bounds eliminate most cases. If there are N admissible neighbors of 0, K≤N+1. Fixing h of the r heavy bits gives an ordinary shortened code of length 15−h with distance at least the smallest d such that d+z min(d,r−h)≥m; hence 16K≤2ʰ A₂(15−h,d). The same short-code table supplies these bounds; cardinal target comparisons are integer-power comparisons, not floating thresholds. The 58 remaining cases have at most 41 neighbors of 0. An exact coloring-bounded clique DFS (at most 30 nodes for any instance) rules out the target in all of them. Thus this entire two-level simplex-coset class is incumbent-optimal. No claim for arbitrary binary code motifs or general three-level, sheared, or real-shifted families. Run the [standalone verifier](verify_two_weight.js) with Node.js; no dependencies.

## Arbitrary translated extra site in 15D equality Hadamard packings

Let C be the 16 rows of a column-shortened Hadamard simplex code, length 15. Rows need not form a group. The periodic packing C+2D₁₅ has min² 8, meets the same density, and admits no additional full translated 2D₁₅ core. In fact its squared covering radius is at most 15/2. Given t, reduce each coordinate to its nearest integer and write b for its parity word, εᵢ for the coordinate error magnitudes, in [0,1/2]. For row c, J={i:bᵢ≠cᵢ}. Distance² to c+2ℤ¹⁵ is Σ εᵢ²+Σ_{i∈J}(1−2εᵢ). Nonempty J lets us realize the required parity by paying at most 4 min_{i∈J} εᵢ. If b equals a row, that row alone gives distance² ≤15εmax²+4(1−εmax)≤23/4.

Otherwise every J is nonempty. On a uniform random row of C the nonconstant columns are pairwise independent fair bits. For any k≤14 coordinates, the probability of agreeing with b on all is ≤1/(k+1)≤(15−k)/15 (Cantelli on their mismatch count, or Bessel); for k=15 it is zero. Apply this bound to the coordinates with εᵢ≤s, integrate in s from 0 to 1/2. Expected parity penalty is ≤(4/15)Σεᵢ. The expected unrestricted squared distance is 15/2−Σεᵢ(1−εᵢ), and εᵢ≤1/2. The stated covering bound follows. This does not forbid multiple simultaneous deletions and insertions.

Further isolated boundary: after deleting one original translate, any new compatible full translate has to come from the deleted binary word, possibly with its opposite parity lift. For b outside C the sum over all 16 squared distances is ≤120, strictly less if any coordinate fractional; an improvement with 15 surviving distances ≥8 would force the deleted distance zero. If b∈C, its row lies within 23/4 and must be the deleted row. Every other row differs on 8 coordinates. In the sum over those 15 rows, Σ d² ≤120+15Σ εᵢ²−12Σ εᵢ, since each coordinate mismatches in eight rows and 4 min_{J} ε≤(1/2)Σ_{J} ε. Thus all errors vanish. The two parity lifts of the deleted word have mutual squared distance 4, so one cannot replace one core by two, even at fractional offsets.
