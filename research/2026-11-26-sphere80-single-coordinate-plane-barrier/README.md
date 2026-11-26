# A complete single-plane barrier in a fixed 11D binary-period construction

**Scope.** Let a scaffold of density `|C|/2^11` use a binary code C and the cubic period 2Z^11, with minimum distance 2. Fix r coordinate directions and one integral residue y in the remaining 11−r coordinates modulo 2. Permit a packing in the r-flat above y, repeated over those transverse periods, and require that the projected scaffold remain at distance at least 2 from that flat (so safety does not rely on tuning axial positions). For r=2,3,5,6 the normalized returns of the A2, A3, D5, E6 guests are respectively 2/sqrt(3), sqrt(2), 2sqrt(2), and 8/sqrt(3). No such single guest, or any no denser guest of this form, improves the known center density 72/2048=9/256. This is not an unrestricted optimality claim for 5D or 6D guests, real transverse shifts, simultaneous interacting fibers, or 11D sphere packing.

**Pullback to all optimal binary codes.** The optimum binary lengths in the r-plane with distance four are A(r,4)=1,1,2,4 respectively (small Plotkin bounds and explicit codes). Every stated guest return lies between A(r,4) and A(r,4)+1. Thus if |C|+return>72, then the integer |C|+A(r,4)≥72. Replace the entire r-dimensional guest by an optimal binary r-code in exactly its integral transverse fiber. All these words are compatible with C by the projection hypothesis, so this gives a binary (11,72,4) code by the known optimum 72. Only the A(r,4) newly implanted words in that code have projected Hamming distance at most 3 from y.

Östergård, Baicheva and Kolev's [complete 96 representative codes](https://users.aalto.fi/~pat/72.html) give the following actual minima over every representative, every choice of r coordinate axes and every integral projected residue y:

| r | A(r,4) | minimum words at projected Hamming distance ≤3 |
|---|---:|---:|
| 2 | 1 | 5 |
| 3 | 1 | 13 |
| 5 | 2 | 38 |
| 6 | 4 | 52 |

Coordinate permutations and binary bit translations preserve this exhaustive minimum. In each case it exceeds A(r,4), contradicting the required pullback. The checker uses exact Hamming distances; the floating-point shapes of the supplied guests are not needed. No new density record.

**Replay.** Download the author's [72_4.txt](https://users.aalto.fi/~pat/opt/72_4.txt) and run `g++ -std=c++17 -O3 verify.cpp -o verify && ./verify 72_4.txt`. The normalized source SHA-256 used in the search is `0f3fe74e631f11d3d3263c1c476d86ad142298ed6e2184f2b06370c771519d23`. The independent check in the round verified all 96 code sizes and minimum distances and reproduced the four minima with a different Hamming-neighborhood accumulation.
