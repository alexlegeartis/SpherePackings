# Contact shell of the 14D triangular packing

The [packing in this directory](README.md) gives a 1620-point kissing configuration centered at a suitable point in a half-integer Witt sheet. It does **not** improve the known unrestricted lower bound of [Ganzhinov, 1932](https://arxiv.org/abs/2207.08266). It has more contacts at that center than the 1422 [minimal vectors of Λ14](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/LAMBDA14.html); density and the largest contact count are different objectives.

To obtain the shell, subtract the center from all its distance-2 neighbors. These vectors have length two and mutually separate by at least two. The exact enumeration `python3 verify_contacts.py` counts every center in the 36-layer period (4800 centers). All 3072 centers in binary Hamming sheets have 1464 neighbors. The maximum over the 1728 centers in Witt sheets is 1620, attained at 16 of them. For example at vertical site (0,0) with binary word 3123 in the bit convention of verify.py:
- 840 in its own sheet: the 24 axial translates and \(16\cdot51\) weight-four differences;
- 768 in the six nearest integer sheets: the 128 words of each at one closest horizontal translate;
- 12 at squared vertical distance three: at each of the six neighboring sheets, one weight-one difference gives two nearest horizontal translates;
- zero at squared vertical distance four. Larger vertical distance cannot contribute.

This also identifies a limitation of keeping the same contact geometry. For a Witt center, all neighbors in any of the six nearest integer layers have twelve horizontal components ±1/2. At a fixed layer their sign vectors have Hamming distance at least four; between adjacent layers of the six-cycle they have distance at least three. Append a 0 to one code and a 1 to the other. Their combined size is at most \(A_2(13,4)=256\), as in [Brouwer's binary code table](https://aeb.win.tue.nl/codes/binary-1.html). Summing around the cycle gives 768, already attained. That argument is independent of linearity or use of the same Hamming cosets. To transfer the layered idea to a better local arrangement requires changing more than the outer codebooks on this ring.
