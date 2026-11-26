# A 22-dimensional carrier inside Bacher's 26-dimensional lattice and its exact phase costs

**Construction component, not a new packing record.** The 27-dimensional lattice target has center density greater than Bacher's (1/√3); the unrestricted comparison is Vardy's (1/√2). The 22+5 gluing below would tie the latter and improve the former if all 64 torus phases worked over the densest 5D lattice. They do not. The 64 exact distance requirements apply to **other** 5D periods and are the reusable output.

## Coordinates and horizontal carrier

Take the top-left 26 by 26 block of [Martinet’s published `bac27a_orig` Gram matrix](https://jamartin.perso.math.cnrs.fr/Lattices/bacher.gp), divided by two, as an integral Gram matrix (G) for (L). It is positive definite, even, has determinant 3 and minimum 4. The supplied script embeds the original matrix and checks the determinant and entire norm-4 shell (117,936 vectors). The parent and comparison come from [Bacher](https://doi.org/10.1007/s002220050180); the unrestricted comparison is [Vardy](https://doi.org/10.1090/S0002-9947-99-02169-8). Center density means (1/\mathrm{covol}) at minimum squared 4.

The four rows `V` in the verifier are twice vectors (z_1,…,z_4) of identical residue modulo the lattice; their Gram matrix is (B=2I_4+J_4). Thus the even-coefficient-sum combinations (K=\{\sum n_i z_i:\sum n_i=0\pmod2\}) form a lattice in (L), of determinant (4\det B=192). This is primitive: an index at least 2 enlargement of this rank-4 lattice would have determinant at most 48 while still of minimum at least 4, below the sharp four-dimensional Hermite threshold 64 (attained by the scaled (D_4)). Put (E=\operatorname{span}K) and (S=L∩E^\perp).

Projection to (E) in the (z_i) coordinates is (B^{-1}(\langle x,z_i\rangle)_i). Modulo (K) the image (H) has 64 elements; the script generates and checks the entire finite group. As (L/(K⊕S)≅H), we obtain (\det S=3·64^2/192=64). Orthogonal projection of (L) onto (E^\perp) lands in (S^*); mod (S) it is injective on these 64 labels by primitivity of (K), so it is the **whole discriminant group** (S^*/S). A coordinate representation of the phase label stores the four (z_i) coefficients times 24 modulo 24 times the even-sum integer lattice. The four representative numerator vectors ((30,18,18,18),(18,6,18,6),(36,0,0,12),(36,12,0,0)) have respective orders (4,4,2,2) and generate the direct product (H).

## Exact required distances

Let (q_c) be the least squared length of the normal-coordinate coset (c) of (K). For nonzero (c), a shortest vector of its (S^*/S) class can be glued to the shortest (K) representative of its matching class. The resulting lattice vector in (L) has even norm at least 4. If the class has a norm-4 lift, its horizontal minimum is (4-q_c). A norm-4 lift with nonzero horizontal component necessarily has shortest normal component: subtracting an appropriate (K) vector otherwise violates minimum 4. The remaining eight nonzero classes have precisely six values (q_c=5/2) and two (q_c=3); none has a norm-4 lift, and each has an explicit norm-6 lift with shortest normal representative. Their horizontal minimum is exactly (6-q_c). The eight witnesses are checked directly against the input Gram matrix.

Therefore the vertical requirement (W_c=4-\min\|S+c\|^2) for any distinct classes is exact. In the following order, (i=a+4b+16c+32d) for the four generators above, the entries are **four times** (W_i); index zero is set to zero, since the vertical period itself must separately have minimum squared at least four.

```
0 3 4 3 2 3 2 7 8 7 4 7 2 7 2 3
4 3 8 3 2 7 2 3 4 7 8 7 2 3 2 7
4 3 8 3 2 7 2 3 4 7 8 7 2 3 2 7
4 7 8 7 6 7 6 3 4 3 8 3 6 3 6 7
```

## Positive target and exact failed specialization

For a 5D period lattice (P), glue all 64 horizontal classes to a group homomorphism (H→\mathbb R^5/P). This is a 27D lattice of covolume (8\operatorname{covol}(P)/64). Separations are sufficient and necessary: (P) has minimum at least 4 and each nonzero vertical class has squared length at least (W_c). A solution with (\det P<192) beats Bacher's lattice benchmark. For (P=√2D_5), the determinant is 128, so a solution would have center density (1/√2), equal to Vardy's unrestricted benchmark. Since (H) has exponent four, any group homomorphism lands on the finite 1/4-period grid, not merely a chosen discretization. Enumerate all four order-two images (two of which must double the order-four generators). Exactly 480 assignments pass their 15 constraints. No choices of order-four lifts pass the 63 constraints. This excludes the proposed full additive gluing on this particular period, not on other 5D lattices or nonlinear motifs. The verifier includes the independent exact integer shell enumeration, eight explicit shell-6 witnesses, and the full D5 homomorphism check.

## Geometric interpretation

The same (E) is a whole flat on which shifting the known Bacher phase cannot improve the 27D layer lattice. Choose (q_B=z_1\bmod L). For any displacement in (E), the lattice sites tested by the odd-coefficient-sum combinations of the four (z_i) form a translate of the checkerboard (D_4) with squared axis lengths (2,2,2,6): a normalized 4 by 4 Hadamard matrix takes the even-sum lattice to itself and diagonalizes (2I+J). Its squared covering radius is 3. Indeed, for nearest-integer coordinate rounding let (t=\min_i w_i(1-2|u_i|)) where (w=(2,2,2,6)) and (|u_i|≤1/2). If parity is wrong, flipping a cheapest coordinate has cost at most (3-t+5t^2/12≤3) for (0≤t≤2); if parity is right, the cost is at most 3 directly. Half corners attain 3 for the checkerboard. This is not a claim that the entire 26D lattice has covering radius (√3).

Run `python3 verify.py` (standard library, one core).
