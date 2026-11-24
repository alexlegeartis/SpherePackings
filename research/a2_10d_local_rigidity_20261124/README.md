# Strict local maximality of the 10D 18-site packing

The 18-site ternary motif from the earlier A2 search is a strict local density maximum, even allowing arbitrary infinitesimal changes to the period lattice and to every site independently. The number of sites in the chosen period cell is fixed. This does not exclude different motifs or perturbations in a larger supercell. No record is claimed.

Write the packing with period Z^10 and metric G, the direct sum of five copies of [[2,-1],[-1,2]]. A label a in {0,1,2} is represented by a(1,2)/3 modulo Z^2. For distinct labels the minimum squared separation in one block is 2/3, with the three nearest displacement vectors equal, up to common sign, to (1,2)/3, (1,-1)/3, (-2,-1)/3. Their mean is zero; their differences span R^2. The attachment has minimum Hamming distance 3, so the intersite minimum squared separation is 2; the period lattice also has squared minimum 2. Its center density is 18 (1/2)^5 / 3^(5/2) = 1/(16 sqrt(3)).

Normalize the squared minimum to 2. Let S be a first-order metric displacement and u_p first-order site motions, after fixing one site to remove common translation. For a feasible direction of nonincreasing period volume, tr(G^(-1) S) <= 0. In each block the native period roots e1, e2, e1+e2 require r^T S_ii r >= 0. Their outer products sum to three times the inverse block Gram matrix. Thus every root derivative is zero, and consequently S_ii = 0.

At any pair p,q of Hamming distance 3, all 27 combinations of nearest displacement vectors in the differing blocks are active contacts. Their derivatives are
    v^T S v + 2 v^T G (u_p - u_q).
They are all nonnegative. Under uniform independent choices of the three lifts in each differing block their mean is zero (the diagonal blocks of S are zero, and each lift has mean zero). Thus all the derivatives vanish. Take the mixed difference over two lift coordinates at a time: this eliminates the linear term and gives 2(v_i-v_i')^T S_ij(v_j-v_j') = 0. Because the lift differences span, S_ij = 0 for each pair of blocks occurring together in some Hamming-3 edge. All ten occur. With S = 0, a single lift difference shows that the block-i parts of u_p and u_q agree along every Hamming-3 edge differing in block i. Each of the five corresponding graphs is connected, so all site motions agree and hence vanish.

For strict locality, any sequence of distinct nearby packings with squared minimum at least 2 and volume no larger would, after dividing the differences in metric and site data by their norm and passing to a subsequence, give a nonzero feasible direction of precisely the kind ruled out. The metric parametrization already quotients ambient rotations, and one site is held fixed for translations.

Run python3 check.py for all finite contact-graph assertions. The productive way to improve this motif is therefore not an arbitrarily small deformation at the same period, but a changed motif or a larger supercell.
