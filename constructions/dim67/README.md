# Dimension 67: 12 translates from orthogonal rank-three and rank-two normals

This explicit periodic nonlattice packing has exact squared minimum **229/30** and center density `(12/120)*(229/120)^(67/2)`. Its occupied fraction is approximately **2.257403109415952e-13**, **1.886142%** above the preceding verified 6-translate 67D packing. It does not change the lattice-only lower bound.

## Ambient input and exact replay

The ambient lattice is the even unimodular Gamma72 of squared minimum 8 in [Nebe (2012)](https://arxiv.org/abs/1008.2862). Use the full 72-by-72 integral Gram matrix G from [Nebe–Sloane's catalogue](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/neb72.html), [archived in this repository](../../research/novel_routes_20261113/higher_parent/neb72.html). Run from this folder:

    node verify.js ../../research/novel_routes_20261113/higher_parent/neb72.html

The checker validates the archived input SHA-256, determinant and positivity, the primitive normal witness, all 66 relative normal distances, the equality witness, the nonlattice test and the exact density gain. It uses the ambient minimum 8 from Nebe's theorem.

## Construction and minimum proof

Use ambient integer coefficient columns x with norm squared xᵀGx. Let A have the following five sparse rows, with zero-based indices and all other entries zero:

    (50:+1, 62:-1)
    ( 3:+1, 62:+1)
    ( 0:+1, 16:+1)
    ( 0:+1,  5:-1, 21:+1)
    ( 0:-1,  1:-1, 19:+1)

Direct multiplication gives M=AGAᵀ=K⊕J where

    K = [[8,-3,-3],[-3,8,-2],[-3,-2,8]]
    J = [[8,-4],[-4,8]]
    det(K)=300, det(J)=48, det(M)=14400.

The A minor on columns (50,62,0,16,5) is +1. Hence the row lattice is primitive and AG sends Z^72 onto Z^5 because G is unimodular. Let H=ker(AG)∩Z^72 and use the orthogonal projection P(x)=x−AᵀM⁻¹AGx. To lift a desired label y, place A_I⁻¹y on the five pivot coordinates I and zero elsewhere to form b_y, and take x_y=G⁻¹b_y. Thus AGx_y=y exactly. The projected normal image is the dual of the row lattice, of covolume 1/120; decomposing an ambient basis into H and normal-image lifts shows H has covolume 120. This is the same projection principle as in the [66D construction](../dim66/README.md).

Use the 12 labels in {0,e₁,e₂,e₃}×{0,e₁,e₂}. The first block has normal squared diameter 1/5, the second 1/6. Thus for distinct labels y,y′ the difference between any centers in H+P(x_y) and H+P(x_y′) comes from a nonzero ambient lattice vector and loses squared length at most 11/30 under projection. Its squared length is at least 8−11/30=229/30. The translates are therefore distinct. Inside one translate the squared minimum is at least 8. The ambient vector with coefficient column e₂₄+e₂₇+e₃₇ (zero-based) has squared norm 8, label (1,0,0,−1,1) in the difference set, and normal squared length 11/30, yielding equality.

Modulo H the coset assignment y↦P(x_y)+H has kernel MZ^5, by primitivity of the row lattice. The double of the first unit label is not congruent to any of the 12 labels modulo MZ^5 (checked exactly); because the packing contains zero, it is not a lattice.

## Comparison

The prior 6-translate packing `gamma72_simplex_rank5.npz` in the [verification ledger](../../research/novel_routes_20261113/higher_parent/gamma72_verified.json) has normal determinant 13122 and certified squared minimum at least 70/9. The squared ratio of the two certified densities is

    (4·13122·687^67)/(14400·700^67) > 1.

The comparison with the older, separately attributed [Gamma72/Leech-block section](../../research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) is +32.920521%. Neither comparison asserts an exhaustive worldwide priority or optimality.
