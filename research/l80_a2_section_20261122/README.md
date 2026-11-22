# An A2 restriction of the 83D Bachoc–Nebe extension

**Occupied sphere-packing density at least 9.17033181889404e−18 in 82D.**
This explicit lattice has minimum squared norm 8 and Gram determinant 3.
Its squared center density is (2^{82}/3), improving the preceding
(2^{86}/55) by a density factor (sqrt{55/48}), or more than
**7.0436048222094%**.

Take the first 82 ordered basis vectors of the [certified 83D L80 + A3
extension](../packing83_85_deep_verified_20261117/README_DEEP.md).
This is an actual sublattice, not an orthogonal projection. The parent
minimum is inherited, and the fourth retained basis vector has squared
norm 8. In the parent block coordinates the three vertical vectors have
Gram (I_3+J_3); the first two have Gram (egin{bmatrix}2&1\\1&2end{bmatrix}),
an A2 basis. The checker recomputes determinants from the actual pinned
full Gram, with exact rational arithmetic and no floating norm search.

The analogous leading 81D prefix has determinant 2 and attained minimum
8, reproducing rather than improving the existing 81D density. This
certificate proves a lower bound and makes no optimality assertion.

Run `python verify.py` with Python 3.9 or later. The only inherited
geometric premise is the previously certified minimum of the 83D
parent; its full replay is documented in the linked package.
