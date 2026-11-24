# Exact local rigidity of the BN80 + D4 84D incumbent

This is a **local optimality result, not a new density bound** or a
claim of global optimality. It applies to the [previously verified
84D lattice](../packing83_85_deep_verified_20261117/README_DEEP.md).
[Proof and scope](PROOF.md): perfection and positive eutaxy rule
out all small lattice strains; fixed BN80 full-fibre periodic motifs
are separately locally rigid for every fixed finite-index supercell.
Perturbations breaking the full fibres or escaping the local
neighbourhood are not excluded.

The replay archive pins the Gram data and literal contact vectors.
Extract into a scratch directory and execute
`OPENBLAS_NUM_THREADS=1 python -B verify.py`, then
`OPENBLAS_NUM_THREADS=1 python -B eutaxy.py`.
The independent replay checks exact membership, local matrix identities,
and positive weights. Archive SHA256:
`702ede0b822fae752d30b6af289bc4d158ee3117fefb79b4c495b7a75bf2b2d7`.
