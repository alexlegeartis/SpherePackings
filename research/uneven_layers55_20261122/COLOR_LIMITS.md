# Limits on using finitely many layer offsets

Suppose a periodic packing in dimension $m+r$ places a full translate of the
same $m$-dimensional packing at each layer position in $\mathbb R^r$. Let $s$
be the number of horizontal translates used, counted modulo the periods of
the horizontal packing. Write $\delta_m$ for the horizontal center density
at its largest packing radius and $\delta_r^{\max}$ for the optimal center
density in dimension $r$. Then
\[
\delta_{m+r}\le s\,\delta_m\,\delta_r^{\max}.
\]
Within one horizontal translate, identical horizontal centers occur at every
layer position of that translate. Their vertical distances therefore must
be at least the common sphere diameter. Each set of vertical positions with
that translate has center density at most $\delta_r^{\max}$. Summing over
the $s$ translates and multiplying by the horizontal density proves the
inequality. It allows irregular layer spacings. Labeling a vertical lattice
by its residue class modulo twice that lattice gives $s\le 2^r$.

## Consequences for the target dimensions

- Starting from Leterrier54, three translates can give center density at
  most $3\delta_{54}/2=90,459.8387\ldots$ in 55D. The printed figure would
  require at least four translates. Starting from the same parent in 56D,
  even matching MW56 would require at least **12** translates, since
  $2\sqrt3\,\delta_{56}/\delta_{54}=11.270789\ldots$.
  [Leterrier54 and the printed 55D comparison](../cross_coset55_20261116/PROOF.md#external-comparison) · [MW56](../cm58_20261116/PROOF.md#published-parent-and-exact-gram-binding).
- Starting from Γ72 and using a label for each residue class modulo twice
  an 8D vertical lattice, any such 80D packing has center density at most
  $2^8\cdot2^{36}/16=2^{40}$. This is **90.20%** of the MW80 benchmark,
  irrespective of the chosen coset depths. [Γ72](https://arxiv.org/abs/1008.2862) · [Optimal 8D center density](https://cohn.mit.edu/sphere-packing/#data) · [MW80](../creative_sections_20261116/PROOF.md#78d-and-79d-antipodes-from-the-reconstructed-mw80).
- Starting from MW80 and using a label for each residue class modulo twice
  a 12D vertical lattice, the 92D center density is at most $2^{12}\delta_{80}$.
  Here the packing fraction bound 1 and $\pi^6/720>1$ give
  $\delta_{12}^{\max}\le1$. The resulting ceiling is **51.09%** of MW92.
  [MW80](../creative_sections_20261116/PROOF.md#78d-and-79d-antipodes-from-the-reconstructed-mw80) · [MW92](../mw92_sections_20261115/PROOF.md#reconstructing-the-published-parent).
