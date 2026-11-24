# Refined quarter-torus depth bounds in 48 and 72 dimensions

Let L be an extremal even unimodular lattice of dimension n=48 or n=72, so its squared minimum is 6 or 8 respectively. Put D=n/12 (4 or 6). For v∈L, write r=v²/2 mod 4, H=dist(v/2,L)², and q=dist(v/4,L)². Then
```
n=48: q≤33/8.
n=72: q≤25/4.
```
Furthermore, the only 48D quarter phase not proved to have depth at most 4 here has r=1, H=5/2 and must have q exactly 33/8, with exactly 3728 contacts. In 72D the only surviving possibilities for q>6 are r=1, H∈{5/2,7/2,9/2}, with q=49/8, and r=2, H∈{3,4,5}, with q=25/4. **No assertion of attainment or new density lower bound is made.** This controls fourth-torsion phases, not the full torus covering radius.

## 1. Exact nearest-half geometry

The half-torus theta bounds give H≤D; its anisotropic case gives ≤D−1/2. A half phase of depth D has an antipodal 2-design contact shell. (Here we use the earlier exact half-torus theorem for these extremal lattices, and not a computation of any particular P48 or Γ72 model.) If H=0 then v/4 is already a half phase; hence assume H>0.

Choose u=v−2w with w∈L and u²=4H (a minimal representative of v+2L). Then v/4=w/2+u/4. Let d0 be the depth of w/2 and z a relative contact there. Since v²/2≡u²/2+2w·u modulo 4, the coset for z·u is Z+c, where c=((r−2H) mod 4)/4 is 0 or 1/2. Antipodality gives
```
q ≤ d0 + H/4 − max_z |z·u|/2.
```
If d0<D, it is at most D−1/2 by half-period length quantization. In that case the largest dot is at least c. If d0=D, the 2-design gives squared dot average D u²/n=H/3. Thus the largest dot is at least the smallest nonnegative s∈Z+c with s²≥H/3. The displayed inequality gives an exact bound in both cases. Finally all squared norms in v/4+L belong to r/8+(1/2)Z, so round down to that progression. For possible H use H≥3/2 in 48D and H≥2 in 72D directly from the inherited minima. Also H∈r/2+Z.

## 2. Quarter modular forms and the signed half shell

Use q-series indeterminate Q=exp(2πiτ), and define the twist f=Σ exp(2πi x·v/4) Q^(x²/2) over x∈L. It is T-invariant and on F=ST⁴S⁻¹ has multiplier i^r. These and −I generate Γ₀(4) by Euclidean reduction of its first column. Under S, Poisson gives the untwisted theta series of v/4+L.

Set X=θ₂(2τ)^4, Y=θ₄(2τ)^4, W=X+Y=θ₃(2τ)^4. For r>0 divide f by P=θ₄(2τ)^r θ₃(2τ)^(4−r), of weight two and the same multiplier; if r=0 use P=1. At cusp zero P has Q-order r/16, and at cusp 1/2 it has Q-order α=(4−r)/4 for r>0 (α=0 when r=0). The quotient is holomorphic at the cusps: cusp zero follows from length quantization, and at 1/2 the lower H≥3/2 suffices. There are no upper-half-plane zeros of P. On passing to z=2τ, Γ₀(4) is conjugate to Γ(2), whose even-weight ring is the homogeneous polynomial ring in X,Y. This follows by distinct leading t-exponents of the monomials and the valence dimension bound.

Suppose q>D. This forces a cusp-zero factor Y^i with i=2D+1 if r=0, and i=2D otherwise. The shift at cusp 1/2 has norm at least H, forcing a distinct factor W^j with j=ceil(H/2−α). Thus
```
f=P Y^i W^j Σ_{u=0}^d c_u X^u Y^(d−u),
d=n/4−(r>0)−i−j.
```
Its Q-coefficients through degree 2 in 48D or 3 in 72D are (1,0,...,0), from the parent minimum. The verifier solves this small linear system exactly. Incompatible rows are excluded. If the system is uniquely solved with c_d>0, the leading coefficient at cusp zero equals 2^(r+4i−n/2)c_d, and its corresponding squared depth is i/2+r/8. Indeed at S, X,Y,W transform respectively to −θ₄(τ/2)^4/4, −θ₂(τ/2)^4/4, and −θ₃(τ/2)^4/4; the common sign in f cancels for n=48,72.

The extra ingredient is the **actual coefficient**, not only the order, at cusp 1/2. Take M=ST^(−2)S⁻¹. First Poisson changes the twist to the coset v/4+L; T^(−2) introduces twist −v/2, up to a constant of modulus one. The next Poisson interchanges the shift and twist. Thus the coefficient at Q^(s/2) of f|M is, up to a common unit complex factor, a sum of unit phases over vectors of norm s in v/2+L. Its absolute value cannot exceed their number.

On z=2τ the same matrix is z↦z/(z+1), which maps X,Y,W to W,−Y,X, respectively. This follows directly from the standard S and T theta transforms. The P factor of weight two is transformed up to a unit to θ₄(2τ)^r θ₂(2τ)^(4−r), giving leading absolute coefficient 2^(4−r). Consequently the absolute coefficient of f|M at Q^(α+j) is
```
2^((4−r if r>0 else 0)+4j) · |Σ (−1)^(d−u)c_u|.
```

## 3. Why six more rows are impossible

In a half coset, for two distinct non-antipodal vectors x,y, both x±y belong to L and have squared norm at least its minimum m. Thus at shell norm m/2 their directions are mutually orthogonal. This bounds their number by 2n: the signed coefficients 304 in 48D at norm 3, and 5840 in 72D at norm 4, exceed those bounds 96 and 144.

In 72D, a half-coset shell of norm 9/2 must have pairwise dot ±1/2 except for equal and antipodal vectors: the minimum-eight bound gives absolute dot ≤1/2 and integrality of the difference gives half-integral dots. Its lines have common angle 1/9. The standard trace-rank bound on their Gram matrix gives at most 72(1−1/81)/(1−72/81)=640 lines, hence at most 1280 vectors. The modular coefficient 10,108 cannot fit either corresponding residue row.

For an anisotropic half phase of maximal depth D−1/2, its untwisted shifted theta series has 108 nearest vectors in 48D or 6272 in 72D. For completeness, divide the half-character theta by θ₃(τ)^2θ₄(τ)^2, then use the Γ(2) homogeneous ring at argument τ. With X₀=θ₂(τ)^4, Y₀=θ₄(τ)^4, W₀=X₀+Y₀ and Z₀=X₀+2Y₀, T-invariance gives the invariant ring in Z₀ and X₀². The S-cusp lower order gives factors Y₀^(D−1) W₀^(D−1), leaving degree d=n/4−1−2(D−1). The coefficients of Z₀^(d−2j)X₀^(2j) are fixed by the same first Fourier coefficients of the parent twist. At S the leading coefficient is 4·16^(D−1) times their sum. This is calculated in the verifier, and is smaller than the respective required signed quarter-cusp coefficients 516 and 29,648.

These exclusions, the short-half geometry and the scalar series comparisons yield the tables and bounds above. The earlier half-torus input uses only level-one harmonic theta theory, Γ(2) and the inherited minima. The present result is independent of a particular choice of extremal lattice.

Replay of the prior elimination stage: `python prior_scalar.py` (SymPy; exact rational and integer arithmetic only).
