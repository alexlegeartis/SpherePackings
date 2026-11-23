# Ternary phase depth bounds for extremal even unimodular 48D and 72D lattices

These are scoped torus bounds, not unrestricted covering-radius or sphere-packing upper bounds. No new packing density is asserted.

Let L be any even unimodular lattice of rank n=48 or72 with minimum squared norm n/12+2. For v∈L set a=v/3, r=v²/2 mod3 in {0,1,2}. Then the squared distance from a to L is at most:

| n | r=0 | r=1 | r=2 |
|---:|---:|---:|---:|
| 48 | 4 | 32/9 | 34/9 |
| 72 | 6 | 50/9 | 52/9 |

Whenever the bound is attained, the corresponding nearest shell has, respectively, 1836, 105, or 450 relative vectors in 48D, and 105084, 6060, or 25632 in 72D. Attainment is not asserted.

## Modular reduction

Write q=exp(2πiτ), ζ=exp(2πi/3), k=n/2. Use the ordinary weight-k slash action and f(τ)=Σ_{x∈L} ζ^(x·v) q^(x²/2). Then f|T=f, since L is even. Poisson summation gives f|S=Σ_{x∈L+a}q^(x²/2), up to replacing a by −a (with the same series). The scalar (−i)^k is 1 in the two ranks under consideration. Because (x+a)²/2 is congruent to v²/18 modulo 1/3, f|S T³=ζ^r f|S. Thus for F=S T³ S⁻¹ we have f|F=ζ^r f.

The matrices T, F and −I generate Γ₀(3). Indeed the first column of a member is (u,3w); left F-power reduces |w| to at most |u|/2, then left T-power reduces |u| to at most 3|w|/2. If w≠0 this decreases |u|, and u cannot become zero since 3 does not divide u. When w=0 the remaining matrix is a signed T-power.

Put U=η(τ)³/η(3τ). This is nowhere zero on the upper half plane, and U|T=U. Dedekind eta's S transformation gives
U|S=−i√3·η(τ)³/η(τ/3).
Under T³ the eta quotient on the right multiplies by exp[(9−1)πi/12]=ζ, hence U|F=ζ U. At the cusp zero the leading exponent of U|S is q^(1/9); at infinity U has constant term1. Therefore f/U^r is holomorphic at both cusps: f|S has only nonnegative exponents in r/9+(1/3)Z, so its lowest possible exponent is at least r/9. Its weight is k−r and character is χ₋₃(d)^r: this holds on T,F and −I and hence the group.

Let A be the theta series of the A₂ root lattice, B=U³, and C=η(3τ)^9/η(τ)^3. The theta transformation of the even determinant-three A₂ and the same eta identities make A,B,C modular on Γ₀(3), with weights 1,3,3 and character χ₋₃. They are holomorphic at both cusps. The nonzero leading coefficients of A,B,C at infinity occur at q-orders 0,0,1. At cusp zero (under S) they occur at orders 0,1/3,0. For A this also follows directly by Poisson from A₂* ≅ A₂/√3; the leading scalar is −i/√3. For B,C the leading scalars at zero are (−i)^3·3^(3/2),(−i)^3·3^(−9/2).

Write k−r=3m+s with 0≤s<3. The m+1 functions A^s B^i C^(m−i), 0≤i≤m, have the required character because s+m ≡r (mod2), and their distinct leading orders at infinity prove independence. By the valence formula for Γ₀(3), of index four, a nonzero form of weight k−r cannot vanish at infinity to order larger than (k−r)/3 (the unitary character does not change the valence formula, and all other contributions are nonnegative). Thus the dimension is at most m+1 and the displayed functions are a basis.

## Exact coefficient and density consequences

In this basis the smallest nonzero i determines the cusp-zero leading order r/9+i/3 of f; the orders are distinct. The known initial coefficients of f are (1,0,0) for rank48 and (1,0,0,0) for rank72, by the ambient minimum. Direct exact truncated-series arithmetic determines the largest possible smallest i:
- rank48: 6,5,5 for r=0,1,2;
- rank72: 9,8,8 for r=0,1,2.

The program verify.py solves the displayed rational systems exactly. The leading cusp-zero coefficient at equality is the first surviving basis coefficient times 3^((r−s+3i−9(m−i))/2); all phases multiply to (−i)^k=1. This gives the shell counts above.

For a 49D or 73D stack with full horizontal parent layers and shifts in (1/3)L/L, at the parent squared minimum (6 or8) any different layers must have vertical spacing at least sqrt2. Layers at the same height cannot be different, because their horizontal distance would be at most 4 or6. Thus the known half-shift lifts at height sqrt2 already meet this upper limit on the density obtainable in this **restricted ternary full-layer family** at the inherited radius. No claim is made about other radii, fractional orders, other parents, or general periodic packings.

Run python3 verify.py (SymPy) for the exact matrix replay.
