# A weighted shear selection improves the 96D lattice bound

We obtain occupied density at least **5.940330816757298e-22**, strictly more than **0.023468%** above the [previous 96D lattice bound](../mw92_eisenstein96_20261128/). The older [Bos–Conway–Sloane lattice comparator](https://neilsloane.com/doc/Me103.pdf) has exact center density 3^48/2^24; the gain over it exceeds 112.5324%. The [stronger unrestricted periodic result](../mw92_96_refined_tail_20261125/) is unchanged. This is an existence proof for a lattice with real horizontal shears, not a numerical 96D Gram matrix and not an upper bound.

Use the [same MW92 parent and audited equal-ball threshold](../mw92_eisenstein96_20261128/) A=11.00098526. In the dichotomy of that proof the generic branch is unchanged. In the fixed-order branch, the parent torus contains a class u with Tu in L, Ju=u modulo L, d(u)^2 at least 34/3, where T=I−J. The parent has squared minimum 16 and covolume 47^7.

Use Eisenstein coordinates x,z in Z[ω] and the tail form

    Q(x,z) = (16/3) N(x+ρ z) + η N(z),
    ρ = (2+ω)(1/3+3583/1000000),    η = 3869092/1000000.

The horizontal shear is x(J)u+z(J)b. The determinant-one block shear removes ρ from the determinant computation, so the fixed branch has squared tail covolume (4η)^2 = 935617056529/3906250000. At z=0 the previous argument is unchanged. Multiplication by units covers all vectors with a given nonzero z orbit. At norms 1, 3 and 4 use representatives 1, 1−ω and 2, respectively. The complete short-shell counts for these choices are 6, 1 and 0. Enumerating both coefficients of x from −4 to 4 is complete: residual coordinates outside this box have absolute value greater than 3, since both coordinates of ρz have absolute value less than 2, and N(i,j)≥(i²+j²)/2. Norm two is impossible and norms at least five have vertical norm at least 5η>16.

For norm one, the six short vectors give just three phase labels, because x(J)u modulo L depends only on i+j modulo three. The largest squared distance demand is D=5834269208999/562500000000. For norm three, the sole short x is −1. Here ρ(1−ω)=1+3d with d=3583/1000000. The required horizontal depth squared is W=16−3η−48d²=274506736333/62500000000.

Choose v with d(v+ju)^2>β for j=0,1,2, where β=10.74136; the exact inequality 3β^46<A^46 gives this from the equal-ball threshold. The five option centers are c=(0,u,−u,v,v+u). Choose arbitrary representatives and write b_i=b+T⁻¹c_i. The fifteen norm-one balls for these five options do not cover the torus since 15D^46<A^46. Thus we can first select b for which all five options pass the norm-one tests.

Put t=61433/100000. Give the first three options equal weights t/3 and the other two equal weights (1−t)/2. If all five failed the strict norm-three test, choose shortest representatives y_i of Tb+c_i−u with squared norm at most W. The weighted sum of torus pair distances squared is at most

    Σ_(i<j) w_i w_j ||y_i−y_j||²
      = Σ_i w_i ||y_i||² − ||Σ_i w_i y_i||² ≤ W.

But the four within-group squared distances are at least 34/3 and all six cross-group distances exceed β. Their weighted sum exceeds

    (34/9)t² + (34/12)(1−t)² + βt(1−t)
      = 4941121710195683/1125000000000000 > W.

Consequently at least one option passes both tests. The generic branch has the unchanged, smaller tail determinant 37137921969123580223529/156250000000000000000. The parent supplies squared minimum exactly 16 and the fixed determinant therefore supplies a worst-case density floor from π^48·2^96/(48!·47^7·sqrt(det)). The gain and the last decimal digits are checked with rational arithmetic and an outward Machin bound on π.

Run `node verify.js` for the new fixed branch and exact density comparison. Run the [previous complete checker](../mw92_eisenstein96_20261128/verify.js) for the inherited generic branch and use its [proof](../mw92_eisenstein96_20261128/) for the parent, covering threshold and dichotomy. No unrestricted table row changes.
