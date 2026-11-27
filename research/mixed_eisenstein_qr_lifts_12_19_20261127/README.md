# Mixed Eisenstein phases in 12D and a real rounding barrier for QR layers in 19D

Sphere 144, 27 November 2026. **No improved record.** Public comparison pinned at [SpherePackings main 221b1c16](https://github.com/alexlegeartis/SpherePackings/blob/221b1c1665aa4b50e006c3bdc64e4dc8b96d11c8/LOWER_BOUNDS_1_96.md). Density throughout this note is **center density**, excluding the unit-ball factor.

## 1. A 12D mixed-prime construction target

Use the orthogonal product of six hexagonal periods. In each plane choose basis $a,b$ of squared lengths 4 and inner product 2. Write a phase as
\[
 \tfrac12(\epsilon_1a+\epsilon_2b)+\tfrac t3(a+b),\qquad
 (\epsilon_1,\epsilon_2,t)\in\mathbb F_2^2\times\mathbb F_3 .
\]
For a phase difference, the local squared torus costs multiplied by 3 are $0,3,4,1$, corresponding respectively to neither component nonzero, binary only, ternary only, and both. The exact 2D calculation uses the form $x^2+xy+y^2$ at integer coordinates with denominator 6 (physical squared distance is one ninth of that form). The box in the verifier is exhaustive: outside it one coordinate exceeds 10 in absolute value, and the form is at least half the Euclidean coordinate square. This construction is a mixed alphabet, not three independent low-dimensional packings.

The full period covolume is $(2\sqrt3)^6=1728$ and its minimum squared length is 4. Thus the actual squared minimum is $D/3$ for an integer $1\le D\le12$. A motif of $M$ phases has center density
\[
 \delta=\frac{M}{1728}(D/12)^6 .
\]
At $D=12$ a winning motif needs only **65** phases. The attributed comparator is the [Coxeter–Todd form](https://doi.org/10.4153/CJM-1953-043-4), [minimum 4 and Gram determinant $3^6$](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/K12.html), center density $1/27$. The literal hexacode $(a,b,c,a+b+c,a+\omega b+\omega^2c,a+\omega^2b+\omega c)$ gives 64 binary phases and attains the comparator at $D=12$. At $D=4$ the same hexacode paired with **every** ternary tag has $64\cdot3^6=46656$ phases, again attaining the comparator: differing binary words contribute at least 1 in four or more positions, and identical binary words with different tags contribute at least 4. These are equality packings, not records or claims of new isometry classes.

### Certificate for the surviving window

Group all possible difference profiles by the counts $(n_{00},n_{10},n_{01},n_{11})$ of their local types (84 profiles). In the Fourier transform group characters by the analogous count of trivial/nontrivial components. The one-coordinate sums of characters over the four character types, with difference types as columns, are
\[
 \begin{pmatrix}
 1&1&1&1\\3&-1&3&-1\\2&2&-1&-1\\6&-2&-3&1
 \end{pmatrix}.
\]
Tensor and sum over permutations by extracting coefficients of the product of the six local polynomials. Call the resulting integer kernel $K(x,y)$ on difference and character profiles. For every motif, the mean oriented difference count $A_x$ per point has $A_0=1$, $\sum A_x=M$, and $\sum_x A_xK(x,y)\ge0$: the latter sum is a sum of absolute squares of character sums divided by $M$.

For each $D$ the provided nonnegative integers $c_y,m$ satisfy $\sum_y c_yK(x,y)\le-m$ for every allowed nonzero profile. Combining these inequalities with the nonnegative Fourier sums gives
\[
 M\le \left\lfloor1+\sum_y c_yK(0,y)/m\right\rfloor.
\]
The supplied dual witness was discovered numerically, then rounded and rescaled; **only** the literal integers and exact inequalities are used here. At $D=12$ it uses nine character profiles and bounds **every arbitrary nonlinear motif** by 69 phases. Thus the constructive range is **65–69**. Other exact integer caps and strict gain requirements are:

| $D$ | Motif upper | Motif needed to improve |
|---:|---:|---:|
| 1 | 2,985,984 | 191,102,977 |
| 2 | 995,328 | 2,985,985 |
| 3 | 199,065 | 262,145 |
| 4 | 46,656 | 46,657 |
| 5 | 12,539 | 12,231 |
| 6 | 5,061 | 4,097 |
| 7 | 1,927 | 1,625 |
| 8 | 1,030 | 730 |
| 9 | 424 | 360 |
| 10 | 184 | 192 |
| 11 | 103 | 108 |
| 12 | 69 | 65 |

No statement is made about arbitrary real phases, unequal scaling, or unrelated periods. At $D=12$, a whole hexacode plus an outside ternary-tag phase cannot work: the most economical such phase has nine short conflicts with the hexacode. These nine correspond, for a three-position ternary support $T$ and binary projection in the hexacode, to the three quaternary lines of weight-four words whose support contains $T$. A finite full-bank exploratory search over all binary rank-four subspaces of this hexacode found no five compatible banks of size 16; this search is not used by the main dual proof. At $D=6$ an instructive positive gate is at least **65 whole 64-word hexacode banks** with coupled ternary tags and quaternary coset offsets. For a three-position tag difference the compatibility is precisely that the coset's unique representative supported on those three positions has a zero entry: any nonzero binary position off that triad already contributes 3. A mass of tags with difference weight two has locked coset offsets; the tested separable repeat-triplet design therefore does not supply the gate. Genuine joint offsets, not mere translate counts, are needed.

## 2. Why parallel complete QR $B^*_{18}$ layers cannot beat 19D

Use the [Bierbrauer–Edel QR pair and Construction $B^*$](https://www.yvesedel.de/Papers/Sphere.pdf), with coordinates divided by two from the integer convention there. This 18D stock $\Gamma$ has minimum square 6 and center density $(3/4)^9$. Let $C,D$ be the parity-extended quadratic residue and nonresidue codes of length 17; they are even dual binary $[18,9,6]$ codes, intersecting in $\{0,\boldsymbol1\}$. Write the banks as $u+2a$ for $u\in C$, $\sum a_i$ even, and $\boldsymbol1/2-v+2b$ for $v\in D$, $\sum b_i$ odd. This is identical to the sign/shift convention in the cited paper since both codes contain $\boldsymbol1$ and are even.

**Theorem.** Every horizontal real translation has a pair with squared separation at most $141/32$:
\[
 d(w,\Gamma-\Gamma)^2\le141/32.
\]
Hence any arrangement of entire copies of this **same orientation** on arbitrary one-dimensional normal sites, with arbitrary translations and no required periodicity, has center density at most
\[
 (3/4)^9\,\frac{(\mu/6)^9\sqrt\mu}{2\sqrt{\mu-141/32}}
\]
for squared minimum $\mu>141/32$. This is nondecreasing when $\mu\ge2679/576$. At $\mu=6$ it gives $(3/4)^9\sqrt{16/17}$, strictly below [Chen–Hu–Li–Wang–Wu's 19D construction](https://arxiv.org/html/2505.02394v1#S3.SS4), $\delta_0=13^{19/2}/(3^9 5^{21/2})=0.088958102129\ldots$. Indeed, normal spacing must be at least $\sqrt{\mu-141/32}$, and the in-layer number density is fixed. Independent unrelated stock orientations are not assumed to share this difference set. Permutations preserving the unordered pair $\{C,D\}$ do not evade it: swapping the codes is exactly translation by $\boldsymbol1/2+2e_0$, as is seen by interchanging the two displayed banks.

**Rounding proof.** Within-bank differences contain all integer points with binary row in $(C\cup D)\setminus\{0,\boldsymbol1\}$, and the even high-digit parity for the two constant rows. Indeed for binary difference $u$ within either linear code, the high parity resulting from a carry varies by the linear functional $u\cdot v$ on that code. It is onto for every nonconstant row because the two codes are dual. Cross-bank differences contain every point $\boldsymbol1/2+p+2c$ with binary $p$ even and $\sum c_i$ odd. In fact $C+D$ is the full even binary space; the carry flag in $-(u+v)$ relative to $p=u+v\bmod2$ is $|u|+|v|+u\cdot v=0\bmod2$, leaving the odd high flag.

For a real point, denote the absolute rounding errors to the nearest integer by $a_i\in[0,1/2]$; let $S=\sum a_i$, $Q=\sum a_i^2$, $M=\max a_i$, and $l=\min a_i$. In the half-integer cross grid the absolute errors become $1/2-a_i$. Two parity flags can be repaired at extra squared cost at most $2-2(M-l)$. For an odd first parity, one coordinate suffices: a move by $+1$ or $-1$ flips the first flag and the two moves give opposite changes of the second flag. Use the cheaper move at maximum half-grid error if it has the needed second flag; otherwise the more expensive move at minimum error has cost at most $1+2\min(1/2-a_i)\le2-2(M-l)$ (the cheap move there is also allowed if correct). For even first parity with incorrect second, use two different coordinates, one cheap at maximum error, one of opposite second flag at another coordinate. If available cheaply, the same bound holds; otherwise use the expensive move at minimum error. Equal errors and rounding ties cause no problem with choosing distinct coordinates. Thus the cross distance is bounded by
\[
 6.5-S+Q-2(M-l)\le 6.5-(8/9)S+MS-2M
\]
using $Q\le MS$ and $l\le S/18$.

For the integer grid, the exact finite QR checker verifies the following for every binary nearest row outside the radius-two balls around the two constant words: either a nonconstant admissible row is within Hamming distance three, or the nearest nonconstant rows are at distance four and each coordinate is included in at least $5/28$ of their difference supports. This gives respective bounds $3+MS$ and $4+(M-5/14)S$ by choosing a nearest support or averaging all those supports. For a nearest binary row at distance one or two to a constant word, the constant word's required high parity can be repaired by reversing a sign at a mismatched coordinate; the extra cost is $4$ times its rounding error, and choosing the smallest gives total distance at most $Q+2$. If the nearest binary row is itself constant then the bound is $Q+4-4M$, by changing one coordinate by 2 if needed. This is at most 4 for $S\le4$, and for $S\ge4$ the cross bound above is at most $6.5-32/9+1<4$.

In the first generic case, maximize the lesser of the two affine bounds in $S$ using their crossing. It is at most $3+(9/8)M(3.5-2M)\le141/32$. In the second, if $M\le5/14$ the first bound is at most 4. Otherwise the crossing gives at most $4+(126/67)(M-5/14)(2.5-2M)\le4+27/67<141/32$. Both functions are increasing on the relevant interval up to $M=1/2$. The QR check enumerates all $2^{18}$ anchor words with a Hamming multi-source breadth-first search, then verifies the exact fractional support bound for every remaining far word; no geometric sample is used in this proof.

## 3. Direct 19D full-parity-row extension is also closed; partial reservoirs remain

One tempting alternative has two parity banks in dimension 19 itself: for each arbitrary binary row $u$ in the integer bank, or $v$ in the half-integer bank, choose **any** high-digit parity flag and fill that entire parity coset of $2\mathbb Z^{19}$. Nonlinear row sets and nonlinear parity flags are allowed. Either or both banks may be present. This family cannot beat $\delta_0$ at any radius. It does not include partial high-digit reservoirs or general positions.

For one binary row, the period is $2D_{19}$ of square minimum 8 and covolume $2^{20}$. Within a bank the exact cross-row minimum is binary Hamming distance: for each mismatched coordinate either sign of a unit difference is available by altering the high digit, so one mismatch can enforce the required parity. Across banks the exact pair minimum is $19/4$ or $27/4$. The all-$\pm1/2$ differences require the flag sum to be $u\cdot v\bmod2$; if parity forbids them, change the closest representative at one coordinate to $\pm3/2$, adding 2 to the squared distance. Thus $\mu>19/4$ implies affine orthogonality of the two row sets: applying the four-pair cancellation removes the two parity functions and shows $(u+u_0)\cdot(v+v_0)=0$. Their affine binary spans total at most 19, so one bank has at most 512 rows.

With $N$ full rows, $\delta=N\mu^{19/2}/2^{39}$. The discrete candidate $\mu$ values and safe row caps are checked by exact squared comparison in the QR script. The only close candidates are $\mu=6$, with at most $1237+512$ rows, $\mu=27/4$, with at most $2\cdot256$, and $\mu=8$ in a single bank, with at most 128. These use [Gijswijt–Mittelmann–Schrijver](https://arxiv.org/abs/1005.4959) for $A(19,6)\le1237$, [Schrijver](https://doi.org/10.1007/s10623-018-0589-8) for $A(20,8)=256$ (parity extension gives $A(19,7)\le256$), and [Östergård](https://aaltodoc.aalto.fi/server/api/core/bitstreams/4bbc45fc-beb8-476c-b125-a1a75221eb0e/content) for $A(19,8)=128$. At $\mu=5$ (and $19/4$) one may use $2A(20,6)\le4558$ by parity extension; at $\mu=4$ a simple puncture and Hamming bound gives at most $2\lfloor 2^{18}/19\rfloor=27594$. At $\mu\le3$ even $2^{20}$ rows lose. At $\mu=7$ only one bank can occur. In particular a one-bank winning configuration at $\mu=8$ would need the equivalent of 129 full rows, so it cannot be obtained by adding new complete rows alone. Partial top-digit patterns or nonbinary phases are necessary for this route.

Replay: run `node verify_eisenstein.js` and `node verify_qr18.js` in this directory. Both use only exact integer inequalities for proof checks.

