# 34D lattice from D6 and four E7 blocks; primitive 33D section

This lattice has squared minimum **8**, covolume **32768**, center density **4** and ordinary density **0.00000318082160059020**. This improves the recent [two-bank nonlattice A6 packing](../e7_a6_twobank34_20261125/) by **32.28756555%**, the earlier [lattice E7 section](../e7_hermitian35_20261123/) by **73.20508076%**, and historical Elkies Q34 (center 3^(33/2)/2^25) by **80.01520114%**. See the [historical Nebe–Sloane catalogue, Q33/Q34](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/density.html). These are lower bounds against the named reference constructions, not optimality or exhaustive-priority claims. A 33D lattice section also reaches 4/√3, matching the [existing unrestricted 33D A5 packing](../e7_a5_twobank33_20261125/) with a lattice, and beating the historical lattice Q33 by 3.93182483%.

## Construction

H is the same 28-dimensional lattice of E7 syndromes as in [the A5 construction](../e7_a5_twobank33_20261125/). All defining data are supplied here. Use E7 simple-root coefficient coordinates, G and the binary-linear parity-mask permutation M from cert.json. Each block parity mask M[t+64h] contributes the syndrome (t,at,h), encoded as t+64(at)+4096h in F64 × F64 × F2. The four multipliers a are 1,14,15,22. F64 uses X^6+X+1. S is the XOR of these four contributions, H_s its fiber at s and H=H_0. S is onto because any two distinct multipliers solve the field equations and the radical sets the remaining parity. For r the E7 coefficient vector of binary mask 82 take T=(0,r/2,r/2,r/2).

D6 consists of integer six-vectors y with even coordinate sum. Its phase x is a six-bit label: the five low bits are the parities of the first five coordinates and the high bit is (sum y_i)/2 modulo 2. This is an additive map with kernel 2D6; the missing coordinate parity is determined by evenness. Define f(x) as the XOR of [88,4613,5130,6396,7485,5014] selected by the six bits of x. Write b(y)=f(x(y)) and u=2212. Take

1. K = {(y,v): y in D6, v in H_{b(y)}};
2. t+K = {(y−(1/2)1,T+v): y in D6, v in H_{u XOR b(y)}}.

The second bank is one disjoint coset. The doubled local displacement −1 has phase 63 and f(63)=4096. The doubled exterior displacement has syndrome S(2T)=4096 (three radical blocks; doubling an integral representative of u contributes zero). Therefore **2t lies in K**, and the union is a lattice containing K with index two. Since det G=2, covol(H)=4·8192=32768. K projects onto D6 with kernel H, hence covol(K)=2·32768 and covol(K∪(t+K))=32768.

## Separation certificate

Reversing the mixed displacement negates all vectors and leaves norms unchanged. The two independent exact checkers give the following complete local/exterior pairs over the 64 phases in each type:

| Difference type | Number of phases | Least local squared norm | Least exterior squared norm |
| --- | ---: | ---: | ---: |
| Within | 1 | 0 | 0 |
| Within | 30 | 2 | 6 |
| Within | 31 | 4 | 4 |
| Within | 2 | 6 | 6 |
| Mixed | 32 | 3/2 | 13/2 |
| Mixed | 32 | 7/2 | 9/2 |

At the zero phase a nonzero local 2D6 vector has squared norm at least eight. In H the checker verifies all nonzero phase assignments; at the all-zero phase nonzero (2E7)^4 periods also have squared norm at least eight. A doubled root in H attains eight. Thus min=8, center density=(√2)^34/32768=4, and occupied density=4π^17/17!.

For the local formula put w=number of odd coordinates, k=(sum y_i)/2 mod 2. For within-bank classes of w>0, the signs of unit entries allow either k without changing w; for w=0,k=1 the minimum is four (one signed 2). In the mixed bank the closest coordinates are half-integers, all of size 1/2, giving w/2 mod 2 for k. If k is the other parity, one entry must move by two, for a least squared penalty of two.

Run `python verify.py` (NumPy) or `node verify_dp.js`. Both verify a positive diagonal factorization of 60G and exactly enumerate all relevant doubled-coordinate E7 phase minima. The Python checker also enumerates the original six-coordinate D6 box and, for each syndrome, independently solves the two finite-field equations and one binary equation by checking 2^15 phase solutions. The Node checker derives all 8192 syndromes simultaneously by min-plus convolution and checks the D6 local formula above. Neither separation proof uses numerical tolerances.

## 33D primitive lattice section

In the 34D lattice let the linear functional be dot product with w=(1,1,1,−1,−1,−1)/2 in its D6 factor, zero in the four E7 factors. Its value is always integral because y in D6 has even coordinate sum and w·1=0. The projection of K onto D6 includes e1+e2, on which the value is 1, so w is primitive. Slicing at w·x=0 yields a rank-33 lattice of covolume 32768√(3/2), with exact squared minimum 8 (a doubled external E7 root remains). Its center density is (√2)^33/(32768√(3/2))=4/√3. The lattice record is distinct from priority for the previously certified unrestricted A5 packing of the same density.