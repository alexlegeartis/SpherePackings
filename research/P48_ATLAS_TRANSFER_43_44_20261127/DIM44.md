# 44D extension: exact canonical optimum of the fixed D4 carrier

We also applied the same true-metric atlas transfer to the primitive normal `√3 D4` on the first four standard E8 coordinates in the same specific P48p/K9 embedding. Its horizontal lattice (L) has Gram determinant (4·3^4=324) and minimum 6. Among all unions of complete canonical translates in `L*/L`, at any common radius, the exact maximum center density is
```
4(17/12)^22 / 18 = 17^22 / (2^43 3^24) ≈ 472.798775937791.
```
This equals the old `T44` density attributed to Conway–Sloane in the [Nebe–Sloane table](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/density.html). It is below the current [Chen–Hu–Li–Wang–Wu comparator](https://arxiv.org/html/2505.02394v1#S3.SSx7), `157^22/(2^22·5^(43/2)·11^23)≈509.61898264`; therefore **no new record**.

In the D4 class order analogous to D5 in the main proof, take `[0,319,322,109]`, whose nonzero doubled standard normal coordinates are `[7,5,5,5]`, `[5,5,5,5]`, `[6,0,0,4]` respectively. The six normal squared distances are all exactly `1/3`, so the antipode guarantee proves their minimum squared separation at least `17/3`. The transferred 39D vectors give matching upper witnesses. The 324 classes have exact squared depths in units 1/60:
```
0:1, 180:3, 220:48, 240:32, 260:72, 280:24, 300:96, 320:24, 340:24.
```
As before, full exact depths use the earlier lower replay, but the sharp density result only uses the literal upper witnesses and the independent normal lower certificate. The transverse enumeration now has four tail coordinates plus the same (j); their squared cost in units 1/60 is `5Σu_i²+2j²`. Normal carry, membership and group arithmetic are unchanged. Anchored cliques and exact density comparison exclude all bigger values. No assertion about other embeddings or real phases. Replay: `node verify_44.js [path-to-39D-archive]`.

