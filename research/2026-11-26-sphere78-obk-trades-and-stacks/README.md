# Exact continuous obstructions around the 11D OBK code packings

These are scoped theorems, not improved density bounds. The input is the 96 published (11,72,4) codes from the [authors' catalogue](https://users.aalto.fi/~pat/72.html), with an exact digest check. They realize center density 9/256 at minimum distance 2.

- [One-to-two real trade](single_trade.md): after deleting one periodic word from any of those 96 cores, there cannot be two arbitrary real replacements in the same 2Z¹¹ cell. In the singleton-hole types the upper intensity conclusion also applies to nonperiodic refilling with that fixed remaining core.
- [Two-to-three real trade](two_trade.md): all 245,376 choices of two periodic deletions also fail to admit three arbitrary real replacements. Unchanged 70-site core and common 2Z¹¹ period only.
- [Cover radius and whole-layer optimum](stacking.md): every core has covering radius √3. Arbitrarily oriented and shifted full copies of the cores, stacked in parallel in 12D at unit radius, need gaps at least 1, so center density ≤9/256. Alternating an explicit code and its complement attains it. This does not bound partially reconstructed, differently scaled or arbitrary 12D packings.

The Python drivers compile exact C++17/OpenMP checks and accept only rational/integer certificates. Floating proposal steps are bounded and cannot make an exclusion pass without the exact integer check. A portable replay can download the pinned original data or use a local copy:
```sh
python verify_single_trade.py
python verify_two_trade.py
python verify_stacking.py
```
