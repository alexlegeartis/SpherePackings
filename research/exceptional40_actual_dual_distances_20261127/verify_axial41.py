#!/usr/bin/env python3
"""Exact endpoint audit for the 41D axial bound."""
from fractions import Fraction as F

ends = [F(0), F(15, 4), F(41, 10), F(47, 10), F(107, 20), F(6)]
caps = [(F(0), 6561), (F(10, 3), 2187), (F(4), 81),
        (F(14, 3), 9), (F(16, 3), 1)]
for i, (t, m) in enumerate(caps):
    for x in ends[i:i + 2]:
        if x == 0:
            assert t == 0
            b = F(0)
        else:
            assert x > t
            b = m**2 * (x / 6)**41 * F(2, 3) / (x - t)
        assert b == 1 if x == 6 else b < 1
assert F(3, 2)**41 / (81**2 * F(2, 3)) == F(3**17, 2**21)**2
print("PASS: 41D axial interval certificates and normalization.")
