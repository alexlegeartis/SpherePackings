#!/usr/bin/env python3
"""Exact incompatible lift certificates for one 28D binary residue code."""

from collections import Counter


GENERATORS = (1520783520, 1010581296, 3284385984, 252645120, 3284352816, 16711920)
COLUMNS = {j: 1 | (j << 1) for j in range(4, 32)}
ADDITIVE_CONFLICT = (4, 8, 10, 20, 22, 24, 26)
PERIODIC_CONFLICT = (
    (6, 2), (6, 3), (8, 2), (8, 3), (10, 2), (10, 5), (11, 3), (11, 10),
    (16, 6), (16, 8), (17, 5), (17, 6), (17, 8), (17, 11), (24, 5), (24, 11), (24, 16),
)


def syndrome(bits):
    answer = 0
    for j, column in COLUMNS.items():
        if (bits >> j) & 1:
            answer ^= column
    return answer


def parity(bits):
    return bits.bit_count() & 1


def verify():
    active = sum(1 << j for j in COLUMNS)
    residues = [(0, 0)]
    for generator in GENERATORS:
        assert 0 < generator and not generator & ~active and syndrome(generator) == 0
        residues += [(low ^ generator, carry ^ (low & generator)) for low, carry in residues]
    weights = Counter(low.bit_count() for low, _ in residues)
    assert weights == {0: 1, 12: 31, 16: 31, 28: 1}
    assert len({low for low, _ in residues}) == 64

    annihilator = {}
    for select, (low, _) in enumerate(residues):
        if low.bit_count() == 12:
            normals = [h for h in range(1, 64) if all(
                not parity(h & column) for j, column in COLUMNS.items() if (low >> j) & 1)]
            assert len(normals) == 1
            annihilator[select] = normals[0]

    # Sum the selected linear equations for the six generator lift syndromes.
    coefficients = [0] * 6
    rhs = 0
    for select in ADDITIVE_CONFLICT:
        h = annihilator[select]
        for j in range(6):
            if (select >> j) & 1:
                coefficients[j] ^= h
        rhs ^= 1 ^ parity(h & syndrome(residues[select][1]))
    assert coefficients == [0] * 6 and rhs == 1

    # Now let each site have its own six-bit syndrome. These pair equations
    # give a contradiction independently of the additive lift restriction.
    coefficients = [0] * 64
    rhs = 0
    for a, b in PERIODIC_CONFLICT:
        h = annihilator[a ^ b]
        coefficients[a] ^= h
        coefficients[b] ^= h
        borrow = residues[b][0] & ~residues[a][0]
        rhs ^= 1 ^ parity(h & syndrome(borrow))
    assert coefficients == [0] * 64 and rhs == 1

    # All pure-even translates are below the target, so each binary residue
    # could occur at most once in any packing with this period.
    attainable = {0}
    for _ in range(3):
        attainable |= {value ^ column for value in attainable for column in COLUMNS.values()}
    assert attainable == set(range(64))
    print("Verified additive and independent-site lift contradictions for the 28D candidate.")


if __name__ == "__main__":
    verify()
