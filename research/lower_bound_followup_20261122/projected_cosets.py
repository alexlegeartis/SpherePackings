#!/usr/bin/env python3
"""Compare wrapped Leech projected-coset packings by exact separation bounds."""

import argparse
from fractions import Fraction
from itertools import combinations, product
import json
from math import isqrt
from pathlib import Path
import sys

import numpy as np
from sympy import Matrix, eye
from sympy.matrices.normalforms import hermite_normal_form


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "lower_bound_search_20261114" / "broad_dimensions_followup"))
from dual_anticodes import maximum_clique


def leech():
    """Return integral generator rows and the complete norm-32 shell."""
    polynomial = sum(1 << j for j in (0, 2, 4, 5, 6, 10, 11))
    binary = [polynomial << i for i in range(12)]
    binary = [row | ((row.bit_count() & 1) << 23) for row in binary]
    masks = [0]
    for row in binary:
        masks += [word ^ row for word in masks]
    words = np.asarray([[(word >> j) & 1 for j in range(24)] for word in masks], dtype=np.int64)
    octads = words[words.sum(axis=1) == 8]
    assert len(octads) == 759 and set(words.sum(axis=1)) == {0, 8, 12, 16, 24}

    units = np.eye(24, dtype=np.int64)
    odd = np.ones(24, dtype=np.int64)
    odd[0] = -3
    generators = np.vstack((4 * (units[:23] - units[23]), 8 * units[23],
                            2 * words[[1 << i for i in range(12)]], odd))
    basis = hermite_normal_form(Matrix(generators.T.tolist()))
    assert abs(basis.det()) == 8**12

    two = np.asarray([a*units[i] + b*units[j] for i, j in combinations(range(24), 2)
                      for a, b in product((-4, 4), repeat=2)], dtype=np.int64)
    signs = np.asarray([[1 - 2*((i >> j) & 1) for j in range(8)]
                        for i in range(256) if i.bit_count() % 2 == 0], dtype=np.int64)
    even_shell = np.zeros((len(octads), 128, 24), dtype=np.int64)
    for i, octad in enumerate(octads):
        even_shell[i][:, octad.nonzero()[0]] = 2 * signs
    sign_words = 1 - 2*words
    odd_shell = np.tile(sign_words, (24, 1, 1))
    for i in range(24):
        odd_shell[i, :, i] *= -3
    shell = np.vstack((two, even_shell.reshape(-1, 24), odd_shell.reshape(-1, 24)))
    assert len(shell) == 196560 and np.all((shell*shell).sum(axis=1) == 32)
    return basis, shell


def separations(normal, gram, basis, shell, extra):
    rank = len(normal)
    k = Matrix(gram.tolist())
    determinant = int(k.det())
    adjugate = np.array(k.adjugate().tolist(), dtype=np.int64)
    assert np.array_equal(normal @ normal.T, 8*gram)
    coordinates = basis.inv() * Matrix(normal.T.tolist())
    assert all(value.q == 1 for value in coordinates)
    pairing = Matrix(normal.tolist()) * basis / 8
    assert all(value.q == 1 for value in pairing)
    assert hermite_normal_form(pairing.applyfunc(int)) == eye(rank)

    def label(values):
        return tuple(map(int, values % determinant))

    zero = (0,) * rank
    classes = {zero}
    for row in adjugate:
        size = 0
        while size != len(classes):
            size = len(classes)
            classes |= {label(row + key) for key in classes}
    assert len(classes) == determinant
    classes.remove(zero)

    marks = shell @ normal.T
    assert np.all(marks % 8 == 0)
    marks = np.unique(marks // 8, axis=0)
    counted = {tuple(map(int, row)) for row in marks}
    assert zero in counted  # The orthogonal period attains minimum 4.
    upper = dict.fromkeys(classes, 4 * determinant)
    for mark in marks:
        address = mark @ adjugate
        key = label(address)
        if key != zero:
            upper[key] = min(upper[key], 4*determinant - int(address @ mark))
    for vector in extra:
        assert all(value.q == 1 for value in basis.inv()*Matrix(vector))
        vector = np.asarray(vector, dtype=np.int64)
        mark = vector @ normal.T
        assert np.all(mark % 8 == 0) and int(vector @ vector) % 8 == 0
        mark //= 8
        address = mark @ adjugate
        key = label(address)
        assert key != zero
        upper[key] = min(upper[key], int(vector @ vector)//8*determinant - int(address @ mark))

    lower = dict.fromkeys(classes, 0)
    # These finite boxes suffice for the pinned cases, as checked by equality
    # with a projected ambient witness for each class.
    limits = [isqrt(3 * int(gram[i, i])) for i in range(rank)]
    for mark in product(*(range(-limit, limit+1) for limit in limits)):
        address = np.asarray(mark, dtype=np.int64) @ adjugate
        norm = int(address @ mark)
        key = label(address)
        if norm > 3*determinant or key == zero:
            continue
        ambient_minimum = 4 if mark in counted else 6
        lower[key] = max(lower[key], ambient_minimum*determinant - norm)
    assert all(0 < lower[key] == upper[key] for key in classes)
    return determinant, lower


def largest_subgroup(allowed, determinant, zero):
    ambient = set(allowed) | {zero}
    start = frozenset({zero})
    known = {start}
    unexplored = [start]
    best = 1
    while unexplored:
        subgroup = unexplored.pop()
        for generator in ambient - subgroup:
            added = set(subgroup)
            offset = generator
            while offset not in subgroup:
                coset = {tuple((a+b) % determinant for a, b in zip(offset, key)) for key in subgroup}
                if not coset <= ambient:
                    break
                added |= coset
                offset = tuple((a+b) % determinant for a, b in zip(offset, generator))
            else:
                extended = frozenset(added)
                if extended not in known:
                    known.add(extended)
                    unexplored.append(extended)
                    best = max(best, len(extended))
    return best


def color_bound(neighbors):
    uncolored = (1 << len(neighbors)) - 1
    bound = 0
    while uncolored:
        bound += 1
        available = uncolored
        while available:
            bit = available & -available
            vertex = bit.bit_length() - 1
            uncolored ^= bit
            available ^= bit
            available &= ~neighbors[vertex]
    return bound


def best_bound(dimension, determinant, separations):
    best_periodic = best_lattice = Fraction(1, determinant)
    thresholds = sorted({min(q, 4*determinant) for q in separations.values()}, reverse=True)
    divisors = [size for size in range(1, determinant+1) if determinant % size == 0]
    for threshold in thresholds:
        allowed = sorted(key for key, norm in separations.items() if norm >= threshold)
        scale = Fraction(threshold, 4*determinant)**dimension / determinant
        if (len(allowed)+1)**2 * scale <= min(best_periodic, best_lattice):
            continue
        neighbors = []
        for key in allowed:
            row = 0
            for i, other in enumerate(allowed):
                difference = tuple((a-b) % determinant for a, b in zip(key, other))
                if separations.get(difference, 0) >= threshold:
                    row |= 1 << i
            neighbors.append(row)
        size_bound = color_bound(neighbors) + 1
        if size_bound**2 * scale > best_periodic:
            clique, _ = maximum_clique(neighbors)
            size_bound = len(clique) + 1
            best_periodic = max(best_periodic, size_bound**2 * scale)
        if max(size for size in divisors if size <= size_bound)**2 * scale > best_lattice:
            zero = (0,)*len(allowed[0])
            lattice_size = largest_subgroup(allowed, determinant, zero)
            best_lattice = max(best_lattice, lattice_size**2 * scale)
    return str(best_periodic), str(best_lattice)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    basis, shell = leech()
    rows = []
    for data in json.loads((HERE / "leech_normals.json").read_text()):
        dimension = data["dimension"]
        normal = np.asarray(data["normal"], dtype=np.int64)
        gram = np.asarray(data["gram"], dtype=np.int64)
        assert len(normal) == 24 - dimension
        determinant, distances = separations(normal, gram, basis, shell, data.get("extra_witnesses", []))
        periodic, lattice = best_bound(dimension, determinant, distances)
        rows.append({"dimension": dimension, "determinant": determinant,
                     "maximum_periodic_center_density_squared": periodic,
                     "maximum_lattice_center_density_squared": lattice})
    result = json.dumps({"rows": rows}, indent=2) + "\n"
    path = HERE / "projected_coset_results.json"
    if args.write:
        path.write_text(result)
    else:
        assert path.read_text() == result
    print("Verified lattice and periodic projected-coset searches for dimensions 19, 20, 21 and 23.")


if __name__ == "__main__":
    main()
