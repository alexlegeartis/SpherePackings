"""Count the 54D coset contacts and verify the 55D local optimality stress."""

from fractions import Fraction
import json
from pathlib import Path
import sys

import numpy as np
from flint import fmpq, fmpq_mat, fmpz_mat, fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx, nmod_mat

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "cm58_20261116"))
from verify_frontier import parents54_55


def closest_contacts(parent, q):
    """Enumerate normalized y polynomials and recover the long-gap contacts."""
    source = HERE.parent / "cm58_20261116" / "frontier"
    raw = json.loads((source / "polynomial_points54.json").read_text())
    assert raw["field_modulus"] == [2, 2, 1, 0, 2, 0, 1]
    field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(3)(raw["field_modulus"]))
    ring = fq_default_poly_ctx(field)
    t = ring.gen()
    twist = t**28
    omega = (t*t+1).roots()[0][0]

    elements = [
        field([(n//3**k) % 3 for k in range(6)]) for n in range(729)
    ]
    subfield = [v for v in elements if v**9 == v]
    assert len(subfield) == 9
    source_points = [(ring([elements[n] for n in row["x"]]),
                      ring([elements[n] for n in row["y"]])) for row in raw["points"]]
    with np.load(source / "local54_55.npz", allow_pickle=False) as data:
        coordinates = fmpz_mat(data["Y"].tolist())

    # The source audit has already bound these points to their lattice.
    reduced, _, rank = coordinates.transpose().rref()
    assert rank == 54
    ids = [next(j for j, value in enumerate(row) if value) for row in reduced.tolist()[:rank]]
    basis = fmpz_mat([coordinates.tolist()[i] for i in ids])
    basis_points = [source_points[i] for i in ids]
    right = (parent*basis.transpose()).inv()
    depth = Fraction(28, 3)

    def pairing(a, b, norm_a):
        x, y = a
        u, v = b
        if x == u:
            assert norm_a == 10 and (y == v or y == -v)
            return norm_a if y == v else -norm_a
        denominator = (x-u)**2
        numerator = (y-v)**2-(x+u)*denominator
        height = max(Fraction(numerator.degree()), depth+denominator.degree())
        height -= numerator.gcd(denominator).degree()
        return (height-norm_a-10)/2

    def rows(points, norm, scale=1):
        return fmpq_mat([[str(scale*pairing(a, b, norm)) for b in basis_points] for a in points])

    def isometry(points):
        matrix = basis.inv()*rows(points, 10)*right
        assert matrix*parent*matrix.transpose() == parent
        assert all(v.denominator == 1 for row in matrix.tolist() for v in row)
        return fmpz_mat([[int(v) for v in row] for row in matrix.tolist()])

    maps = [isometry([(x+omega, y) for x, y in basis_points])]
    for j in range(6):
        a = field.gen()**j
        shift = t+a
        c = (t**3+t-a**28).roots()[0][0]
        correction = c+a.frobenius(3)*t-a.frobenius(4)*t**3+a.frobenius(5)*t**9
        assert correction**3+correction == twist.compose(shift)-twist
        maps.append(isometry([(x.compose(shift)+correction, y.compose(shift)) for x, y in basis_points]))
    q3 = fmpz_mat([[int(value) for value in (3*q).tolist()[0]]])
    assert rows([(ring.zero(), t**14)], depth, 3)*right == q3

    powers = []
    for matrix in maps:
        identity = matrix**0
        assert matrix**3 == identity
        assert all(int(v) % 3 == 0 for v in (q3*matrix-q3).tolist()[0])
        powers.append((identity, matrix, matrix**2))

    # Use the inverse action on Q to test contacts without materializing the
    # entire coset shell. The words run through rho and the six translations.
    orbit = [(q3, ())]
    for block in powers:
        orbit = [(vector*matrix, word+(j,))
                 for vector, word in orbit for j, matrix in enumerate(block)]
    assert len({tuple(map(int, vector.tolist()[0])) for vector, _ in orbit}) == 2187
    orbit_dual = fmpz_mat([vector.tolist()[0] for vector, _ in orbit])*parent

    def parameters():
        for u in elements:
            for v in subfield:
                if u.is_zero():
                    yield u, v, field.zero(), field.zero()
                else:
                    quadratic = t*t+(v**6/u**3)*t+u**15+u
                    for cube, _ in quadratic.roots():
                        w = cube.frobenius(5)
                        yield u, v, w, -v*w/u

    def normalized_points():
        for u, v, w, z in parameters():
            y = t**14+u*t**10+v*t**7+u**2*t**6+w*t**4-u*v*t**3-u**3*t**2+z*t+v*v-u*w
            target = y*y-twist
            x = [field.zero() for _ in range(10)]
            for j in range(1, 28):
                coefficient = target[j]
                if j % 3 == 0:
                    coefficient -= x[j//3]**3
                if j <= 9:
                    x[j] = coefficient
                elif not coefficient.is_zero():
                    break
            else:
                constants = (t**3+t-target[0]).roots()
                if constants:
                    assert len(constants) == 3
                    x[0] = constants[0][0]
                    assert y*y == ring(x)**3+ring(x)+twist
                    yield ring(x), y

    points = list(normalized_points())
    assert len(points) == 1143 and len({str(y) for _, y in points}) == len(points)
    normalized = rows(points, depth, 3)*right
    assert all(v.denominator == 1 for row in normalized.tolist() for v in row)
    seeds = fmpz_mat([[int(v) for v in row] for row in normalized.tolist()])
    seed_rows = seeds.tolist()
    for row in seed_rows:
        vector = fmpz_mat([row])
        assert (vector*parent*vector.transpose())[0, 0] == 84
        assert all(int(v) % 3 == 0 for v in (vector-q3).tolist()[0])

    dots = np.array((seeds*orbit_dual.transpose()).tolist(), dtype=np.int64)
    long_contacts = set()
    for i, j in zip(*np.where(dots == -42)):
        vector = fmpz_mat([seed_rows[int(i)]])
        word = orbit[int(j)][1]
        for exponent, block in zip(reversed(word), reversed(powers)):
            vector = vector*block[-exponent % 3]
        assert (vector*parent*q3.transpose())[0, 0] == -42
        long_contacts.add(tuple(map(int, vector.tolist()[0])))
    assert len(long_contacts) == np.count_nonzero(dots == -42) == 584
    m3 = fmpz_mat(sorted(long_contacts))
    directions = m3*2+fmpz_mat(m3.nrows(), 1, [1]*m3.nrows())*q3
    assert nmod_mat(directions.tolist(), 101).rank() == 53
    assert all(tuple(-a-b for a, b in zip(row, q3.tolist()[0])) in long_contacts
               for row in long_contacts)
    return len(points)*len(orbit), fmpq_mat(directions)/6


def verify_stress(parent, q, roots, directions):
    """Construct positive parent-root weights and check their tensor identity."""
    period_weight = fmpq(3, 76)
    short_isotropic_weight = period_weight*fmpq(8, 5)*fmpq(28, 3)/54
    target = (fmpq(1, 2)-short_isotropic_weight)*parent.inv()
    target -= period_weight*(fmpq(22, 15)*directions.transpose()*directions/directions.nrows()
                            +fmpq(3, 10)*q.transpose()*q)

    root_rows = np.asarray(roots.tolist(), dtype=np.int64)
    left, right = np.triu_indices(54)
    features = fmpz_mat((root_rows[:, left]*root_rows[:, right]).T.tolist())
    target_column = fmpq_mat([[target[int(i), int(j)]] for i, j in zip(left, right)])
    config = json.loads((HERE / "stress_seed.json").read_text())
    seed = fmpz_mat(roots.nrows(), 1)
    for i, value in config["root_weights"]:
        assert 0 <= i < roots.nrows() and seed[i, 0] == 0 and value > 0
        seed[i, 0] = value
    denominator = config["denominator"]
    rhs = target_column - fmpq_mat(features*seed)/denominator
    rhs_num, scale = rhs.numer_denom()
    transpose = features.transpose()
    coefficients = (features*transpose).solve(rhs_num)/scale
    weights = transpose*coefficients+fmpq_mat(seed)/denominator
    assert features*weights == target_column
    assert all(row[0] > 0 for row in weights.tolist())


def main():
    _, (parent, q, roots, closest) = parents54_55(return_data=True)
    assert closest.nrows() == 2187 and roots.nrows() == 8748
    count, directions = closest_contacts(parent, q)
    assert count == 2499741
    print(f"Verified {count:,} closest 54D coset points and {directions.nrows()} long-gap contacts.",
          flush=True)
    verify_stress(parent, q, roots, directions)
    print("Verified positive contact stresses for the strict 55D local maximum.")


if __name__ == "__main__":
    main()
