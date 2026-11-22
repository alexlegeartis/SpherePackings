"""Verify the characteristic-three 80D lattice and the small-index obstruction."""

from fractions import Fraction
import json
from math import lcm
from pathlib import Path

import numpy as np
from flint import fmpq_mat, fmpz_mat, fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx, nmod_mat

if not __debug__:
    raise RuntimeError("Run without -O.")


class Curve:
    """Points are represented as (x/z^2, y/z^3), with polynomial entries."""

    def __init__(self, modulus):
        self.field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(3)(modulus))
        self.ring = fq_default_poly_ctx(self.field)
        self.t = self.ring.gen()
        self.twist = self.t**41
        self.one = self.ring.one()
        self.field_degree = len(modulus)-1

    def polynomial(self, encoded):
        def decode(n):
            return self.field([(n//3**i) % 3 for i in range(self.field_degree)])
        return self.ring([decode(n) for n in encoded])

    def point(self, row):
        result = self.polynomial(row[0]), self.polynomial(row[1]), self.one
        self.check(result)
        return result

    def check(self, point):
        x, y, z = point
        assert y*y == x**3+x*z**4+self.twist*z**6
        assert x.gcd(z) == self.one and z.leading_coefficient() == self.field.one()

    def cube_root(self, polynomial):
        coefficients = polynomial.coeffs()
        assert all(c.is_zero() for i, c in enumerate(coefficients) if i % 3)
        return self.ring([c.frobenius(self.field_degree-1) for c in coefficients[::3]])

    def square_root(self, polynomial):
        degree = polynomial.degree()
        assert degree >= 0 and degree % 2 == 0
        assert polynomial.leading_coefficient() == self.field.one()
        size = degree//2
        values = [self.field.zero()]*size+[self.field.one()]
        for i in range(size-1, -1, -1):
            known = sum((values[j]*values[size+i-j] for j in range(i+1, size)),
                        self.field.zero())
            values[i] = 2*(polynomial[size+i]-known)
        result = self.ring(values)
        assert result*result == polynomial
        return result

    def inverse_phi(self, point, tau):
        x, y, z = point
        result = (self.cube_root(-x-(self.twist+tau)*z*z),
                  self.cube_root(y/tau), self.cube_root(z))
        self.check(result)
        return result

    def add_raw(self, first, second):
        if first is None:
            return second
        if second is None:
            return first
        x, y, z = first
        X, Y, Z = second
        dx = X*z*z-x*Z*Z
        dy = Y*z**3-y*Z**3
        if dx.is_zero():
            if not dy.is_zero() or y.is_zero():
                return None
            a, b = 2*y, z**3
            xa = b*b*z*z-2*x*a*a
            return xa, b*z*(x*a*a-xa)-y*a**3, a*z
        xa = dy*dy-dx*dx*(x*Z*Z+X*z*z)
        ya = dy*(x*dx*dx*Z*Z-xa)-y*(dx*Z)**3
        return xa, ya, dx*z*Z

    def add(self, first, second):
        raw = self.add_raw(first, second)
        if raw is None:
            return None
        x, y, z = raw
        factor = self.square_root(x.gcd(z*z))
        x, y, z = x//(factor*factor), y//(factor**3), z//factor
        scale = 1/z.leading_coefficient()
        result = x*scale**2, y*scale**3, z*scale
        self.check(result)
        return result

    @staticmethod
    def degree_y(point):
        if point is None:
            return 0
        _, y, z = point
        denominator = z**3
        return max(y.degree(), denominator.degree())-y.gcd(denominator).degree()

    def pairing(self, first, second):
        return Fraction(self.degree_y(self.add_raw(first, second))
                        -self.degree_y(first)-self.degree_y(second), 3)


def integer_matrix(matrix):
    rows = matrix.tolist()
    assert all(value.denominator == 1 for row in rows for value in row)
    return fmpz_mat([[int(value) for value in row] for row in rows])


def lattice_basis(rows, dimension):
    order = lcm(*(int(value.denominator) for row in rows for value in row))
    numerators = fmpz_mat([[int(value*order) for value in row] for row in rows])
    return fmpq_mat(fmpz_mat(numerators.hnf().tolist()[:dimension]))/order


def audit():
    data = json.loads(Path(__file__).with_name("certificate.json").read_text())
    curve = Curve(data["field_modulus"])
    assert curve.field_degree == 16
    points = [curve.point(row) for row in data["points"]]
    tau = (points[1][0]-points[0][0])[0]
    assert tau*tau == -curve.field.one()

    coordinates = fmpz_mat(data["Y"])
    dim = coordinates.ncols()
    assert dim == 80 and len(points) == coordinates.nrows() == 134
    reduced, _, rank = coordinates.transpose().rref()
    assert rank == dim
    indices = [next(j for j, value in enumerate(row) if value) for row in reduced.tolist()[:dim]]
    initial_rows = coordinates.tolist()
    frame = fmpz_mat([initial_rows[i] for i in indices])
    frame_points = [points[i] for i in indices]

    def cross(left, right):
        return fmpq_mat([[str(curve.pairing(p, q)) for q in right] for p in left])

    frame_inverse = frame.inv()
    gram = frame_inverse*cross(frame_points, frame_points)*frame_inverse.transpose()
    assert coordinates*gram*coordinates.transpose() == cross(points, points)
    assert gram.det() == 3**80
    assert lattice_basis(fmpq_mat(coordinates).tolist(), dim) == frame**0
    pairing_inverse = (gram*frame.transpose()).inv()

    def morphism(fn):
        transformed = [fn(*point) for point in frame_points]
        for point in transformed:
            curve.check(point)
        matrix = frame_inverse*cross(transformed, frame_points)*pairing_inverse
        assert matrix*gram*matrix.transpose() == gram
        return matrix

    rho = morphism(lambda x, y, z: (x+tau*z*z, y, z))
    identity = rho**0
    assert rho*rho+rho+identity == 0*identity
    for i in range(0, len(points), 2):
        x, y, z = points[i]
        assert points[i+1] == (x+tau, y, z)
        assert fmpq_mat([initial_rows[i]])*rho == fmpq_mat([initial_rows[i+1]])

    phi_inverse = (rho-identity).inv()
    physical_points = list(points)
    rows = fmpq_mat(coordinates).tolist()
    for expression in data["inverse_phi_sums"]:
        amount = fmpq_mat(1, dim)
        point = None
        for i, sign in expression:
            assert 0 <= i < len(physical_points) and sign in [-1, 1]
            x, y, z = physical_points[i]
            point = curve.add(point, (x, sign*y, z))
            amount += sign*fmpq_mat([rows[i]])
        point = curve.inverse_phi(point, tau)
        amount = amount*phi_inverse
        assert cross([point], frame_points)*pairing_inverse == amount
        assert 3*(amount*gram*amount.transpose())[0, 0] == 2*curve.degree_y(point)
        physical_points.append(point)
        rows.append(amount.tolist()[0])

    basis = lattice_basis(rows, dim)
    basis_inverse = basis.inv()
    final_gram = integer_matrix(basis*gram*basis.transpose())
    assert final_gram.det() == 3**48

    roots = [curve.point(row) for row in data["norm14_points"]]
    root_coordinates = integer_matrix(cross(roots, frame_points)*pairing_inverse*basis_inverse)
    for point, row in zip(roots, root_coordinates.tolist()):
        vector = fmpz_mat([row])
        assert curve.degree_y(point) == 21
        assert (vector*final_gram*vector.transpose())[0, 0] == 14

    t = curve.t
    zeta = next(root for root, _ in (t**41-1).roots() if root != curve.field.one())
    rotate = morphism(lambda x, y, z: tuple(f.compose(zeta*t) for f in (x, y, z)))
    conjugate = morphism(lambda x, y, z: tuple(curve.ring([c**3 for c in f.coeffs()])
                                               for f in (x, y, z)))
    unit = morphism(lambda x, y, z: (-x.compose(-t), tau*y.compose(-t), z.compose(-t)))

    def integral_action(matrix):
        transform = integer_matrix(basis*matrix*basis_inverse)
        assert transform*final_gram*transform.transpose() == final_gram
        return transform

    rho_final = integral_action(rho)
    rho_minus_one = rho_final-rho_final**0
    assert rho_minus_one**2 == -3*rho_final
    kernel, nullity = nmod_mat(rho_minus_one.tolist(), 3).nullspace()
    assert nullity == 40
    quotient = np.asarray(kernel.tolist(), dtype=np.int64)[:, :nullity]
    rotation, frob, quarter_turn = [
        np.array([[int(entry) % 3 for entry in row] for row in integral_action(matrix).tolist()],
                 dtype=np.int64)
        for matrix in (rotate, conjugate, unit)
    ]

    residues = set()
    for start in root_coordinates.tolist():
        a = np.array(start, dtype=np.int64) % 3
        for row in (a, a@quarter_turn % 3):
            for _ in range(16):
                v = row
                for _ in range(41):
                    label = v@quotient % 3
                    leading = label[np.flatnonzero(label)[0]]
                    residues.add(tuple(map(int, leading*label % 3)))
                    v = v@rotation % 3
                row = row@frob % 3
    assert len(residues) == 1312
    labels = np.array(sorted(residues), dtype=np.int64)
    assert nmod_mat(labels.tolist(), 3).rank() == 40

    i, j = np.triu_indices(40)
    dyads = labels[:, i]*labels[:, j] % 3
    augmented = np.column_stack((dyads, np.ones(len(labels), dtype=int)))
    assert nmod_mat(dyads.tolist(), 3).rank() == 808
    assert nmod_mat(augmented.tolist(), 3).rank() == 809

    reference = Fraction(7, 2)**80 / 41**12
    assert Fraction(7, 2)**80 / 3**48 < reference
    assert Fraction(4)**80 / 3**50 > reference > Fraction(4)**80 / (16*3**48)
    print("Verified the 80D determinant and norm-14 witnesses.")
    print("Verified 1,312 short-vector residue rays and the index-two and index-three obstructions.")


if __name__ == "__main__":
    audit()
