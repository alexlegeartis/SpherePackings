"""Rebuild the Gamma72 gluing and prove all required half-coset exclusions."""

import argparse
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import isqrt
from pathlib import Path
import re
import subprocess
import sys
from tempfile import TemporaryDirectory

from flint import fmpq_mat, fmpz_mat

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "mw92_extensions_20261122"))
from verify_mw92 import display_floor, gain_floor, pi_interval, unit_ball_squared

SOURCE = HERE.parent / "novel_routes_20261113" / "higher_parent" / "neb72.html"
BASELINE = HERE.parent / "range48_96_20261115" / "results.json"
COMPARISONS = HERE.parent / "covering_and_deep_cosets_20261115" / "literature_comparisons.json"
BASELINE_COMMIT = "0974b5afeac553a218a66a6f5340d926ecc8b01f"
FLAGS = ["-O3", "-std=c++17"]


def identity(n):
    return fmpz_mat([[int(i == j) for j in range(n)] for i in range(n)])


def integral(matrix):
    assert all(x.denominator == 1 for row in matrix.tolist() for x in row)
    return fmpz_mat([[int(x) for x in row] for row in matrix.tolist()])


def block(matrix, r, c, n):
    return fmpz_mat([[matrix[r+i, c+j] for j in range(n)] for i in range(n)])


def blocks(rows):
    return fmpz_mat([sum(part, []) for row in rows
                     for part in zip(*(matrix.tolist() for matrix in row))])


def positive(matrix):
    data = matrix.tolist()
    cls = type(matrix)
    for k in range(1, matrix.nrows()+1):
        assert cls([row[:k] for row in data[:k]]).det() > 0


def catalogue_gram():
    body = SOURCE.read_text().split('NAME="GRAM"', 1)[1].split('NAME="DET"', 1)[0]
    dim, *values = map(int, re.findall(r"-?\d+", re.sub("<[^>]*>", " ", body)))
    assert dim == 72 and len(values) == dim*dim
    gram = fmpz_mat(dim, dim, values)
    assert gram == gram.transpose() and gram.det() == 1
    assert all(gram[i, i] % 2 == 0 for i in range(dim))
    positive(gram)
    return gram


def parity(row):
    return tuple(int(x) % 2 for x in row.tolist()[0])


def glue(gram, change):
    # Recover the orthogonal coordinates from the published Gram.
    ident = identity(24)
    zero = 0*ident
    original = integral(fmpq_mat(block(gram, 0, 0, 24))/2)
    a_original = integral(block(gram, 0, 48, 24)*original.inv())
    b_original = ident-a_original
    assert abs(change.det()) == 1
    inverse = integral(change.inv())
    leech = change*original*change.transpose()
    a = change*a_original*inverse
    b = ident-a
    assert a*a-a+2*ident == zero
    assert a*leech*a.transpose() == 2*leech
    transform = blocks([[ident, ident, a_original],
                        [zero, b_original, b_original],
                        [zero, zero, 2*ident]])
    transform *= blocks([[inverse, zero, zero], [zero, inverse, zero], [zero, zero, inverse]])
    ambient = blocks([[leech, zero, zero], [zero, leech, zero], [zero, zero, leech]])
    assert transform*ambient*transform.transpose() == 2*gram
    assert abs(transform.det()) == 2**36
    integral(2*transform.inv())

    # Both idempotents have rank 12 mod 2, checked when making their masks.
    # The index then proves these congruences describe the whole image.
    for row in transform.tolist():
        u, v, w = [fmpz_mat([row[i:i+24]]) for i in (0, 24, 48)]
        assert parity((u-v)*a) == parity((u-w)*a) == (0,)*24
        assert parity((u+v+w)*b) == (0,)*24
    return leech, a, b, transform


def masks(matrix):
    """Select independent output columns of a binary projection."""
    pivots = {}
    chosen = []
    for j in range(24):
        col = sum((int(matrix[i, j]) % 2) << i for i in range(24))
        while col:
            p = col.bit_length()-1
            if p not in pivots:
                pivots[p] = col
                chosen.append(j)
                break
            col ^= pivots[p]
    assert len(chosen) == 12
    return [sum((int(matrix[i, j]) % 2) << k for k, j in enumerate(chosen))
            for i in range(24)]


def frame_residues(to_frame):
    """Generate the complete image of the Leech basis modulo four."""
    def planes(row):
        return tuple(sum(((int(v) >> b) & 1) << i for i, v in enumerate(row))
                     for b in range(2))

    def add(u, v):
        return u[0] ^ v[0], u[1] ^ v[1] ^ (u[0] & v[0])

    image = {(0, 0)}
    for row in to_frame.tolist():
        word = planes(row)
        image = {add(u, shift) for u in image
                 for shift in ((0, 0), word, add(word, word), add(add(word, word), word))}
    assert len(image) == 8192
    return [[(lo >> i & 1) + 2*(hi >> i & 1) for i in range(24)] for lo, hi in sorted(image)]


def line(values):
    return " ".join(map(str, values)) + "\n"


def decoder_input(leech, a, b, frame):
    assert frame.nrows() == frame.ncols() == 24
    assert frame*leech*frame.transpose() == 8*identity(24)
    forward = leech*frame.transpose()
    assert forward*frame == 8*identity(24)
    image = frame_residues(forward)
    return ("".join(map(line, frame.tolist() + forward.tolist() + [masks(a), masks(b)]))
            + line([len(image)]) + "".join(map(line, image)))


def required_cutoffs(gram, numerators, vertical, scale, orbit_representative):
    r = vertical.nrows()
    assert vertical == vertical.transpose()
    positive(vertical)
    inverse = vertical.inv()
    # If zGz^T < 8, Cauchy-Schwarz gives z_i^2 < 8*(G^-1)_ii.
    limits = [isqrt(int(8*scale*inverse[i, i])) for i in range(r)]
    basis = fmpz_mat(numerators.tolist()[:r])
    data = vertical.tolist()
    color_minima = {}
    for z in product(*(range(-a, a+1) for a in limits)):
        if not any(z):
            continue
        norm = sum(z[i]*data[i][j]*z[j] for i in range(r) for j in range(r))
        if norm < 8*scale:
            color = tuple(x % 2 for x in z)
            assert any(color)
            color_minima[color] = min(color_minima.get(color, 8*scale), norm)

    cutoffs = {}
    seen = set()
    for color in product((0, 1), repeat=r):
        if not any(color):
            continue
        row = fmpz_mat([color])*basis
        label = parity(row)
        assert any(label) and label not in seen
        seen.add(label)
        norm = color_minima.get(color, 8*scale)
        residue = int((row*gram*row.transpose())[0, 0]) % 4
        cut = max((k for k in range(8, 32)
                   if k % 4 == residue and k*scale < 4*(8*scale-norm)), default=0)
        if cut:
            key = orbit_representative(label)
            cutoffs[key] = max(cutoffs.get(key, 0), cut)
    return cutoffs


def density_rows(parent, B, vertical, scale):
    previous = {row["dimension"]: Fraction(row["best_checked"]["density_squared"])
                for row in json.loads(BASELINE.read_text())["rows"]}
    external = {row["dimension"]: Fraction(row["comparison_density_squared"])
                for row in json.loads(COMPARISONS.read_text())["rows"]
                if 73 <= row["dimension"] <= 77}
    reference78 = HERE.parent / "packing78_verified_20261117" / "packing78_independent_20261117" / "result.json"
    comparator = json.loads(reference78.read_text())
    previous[78] = Fraction(comparator["center_density_squared"])
    external[78] = Fraction(comparator["external_center_density_squared"])
    pi_lower, _ = pi_interval()
    rows = []
    P = fmpq_mat(parent)
    for r in range(1, vertical.nrows()+1):
        n = 72+r
        V = fmpq_mat(block(vertical, 0, 0, r))/scale
        numerators = fmpq_mat(B.tolist()[:r])
        cross = P*numerators.transpose()/2
        bottom = numerators*P*numerators.transpose()/4+V
        full = fmpq_mat([upper+right for upper, right in zip(P.tolist(), cross.tolist())]
                         + [left+lower for left, lower in zip(cross.transpose().tolist(), bottom.tolist())])
        det = Fraction(str(V.det()))
        assert full.det() == V.det()
        squared = Fraction(2**n)/det
        assert squared > external[n]
        if n != 78:
            assert squared > previous[n]
        row = {
            "dimension": n,
            "kind": "explicit lattice",
            "minimum_squared": 8,
            "gram_determinant": str(det),
            "center_density_squared": str(squared),
            "occupied_density_floor": display_floor(squared*unit_ball_squared(n, pi_lower)),
            "gain_over_repository_percent_floor": gain_floor(squared, previous[n]),
            "named_external_comparison_center_density_squared": str(external[n]),
            "gain_over_named_external_comparison_percent_floor": gain_floor(squared, external[n], 6),
        }
        rows.append(row)
    return rows


def decode(geometry, transform, tests):
    with TemporaryDirectory() as directory:
        p = Path(directory)
        (p/"geometry.txt").write_text(geometry)
        anchors = []
        for label, cutoff, expected in tests:
            images = (fmpz_mat([label])*transform).tolist()[0]
            anchors.append(line([2*cutoff, *(int(x) % 4 for x in images)]))
        (p/"anchors.txt").write_text("".join(anchors))
        subprocess.run(["g++", *FLAGS, str(HERE/"decode.cpp"), "-o", str(p/"decode")], check=True)
        readings = subprocess.check_output([str(p/"decode"), str(p/"geometry.txt"), str(p/"anchors.txt")],
                                           text=True).splitlines()
    assert len(readings) == len(tests)
    reports = []
    for (label, cutoff, expected), report in zip(tests, readings):
        status, best, pairs, *blocks = report.split()
        assert status == expected, (cutoff, report)
        assert int(best) == (2*cutoff if status == "FOUND" else 2*cutoff+1)
        assert len(blocks) == 3
        reports.append({
            "class": "".join(map(str, label)),
            "gamma_norm_cutoff": cutoff,
            "status": status,
            "pairs_checked": int(pairs),
            "block_norm_vectors": [int(block.split(":")[1]) for block in blocks],
        })
    return reports


def verify():
    saved = json.loads((HERE/"certificate.json").read_text())
    parent = catalogue_gram()
    change = fmpz_mat(saved["leech_basis_change"])
    leech, a, b, transform = glue(parent, change)
    frame = fmpz_mat(saved["leech_frame"])
    geometry = decoder_input(leech, a, b, frame)
    B = fmpz_mat(saved["half_coset_numerators_in_parent_basis"])
    assert all(set(row) <= {0, 1} for row in B.tolist())
    scale = saved["vertical_gram_scale"]
    assert isinstance(scale, int) and scale > 0
    vertical = fmpz_mat(saved["vertical_gram_integer"])
    assert B.nrows() == vertical.nrows() == vertical.ncols()
    ident, zero = identity(24), fmpz_mat(24, 24)
    s = fmpz_mat(saved["order_three_block"])
    S = blocks([[s, zero, zero], [zero, s, zero], [zero, zero, s]])
    assert s*s*s == ident and S*parent*S.transpose() == parent
    w = fmpz_mat([B.tolist()[0]])
    iterates = (identity(72), S, S*S)
    orbit = tuple(parity(w*operator) for operator in iterates)
    assert len(set(orbit)) == 3 and parity(fmpz_mat([B.tolist()[1]])) == orbit[1]
    assert all((x+y) % 2 == z for x, y, z in zip(*orbit))
    assert (w*parent*w.transpose())[0, 0] % 4 == 2
    for entry in B.tolist()[2:4]:
        radical = fmpz_mat([entry])-w
        assert parity(radical*S) == parity(radical)

    def orbit_representative(label):
        vector = fmpz_mat([label])
        return min(parity(vector*operator) for operator in iterates)

    cutoffs = {}
    for k in range(1, B.nrows()+1):
        local = required_cutoffs(parent, B, block(vertical, 0, 0, k), scale, orbit_representative)
        for label, cut in local.items():
            cutoffs[label] = max(cutoffs.get(label, 0), cut)
    assert cutoffs[min(orbit)] == 18
    assert sum(k == 18 for k in cutoffs.values()) == 4
    control = (1,)+(0,)*71
    assert parent[0, 0] == 8
    cases = [(control, 7, "EXCLUDED"), (control, 8, "FOUND")]
    cases.extend((label, cutoff, "EXCLUDED") for label, cutoff in sorted(cutoffs.items()))

    return {
        "comparison_base_commit": BASELINE_COMMIT,
        "hashes": {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                   for path in (SOURCE, HERE/"certificate.json", HERE/"decode.cpp")},
        "enumeration": decode(geometry, transform, cases),
        "rows": density_rows(parent, B, vertical, scale),
    }


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run without Python -O.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    result = verify()
    path = HERE/"results.json"
    if args.write_results:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        assert result == json.loads(path.read_text())
    for row in result["rows"]:
        print(row["dimension"], row["occupied_density_floor"],
              row["gain_over_repository_percent_floor"] + "%")
