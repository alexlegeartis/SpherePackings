"""Verify the P48n quaternionic gluing and replay its coset exclusions."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
from fractions import Fraction
from itertools import product
import json
from pathlib import Path
import re
import subprocess
import sys
from tempfile import TemporaryDirectory

from flint import fmpq_mat, fmpz_mat
import sympy

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "mw92_extensions_20261122"))
from verify_mw92 import display_floor, gain_floor, pi_interval, unit_ball_squared

CATALOGUE = (HERE.parent / "density_breakthroughs_20261113"
             / "seed_deformations" / "sources" / "CQ48a.html")
ENUMERATOR = (HERE.parent / "covering_and_deep_cosets_20261115"
              / "explicit_p48" / "interval_enum.cpp")
FLAGS = ["-O3", "-std=c++17", "-fno-fast-math", "-ffp-contract=off", "-DCUTOFF=12"]
BASE_COMMIT = "9a9851887da20ca75825342bfb09e8ab031f4153"
ROOT_GRAM = [
    [2, 1, 0, 0, -1],
    [1, 2, 1, 1, 0],
    [0, 1, 2, 0, 0],
    [0, 1, 0, 2, 0],
    [-1, 0, 0, 0, 2],
]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def catalogue_gram():
    body = CATALOGUE.read_text().split('NAME="GRAM"')[1].split("<p><li>")[0]
    header, flag, *entries = map(int, re.findall(r"-?\d+", re.sub("<[^>]*>", " ", body)))
    assert (header, flag) == (48, 0) and len(entries) == 48*49//2
    lower = iter(entries)
    gram = fmpz_mat(48, 48)
    for i in range(48):
        for j in range(i+1):
            gram[i, j] = gram[j, i] = next(lower)
    assert gram.det() == 1
    assert all(gram[i, i] % 2 == 0 for i in range(48))
    assert all(fmpz_mat([row[:k] for row in gram.tolist()[:k]]).det() > 0
               for k in range(1, 49))
    return gram


def interval_inputs(gram, residues):
    matrix = sympy.Matrix(gram.tolist())
    lower, diagonal = matrix.LDLdecomposition(hermitian=False)
    assert lower * diagonal * lower.T == matrix
    ldl = []
    for i in range(48):
        assert diagonal[i, i] > 0
        for value in [diagonal[i, i], *(lower[j, i] for j in range(48))]:
            a, b = sympy.floor(value*2**40), sympy.ceiling(value*2**40)
            assert -2**53 < a <= value*2**40 <= b < 2**53
            ldl.append(f"{a} {b}\n")
    first_minimum = next(i for i in range(48) if gram[i, i] == 6)
    control = [int(i == first_minimum) for i in range(48)]
    line = lambda values: " ".join(map(str, values)) + "\n"
    return {
        "ldl.txt": "".join(ldl),
        "gram.txt": "".join(map(line, gram.tolist())),
        **{f"{name}.txt": line(residue) for name, residue in residues.items()},
        "control.txt": line(control),
    }


def roots():
    gram = fmpz_mat(ROOT_GRAM)
    assert gram.det() == 4
    assert all(fmpz_mat([row[:k] for row in ROOT_GRAM[:k]]).det() > 0 for k in range(1, 6))
    inverse = gram.inv()
    # A coefficient outside [-2,2] violates z_i^2 <= 2*(Gram^-1)_ii for roots.
    assert all(2*inverse[i, i] < 9 for i in range(5))
    out = []
    for z in product(range(-2, 3), repeat=5):
        row = fmpz_mat([z])
        if (row * gram * row.transpose())[0, 0] == 2:
            out.append(z)
    assert len(out) == 40 and sum(z[4] == 0 for z in out) == 24
    return out


def density_rows(parent, lifts):
    old = {row["dimension"]: Fraction(row["density_squared"])
           for name in ("p48_results.json", "section53_results.json")
           for row in json.loads((HERE.parent / "adaptive_layers_20261115" / name).read_text())["rows"]}
    pi_lower, _ = pi_interval()
    rows = []
    for r in (3, 4, 5):
        n = 48 + r
        vertical = fmpq_mat([row[:r] for row in ROOT_GRAM[:r]])
        B = fmpq_mat(lifts.tolist()[:r])
        G = fmpq_mat(parent)
        cross = G * B.transpose() / 2
        bottom = B * G * B.transpose() / 4 + vertical
        full = fmpq_mat([upper + right for upper, right in zip(G.tolist(), cross.tolist())]
                        + [left + row for left, row in zip(cross.transpose().tolist(), bottom.tolist())])
        assert full.det() == vertical.det() == 4
        center_squared = Fraction(3, 2)**n / 4
        assert center_squared > old[n]
        rows.append({
            "dimension": n,
            "kind": "explicit lattice",
            "squared_minimum": 6,
            "gram_determinant": 4,
            "vertical_gram": [[str(x) for x in row] for row in vertical.tolist()],
            "half_coset_numerators_in_catalogue_basis": [[int(x) for x in row] for row in B.tolist()],
            "center_density_squared": str(center_squared),
            "occupied_density_floor": display_floor(center_squared * unit_ball_squared(n, pi_lower)),
            "prior_repository_center_density_squared": str(old[n]),
            "gain_over_repository_percent_floor": gain_floor(center_squared, old[n]),
        })
    return rows


def setup():
    saved = json.loads((HERE / "certificate.json").read_text())
    parent = catalogue_gram()
    ident = fmpz_mat([[int(i == j) for j in range(48)] for i in range(48)])
    U, X, Y = [fmpz_mat(saved[name]) for name in ("U", "X", "Y")]
    assert X*X == Y*Y == -ident and X*Y == -Y*X
    assert 2*U == ident + X + Y + X*Y
    operators = [ident, U, X, Y]
    root_list = roots()
    for z in root_list:
        if z[4] == 0:
            operator = sum((k*A for k, A in zip(z, operators)), fmpz_mat(48, 48))
            assert operator * parent * operator.transpose() == parent

    change = fmpz_mat(saved["unimodular_basis_change"])
    assert abs(change.det()) == 1
    reduced = change * parent * change.transpose()
    residue, fifth = (saved[f"{name}_residue_in_reduced_basis"] for name in ("seed", "fifth"))
    for bits in [residue, fifth]:
        assert set(bits) <= {0, 1} and len(bits) == 48 and any(bits)
    seed = fmpz_mat([residue]) * change
    extra = fmpz_mat([fifth]) * change

    def parity(row):
        return tuple(int(x % 2) for x in row.tolist()[0])

    lifts = fmpz_mat([parity(seed*A) for A in operators] + [parity(extra)])
    seen = set()
    for color in product((0, 1), repeat=5):
        if any(color):
            row = fmpz_mat([color]) * lifts
            label = parity(row)
            assert any(label) and label not in seen
            seen.add(label)
            # A nonroot class needs depth 2; squared norms are 0 mod 4.
            assert (row * parent * row.transpose())[0, 0] % 4 == 0
    assert len(seen) == 31

    new_labels = {parity(fmpz_mat([z]) * lifts) for z in root_list if z[4]}
    assert len(new_labels) == 8
    # The signed quaternion group acts through I, X, Y, XY on binary classes.
    group = [ident, X, Y, X*Y]
    assert all(A*parent*A.transpose() == parent for A in group)
    change_inv = change.inv()
    to_reduced = fmpz_mat([[int(x) for x in row] for row in change_inv.tolist()])
    assert to_reduced * change == ident
    residues = {"seed": residue}
    remaining = set(new_labels)
    while remaining:
        label = fmpz_mat([min(remaining)])
        residues[f"fifth_{len(residues)}"] = parity(label * to_reduced)
        remaining.difference_update(parity(label * A) for A in group)
    assert len(residues) == 6
    return interval_inputs(reduced, residues), reduced, density_rows(parent, lifts)


def check_reports(reports, reduced, inputs):
    control = reports["control"]
    assert set(reports) == {name[:-4] for name in inputs if name not in ("ldl.txt", "gram.txt")}
    for name, a in reports.items():
        assert a["cutoff"] == 12 and a["shard_id"] == 0 and a["shard_count"] == 1
        assert a["prefixes"] == a["assigned_prefixes"] > 0
        if name != "control":
            assert a["status"] == "COMPLETE" and not a["found"]
            assert a["norm"] == -1 and not a["witness"]
    assert control["status"] == "EXACT_REJECT" and control["found"] and control["norm"] == 6
    z = fmpz_mat([control["witness"]])
    assert (z * reduced * z.transpose())[0, 0] == 6
    assert [int(x % 2) for x in z.tolist()[0]] == list(map(int, inputs["control.txt"].split()))


def replay(inputs):
    with TemporaryDirectory() as directory:
        path = Path(directory)
        for name, text in inputs.items():
            (path / name).write_text(text)
        binary = path / "exclude"
        subprocess.run(["g++", *FLAGS, str(ENUMERATOR), "-o", str(binary)], check=True)
        def run(name):
            expected_code = int(name == "control")
            process = subprocess.run([str(binary), "ldl.txt", "gram.txt", name + ".txt", "0", "1"],
                                     cwd=path, capture_output=True, text=True)
            if process.returncode != expected_code:
                raise RuntimeError(process.stdout + process.stderr)
            report = json.loads(process.stdout)
            del report["seconds"]
            print(name, report["status"], flush=True)
            return name, report

        names = [name[:-4] for name in inputs if name not in ("ldl.txt", "gram.txt")]
        with ThreadPoolExecutor(max_workers=4) as pool:
            return dict(pool.map(run, names))


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run without Python -O.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay", action="store_true")
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    if args.write_results and not args.replay:
        parser.error("--write-results requires a fresh --replay")

    inputs, reduced, rows = setup()
    report_path = HERE / "results.json"
    reports = replay(inputs) if args.replay else json.loads(report_path.read_text())["enumeration"]
    check_reports(reports, reduced, inputs)
    result = {
        "comparison_base_commit": BASE_COMMIT,
        "catalogue_sha256": sha256(CATALOGUE),
        "enumerator_sha256": sha256(ENUMERATOR),
        "certificate_sha256": sha256(HERE / "certificate.json"),
        "compiler_flags": FLAGS,
        "enumeration_input_sha256": {name: hashlib.sha256(body.encode()).hexdigest()
                                     for name, body in inputs.items()},
        "enumeration": reports,
        "rows": rows,
    }
    if args.write_results:
        report_path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        assert result == json.loads(report_path.read_text())
    for row in rows:
        print(row["dimension"], row["occupied_density_floor"],
              row["gain_over_repository_percent_floor"] + "%")
