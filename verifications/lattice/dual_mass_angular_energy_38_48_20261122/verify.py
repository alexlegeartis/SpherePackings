"""Replay the pinned parent proofs and verify the 38–48D improvements."""

import argparse
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from math import ceil, floor
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

from flint import arb, ctx, fmpq, fmpq_mat, fmpq_poly


HERE = Path(__file__).resolve().parent
RANKS = range(38, 49)
ARCHIVE_SHA = "bd9490675ecc38e96c4dddbc484ee7a99bae0a231902d5f3ac974eae9d2143b6"
MANIFEST_SHA = "c757fa3830618af1d0ea1c19762f0bd54c0bbfe4c4aae866006989bcccea191d"
NINA_MANIFEST_SHA = "a5cb5e46a7e72de2090d958774f6ec0ca44f9f21ca8d1ca12f4fb11f5865c0db"
DENSITY_SCALE = 10**34


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def interval(value):
    return arb(fmpq(str(value)))


@lru_cache(maxsize=None)
def ball_volume(n):
    return arb.pi() ** interval(Fraction(n, 2)) / interval(Fraction(n + 2, 2)).gamma()


def parent(archive, n):
    folder = archive / "p"
    row = read(folder / "BOUNDS.json")[str(n)]
    flag = None
    if row["mode"] == "flag":
        flag = read(folder / row["certificate"])
        row = row["parent"]
    if row["mode"] == "baseline":
        folder /= "baseline"
        row = read(folder / "BOUNDS.json")[str(n)]
    assert row.get("mode") != "reciprocal-plane"
    name = row["certificate"]
    if "-rigorous" in name:
        name = name.split("-rigorous", 1)[0] + ".json"
    source = read(folder / name)
    assert source["n"] == n and "d" in source and not source.get("dual_gap")
    return folder, source, flag


def flag_terms(flag):
    if flag is None:
        return Fraction(0), Fraction(0)
    return Fraction(flag["geometry"]["target"]), Fraction(2 * flag["lower"], flag["unit"])


def verify_dual_mass(n, bounds, source, flag, witness):
    from audit_angular import nonnegative
    from independent_fourier import transform

    assert flag is None
    old = Fraction(bounds[str(n)]["upper"])
    target = Fraction(999, 1000) * old
    assert 0 < target < old
    assert target == Fraction(witness["threshold"])
    scale = 4 * arb.pi()**2 / interval(source["tau"])
    covolume = scale ** interval(Fraction(n, 2)) * 2**n * interval(old) / ball_volume(n)
    hermite = 4 * (interval(old) / ball_volume(n)) ** interval(Fraction(2, n))
    predecessor = Fraction(bounds[str(n - 1)]["upper"])
    section_ratio = 2 * ball_volume(n - 1) * interval(target / predecessor) / ball_volume(n)
    a_low = scale * section_ratio**2
    a_high = hermite * covolume ** interval(Fraction(2, n))
    lo, hi = map(Fraction, witness["squared_length_interval"])
    assert 0 < lo < hi
    assert interval(lo) <= a_low and interval(hi) >= a_high

    _, q = transform(n, source["d"])
    denominator = Fraction(str(q(0)))
    payment = Fraction(witness["payment"])
    assert denominator > 0 and payment > 0
    breaks = [lo, *map(Fraction, range(floor(lo) + 1, ceil(hi))), hi]
    for left, right in zip(breaks, breaks[1:]):
        exponential_upper = (interval(right) / 2).exp().upper().fmpq()
        assert nonnegative(q - fmpq(str(payment)) * exponential_upper, left, right)
    return max(target, old * denominator / (denominator + 2 * payment))


def principal_minors_nonnegative(matrix):
    assert matrix == matrix.transpose()
    m = matrix.nrows()
    return all(fmpq_mat([[matrix[i, j] for j in ids] for i in ids]).det() >= 0
               for size in range(1, m + 1) for ids in combinations(range(m), size))


def verify_angular_energy(archive, n, bound, folder, source, flag, witness):
    from audit_angular import nonnegative
    from portable_local import restriction, zonals

    index = witness["atom_index"]
    atom = source["kernels"][index]
    if "local_kernel" in atom:
        path = folder / atom["local_kernel"]
    else:
        path = next(p for p in (folder / "upstream").rglob("frozen.json") if digest(p) == atom["source_sha"])
    assert digest(path) == atom["source_sha"]
    data = read(path)
    origin_divisor = fmpq(data["q"])
    if bound["method"] == "origin-reallocation":
        evidence = read(archive / bound["evidence"])
        if index == evidence["atom_index"]:
            path = archive / "kernels" / evidence["kernel"]
            assert digest(path) == evidence["kernel_sha"]
            data = read(path)
    assert not fmpq(str(data.get("angular_debit", 0)))
    assert digest(path) == witness["kernel_sha"]
    radius = Fraction(data.get("radius", "19/13"))
    dilation = Fraction(atom["dilation"])
    cutoff = Fraction(atom.get("cutoff", str(radius / dilation)))
    assert 1 <= dilation and 1 < cutoff <= min(Fraction(2), radius / dilation)
    angular_cap = 1 - 1 / (2 * cutoff)

    powers = data["powers"]
    m = len(powers)
    step = powers[1] - powers[0]
    assert m > 1 and step > 0 and powers[0] > 0
    assert powers == [powers[0] + j * step for j in range(m)]
    v = fmpq_mat(m, 1, list(map(fmpq, witness["direction"])))
    ell = fmpq_mat(m, 1, list(map(fmpq, witness["radial_lower"])))
    h = fmpq_poly([0, 1])
    v_shape = fmpq_poly([v[j, 0] for j in range(m)])(h**step)
    ell_shape = fmpq_poly([ell[j, 0] for j in range(m)])(h**step)
    assert nonnegative(v_shape, 0, 1)
    assert nonnegative(h**powers[0] * v_shape**2 - ell_shape, 0, 1)

    alpha = fmpq(atom["alpha"])
    assert alpha > 0
    blocks = [fmpq_mat([list(map(fmpq, row)) for row in matrix]) * alpha * fmpq(weight) / fmpq(data["q"])**2
              for matrix, weight in zip(data["B"], data["co"])]
    weights = restriction(n, data["n"], len(blocks) - 1)
    restricted = [fmpq_mat(m, m) for _ in blocks]
    for matrix, coefficients in zip(blocks, weights):
        for j, coefficient in enumerate(coefficients):
            restricted[j] += matrix * coefficient
    coefficients = list(map(fmpq, witness["angular_coefficients"]))
    assert len(coefficients) == len(restricted) and min(coefficients) >= 0
    angular = sum((c * p for c, p in zip(coefficients, zonals(n, len(coefficients) - 1))), fmpq_poly())
    if witness.get("angular_square"):
        constant, factors = angular.factor()
        assert constant > 0 and all(multiplicity % 2 == 0 for _, multiplicity in factors)
    else:
        assert nonnegative(angular, 0, angular_cap)
    outer = v * v.transpose()
    rest = [matrix - c * outer for matrix, c in zip(restricted, coefficients)]
    assert all(principal_minors_nonnegative(matrix) for matrix in rest)
    assert rest[0].det() > 0 and restricted[0].det() > 0

    z = fmpq_mat(m, 1, [fmpq(x) / origin_divisor for x in atom["z"]])
    shifted = z + sum(coefficients) * ell
    old_allowance = (z.transpose() * restricted[0].solve(z))[0, 0]
    new_allowance = (shifted.transpose() * rest[0].solve(shifted))[0, 0]
    saving = old_allowance - new_allowance
    assert saving == fmpq(witness["origin_saving"]) and saving > 0
    threshold, flag_mass = flag_terms(flag)
    denominator = sum(Fraction(x) * (-1)**j for j, x in enumerate(source["d"])) + flag_mass
    assert denominator > 0
    factor = (interval(source["tau"]) / 8) ** interval(Fraction(n, 2)) / interval(Fraction(n + 2, 2)).gamma()
    debit = (factor * arb(saving) / interval(denominator) * DENSITY_SCALE).lower().floor().unique_fmpz()
    assert debit is not None and debit > 0
    return max(threshold, Fraction(bound["upper"]) - Fraction(int(debit), DENSITY_SCALE))


def verify(archive, dims, nina_manifest):
    assert digest(archive / "MANIFEST.json") == MANIFEST_SHA
    rows = read(HERE / "certificates.json")
    assert sorted(map(int, rows)) == list(RANKS)
    if nina_manifest is not None:
        assert digest(nina_manifest) == NINA_MANIFEST_SHA
        nina = {r["n"]: Fraction(r["claims"]["density"]) for r in read(nina_manifest)["certificates"]}
        assert all(Fraction(rows[str(n)]["nina_upper"]) == nina[n] for n in RANKS)

    needed = set(dims)
    for n in dims:
        if n in (42, 48):
            needed.update({n - 1, n - 2})
        elif rows[str(n)]["method"] == "dual-mass":
            needed.add(n - 1)
    subprocess.run([sys.executable, "-B", "-u", str(archive / "verify.py"), *map(str, sorted(needed))],
                   cwd=archive, check=True)
    sys.path.insert(0, str(archive / "p"))
    ctx.prec = 512
    bounds = read(archive / "BOUNDS.json")
    for n in dims:
        row = rows[str(n)]
        folder, source, flag = parent(archive, n)
        witness = row["witness"]
        if row["method"] == "dual-mass":
            value = verify_dual_mass(n, bounds, source, flag, witness)
        else:
            assert row["method"] == "angular-energy"
            value = verify_angular_energy(archive, n, bounds[str(n)], folder, source, flag, witness)
        old = Fraction(row["benchmark_upper"])
        assert old == Fraction(bounds[str(n)]["upper"])
        assert 0 < value == Fraction(row["upper"]) <= Fraction(row["display"]) < old
        assert value <= Fraction(99, 100) * Fraction(row["nina_upper"])
        print(f"Verified new rank {n} upper {row['display']} ({row['method']})", flush=True)


def main():
    if not __debug__:
        raise RuntimeError("Verification requires assertions.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("dimensions", nargs="*", type=int)
    parser.add_argument("--nina-manifest", type=Path)
    args = parser.parse_args()
    dims = args.dimensions or list(RANKS)
    assert all(n in RANKS for n in dims)
    archive = args.archive.resolve()
    if archive.is_dir():
        verify(archive, dims, args.nina_manifest)
    else:
        assert digest(archive) == ARCHIVE_SHA
        with TemporaryDirectory() as temp:
            with ZipFile(archive) as proof:
                proof.extractall(temp)
            verify(Path(temp), dims, args.nina_manifest)


if __name__ == "__main__":
    main()
