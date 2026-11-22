"""Exact reciprocal and radial-moment refinements of the pinned lattice bounds."""

import argparse
from fractions import Fraction as F
from functools import lru_cache
from hashlib import file_digest
import json
from math import comb
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

from flint import arb, ctx, fmpq, fmpq_mat, fmpq_poly


HERE = Path(__file__).resolve().parent
ARCHIVE_SHA = {
    "native": "4cae881e275464d6555a9beff843c887141830da8978fd4ee8b75f0b4d968b8f",
    "high": "bd9490675ecc38e96c4dddbc484ee7a99bae0a231902d5f3ac974eae9d2143b6",
}
V13_UPPER = F("0.042023513562578")
V14_UPPER = F(169984398087502876851, 5000000000000000000000)
DENSITY_UNIT = 10**27
MANIFEST_SHA = {
    "native": "fadef3b64904ab5d89b477fba2482a5b0a2e2e8669f57c3cac44a20205309d32",
    "high": "c757fa3830618af1d0ea1c19762f0bd54c0bbfe4c4aae866006989bcccea191d",
}
DIM14_SHA = {
    "replay14.py": "bee546244a8a24662391270aad69412bf2edecd1446fd52c1d8213a5a866686e",
    "n14/geometry/sharp-ref-fourier.json": "4b857c09f74b0c4f47ce618739ecbc9e35dc9717fca93e7c9d287c7ce3b81968",
    "n14/tightmatching14-final-packing-proof.json": "0b85c72fc27c26e862d5dcea6988dcfd78ad6f3d0c6b4dc985b254fde701f1fe",
}


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    with path.open("rb") as stream:
        return file_digest(stream, "sha256").hexdigest()


def rational(value):
    return fmpq(str(value))


def ball(value):
    return arb(rational(value))


@lru_cache(maxsize=None)
def volume(n):
    return arb.pi() ** ball(F(n, 2)) / ball(F(n + 2, 2)).gamma()


def outward(a, b, units=10**6):
    return (
        F(int((a * units).lower().floor().unique_fmpz()), units),
        F(int((b * units).upper().ceil().unique_fmpz()), units),
    )


def line_domain(n, tau, upper, floor, prior):
    scale = 4 * arb.pi() ** 2 / ball(tau)
    lower = scale * (2 * volume(n - 1) * ball(floor / prior) / volume(n)) ** 2
    higher = scale * (2**n * ball(upper) / volume(n)) ** ball(F(4, n))
    return outward(lower, higher)


def density_ceiling(value):
    return F(-((-value.numerator * DENSITY_UNIT) // value.denominator), DENSITY_UNIT)


def signed_lower_mesh(q, subdivisions=32, end=128, units=10**24):
    from independent_fourier import gaussian_bins
    from polynomial_tools import signed_uppers

    edges = [F(i, subdivisions) for i in range(end * subdivisions + 1)]
    signed = signed_uppers(-q, edges, 1, units)
    checked = gaussian_bins(-q, 1, edges, signed, units)
    values = [max(0, -v) for v in signed]

    @lru_cache(maxsize=65536)
    def minimum(first, last):
        return min(values[first : last + 1])

    def lower(a, b):
        first, last = int(a * subdivisions), int(b * subdivisions)
        assert 0 <= first <= last
        if last >= len(values):
            return 0
        return minimum(first, last)

    return lower, units, checked


def cover_ray(q, ranges, payment):
    from independent_fourier import nonnegative_ray

    left, right = ranges
    assert 0 < left < right and payment > 0
    assert nonnegative_ray(q, 0)
    mesh, units, analytic = signed_lower_mesh(q)
    needed = payment * units / 2
    pending = [(left, right, 0)]
    visited = 0
    while pending:
        a, b, depth = pending.pop()
        visited += 1
        if sum(mesh(j * j * a, j * j * b) for j in range(1, 4)) >= needed:
            continue
        assert depth < 36 and b - a > F(1, 10**5), (a, b, payment)
        m = (a + b) / 2
        pending.extend([(a, m, depth + 1), (m, b, depth + 1)])
    return visited, analytic


def bernstein_ray(q, span, payment):
    """Check one reciprocal pair against an upper exponential on each bin."""
    left, right = span
    degree = q.degree()
    checks = 0
    while left < right:
        end = min(right, left + F(1, 8))
        shift = q(fmpq_poly([rational(left), rational(end - left)]))
        shift -= rational(payment / 2) * (ball(end) / 2).exp().upper().fmpq()
        coefficients = [shift[j] / comb(degree, j) for j in range(degree + 1)]
        assert all(
            sum(comb(i, j) * coefficients[j] for j in range(i + 1)) >= 0
            for i in range(degree + 1)
        ), (left, end, payment)
        left = end
        checks += 1
    return (checks,)


def output(n, old, floor, denominator, old_payment, new_payment, budget):
    assert denominator > 0 and new_payment > old_payment >= 0
    result = max(floor, old * (denominator + old_payment) / (denominator + new_payment))
    upper = density_ceiling(result)
    assert upper < old
    row = {
        "n": n,
        "old": str(old),
        "upper": str(upper),
        "floor": str(floor),
        "payment": str(new_payment),
        "old_payment": str(old_payment),
        "checked_boxes": list(budget),
    }
    print(json.dumps(row), flush=True)
    return row


def dim14(benchmark, payment):
    for name, expected in DIM14_SHA.items():
        assert digest(benchmark / "sections" / name) == expected, name
    mark = benchmark / "sections" / "n14"
    proof = read(mark / "tightmatching14-final-packing-proof.json")
    assert F(proof["upper"]) == V14_UPPER
    subprocess.run(
        [sys.executable, "-B", str(benchmark / "sections" / "replay14.py")], check=True
    )
    source = read(mark / "geometry" / "sharp-ref-fourier.json")
    tau = F(source["tau0"])
    # This parent stores unnormalized generalized Laguerre coefficients.
    coefficients = list(map(rational, source["c"]))
    terms = []
    factorial = 1
    for k in range(len(coefficients)):
        if k:
            factorial *= k
        terms.append(
            sum(
                (
                    (-1) ** (j + k) * comb(j + 6, j - k) * coefficients[j]
                    for j in range(k, len(coefficients))
                ),
                fmpq(0),
            )
            / factorial
        )
    q = fmpq_poly(terms)
    floor = F("0.033818453867395")
    assert volume(14) / arb(314).sqrt() < ball(floor)
    domain = line_domain(14, tau, V14_UPPER, floor, V13_UPPER)
    checked = bernstein_ray(q, domain, payment)
    return output(14, V14_UPPER, floor, F(str(q(0))), F(0), payment, checked)


def dim15(path, payment):
    base = path / "native"
    sys.path.insert(0, str(base))
    from independent_fourier import transform
    from review_native_final import replay

    name = "minfix-lift15a-final-packing-proof.json"
    record = read(base / name)
    old = F(replay(name)["upper"])
    assert old == F(record["upper"])
    source = read(base / "lift15a.json")
    _, q = transform(15, source["d"])
    floor = max(F(record["profile_threshold"]), F(989, 1000) * old)
    span = line_domain(15, F(source["tau"]), old, floor, V14_UPPER)
    return output(
        15, old, floor, F(str(q(0))), F(0), payment, cover_ray(q, span, payment)
    )


def dim36(path):
    subprocess.run(
        [sys.executable, "-B", str(path / "verify.py"), "36"], cwd=path, check=True
    )
    sys.path.insert(0, str(path / "p"))
    from independent_fourier import transform
    from independent_compact import compact_data
    from portable_local import matrices

    old = F(read(path / "BOUNDS.json")["36"]["upper"])
    source = read(path / "p" / "pd36b.json")
    row = read(path / "p" / "BOUNDS.json")["36"]
    assert row["mode"] == "reciprocal-plane" and row["certificate"] == "pd36b.json"
    assert F(row["upper"]) == old and len(source["kernels"]) == 1
    (atom,) = source["kernels"]
    file = path / "p" / atom["local_kernel"]
    data, R, powers, qscale, _, B = matrices(file, atom["source_sha"], 36)
    assert (
        powers == sorted(powers)
        and len(powers) == 3
        and 1 < R / rational(atom["dilation"]) <= 2
    )
    assert not F(data.get("angular_debit", 0))
    B = B * rational(atom["alpha"]) / qscale**2
    z = fmpq_mat(3, 1, [rational(c) / qscale for c in atom["z"]])
    w = fmpq_mat(3, 1, [0, -1, 1])
    u = B.solve(-z)
    gap = (w.transpose() * u)[0, 0]
    assert gap > 0
    drop = gap**2 / (w.transpose() * B.solve(w))[0, 0]
    assert drop > 0
    _, q = transform(36, source["d"])
    reciprocal, *_ = compact_data(dict(n=36, kernels=source["reciprocal_kernels"]))
    assert 0 <= reciprocal < q(0)
    factor = (ball(source["tau"]) / 8) ** 18 / ball(19).gamma()
    debit = F(
        int((factor * arb(drop / q(0)) * DENSITY_UNIT).lower().floor().unique_fmpz()),
        DENSITY_UNIT,
    )
    floor = F(source["dual_target"])
    new = max(floor, old - debit)
    assert new < old
    row = {
        "n": 36,
        "old": str(old),
        "upper": str(new),
        "floor": str(floor),
        "density_debit": str(debit),
    }
    print(json.dumps(row), flush=True)
    return row


def check_manifest(path, filename):
    manifest = read(path / filename)
    for name, hash_ in manifest.items():
        f = (path / name).resolve()
        assert f.is_relative_to(path) and digest(f) == hash_, name


def run(n, argument, payment):
    ctx.prec = 512
    if n == 14:
        return dim14(argument, payment)
    archive = "native" if n == 15 else "high"
    if argument.is_dir():
        assert digest(argument / "MANIFEST.json") == MANIFEST_SHA[archive]
        check_manifest(argument, "MANIFEST.json")
        return dim15(argument, payment) if n == 15 else dim36(argument)
    assert digest(argument) == ARCHIVE_SHA[archive]
    with TemporaryDirectory() as temp:
        with ZipFile(argument) as proof:
            proof.extractall(temp)
        path = Path(temp)
        return dim15(path, payment) if n == 15 else dim36(path)


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Verification requires assertions.")
    certificates = read(HERE / "certificates.json")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dimension", type=int, choices=sorted(map(int, certificates)))
    parser.add_argument("source", type=Path)
    args = parser.parse_args()
    expected = certificates[str(args.dimension)]
    payment = F(expected["payment"]) if "payment" in expected else None
    row = run(args.dimension, args.source.resolve(), payment)
    assert row == expected
