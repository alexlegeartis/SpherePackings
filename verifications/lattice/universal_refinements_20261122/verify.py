"""Exact reciprocal and radial-moment refinements of the pinned lattice bounds."""

import argparse
from fractions import Fraction as F
from functools import lru_cache
from hashlib import file_digest
from itertools import product
import json
from math import comb
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

from flint import arb, arb_poly, ctx, fmpq, fmpq_mat, fmpq_poly


HERE = Path(__file__).resolve().parent
ARCHIVE_SHA = {
    "native": "4cae881e275464d6555a9beff843c887141830da8978fd4ee8b75f0b4d968b8f",
    "high": "bd9490675ecc38e96c4dddbc484ee7a99bae0a231902d5f3ac974eae9d2143b6",
}
NATIVE_PROOFS = {15: "minfix-lift15a", 19: "nativeflag19v2", 23: "flag23a"}
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


def stationary_lower_bounds(q, units, end=128):
    """Enclose the minimum of the Gaussian polynomial on a closed interval."""
    derivative = 2 * q.derivative() - q
    roots = derivative.complex_roots()
    assert sum(m for _, m in roots) == derivative.degree()
    pol = arb_poly([arb(r) for r in q.coeffs()])

    def value(point):
        return pol(point) * (-point / 2).exp()

    real = [
        root.real
        for root, _ in roots
        if root.imag.contains(0) and not (root.real < 0 or root.real > end)
    ]
    interior = [(point, value(point)) for point in real]

    @lru_cache(maxsize=16384)
    def endpoint(point):
        return arb(q(rational(point))) * (-ball(point) / 2).exp()

    def lower(a, b):
        assert 0 <= a <= b
        if b > end:
            return 0
        left, right = ball(a), ball(b)
        candidates = [endpoint(a), endpoint(b)]
        candidates.extend(
            val for point, val in interior if not (point < left or point > right)
        )
        return max(
            0, min(int((v * units).lower().floor().unique_fmpz()) for v in candidates)
        )

    return lower


def cover_plane(n, q, g, payment):
    """Cover the parent reduced-plane domain using exact signed bounds."""
    # It suffices to take one representative of each opposite pair.
    from independent_fourier import nonnegative_ray

    mesh, units, analytic = signed_lower_mesh(q)
    assert nonnegative_ray(q, F(g["gap"]))
    local = stationary_lower_bounds(q, units)
    needed = payment * units / 2
    G, L = g["R"], g["L"]
    index = sorted(
        (
            (i, j)
            for i, j in product(range(-4, 5), range(5))
            if j > 0 or (j == 0 and i > 0)
        ),
        key=lambda pair: sum(t * t for t in pair),
    )
    assert len(index) == len(set(index))
    granularity = 16
    G *= granularity
    L *= granularity
    area = g["area"] * granularity**2
    projection = int(g["projection_power_upper"]) * granularity**n
    pending = [tuple(granularity * v for v in g["box"])]
    checked = discarded = 0
    while pending:
        a, A, b, B, l, h = pending.pop()
        # These are the parent area, projection and reduction inequalities.
        if (
            A * B < area
            or a * b ** (n - 1) > projection
            or a * (L * L - h * h) > B * L * L
        ):
            discarded += 1
            continue
        amount = 0
        lengths = []
        for i, j in index:
            lo, hi = i * L + j * l, i * L + j * h
            s = max(0, lo, -hi)
            t = max(abs(lo), abs(hi))
            span = (
                F(a * s * s + b * j * j * L * L, G * L * L),
                F(A * t * t + B * j * j * L * L, G * L * L),
            )
            amount += mesh(*span)
            lengths.append(span)
            if amount >= needed:
                break
        if amount >= needed:
            checked += 1
            continue
        # A small box can be resolved using the stationary-point enclosures.
        if max(F(A - a, G), F(B - b, G), F((h - l) * max(A, G), L * G)) <= F(1, 40):
            refined = 0
            for span in lengths:
                refined += local(*span)
                if refined >= needed:
                    break
            if refined >= needed:
                checked += 1
                continue
        axes = [i for i, v in enumerate((A - a, B - b, h - l)) if v >= 2]
        assert axes, (a, A, b, B, l, h, payment)
        axis = max(
            axes,
            key=lambda k: (F(A - a, G), F(B - b, G), F((h - l) * max(A, G), L * G))[k],
        )
        node = [a, A, b, B, l, h]
        left, right = 2 * axis, 2 * axis + 1
        middle = (node[left] + node[right]) // 2
        assert node[left] < middle < node[right], (node, payment, F(2 * amount, units))
        one, two = node.copy(), node.copy()
        one[right] = two[left] = middle
        pending.extend([tuple(one), tuple(two)])
    assert checked > 0
    return checked, analytic, discarded


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


def native(n, path, payment):
    base = path / "native"
    sys.path.insert(0, str(base))
    from independent_fourier import transform
    from review_native_final import replay

    name = NATIVE_PROOFS[n] + "-final-packing-proof.json"
    record = read(base / name)
    old = F(replay(name)["upper"])
    assert old == F(record["upper"]) and record["n"] == n
    if n == 15:
        source = read(base / "lift15a.json")
        _, q = transform(n, source["d"])
        floor = max(F(record["profile_threshold"]), F(989, 1000) * old)
        span = line_domain(n, F(source["tau"]), old, floor, V14_UPPER)
        return output(
            n, old, floor, F(str(q(0))), F(0), payment, cover_ray(q, span, payment)
        )

    if n == 19:
        source_path = base / "nativeflag19v2.json"
        source = read(source_path)
        plane = read(base / "nativeflag19v2-reciprocal-flag.json")
        geometry = source["dual_flag"]
        floor = max(F(record["profile_threshold"]), F(geometry["target"]))
    else:
        plane = read(base / "flag23a-dual-plane.json")
        geometry = plane["geometry"]
        source_path = base / plane["prior_fourier"]
        source = read(source_path)
        floor = F(geometry["target"])
    assert digest(source_path) == plane["prior_sha"]
    assert geometry["n"] == source["n"] == n
    previous = F(2 * plane["lower"], plane["unit"])
    _, q = transform(n, source["d"])
    checked = cover_plane(n, q, geometry, payment)
    return output(n, old, floor, F(str(q(0))), previous, payment, checked)


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
    archive = "native" if n in NATIVE_PROOFS else "high"
    if argument.is_dir():
        assert digest(argument / "MANIFEST.json") == MANIFEST_SHA[archive]
        check_manifest(argument, "MANIFEST.json")
        return native(n, argument, payment) if archive == "native" else dim36(argument)
    assert digest(argument) == ARCHIVE_SHA[archive]
    with TemporaryDirectory() as temp:
        with ZipFile(argument) as proof:
            proof.extractall(temp)
        path = Path(temp)
        return native(n, path, payment) if archive == "native" else dim36(path)


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
