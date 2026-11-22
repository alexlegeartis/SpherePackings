"""Replay the native packing bounds and check the superlattice deductions."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

from flint import arb, ctx, fmpq, fmpq_poly


HERE = Path(__file__).resolve().parent
ARCHIVE_SHA = "4cae881e275464d6555a9beff843c887141830da8978fd4ee8b75f0b4d968b8f"
TAGS = {17: "oracle17a", 18: "oracle18a", 20: "mix20a", 21: "mix21a", 22: "mix22a"}


def read(path):
    return json.loads(path.read_bytes())


def rational(x):
    return fmpq(str(x))


def volume(n):
    return arb.pi()**arb(fmpq(n, 2)) / arb(fmpq(n + 2, 2)).gamma()


def check_spatial(polynomial, scales, edges, prices, count, short, long):
    first = [t * (polynomial.derivative() - polynomial / 2) for t in scales]
    second = [t * (p.derivative() - p / 2) for t, p in zip(scales, first)]
    for a, b, price in zip(edges, edges[1:], prices):
        ceiling = arb(rational(price + (short if b <= 2 else long))) / count
        pending = [(rational(a), rational(b), 0)]
        while pending:
            left, right, depth = pending.pop()
            middle = (left + right) / 2
            width = (right - left) / 2
            value, slope, curvature = arb(0), arb(0), arb(0)
            for t, d1, d2 in zip(scales, first, second):
                exponential = arb(-t * middle / 2).exp()
                value += arb(polynomial(t * middle)) * exponential
                slope += arb(d1(t * middle)) * exponential
                local = d2(fmpq_poly([t * middle, t * width]))
                curvature += arb(sum(map(abs, local.coeffs()), fmpq(0))) * arb(-t * left / 2).exp()
            if value + abs(slope) * width + curvature * (width**2 / 2) <= ceiling:
                continue
            assert value <= ceiling and depth < 45
            pending.extend([(left, middle, depth + 1), (middle, right, depth + 1)])


def verify_refit(native, old, saturation):
    from independent_fourier import nonnegative_ray, transform
    from native_lift_audit import exact_fields
    from review_native_final import replay

    n = 18
    prior_rank = Fraction(replay("minfix-oracle17a-final-packing-proof.json")["upper"])
    refit = read(HERE / "refit-18.json")
    source = read(native / "oracle18a.json")
    assert refit["n"] == n and Fraction(refit["threshold"]) == Fraction(97, 100) * old
    weights = list(map(Fraction, refit["prices"]))
    assert len(weights) == len(source["prices"]) and min(weights) > 0
    source["prices"] = [[kind, str(weight)] for (kind, _), weight in zip(source["prices"], weights)]
    edges = list(map(Fraction, source["E"]))
    prices, cost, threshold, _ = exact_fields(n, source, edges)
    assert not threshold and len(prices) == len(edges) - 1 and 2 in edges

    p, q = transform(n, refit["d"])
    tau = rational(refit["tau"])
    target = Fraction(refit["threshold"])
    gap = rational(refit["gap"])
    lower = 4 * arb.pi()**2 / arb(tau) * (2 * volume(n - 1) * arb(rational(target / prior_rank)) / volume(n))**2
    assert 0 < arb(gap) <= lower and q(0) > 0
    assert nonnegative_ray(q, gap) and nonnegative_ray(-p, tau * rational(edges[-1]))
    short, long = map(Fraction, (refit["short"], refit["long"]))
    assert min(short, long) >= 0
    classes = 2**n - 1
    scales = [tau * j*j for j in range(1, source["multiples"] + 1)]
    check_spatial(p, scales, edges, prices, 2 * classes, short, long)

    z = (1 - 2 / edges[-1])**2
    mean = Fraction(3, n * (n + 2)) - z / n
    assert 2 < edges[-1] < 4 and mean > 0 and z <= Fraction(6, n + 4)
    capacity = (1 - z) // mean
    charge = capacity * long + max(Fraction(0), short - long)
    ratio = (Fraction(str(p(0))) + cost + charge) / Fraction(str(q(0)))
    density = (arb(tau) / 8)**9 * arb(rational(ratio)) / arb(10).gamma()
    upper = Fraction(refit["upper"])
    assert density <= arb(rational(upper)) and target <= upper < saturation < old
    assert upper <= Fraction(refit["display"])
    print(f"Verified dimension 18 Fourier refit upper {refit['display']}", flush=True)


def verify(extracted, dimensions):
    native = extracted / "native"
    sys.path[:0] = [str(native), str(extracted / "support")]
    from independent_fourier import transform
    from native_lift_audit import exact_fields
    from review_native_final import replay

    rows = read(HERE / "certificates.json")
    assert sorted(map(int, rows)) == sorted(TAGS)
    ctx.prec = 512
    for n in dimensions:
        tag = TAGS[n]
        name = f"minfix-{tag}-final-packing-proof.json"
        prior = read(native / name)
        assert replay(name)["upper"] == prior["upper"]
        row = rows[str(n)]
        old = Fraction(prior["upper"])
        assert n == prior["n"] == row["n"] and old == Fraction(row["old"])
        source = read(native / f"{tag}.json")
        heights = read(native / f"{tag}-enclosures.json")
        edges = list(map(Fraction, heights["E"]))
        prices, _, threshold, _ = exact_fields(n, source, edges)
        assert threshold == Fraction(prior["profile_threshold"])
        residual = [Fraction(sum(v), heights["unit"]) - price
                    for v, price in zip(zip(*heights["rows"]), prices)]
        assert len(residual) == len(edges) - 1 and max(residual) <= 0
        assert Fraction(prior["ordinary_error"]) == 0 and Fraction(prior["minimum_class_error"]) <= 0

        count = prior["selected_exact_minimum_pairs"]
        classes = prior["nonzero_classes"]
        assert count == n * (n + 1) // 2 and classes == 2**n - 1
        saving = Fraction(0)
        for codimension, endpoint, mandatory, margin in row["schedule"]:
            assert 0 < codimension < n and count < mandatory == 2**(codimension + 1) - 1
            end, margin = Fraction(endpoint), Fraction(margin)
            assert end**n > 2**(2 * codimension) and end in edges
            j = edges.index(end)
            assert 0 < margin <= -max(residual[:j])
            saving += Fraction(mandatory - count, classes) * margin
            count = mandatory
        _, fourier = transform(n, source["d"])
        denominator = Fraction(str(fourier(0)))
        assert saving > 0 and denominator > 0
        ratio = Fraction(prior["rho"]) - saving / denominator
        assert 0 < ratio == Fraction(row["ratio"])
        factor = (arb(fmpq(source["tau"])) / 8)**arb(fmpq(n, 2)) / arb(fmpq(n + 2, 2)).gamma()
        upper = Fraction(row["upper"])
        assert threshold <= upper < old
        assert factor * arb(fmpq(str(ratio))) <= arb(fmpq(str(upper)))
        assert upper <= Fraction(row["display"]) < old
        print(f"Verified dimension {n} upper {row['display']} with {count} forced parity classes", flush=True)
        if n == 18:
            verify_refit(native, old, upper)


def main():
    if not __debug__:
        raise RuntimeError("Verification requires assertions.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("dimensions", type=int, nargs="*")
    args = parser.parse_args()
    dimensions = args.dimensions or sorted(TAGS)
    assert all(n in TAGS for n in dimensions)
    assert sha256(args.archive.read_bytes()).hexdigest() == ARCHIVE_SHA
    with TemporaryDirectory(prefix="lattice-saturation-") as tmp:
        root = Path(tmp)
        with ZipFile(args.archive) as archive:
            manifest = json.loads(archive.read("MANIFEST.json"))
            assert set(archive.namelist()) == set(manifest) | {"MANIFEST.json"}
            for name, digest in manifest.items():
                path = (root / name).resolve()
                assert path.is_relative_to(root) and path != root
                content = archive.read(name)
                assert sha256(content).hexdigest() == digest, name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(content)
        verify(root, dimensions)


if __name__ == "__main__":
    main()
