"""Replay the lattice density bounds from the pinned benchmark archive."""
import argparse
from bisect import bisect_left, bisect_right
from fractions import Fraction
from hashlib import sha256
from math import floor
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

HERE = Path(__file__).resolve().parent
BENCHMARK_COMMIT = "d7ffd0df47dfeb11814631858205f05ed73ec8ff"
BENCHMARK_ARCHIVE_SHA256 = "bd9490675ecc38e96c4dddbc484ee7a99bae0a231902d5f3ac974eae9d2143b6"
BENCHMARK_MANIFEST_SHA256 = "c757fa3830618af1d0ea1c19762f0bd54c0bbfe4c4aae866006989bcccea191d"
REPRICE_DIMENSIONS = tuple(range(25, 36))
ORIGIN_DIMENSIONS = tuple(map(int, json.loads((HERE / "origin-certificates.json").read_bytes())))
DIMENSIONS = tuple(sorted(REPRICE_DIMENSIONS + ORIGIN_DIMENSIONS))
DISPLAY_SCALE = 10**22


def read_json(path):
    return json.loads(path.read_bytes())


def benchmark_root(archive, dimension):
    if dimension == 31:
        return archive / "p2"
    if dimension == 32:
        return archive / "p3"
    if dimension <= 30:
        return archive / "p"
    return archive / "p" / "baseline"


def lift_weights(edges, weights, target_edges):
    """A target bin may cross a source edge. Use the smallest valid coefficient."""
    assert len(edges) == len(weights) + 1
    assert edges[0] <= target_edges[0] < target_edges[-1] <= edges[-1]
    return [
        min(weights[bisect_right(edges, lo) - 1:bisect_left(edges, hi)])
        for lo, hi in zip(target_edges, target_edges[1:])
    ]


def upstream_angular_cut(dimension, cut, edges, tangent):
    from inherited_angular import verify

    verify(dimension, edges, cut)
    end = Fraction(cut["source_interval"][1])
    diagonal = Fraction(cut["f1"])
    scale = Fraction(cut["constant"])
    weights = []
    for upper, cap in zip(edges[1:], cut["caps"]):
        count = int(upper <= end)
        payment = Fraction(cap["M"]) * cap["C"]
        weights.append(((2*tangent + 1 - diagonal)*count - payment) / scale)
    return weights, Fraction(tangent*(tangent + 1), scale)


def angular_bank_cut(dimension, cut, tangent):
    from audit_angular import check

    original_weights, original_budget = check(dimension, cut)
    edges = list(map(Fraction, cut["grid"]))
    old_tangent = int(cut["tangent"])
    if not old_tangent:
        assert tangent == 0
        return edges, original_weights, original_budget * ((1 << dimension) - 1)
    start, stop = map(Fraction, (cut["A"], cut["B"]))
    epsilon = Fraction(cut["eps"])
    leading = 1 - epsilon
    diagonal = sum(map(Fraction, cut["coeffs"])) - epsilon
    scale = 2*old_tangent + 1
    weights = []
    for lo, hi, old in zip(edges, edges[1:], original_weights):
        count = int(start <= lo < hi <= stop)
        payment = ((scale * leading - diagonal)*count - scale*old)
        assert payment >= 0
        weights.append((((2*tangent + 1)*leading - diagonal)*count - payment)/scale)
    return edges, weights, leading*Fraction(tangent*(tangent + 1), scale)


def verify_dimension(archive, dimension, record):
    root = benchmark_root(archive, dimension)
    sys.path[:0] = [str(root), str(archive)]
    from audit_geometry import check_half, check_prefix
    from audit_price import check_density
    from external_profiles import cert as inherited_source
    from independent_master import integer_census, pooled, verify as verify_parent
    from flint import arb, fmpq

    registry = read_json(root / "BOUNDS.json")
    proof_path = root / registry[str(dimension)]["certificate"]
    parent = verify_parent(proof_path, details=True, portable=True)
    assert parent["report"]["status"] == "PASS"
    assert parent["report"]["n"] == dimension
    assert parent["report"]["unconditional_fourier_sign"]
    assert not parent["c"].get("dual_gap")

    proof = read_json(proof_path)
    source = parent["c"]
    heights = read_json(root / proof["height_file"])
    unit = heights["unit"]
    edges = parent["E"]
    short = edges[:edges.index(2) + 1]
    classes = (1 << dimension) - 1
    original_tag = proof_path.stem.split("-rigorous", 1)[0]
    bank_path = root / (original_tag.rsplit("-r", 1)[0] + "-windows.json")
    bank = read_json(bank_path) if bank_path.exists() else []
    inherited, _ = inherited_source(dimension)
    inherited_edges = list(map(Fraction, inherited["E"]))
    prefix = read_json(root / f"prefix-{dimension}-degree7.json")

    # Preserve the independently checked compact origin and the fixed nonangular rows.
    origin = parent["bonus"]
    costs = [Fraction(0)] * len(heights["heights"])
    for old_cut, multiplier in source["prices"]:
        price = Fraction(multiplier)
        assert price >= 0
        kind = old_cut["kind"]
        if kind == "upstream-angular":
            from inherited_angular import verify
            _, previous_budget = verify(dimension, inherited_edges, inherited["angular_cuts"][old_cut["index"]])
            origin -= price * previous_budget
        elif kind == "angular-bank":
            from audit_angular import check
            _, previous_budget = check(dimension, bank[old_cut["index"]])
            origin -= price * previous_budget
        else:
            end = Fraction(old_cut["end"])
            if kind == "prefix":
                count = check_prefix(dimension, prefix[str(end)])
            else:
                assert kind == "half-sum"
                table = read_json(root / f"half-sum-counts-{dimension}-16.json")
                count = check_half(next(row for row in table if Fraction(row["t"]) == end))
            assert count == old_cut["cap"]
            costs = [value + price * (upper <= end) for value, upper in zip(costs, edges[1:])]
    assert origin >= 0

    for witness in record["rows"]:
        price = Fraction(witness["price"])
        index = witness["i"]
        tangent = witness["k"]
        assert price > 0 and type(index) is int and index >= 0
        assert type(tangent) is int and tangent >= 0
        kind = witness["kind"]
        if kind == "upstream-angular":
            weights, allowance = upstream_angular_cut(
                dimension, inherited["angular_cuts"][index], inherited_edges, tangent
            )
            grid = inherited_edges[:inherited_edges.index(2) + 1]
        else:
            assert kind == "angular-bank"
            grid, weights, allowance = angular_bank_cut(dimension, bank[index], tangent)
        for i, value in enumerate(lift_weights(grid, weights, short)):
            costs[i] += price * value
        origin += price * allowance / classes

    if dimension == 31:
        # The reallocation keeps the radial diagonal exactly fixed.
        from verify_origin_saving import verify
        benchmark = read_json(archive / "BOUNDS.json")["31"]
        evidence = read_json(archive / benchmark["evidence"])
        assert verify(root, evidence, archive / "kernels" / evidence["kernel"])["status"] == "PASS"
        origin -= Fraction(evidence["saving"])

    scores = [height - floor(unit * price) for height, price in zip(heights["heights"], costs)]
    census = integer_census(dimension, edges, scores, unit, proof)
    numerator = Fraction(str(parent["p"](0))) + origin + pooled(dimension, edges, census, parent["prefix"])
    denominator = Fraction(str(parent["q"](0)))
    assert numerator > 0 and denominator > 0
    ratio = numerator / denominator
    tau = parent["tau"]

    value = (arb(fmpq(str(tau))) / 8)**arb(fmpq(dimension, 2))
    value *= arb(fmpq(str(ratio))) / arb(fmpq(dimension + 2, 2)).gamma()
    digits = int((value * DISPLAY_SCALE).upper().ceil().unique_fmpz())
    upper = Fraction(digits, DISPLAY_SCALE)
    check_density(dimension, tau, ratio, upper)
    previous = Fraction(read_json(archive / "BOUNDS.json")[str(dimension)]["upper"])
    assert upper < previous
    assert Fraction(record["upper"]) == upper
    assert Fraction(record["benchmark_upper"]) == previous
    print(json.dumps({
        "dimension": dimension,
        "upper": f"{digits // DISPLAY_SCALE}.{digits % DISPLAY_SCALE:022d}",
        "benchmark_upper": str(previous),
        "reduction_percent": float(100 * (previous-upper) / previous),
    }), flush=True)


def verify_origin_case(archive, dimension, record):
    root = archive / "p"
    subprocess.run(
        [sys.executable, "-B", str(root / "verify.py"), str(dimension)],
        cwd=root,
        check=True,
        stdout=sys.stderr,
    )
    sys.path[:0] = [str(root), str(archive)]
    from verify_origin_saving import verify

    assert record["n"] == dimension
    kernel = (HERE / record["kernel"]).resolve()
    assert kernel.is_relative_to(HERE / "kernels")
    result = verify(root, record, kernel)
    assert result["status"] == "PASS"
    assert Fraction(result["upper"]) == Fraction(record["upper"])
    upper = Fraction(record["display"])
    previous = Fraction(read_json(archive / "BOUNDS.json")[str(dimension)]["upper"])
    assert Fraction(record["benchmark_upper"]) == previous
    assert Fraction(result["upper"]) <= upper < previous
    print(json.dumps({
        "dimension": dimension,
        "upper": record["display"],
        "benchmark_upper": str(previous),
        "reduction_percent": float(100 * (previous-upper) / previous),
    }), flush=True)


def check_manifest(archive):
    manifest_file = archive / "MANIFEST.json"
    assert sha256(manifest_file.read_bytes()).hexdigest() == BENCHMARK_MANIFEST_SHA256
    manifest = read_json(manifest_file)
    for name, expected in manifest.items():
        path = (archive / name).resolve()
        assert path.is_relative_to(archive)
        assert sha256(path.read_bytes()).hexdigest() == expected, name


def run(archive, dimensions):
    check_manifest(archive)
    for dimension in dimensions:
        subprocess.run(
            [sys.executable, "-B", str(HERE / "verify.py"), str(archive), "--worker", str(dimension)],
            check=True,
        )


def main():
    if not __debug__:
        raise RuntimeError("Run the checker with Python assertions enabled.")
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("archive", type=Path, help="Original n25-48-proofs-r9.zip or its extracted directory")
    cli.add_argument("dimensions", type=int, nargs="*", choices=DIMENSIONS)
    cli.add_argument("--worker", type=int, choices=DIMENSIONS, help=argparse.SUPPRESS)
    args = cli.parse_args()
    archive = args.archive.resolve(strict=True)
    if args.worker:
        if args.worker in ORIGIN_DIMENSIONS:
            record = read_json(HERE / "origin-certificates.json")[str(args.worker)]
            verify_origin_case(archive, args.worker, record)
        else:
            record = read_json(HERE / "certificates.json")[str(args.worker)]
            verify_dimension(archive, args.worker, record)
        return
    dimensions = args.dimensions or DIMENSIONS
    if archive.is_dir():
        run(archive, dimensions)
    else:
        assert sha256(archive.read_bytes()).hexdigest() == BENCHMARK_ARCHIVE_SHA256
        with tempfile.TemporaryDirectory(prefix="lattice-certificate-") as folder:
            with zipfile.ZipFile(archive) as package:
                package.extractall(folder)
            run(Path(folder).resolve(), dimensions)


if __name__ == "__main__":
    main()
