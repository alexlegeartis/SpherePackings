"""Replay the parent bounds and verify exact compact-plane refit certificates."""

import argparse
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

from flint import arb, ctx, fmpq, fmpq_mat, fmpq_poly

HERE = Path(__file__).resolve().parent
ARCHIVE_SHA = "bd9490675ecc38e96c4dddbc484ee7a99bae0a231902d5f3ac974eae9d2143b6"
MANIFEST_SHA = "c757fa3830618af1d0ea1c19762f0bd54c0bbfe4c4aae866006989bcccea191d"


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


@lru_cache(maxsize=None)
def volume(n):
    return arb.pi() ** arb(fmpq(n, 2)) / arb(fmpq(n + 2, 2)).gamma()


def verify(row, parent):
    from free_plane_geometry import signs
    from independent_compact import differential_check
    from independent_fourier import transform
    from polynomial_tools import sign_half

    n = row["n"]
    assert type(n) is int and n in range(36, 49)
    tau = fmpq(row["tau"])
    radius = fmpq(row["radius"])
    gap = fmpq(row["gap"])
    target = fmpq(row["target"])
    benchmark = fmpq(parent[str(n)]["upper"])
    assert tau > 0 and 1 < radius and radius**2 < 3
    assert 0 <= target < benchmark and gap >= 0
    assert (target == 0) == (gap == 0)
    if gap:
        previous = fmpq(parent[str(n - 1)]["upper"])
        reciprocal_floor = 4 * arb.pi() ** 2 / arb(tau) * (
            2 * volume(n - 1) * arb(target / previous) / volume(n)
        ) ** 2
        assert arb(gap) < reciprocal_floor

    powers = row["powers"]
    m = len(powers)
    assert 2 <= m <= 5 and all(type(exponent) is int for exponent in powers)
    step = powers[1] - powers[0]
    assert powers[0] > 0 and step > 0
    assert powers == [powers[0] + k * step for k in range(m)]

    blocks = [fmpq_mat([list(map(fmpq, data)) for data in block]) for block in row["A"]]
    assert blocks
    for block in blocks:
        assert block.nrows() == block.ncols() == m
        assert block == block.transpose()
        for size in range(1, m + 1):
            for indices in combinations(range(m), size):
                assert fmpq_mat([[block[i, j] for j in indices] for i in indices]).det() >= 0
    assert blocks[0].det() > 0
    assert len(row["z"]) == m
    z = fmpq_mat(m, 1, list(map(fmpq, row["z"])))
    origin = (z.transpose() * blocks[0].solve(z))[0, 0]

    p, q = transform(n, row["d"])
    assert q(0) > 0
    for label, polynomial, start in [("reciprocal", q, gap), ("spatial", -p, tau * radius)]:
        result = sign_half(polynomial, start)
        assert result["pass_sign"], (n, label, result)

    h = fmpq_poly([radius / (radius - 1), -1 / (radius - 1)])
    features = [h**power for power in powers]
    total = sum(blocks, fmpq_mat(m, m))
    radial = fmpq_poly()
    for i in range(m):
        radial += 2 * z[i, 0] * features[i]
        for j in range(m):
            radial += 2 * total[i, j] * features[i] * features[j]
    radial_nodes, radial_depth = differential_check(p, tau, [(radius, radial)])

    block_count = len(blocks)
    kernel = {
        "n": n,
        "radius": str(radius),
        "powers": powers,
        "co": ["1"] * block_count,
        "q": str(block_count),
        "B": [
            [[str(block[i, j] * block_count**2) for j in range(m)] for i in range(m)]
            for block in blocks
        ],
    }
    plane_nodes = []
    for kind in ("cap", "triangle"):
        result = signs(kernel, kind, limit=8000)
        assert result["status"] == "PASS", (n, result)
        plane_nodes.append(result["nodes"])

    numerator = p(0) + origin
    assert numerator > 0
    bound = (arb(tau) / 8) ** arb(fmpq(n, 2)) / arb(fmpq(n + 2, 2)).gamma()
    bound *= arb(numerator / q(0))
    upper = fmpq(str(Fraction(row["upper"])))
    assert 0 < upper < benchmark and target <= upper
    assert bound <= arb(upper)
    reduction = 100 * (1 - float(Fraction(str(upper / benchmark))))
    print(f"Verified rank {n} upper {row['upper']} ({reduction:.8f}% reduction), "
          f"plane nodes {plane_nodes}, radial nodes {radial_nodes}, radial depth {radial_depth}",
          flush=True)


def replay(archive, dimensions):
    assert digest(archive / "MANIFEST.json") == MANIFEST_SHA
    records = read(HERE / "certificates.json")
    assert set(dimensions) <= set(map(int, records))
    needed = set(dimensions) | {n - 1 for n in dimensions if fmpq(records[str(n)]["gap"])}
    if 42 in needed:
        needed.add(40)
    if 48 in needed:
        needed.update({46, 47})
    subprocess.run(
        [sys.executable, "-B", "-u", str(archive / "verify.py"), *map(str, sorted(needed))],
        cwd=archive,
        check=True,
    )
    sys.path.insert(0, str(archive / "p"))
    ctx.prec = 512
    parent = read(archive / "BOUNDS.json")
    for n in dimensions:
        row = records[str(n)]
        assert row["n"] == n
        verify(row, parent)


def main():
    if not __debug__:
        raise RuntimeError("Verification requires assertions.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("dimensions", nargs="*", type=int)
    args = parser.parse_args()
    dimensions = args.dimensions or sorted(map(int, read(HERE / "certificates.json")))
    archive = args.archive.resolve()
    if archive.is_dir():
        replay(archive, dimensions)
        return
    assert digest(archive) == ARCHIVE_SHA
    with TemporaryDirectory() as temp:
        with ZipFile(archive) as proof:
            proof.extractall(temp)
        replay(Path(temp), dimensions)


if __name__ == "__main__":
    main()
