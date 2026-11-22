"""Rebuild the source-attributed 1–96 packing lower-bound table."""

import argparse
from collections import defaultdict
from fractions import Fraction
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / "mw92_extensions_20261122"))
from verify_mw92 import display_floor, pi_interval, unit_ball_squared


def data(relative_path):
    return json.loads((ROOT / relative_path).read_text())


def source_link(evidence):
    if not evidence.startswith("https://"):
        assert (ROOT / evidence.split("#")[0]).exists(), evidence
    return evidence


def build_table():
    candidates = defaultdict(list)

    def add(n, squared_density, label, evidence):
        candidates[n].append((Fraction(squared_density), label, source_link(evidence)))

    cohn = json.loads((HERE / "cohn_inputs.json").read_text())
    assert [row["dimension"] for row in cohn["rows"]] == list(range(1, 49))
    for row in cohn["rows"]:
        n = row["dimension"]
        center_squared = Fraction(1)
        for base, raw in row["center_density_factors"]:
            exponent = 2*Fraction(raw)
            assert exponent.denominator == 1
            center_squared *= Fraction(base)**int(exponent)
        label = "Cohn’s literature table"
        if n == 36:
            audit = data("constructions/dim36/result.json")
            label = "Contact-preserving metric deformation"
            source = "constructions/dim36/README.md"
            assert center_squared == Fraction(2**24, 3**3*7**6)**2
            assert audit["dimension"] == 36
            assert center_squared == Fraction(audit["center_density"])**2
        else:
            source = cohn["source"] + "#data"
        add(n, center_squared, label, source)

    history = data("research/range48_96_20261115/results.json")
    for row in history["rows"]:
        best = row["best_checked"]
        add(row["dimension"], best["density_squared"], best["name"], best["evidence"])

    package = "research/packing_binary57_58_20261117"
    for row in data(f"{package}/RESULTS.json")["rows"]:
        add(row["dimension"], row["center_density_squared"], row["method"], f"{package}/README.md")

    for n, package, item, label in [
        (78, "packing78_verified", "packing78_independent_20261117/result.json", "Three-fiber MW80 projection"),
        (81, "packing81_verified", "l80_case101_independent_20261117/extension_result.json", "Bachoc–Nebe deep half-coset lattice"),
        (82, "packing82_verified", "l80_case108203_independent_20261117/extension_result.json", "Bachoc–Nebe mixed-depth lattice"),
    ]:
        package = f"research/{package}_20261117"
        row = data(f"{package}/{item}")
        assert row["dimension"] == n
        add(n, row["center_density_squared"], label, f"{package}/README.md")

    package = "research/packing83_85_deep_verified_20261117"
    for n in (83, 84, 85):
        row = data(f"{package}/deep_extension/verified{n}.json")
        assert row["dimension"] == n
        add(n, row["center_density_squared"], "Bachoc–Nebe parity extension", f"{package}/README_DEEP.md")

    package = "research/packing_release8590_20261117"
    for row in data(f"{package}/RESULTS.json")["rows"]:
        file = row["file"].replace(".json", "_verified.json")
        audit = data(f"{package}/packing_resume_20261117/{file}")
        n = row["dimension"]
        assert audit["dimension"] == n
        add(n, audit["density_squared"], "Correlated MW92 antipode", f"{package}/README.md")

    package = "research/mw92_flags9495_verified_20261117"
    for filename in ("RESULT.json", "RESULT95.json"):
        row = data(f"{package}/{filename}")
        add(row["dimension"], row["center_density_squared"], "MW92 torus offset flags", f"{package}/README.md")

    for package, label in [
        ("p48n_hurwitz_20261122", "Explicit P48n quaternionic lattice"),
        ("gamma72_leech_layers_20261122", "Explicit Γ72 lattice with Leech-block verification"),
        ("mw92_extensions_20261122", "MW92 extension"),
    ]:
        path = f"research/{package}"
        for row in data(f"{path}/results.json")["rows"]:
            add(row["dimension"], row["center_density_squared"],
                "MW92 nine-translate periodic packing" if row["kind"] == "periodic existence" else label,
                f"{path}/README.md")

    assert set(candidates) == set(range(1, 97))
    pi_lower, _ = pi_interval()
    lines = [
        "# Sphere-packing lower bounds in dimensions 1–96", "",
        "The table gives occupied-volume density floors. It uses Cohn’s literature table for",
        "1–48D and the strongest checked constructions in the repository’s 48–96D audit",
        "and later certificates. The source coverage reflects the 22 November 2026 review.",
        "The bounds in 1–3D, 8D and 24D are exact optima. The 55D entry is the proved",
        "lower bound documented by its linked construction. The printed Leterrier55",
        "challenge value is discussed in that proof.", "",
        "| Dimension | Occupied density at least | Construction or source |",
        "|---:|---:|---|",
    ]
    for n, choices in sorted(candidates.items()):
        center_squared, title, source = max(choices, key=lambda row: row[0])
        density = display_floor(center_squared*unit_ball_squared(n, pi_lower)).lower()
        if n == 1:
            density = "1"
        lines.append(f"| {n} | {density} | [{title}]({source}) |")

    lines.extend([
        "", "Rebuild the table with",
        "`python research/lower_bounds_1_96_20261122/verify_table.py --write`.",
        "Run the same command without `--write` to check the committed table.",
        "The inputs retain [Cohn’s centre-density factors](research/lower_bounds_1_96_20261122/cohn_inputs.json).",
        "The [previous literature review](research/covering_and_deep_cosets_20261115/PRIORITY.md)",
        "and [higher-dimensional comparison audit](research/range48_96_20261115/) document",
        "the scope of the underlying searches.", "",
    ])
    return "\n".join(lines)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    path = ROOT / "LOWER_BOUNDS_1_96.md"
    table = build_table()
    if args.write:
        path.write_text(table)
    else:
        assert table == path.read_text()
    print("Verified density floors and sources for 96 dimensions.")
