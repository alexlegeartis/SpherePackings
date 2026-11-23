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

    def add(n, squared_density, label, evidence, digits=None):
        candidates[n].append((Fraction(squared_density), label, source_link(evidence), digits))

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

    package = "research/vardy28_stack29_20261122"
    for row in data(f"{package}/results.json")["rows"]:
        assert row["dimension"] == 29 and row["kind"] == "explicit nonlattice periodic"
        add(row["dimension"], row["center_density_squared"],
            "Alternating Vardy-type layers", f"{package}/README.md")

    package = "research/vardy30_shear31_20261123"
    row = data(f"{package}/verification.log")
    assert row["dimension"] == 31 and row["status"] == "PASS"
    add(31, row["squared_center_density"], "Order-four Vardy-type shear", f"{package}/README.md")

    package = "research/e7_hermitian35_20261123"
    proof = data(f"{package}/results.json")
    assert proof["status"] == "PASS"
    for row in proof["rows"]:
        n = row["dimension"]
        assert n in (34, 35, 36, 37)
        new = Fraction(row["center_density_squared"])
        assert new == Fraction(row["minimum_squared"], 4)**n / int(row["gram_determinant"])
        assert new > Fraction(row["preceding_center_density_squared"]) == candidates[n][0][0]
        add(n, new, row["method"], f"{package}/README.md")

    package = "research/e7_37_section_20261123"
    sec37 = data(f"{package}/results.json")
    assert sec37["status"] == "PASS" and sec37["kind"] == "explicit lattice"
    assert sec37["dimension"] == 37
    new = Fraction(sec37["center_density_squared"])
    assert new == Fraction(sec37["minimum_squared"], 4)**37 / int(sec37["gram_determinant"])
    assert new > Fraction(sec37["preceding_center_density_squared"]) == max(x[0] for x in candidates[37])
    add(37, new, sec37["method"], f"{package}/README.md")

    package = "research/e7_periodic38_20261123"
    row = data(f"{package}/results.json")
    assert row["dimension"] == 38 and row["status"] == "PASS"
    assert row["kind"] == "explicit nonlattice periodic"
    new = Fraction(row["center_density_squared"])
    assert new == Fraction(row["minimum_squared"], 4)**38 * row["motif_size"]**2 / int(row["period_gram_determinant"])
    assert new > Fraction(row["preceding_center_density_squared"]) == candidates[38][0][0]
    package = "research/e7_lattice38_20261123"
    lattice = data(f"{package}/results.json")
    assert lattice["status"] == "PASS" and lattice["kind"] == "explicit lattice"
    assert lattice["dimension"] == 38
    assert new == Fraction(lattice["center_density_squared"])
    assert new == Fraction(lattice["minimum_squared"], 4)**38 / int(lattice["gram_determinant"])
    add(38, new, lattice["method"], f"{package}/README.md")

    history = data("research/range48_96_20261115/results.json")
    for row in history["rows"]:
        best = row["best_checked"]
        add(row["dimension"], best["density_squared"], best["name"], best["evidence"])

    package = "research/ne64_dual_matching59_20261123"
    row = data(f"{package}/data/results.json")
    assert row["dimension"] == 59 and row["status"] == "PASS" and row["periodic_not_lattice"]
    add(59, row["center_density_squared_floor"], "Ne64 antipode with matched normal", f"{package}/README.md")

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

    add(81, Fraction(7, 2)**81 / (4*41**12), "MW80 half-coset lift",
        "research/mw80_half_coset_81_20261123/README.md")

    package = "research/mw92_small_simplexes_20261122"
    for row in data(f"{package}/verified_results.json"):
        assert row["dimension"] in (86, 87)
        add(row["dimension"], row["center_density_squared"],
            "MW92 projected-simplex packing", f"{package}/README.md", row["occupied_density_floor"])

    package = "research/l80_a2_section_20261122"
    row = data(f"{package}/results.json")
    assert row["dimension"] == 82
    add(row["dimension"], row["center_density_squared"],
        "A2 section of the 83D Bachoc–Nebe extension", f"{package}/README.md")

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
        ("uneven_layers55_20261122", "Three unevenly spaced Leterrier54 layers"),
        ("gamma72_leech_layers_20261122", "Explicit Γ72 lattice with Leech-block verification"),
        ("gamma72_depth6_20261122", "Explicit Γ72 half-coset lattice"),
        ("gamma72_tilted_d4_20261122", "Tilted Γ72–D4 lattice"),
        ("gamma72_root_tilt_20261122", "Half-root Γ72–D4 lattice"),
        ("mw92_extensions_20261122", "MW92 extension"),
        ("mw92_order3_dichotomy_20261122", "MW92 order-three lattice dichotomy"),
    ]:
        path = f"research/{package}"
        for row in data(f"{path}/results.json")["rows"]:
            add(row["dimension"], row["center_density_squared"],
                "MW92 nine-translate periodic packing" if row["kind"] == "periodic existence" else label,
                f"{path}/README.md")

    package = "mw92_cfr_taylor13_20261123"
    path = f"research/{package}"
    for row in data(f"{path}/result.json")["rows"]:
        assert row["dimension"] in (93, 94, 95, 96)
        add(row["dimension"], row["center_density_squared"], "MW92 degree-13 covering refinement",
            f"{path}/README.md", row["occupied_density_floor"])

    for row in data("research/lattice_lower_49_96_20261122/mw92_bare_results.json")["rows"]:
        assert row["kind"] == "audited lattice section"
        add(row["dimension"], row["center_density_squared"], "Primitive MW92 bare lattice section",
            "research/lattice_lower_49_96_20261122/SOURCE_NOTES.md#8691d-mw92-bare-kernels")

    assert set(candidates) == set(range(1, 97))
    pi_lower, _ = pi_interval()
    lines = [
        "# Sphere-packing lower bounds in dimensions 1–96", "",
        "The table gives occupied-volume density floors. It uses Cohn’s literature table for",
        "1–48D and replaces its bounds where this repository has verified improvements.",
        "For 49–96D it uses the strongest checked constructions in the repository’s audit",
        "and later certificates. The source coverage reflects the 22 November 2026 review.",
        "The bounds in 1–3D, 8D and 24D are exact optima. The 55D entry is the proved",
        "lower bound documented by its linked construction. The printed Leterrier55",
        "challenge value is discussed in that proof.", "",
        "| Dimension | Occupied density at least | Construction or source |",
        "|---:|---:|---|",
    ]
    for n, choices in sorted(candidates.items()):
        center_squared, title, source, exact_floor = max(choices, key=lambda row: row[0])
        occupied_squared = center_squared * unit_ball_squared(n, pi_lower)
        if exact_floor is not None:
            assert Fraction(exact_floor)**2 <= occupied_squared
        density = (exact_floor or display_floor(occupied_squared)).lower()
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
