#!/usr/bin/env python3
"""Check vertical lattice extensions and count their shortest vectors exactly."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from functools import cache
from itertools import combinations
import json
from math import isqrt
from pathlib import Path
import subprocess
import tempfile

from replay import HERE, compile_screen, seeds


def scan(binary, seed, height_squared, offset=0):
    active = sum(1 << j for j in seed.columns)
    height_times_16 = int(16 * height_squared)
    header = f"{active} {seed.double_weight} {offset} {64-height_times_16} {len(seed.rows)} {64-4*height_times_16}"
    rows = [header] + [f"{low} {high}" for low, high in seed.rows]
    output = subprocess.run([binary, "--stack"], input="\n".join(rows), text=True,
                            capture_output=True, check=True).stdout.splitlines()
    checked, count, index_two, first_low, first_high, *_ = map(int, output[0].split())
    assert checked == 1 << (len(seed.columns) - len(seed.rows))
    assert index_two == 0
    if offset:
        candidate = (first_low, first_high) if count else None
        return 64 * checked, count, candidate
    candidates = [tuple(map(int, row.split())) for row in output[1:]]
    assert count == len(candidates) == len(set(candidates))
    return 64 * checked, candidates


def contact_counter(seed):
    residues = seed.residues()

    @cache
    def count_residue(low, high, target):
        counts = {(0, 0): 1}
        for j, column in seed.columns.items():
            scale = 4 * seed.weight(j)
            bound = isqrt(target // scale)
            parity = (low >> j) & 1
            values = [(column if (v >> 1) & 1 else 0, scale * v * v)
                      for v in range(-bound, bound + 1) if v & 1 == parity]
            next_counts = {}
            for (h, cost), count in counts.items():
                for shift, extra_cost in values:
                    total = cost + extra_cost
                    if total <= target:
                        key = h ^ shift, total
                        next_counts[key] = next_counts.get(key, 0) + count
            counts = next_counts
        return counts.get((high, target), 0)

    def count_coset(low, high, target):
        return sum(count_residue(low ^ bits, high ^ seed.syndrome(upper ^ (low & bits)), target)
                   for bits, upper in residues.items())

    return count_coset


def certify_distances(seed, height_squared, translate, offset=0):
    low, high = translate
    h = int(16 * height_squared)
    double = seed.syndrome(low)
    residues = seed.residues()

    def nearest(bits, syndrome, shifted=0):
        return min(seed.min_cost(bits ^ r, syndrome ^ seed.syndrome(s ^ (bits & r)), shifted)
                   for r, s in residues.items())

    first_nearest = nearest(low, high, offset)
    second_nearest = nearest(offset, double)
    assert first_nearest + h >= 64 and second_nearest + 4*h >= 64
    # For a quarter-grid phase, -3t differs from t by an integer shift on the
    # offset positions. For a half-grid phase, 3t and -t are the same coset.
    third_nearest = nearest(low, high ^ seed.syndrome(offset), offset) if offset else first_nearest
    assert third_nearest + 9*h >= 64 and 16*h >= 64
    covolume_squared = (2 ** (2 * (6 - len(seed.rows)))) * (2 if seed.double_weight else 1)
    return {"low_bits": low, "high_syndrome": high,
            "horizontal_squared_distances_times_16": [first_nearest, second_nearest],
            "center_density_squared": str(1 / (covolume_squared * height_squared))}


def certify_stack(seed, height_squared, translate, count_contacts):
    proof = certify_distances(seed, height_squared, translate)
    low, high = translate
    h = int(16 * height_squared)
    first = count_contacts(low, high, 64-h)
    second = count_contacts(0, seed.syndrome(low), 64-4*h)
    fourth = int(16*h == 64)
    proof["kissing_number"] = count_contacts(0, 0, 64) + 2 * (first + second + fourth)
    return proof


def quarter_grid(binary, n29, n30, jobs):
    cases = [(seed, height, positions)
             for seed, height in [(n29, Fraction(7, 16)), (n30, Fraction(5, 8))]
             for size in (1, 2) for positions in combinations(seed.columns, size)]

    def run_case(case):
        seed, height, positions = case
        offset = sum(1 << j for j in positions)
        checked, count, candidate = scan(binary, seed, height, offset)
        proof = certify_distances(seed, height, candidate, offset) if candidate else None
        return {"horizontal_dimension": len(seed.columns), "height_squared": str(height),
                "offset_positions": positions, "cosets_checked": checked,
                "admissible_phases": count, "first_certificate": proof}

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        rows = list(pool.map(run_case, cases))
    return {"cases_checked": len(rows), "cosets_checked": sum(row["cosets_checked"] for row in rows),
            "height_squared_by_horizontal_dimension": {"29": "7/16", "30": "5/8"},
            "cases_with_admissible_phases": [row for row in rows if row["admissible_phases"]]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quarter-grid", action="store_true")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    n29, n30, *_ = list(seeds())
    with tempfile.TemporaryDirectory(prefix="quaternary-stacks-") as path:
        binary = Path(path) / "screen"
        compile_screen(binary)
        if args.quarter_grid:
            report = quarter_grid(binary, n29, n30, args.jobs)
            name = "quarter_grid_stack_results.json"
            print(f"Checked {report['cases_checked']} quarter-grid stacking cases. "
                  f"{len(report['cases_with_admissible_phases'])} had admissible phases.")
        else:
            report = half_grid(binary, n29, n30)
            name = "stack_results.json"
            classes = sorted({certificate["kissing_number"] for trial in report["trials"]
                              for certificate in trial["certificates"] if trial["horizontal_dimension"] == 29})
            print(f"30D density ties with distinct kissing numbers {classes}.")
    contents = json.dumps(report, indent=2) + "\n"
    path = HERE / name
    if args.write:
        path.write_text(contents)
    else:
        assert path.read_text() == contents


def half_grid(binary, n29, n30):
    count29, count30 = contact_counter(n29), contact_counter(n30)
    # The vector supported on {2,4,8,14} projects to half a unit on the folded
    # axis and one unit at positions 4, 8 and 14 in the original 30D stack.
    original_high = (1 << 4) | (1 << 8) | (1 << 14)
    original = certify_stack(n29, Fraction(1, 2), (1 << 2, n29.syndrome(original_high)), count29)
    assert original["kissing_number"] == count30(0, 0, 64)
    trials = []
    for seed, count in [(n29, count29), (n30, count30)]:
        for height in (Fraction(1, 4), Fraction(1, 2)):
            checked, candidates = scan(binary, seed, height)
            trials.append({"horizontal_dimension": len(seed.columns), "height_squared": str(height),
                           "cosets_checked": checked, "certificates": [
                               certify_stack(seed, height, candidate, count) for candidate in candidates
                           ]})
    return {"original_30d_stack": original, "horizontal_29d_kissing_number": count29(0, 0, 64),
            "trials": trials}


if __name__ == "__main__":
    main()
