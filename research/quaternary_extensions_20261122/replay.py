#!/usr/bin/env python3
"""Replay finite translate searches around the 29D and 30D quaternary lattices."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, replace
from itertools import combinations
import json
from pathlib import Path
import subprocess
import tempfile


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "quaternary_29_30_20261122" / "witness.json"


@dataclass(frozen=True)
class Seed:
    columns: dict[int, int]
    double_weight: int
    rows: tuple[tuple[int, int], ...]
    restored_position: int | None = None
    restored_lifts: int = 0

    def syndrome(self, bits):
        answer = 0
        for j, column in self.columns.items():
            if bits & (1 << j):
                answer ^= column
        return answer

    def residues(self):
        residues = {0: 0}
        for low, high in self.rows:
            assert self.syndrome(low) == 0
            extension = {low ^ bits: high ^ upper ^ (low & bits)
                         for bits, upper in residues.items()}
            assert not residues.keys() & extension.keys()
            residues.update(extension)
        return residues

    def weight(self, j):
        return 1 + ((self.double_weight >> j) & 1)

    def min_cost(self, low, high_syndrome, offset):
        """Integer shortest-path replay, in units of 1/16 of final squared norm."""
        distances = [0] + [10**6] * 63
        for j, column in self.columns.items():
            costs = (1, 9, 9, 1) if (offset >> j) & 1 else (0, 4, 16, 4)
            symbol = (low >> j) & 1
            a = costs[symbol] * self.weight(j)
            b = costs[symbol + 2] * self.weight(j)
            distances = [min(distances[h] + a, distances[h ^ column] + b)
                         for h in range(64)]
        return distances[high_syndrome]

    def verify_minimum(self):
        span = {0}
        assert len(set(self.columns.values())) == len(self.columns)
        assert all(self.columns.values())
        for column in self.columns.values():
            span |= {h ^ column for h in span}
        assert len(span) == 64
        # Nonzero binary torsion words must have weighted Hamming weight >= 4.
        for size in range(1, 4):
            for support in combinations(self.columns, size):
                if sum(self.weight(j) for j in support) < 4:
                    assert self.syndrome(sum(1 << j for j in support)) != 0
        assert any(self.weight(j) == 1 for j in self.columns)
        return all(self.min_cost(low, self.syndrome(upper), 0) >= 64
                   for low, upper in self.residues().items() if low)


def seeds():
    data = json.loads(SOURCE.read_text())
    rows = tuple(zip(data["residue_masks"], data["lift_masks"], strict=True))
    n30_columns = {j: 1 | (j << 1) for j in range(32)
                   if j not in data["omitted_coordinates"]}
    assert len(n30_columns) == 30
    # Fold equal coordinates 2 and 3 in the section. The combined parity column
    # is their XOR and the combined squared-length weight is two.
    assert all(((mask >> 2) ^ (mask >> 3)) & 1 == 0 for row in rows for mask in row)
    n29_columns = n30_columns.copy()
    n29_columns[2] ^= n29_columns.pop(3)
    n29 = Seed(n29_columns, 1 << 2,
               tuple((low & ~(1 << 3), high & ~(1 << 3)) for low, high in rows))
    n30 = Seed(n30_columns, 0, rows)
    assert n29.verify_minimum() and n30.verify_minimum()
    yield n29
    yield n30
    for j in data["omitted_coordinates"]:
        columns = n30_columns.copy()
        columns[j] = 1 | (j << 1)
        # Keep the original six residue words. Try both choices of a high bit
        # in each lifted generator on the restored coordinate.
        for lifts in range(1 << len(rows)):
            lifted = tuple((low, high | (((lifts >> i) & 1) << j))
                           for i, (low, high) in enumerate(rows))
            seed = replace(n30, columns=columns, rows=lifted,
                           restored_position=j, restored_lifts=lifts)
            if seed.verify_minimum():
                yield seed


def cases():
    for seed in seeds():
        yield seed, (), 64
        if seed.restored_lifts == 0:
            for j in sorted(seed.columns):
                yield seed, (j,), 64
        if seed.restored_position is None:
            yield seed, (), 48


def hyperplanes(seed):
    """Take the inverse image of every binary residue hyperplane in the seed."""
    for normal in range(1, 1 << len(seed.rows)):
        pivot = (normal & -normal).bit_length() - 1
        pivot_low, pivot_high = seed.rows[pivot]
        rows = []
        for i, (low, high) in enumerate(seed.rows):
            if i == pivot:
                continue
            if (normal >> i) & 1:
                low, high = low ^ pivot_low, high ^ pivot_high ^ (low & pivot_low)
            rows.append((low, high))
        yield normal, replace(seed, rows=tuple(rows))


def replacement_cases():
    for seed in seeds():
        if seed.restored_lifts == 0:
            for normal, sublattice in hyperplanes(seed):
                yield sublattice, (), 64, normal


def near_contact_cases():
    for seed in seeds():
        if seed.restored_lifts == 0:
            if seed.double_weight:
                yield seed, (2,), 62
            for offset in combinations(sorted(seed.columns), 2):
                yield seed, offset, 62


def run(binary, case):
    seed, offset_positions, minimum, *normal = case
    offset = sum(1 << j for j in offset_positions)
    active = sum(1 << j for j in seed.columns)
    rows = [f"{active} {seed.double_weight} {offset} {minimum} {len(seed.rows)}"]
    rows.extend(f"{low} {high}" for low, high in seed.rows)
    command = [binary, "--pairs"] if normal else [binary]
    raw = subprocess.run(command, input="\n".join(rows), capture_output=True, text=True, check=True)
    lines = raw.stdout.splitlines()
    checked, accepted, index_two, low, high, lattice_low, lattice_high = map(int, lines[0].split())
    assert checked == 1 << (len(seed.columns) - len(seed.rows))
    residues = seed.residues()

    def witness(low, high):
        distance = min(seed.min_cost(low ^ bits, high ^ seed.syndrome(upper ^ (low & bits)), offset)
                       for bits, upper in residues.items())
        assert distance >= minimum
        doubles_into_seed = offset in residues and seed.syndrome(low ^ residues[offset]) == 0
        return {"low_bits": low, "high_syndrome": high,
                "squared_distance_times_16": distance, "doubles_into_seed": doubles_into_seed}

    row = {"dimension": len(seed.columns), "restored_position": seed.restored_position,
           "restored_lifts": seed.restored_lifts, "shifted_positions": offset_positions,
           "minimum_squared_times_16": minimum, "cosets_checked": 64 * checked,
           "admissible_translates": accepted, "index_two_extensions": index_two}
    assert 0 <= index_two <= accepted <= row["cosets_checked"]
    if accepted:
        row["first_translate"] = witness(low, high)
    if index_two:
        row["first_index_two_translate"] = witness(lattice_low, lattice_high)
        assert row["first_index_two_translate"]["doubles_into_seed"]
    if normal:
        row["binary_hyperplane_normal"] = normal[0]
        edges, order_two_edges, cyclic_generators, *coordinates = map(int, lines[1].split())
        row.update(compatible_pairs=edges, compatible_order_two_pairs=order_two_edges,
                   cyclic_generators=cyclic_generators)

        def pair(coordinates):
            x, h, y, k = coordinates
            return [witness(x, h), witness(y, k),
                    witness(x ^ y, h ^ k ^ seed.syndrome(y & ~x))]

        if edges:
            row["first_compatible_pair_and_difference"] = pair(coordinates[:4])
        if order_two_edges:
            proofs = pair(coordinates[4:8])
            assert all(proof["doubles_into_seed"] for proof in proofs)
            row["first_index_four_pair_and_difference"] = proofs
        if cyclic_generators:
            x, h = coordinates[8:]
            row["first_index_four_generator_and_double"] = [witness(x, h), witness(0, seed.syndrome(x))]
    return row


def compile_screen(binary):
    subprocess.run(["c++", "-O3", "-std=c++20", "-Wall", "-Wextra", "-Wpedantic",
                    HERE / "screen.cpp", "-o", binary], check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=1)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--repack", action="store_true", help="replace one coset of each index-two sublattice")
    mode.add_argument("--near-contact", action="store_true", help="check two shifted axes at squared separation 31/8")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="quaternary-extensions-") as path:
        binary = Path(path) / "screen"
        compile_screen(binary)
        if args.repack:
            trial_cases = replacement_cases()
        elif args.near_contact:
            trial_cases = near_contact_cases()
        else:
            trial_cases = cases()
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            rows = list(pool.map(lambda case: run(binary, case), trial_cases))
    assert all(row["admissible_translates"] and row["index_two_extensions"]
               for row in rows if row["minimum_squared_times_16"] == 48)
    results = {"rows": rows}
    filename = "results.json"
    if args.repack:
        filename = "replacement_results.json"
    elif args.near_contact:
        filename = "near_contact_results.json"
        results = {"cases_checked": len(rows), "cosets_checked": sum(row["cosets_checked"] for row in rows),
                   "minimum_squared_times_16": 62,
                   "cases_with_admissible_translates": [row for row in rows if row["admissible_translates"]]}
    contents = json.dumps(results, indent=2) + "\n"
    output = HERE / filename
    if args.write:
        output.write_text(contents)
    else:
        assert output.read_text() == contents
    minimum = 62 if args.near_contact else 64
    targets = [row for row in rows if row["minimum_squared_times_16"] == minimum]
    print(f"Checked {len(targets)} offset/seed cases at squared separation {minimum / 16:g}.")
    for n in (29, 30, 31):
        selected = [row for row in targets if row["dimension"] == n]
        print(f"{n}D: {sum(row['admissible_translates'] for row in selected):,} admissible translates, "
              f"{sum(row['index_two_extensions'] for row in selected):,} index-two extensions")
        if args.repack:
            print(f"  {sum(row['compatible_pairs'] for row in selected):,} compatible pairs, "
                  f"{sum(row['compatible_order_two_pairs'] for row in selected):,} with order two, "
                  f"{sum(row['cyclic_generators'] for row in selected):,} with order four")


if __name__ == "__main__":
    main()
