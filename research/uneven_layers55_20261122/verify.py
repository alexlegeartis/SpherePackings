"""Verify the minimum, period and density of the three-layer 55D packing."""

import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys

from flint import fmpq, fmpq_mat

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent / "cm58_20261116"),
                str(HERE.parent / "mw92_extensions_20261122")]
from verify_frontier import parents54_55
from verify_mw92 import display_floor, gain_floor, pi_interval, unit_ball_squared


def horizontal_minimum(displacement, depth):
    """Use the minimum in each of the three order-three cosets."""
    displacement %= 3
    if displacement <= 1 or displacement >= 2:
        return depth * min(displacement, 3 - displacement)**2
    alpha = displacement - 1
    return depth * (1 - alpha * (1 - alpha))


def audit():
    _, (parent, q, _, closest) = parents54_55(return_data=True)
    gram = fmpq_mat(parent)
    depth = F(str((q * gram * q.transpose())[0, 0]))
    assert depth == F(28, 3) and parent.det() == 3**25
    assert all(x.denominator == 1 for x in (3 * q).tolist()[0])
    # One simultaneous nearest vector in the two nonzero cosets attains
    # the interpolation formula throughout the middle segment.
    dots = (closest * gram * q.transpose()).tolist()
    index = dots.index([fmpq(str(-depth / 2))])
    common = fmpq_mat([closest.tolist()[index]])
    assert (common * gram * common.transpose())[0, 0] == fmpq(str(depth))
    assert ((common + q) * gram * (common + q).transpose())[0, 0] == fmpq(str(depth))

    # Heights are in units of sqrt(2/3). Horizontal shifts are in units of Q.
    steps = [(F(1), F(1)), (F(1), F(1)), (F(6, 5), F(9, 5))]
    unit_squared = F(2, 3)
    separation = []
    for length in range(1, 4):
        row = []
        for start in range(3):
            segment = [steps[(start + i) % 3] for i in range(length)]
            offset, height = map(sum, zip(*segment))
            squared = horizontal_minimum(offset, depth) + unit_squared * height**2
            assert squared >= 10
            row.append(str(squared))
        separation.append(row)
    assert min(height for _, height in steps) == 1
    assert 16 * unit_squared > 10

    # Remove 3Q from the period shift using the embedded parent lattice.
    offset, height = map(sum, zip(*steps))
    shift = q * fmpq(str(offset % 3))
    cross = gram * shift.transpose()
    height_squared = unit_squared * height**2
    bottom = F(str((shift * cross)[0, 0])) + height_squared
    assert bottom == 10  # Attainment by the period vector.
    period = fmpq_mat([a + b for a, b in zip(gram.tolist(), cross.tolist())]
                      + [cross.transpose().tolist()[0] + [str(bottom)]])
    determinant = F(str(period.det()))
    assert determinant == 3**25 * height_squared

    motif = []
    for level in range(3):
        vertical_coefficient = F(level, 1) / height
        horizontal = level*q - shift*fmpq(str(vertical_coefficient))
        motif.append([str(x) for x in horizontal.tolist()[0]] + [str(vertical_coefficient)])
    squared_density = 9 * F(5, 2)**55 / determinant
    earlier = json.loads((HERE.parent / "cross_coset55_20261116" / "results.json").read_text())
    previous = F(earlier["density_squared"])
    printed = F(5, 2)**55 / (2 * 3**24)
    assert squared_density > previous
    assert squared_density / printed == F(225, 361)

    # Constants in the radius-dependent upper bound within this family.
    cutoff = F(149, 150)
    assert depth * height_squared / 75 * cutoff**53 < 1
    assert F(14, 5)**2 + 14 * F(2, 5)**2 < F(149, 10)
    assert F(19, 5)**2 + F(14, 25) == 15
    assert 9 * F(1, 10) < F(114, 5) - 1
    assert 15 * cutoff**14 > F(225, 19)

    pi_lower, _ = pi_interval()
    return {
        "rows": [{
            "dimension": 55,
            "kind": "explicit three-translate periodic packing",
            "squared_minimum": 10,
            "center_density_squared": str(squared_density),
            "occupied_density_floor": display_floor(squared_density * unit_ball_squared(55, pi_lower)),
            "preceding_center_density_squared": str(previous),
            "gain_over_repository_percent_floor": gain_floor(squared_density, previous),
            "ratio_to_printed_challenge_center_density": "15/19",
        }],
        "period_gram_determinant": str(determinant),
        "period_gram": [[str(x) for x in row] for row in period.tolist()],
        "motif_in_period_basis": motif,
        "minimum_distance_squared_by_gap_count_and_start": separation,
        "simultaneous_nearest_vector_for_middle_segment": [str(x) for x in common.tolist()[0]],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    path = HERE / "results.json"
    report = audit()
    if args.write_results:
        path.write_text(json.dumps(report, indent=2) + "\n")
    else:
        assert report == json.loads(path.read_text())
    print("Verified the 55D minimum, period and density.")
