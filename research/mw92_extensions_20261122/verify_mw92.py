"""Verify the MW92 extension certificates with exact rational arithmetic."""

import argparse
import json
from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import product
from math import factorial, isqrt, prod
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "cross_dimension_lower_20261115"))
from verify_counting import determinant
from verify_to96 import pi_interval

PARENT_MINIMUM = 16
PARENT_COVOLUME = 47**7
BASE_COMMIT = "9a9851887da20ca75825342bfb09e8ab031f4153"
DEPTHS = {
    1: Fraction("11.000980284"),
    4: Fraction("10.67439146"),
    7: Fraction("10.54531822"),
    8: Fraction("10.51475106"),
    18: Fraction("10.33101170"),
}


def covering_bound():
    n = 92
    mean = Fraction(2, n + 1)
    second = Fraction(24*n + 4*n*(n-1), prod(range(n, n+4)))
    third = Fraction(720*n + 144*n*(n-1) + 8*n*(n-1)*(n-2),
                     prod(range(n, n+6)))
    variance = second - mean**2
    centered_third = third - 3*mean*second + 2*mean**3
    a = n // 2
    correction = (1 + a*(a+1)*variance / (2*(1+mean)**2)
                  - a*(a+1)*(a+2)*centered_third / (6*(1+mean)**3))
    return (n+1) * Fraction(n, n+3)**a * correction


def short_vectors(gram):
    dim = len(gram)
    assert all(gram[i][j] == gram[j][i] for i in range(dim) for j in range(dim))
    assert all(determinant([row[:i] for row in gram[:i]]) > 0
               for i in range(1, dim + 1))
    det = determinant(gram)
    # Cauchy-Schwarz gives z_i^2 <= (G^-1)_ii z^T G z.
    limits = []
    for i in range(dim):
        minor = [[x for j, x in enumerate(row) if j != i]
                 for k, row in enumerate(gram) if k != i]
        diagonal = determinant(minor) / det
        limits.append(isqrt(int(PARENT_MINIMUM * diagonal)))
    for z in product(*(range(-b, b+1) for b in limits)):
        if any(z):
            norm = sum(gram[i][j]*z[i]*z[j] for i in range(dim) for j in range(dim))
            if norm < PARENT_MINIMUM:
                yield z, norm


def character(z, coefficients):
    return sum(x*y for x, y in zip(z, coefficients)) % 3


def periodic():
    a = Fraction(16, 3)
    b = 16 - DEPTHS[7]
    gram = [
        [a, -a/2, 0, 0],
        [-a/2, b, a/2-b, -b/2],
        [0, a/2-b, b, (b-a)/2],
        [0, -b/2, (b-a)/2, b],
    ]
    first, second = (1, 0, 1, 1), (0, 1, 2, 1)
    assert len({(character(z, first), character(z, second))
                for z in product(range(3), repeat=4)}) == 9
    minima = {}
    for z, norm in short_vectors(gram):
        color = character(z, first), character(z, second)
        assert color != (0, 0)
        depth = DEPTHS[7] if color[1] else DEPTHS[1]
        assert norm + depth >= PARENT_MINIMUM
        minima[color] = min(minima.get(color, norm), norm)
    assert len(minima) == 8
    return gram, {
        "translate_count": 9,
        "characters_mod_3": [list(first), list(second)],
        "short_vertical_minima_by_color": {
            str(color): str(norm) for color, norm in sorted(minima.items())
        },
    }


def shear(dim):
    counts = (1, 4, 8, 18)[:dim]
    heights = [16 - DEPTHS[count] for count in counts]
    gram = [[heights[i] if i == j else heights[min(i, j)]/2
             for j in range(dim)] for i in range(dim)]
    forbidden = [set() for _ in range(dim)]
    for z, norm in short_vectors(gram):
        j = max(i for i, x in enumerate(z) if x)
        assert abs(z[j]) == 1
        if z[j] > 0:
            forbidden[j].add(z[:j])
            assert norm + DEPTHS[counts[j]] >= PARENT_MINIMUM
    assert tuple(map(len, forbidden)) == counts
    return gram, {"forbidden_centers_by_step": list(counts)}


def unit_ball_squared(n, pi_bound):
    if n % 2 == 0:
        return pi_bound**n / factorial(n//2)**2
    return 2**(n+1) * pi_bound**(n-1) / prod(range(1, n+1, 2))**2


def display_floor(value_squared):
    with localcontext() as ctx:
        ctx.prec = 60
        value = (Decimal(value_squared.numerator) / Decimal(value_squared.denominator)).sqrt()
        unit = Decimal(1).scaleb(value.adjusted() - 14)
        scale = Fraction(unit)
    ratio = value_squared / scale**2
    k = isqrt(ratio.numerator // ratio.denominator)
    assert (k*scale)**2 <= value_squared < ((k+1)*scale)**2
    return f"{k*unit:.14E}"


def gain_floor(value_squared, prior_squared, decimals=4):
    divisor = 10**decimals
    scaled = value_squared * (100*divisor)**2 / prior_squared
    k = isqrt(scaled.numerator // scaled.denominator) - 100*divisor
    sign = "-" if k < 0 else ""
    return f"{sign}{abs(k)//divisor}.{abs(k)%divisor:0{decimals}d}"


def comparison_squares():
    older = HERE.parent / "adaptive_layers_20261115" / "covering_results.json"
    previous = {row["dimension"]: Fraction(row["density_squared"])
                for row in json.loads(older.read_text())["rows"]}
    for filename in ["RESULT.json", "RESULT95.json"]:
        prior = HERE.parent / "mw92_flags9495_verified_20261117" / filename
        row = json.loads(prior.read_text())
        previous[row["dimension"]] = Fraction(row["center_density_squared"])

    # MW92 multiplied by minimum-matched A1, A2, A3. BCS Table 1 in 96D.
    external = {92 + r: Fraction(2**184, 47**14 * 2**r * (r+1))
                for r in range(1, 4)}
    external[96] = Fraction(3**96, 2**48)
    return previous, external


def certify():
    pi_lower, pi_upper = pi_interval()
    covering = covering_bound()
    for count, depth in DEPTHS.items():
        assert 0 < depth < PARENT_MINIMUM
        assert count * pi_upper**46 * depth**46 < factorial(46) * PARENT_COVOLUME * covering

    previous, external = comparison_squares()
    rows = []
    constructions = [(92+r, "lattice existence", *shear(r)) for r in range(1, 5)]
    constructions.append((96, "periodic existence", *periodic()))
    for n, kind, gram, checks in constructions:
        covolume_squared = PARENT_COVOLUME**2 * determinant(gram)
        center_squared = Fraction(4**n) / covolume_squared
        assert center_squared > external[n]
        rows.append({
            "dimension": n,
            "kind": kind,
            "squared_minimum": PARENT_MINIMUM,
            "vertical_gram": [[str(x) for x in row] for row in gram],
            "covolume_per_center_squared": str(covolume_squared),
            "center_density_squared": str(center_squared),
            "occupied_density_floor": display_floor(center_squared * unit_ball_squared(n, pi_lower)),
            "external_lattice_comparison_center_density_squared": str(external[n]),
            "gain_over_external_lattice_comparison_percent_floor": gain_floor(center_squared, external[n]),
            "repository_unrestricted_comparison_center_density_squared": str(previous[n]),
            "gain_over_repository_unrestricted_comparison_percent_floor": gain_floor(center_squared, previous[n]),
            **checks,
        })
    assert Fraction(rows[0]["center_density_squared"]) == previous[93]
    assert Fraction(rows[-1]["center_density_squared"]) > previous[96]
    return {
        "comparison_base_commit": BASE_COMMIT,
        "covering_density_lower_bound": str(covering),
        "strict_torus_depths_by_ball_count": {str(k): str(v) for k, v in DEPTHS.items()},
        "rows": rows,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    result = certify()
    path = HERE / "results.json"
    if args.write_results:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        assert result == json.loads(path.read_text())
    for row in result["rows"]:
        print(row["dimension"], row["kind"], row["occupied_density_floor"],
              row["gain_over_repository_unrestricted_comparison_percent_floor"] + "%")
