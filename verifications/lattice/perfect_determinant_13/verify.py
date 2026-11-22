"""Compute an explicit strict upper bound below the rank-13 Mordell value."""

from fractions import Fraction
from math import comb
import json

from flint import arb, ctx, fmpq, fmpz


def farey_neighbors(target, max_denominator):
    """Consecutive fractions of bounded denominator bracketing a rational target."""
    num, den = target.as_integer_ratio()
    p0, q0, p1, q1 = 0, 1, 1, 0
    while den:
        digit, remainder = divmod(num, den)
        if q0 + digit * q1 > max_denominator:
            break
        p0, q0, p1, q1 = p1, q1, p0 + digit * p1, q0 + digit * q1
        num, den = den, remainder
    assert den and q1 > 0
    step = (max_denominator - q0) // q1
    return sorted((Fraction(p1, q1), Fraction(p0 + step*p1, q0 + step*q1)))


def verify():
    n = 13
    rank = n * (n + 1) // 2
    mordell_power = Fraction(729, 4)**13
    # D > 469 and b * 2^39 * D^7 <= 4^91 by the dyad determinant bound.
    assert 469**11 < mordell_power
    b_max = 4**rank // (2**(comb(n, 2)//2) * 469**7)
    max_denominator = b_max**n

    # Use an integer root only to find the neighbors. The following rational
    # comparisons check them against the exact irrational determinant bound.
    bits = 2*max_denominator.bit_length() + 96
    scale = 1 << bits
    scaled = (mordell_power.numerator * scale**11) // mordell_power.denominator
    root = int(fmpz(scaled).root(11))
    left, right = farey_neighbors(Fraction(root, scale), max_denominator)
    assert left**11 < mordell_power < right**11
    assert right.numerator*left.denominator - left.numerator*right.denominator == 1
    assert max(left.denominator, right.denominator) <= max_denominator
    assert left.denominator + right.denominator > max_denominator

    ctx.prec = bits + 256
    rational = lambda x: fmpq(x.numerator, x.denominator)
    volume = 128*arb.pi()**6/135135
    benchmark = volume/(arb(rational(mordell_power))**fmpq(1,22))
    improved = volume/arb(rational(right)).sqrt()

    def ceiling(x):
        return Fraction(int((x*scale).upper().ceil().unique_fmpz()), scale)

    upper = ceiling(improved)
    assert improved <= arb(rational(upper)) < benchmark
    gain = 100*(1-arb(rational(upper))/benchmark)
    assert gain > 0
    return {
        "dimension": n,
        "contact_determinant_limit": b_max,
        "determinant_denominator_limit": str(max_denominator),
        "farey_lower": str(left),
        "determinant_floor": str(right),
        "density_upper": str(upper),
        "benchmark_rational_ceiling": str(ceiling(benchmark)),
        "percent_below_exact_mordell": gain.str(12),
    }


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    print(json.dumps(verify(), indent=2))
