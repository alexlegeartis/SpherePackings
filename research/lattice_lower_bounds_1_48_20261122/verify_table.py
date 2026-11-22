"""Rebuild or check the 1–48D lattice-only lower-bound catalog."""
import argparse
from decimal import Decimal, ROUND_FLOOR, localcontext
from fractions import Fraction as F
from math import factorial, prod
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEST = ROOT / "LATTICE_LOWER_BOUNDS_1_48.md"
SLIP = "research/bacher_slip_29_30_20261122"
COHN = "https://cohn.mit.edu/sphere-packing/"
LABELS = {29: "4/sqrt(33)", 30: "8/sqrt(105)"}


def atan_interval(q):
    n = 45
    s = sum(F((-1)**k, (2*k + 1)*q**(2*k + 1)) for k in range(n))
    return s - F(1, (2*n + 1)*q**(2*n + 1)), s


def pi_interval():
    a, b = atan_interval(5)
    c, d = atan_interval(239)
    return 16*a - 4*d, 16*b - 4*c


def ball_squared(n, pi):
    if n % 2 == 0:
        return pi**n / factorial(n//2)**2
    return 2**(n+1) * pi**(n-1) / prod(range(1, n+1, 2))**2


def density_squared(factors):
    value = F(1)
    for a, b in factors:
        e = 2 * F(b)
        assert a > 0 and e.denominator == 1
        value *= F(a)**e.numerator
    return value


def center_formula(factors):
    return " * ".join(
        str(a) if b == "1" else f"{a}^{b}" if F(b).denominator == 1 else f"{a}^({b})"
        for a, b in factors
    ) or "1"


def floor_sqrt(lo, hi, digits=14):
    """Return a decimal floor, checked by exact squaring at both endpoints."""
    assert 0 < lo <= hi
    if lo == hi == 1:
        return "1"
    with localcontext() as c:
        c.prec = 100
        s = (Decimal(lo.numerator) / Decimal(lo.denominator)).sqrt()
        unit = Decimal(10)**(s.adjusted() + 1 - digits)
        floor = s.quantize(unit, rounding=ROUND_FLOOR)
    assert F(floor)**2 <= lo <= hi < F(floor + unit)**2
    return f"{floor:.{digits-1}e}"


def render():
    inp = json.loads((HERE / "cohn_inputs.json").read_text())
    rows = inp["rows"]
    assert [r[0] for r in rows] == list(range(1, 49))
    c36_sq = F(16777216, 3176523)**2
    assert density_squared(rows[35][1]) == c36_sq

    proofs = json.loads((ROOT / SLIP / "results.json").read_text())["rows"]
    improvements = {r["dimension"]: r for r in proofs}
    assert set(improvements) == LABELS.keys()
    pi_lo, pi_hi = pi_interval()
    t = ["""# Lattice sphere-packing lower bounds, dimensions 1–48

These are densities attained by **lattice packings**. A genuinely nonlattice
union of translates does not qualify. The baseline is the lattice-only table on
[Henry Cohn’s sphere-packing page](https://cohn.mit.edu/sphere-packing/#lattice-data),
checked on 22 November 2026. The [unrestricted lower-bound catalog](LOWER_BOUNDS_1_96.md)
is separate. This table does not claim optimality or exhaustive priority.

## Updates for Cohn’s lattice table

| Dimension | Listed center density | New exact center density | Occupied density (lower bound) | Gain vs. listed lattice | Construction |
|---:|---:|---:|---:|---:|---|
""".rstrip()]
    for n, factors, ref in rows:
        if n not in improvements:
            continue
        r = improvements[n]
        old = density_squared(factors)
        new = F(r["center_density_squared"])
        assert old == F(r["previous_lattice_center_density_squared"]) < new
        assert new == F(r["minimum_squared"], 4)**n / F(r["gram_determinant"])
        occupied = floor_sqrt(new*ball_squared(n, pi_lo), new*ball_squared(n, pi_hi))
        assert F(r["occupied_density_floor"])**2 <= new*ball_squared(n, pi_lo)
        gain = F(r["gain_percent_floor"])
        assert (1+gain/100)**2 <= new/old < (1+(gain+F(1, 10**11))/100)**2
        gain_display = f"{int(gain*100)//100}.{int(gain*100)%100:02d}%"
        t.append(f"| {n} | `{center_formula(factors)}` | **`{LABELS[n]}`** | **{occupied}** | **{gain_display}** | [Shifted Bacher lattice]({SLIP}/) |")
    t.append("""
The 29D and 30D packings are strictly denser than the lattice entries B29 and
Q30, respectively, but do not beat the unrestricted B29* and T30 entries.
The repository’s [36D lattice](constructions/dim36/) improved on KP36 by
18.97%; Cohn already lists that result in both tables, so it is not a new
update for his table.

## Full 1–48D catalog

For a lattice with shortest nonzero vector of length `2r` and Gram matrix
`G`, the center density is `r^n/sqrt(det(G))`. The occupied density is the
center density times the unit-ball volume `pi^(n/2)/Gamma(1+n/2)`. Decimal
occupied densities below are rounded downward. Source numbers point to the
construction references in the pinned Cohn table, as listed on the review date.

| Dimension | Exact center density | Occupied density (lower bound) | Construction or source |
|---:|---:|---:|---|
""".rstrip())
    for n, factors, ref in rows:
        sq, expression = density_squared(factors), center_formula(factors)
        source = f"[Cohn, ref. {ref}]({COHN}#{ref})" if ref is not None else "Integer lattice"
        if n in improvements:
            sq = F(improvements[n]["center_density_squared"])
            expression = f"**`{LABELS[n]}`**"
            source = f"[Shifted Bacher lattice]({SLIP}/) (improves Cohn’s listed lattice)"
        elif n == 36:
            expression = "`16777216/3176523`"
            source = "[Metric deformation](constructions/dim36/) (already listed by Cohn)"
        else:
            expression = f"`{expression}`"
        floor = floor_sqrt(sq*ball_squared(n, pi_lo), sq*ball_squared(n, pi_hi))
        t.append(f"| {n} | {expression} | {floor} | {source} |")
    t.append(f"""
## Reproduce

[Source factors and reference numbers](research/lattice_lower_bounds_1_48_20261122/cohn_inputs.json)
record the lattice-only table as read on 22 November. Run
`python3 research/lattice_lower_bounds_1_48_20261122/verify_table.py` to check
this catalog (add `--write` to rebuild it). This checks density arithmetic
with rational enclosures of π; it does not re-prove the source packings.
For the new lattices, run `python3 {SLIP}/verify.py` for
the separate complete minimum and determinant proof.
""".rstrip())
    return "\n".join(t) + "\n"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--write", action="store_true")
    a = p.parse_args()
    result = render()
    if a.write:
        DEST.write_text(result, encoding="utf-8")
    elif DEST.read_text(encoding="utf-8") != result:
        raise SystemExit("Lattice lower-bound catalog is out of date")
    print("PASS: 48 lattice rows, two gains against Cohn, established 36D tie")


if __name__ == "__main__":
    main()
