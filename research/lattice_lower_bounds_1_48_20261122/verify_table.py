"""Rebuild or check the 1–48D lattice-only lower-bound catalog."""
import argparse
from decimal import Decimal, ROUND_FLOOR, localcontext
from fractions import Fraction as F
from math import factorial, isqrt, prod
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEST = ROOT / "LATTICE_LOWER_BOUNDS_1_48.md"
PACK = "research/quaternary_29_30_20261122"
SLIP = "research/bacher_slip_29_30_20261122"
E7 = "research/e7_hermitian35_20261123"
E7L38 = "research/e7_lattice38_20261123"
E7L37 = "research/e7_37_section_20261123"
COHN = "https://cohn.mit.edu/sphere-packing/"
LABELS = {29: "1/sqrt(2)", 30: "1"}
METHODS = {29: "Section of quaternary-code lattice", 30: "Quaternary-code lattice"}


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

    proofs = json.loads((ROOT / PACK / "results.json").read_text())["rows"]
    predecessor = {r["dimension"]: F(r["center_density_squared"])
                   for r in json.loads((ROOT / SLIP / "results.json").read_text())["rows"]}
    improvements = {r["dimension"]: r for r in proofs}
    assert set(improvements) == LABELS.keys()
    hermitian = json.loads((ROOT / E7 / "results.json").read_text())
    assert hermitian["status"] == "PASS"
    e7_rows = {r["dimension"]: r for r in hermitian["rows"]}
    assert set(e7_rows) == set(range(34, 38))
    lift = json.loads((ROOT / E7L38 / "results.json").read_text())
    assert lift["status"] == "PASS" and lift["kind"] == "explicit lattice"
    assert lift["dimension"] == 38
    sec37 = json.loads((ROOT / E7L37 / "results.json").read_text())
    assert sec37["status"] == "PASS" and sec37["kind"] == "explicit lattice"
    assert sec37["dimension"] == 37 and F(sec37["center_density_squared"]) > F(sec37["preceding_center_density_squared"])
    assert F(sec37["preceding_center_density_squared"]) == F(e7_rows[37]["center_density_squared"])
    assert F(sec37["cohn_center_density_squared"]) == F(e7_rows[37]["preceding_center_density_squared"])
    e7_rows[37] = sec37
    e7_rows[38] = lift
    improvements.update(e7_rows)
    e7_path = {n: {37:E7L37,38:E7L38}.get(n,E7) for n in e7_rows}
    pi_lo, pi_hi = pi_interval()
    t = ["""# Lattice sphere-packing lower bounds, dimensions 1–48

[Companion lattice-only catalog for 49D–96D](LATTICE_LOWER_BOUNDS_49_96.md).

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
        if n in e7_rows:
            assert old == F(r.get("cohn_center_density_squared", r["preceding_center_density_squared"])) < new
            method, expression, package = r["method"], r["center_density"], e7_path[n]
        else:
            assert old == F(r["cohn_lattice_center_density_squared"])
            assert old < predecessor[n] == F(r["preceding_repository_center_density_squared"]) < new
            assert new == F(r["cohn_unrestricted_center_density_squared"])
            method, expression, package = METHODS[n], LABELS[n], PACK
        assert new == F(r["minimum_squared"], 4)**n / F(r["gram_determinant"])
        occupied = floor_sqrt(new*ball_squared(n, pi_lo), new*ball_squared(n, pi_hi))
        q = 10**8 * new/old
        pct = isqrt(q.numerator//q.denominator)-10000
        if F((2*pct+20001)**2, 4*10**8) <= new/old:
            pct += 1
        gain_display = f"{pct//100}.{pct%100:02d}%"
        t.append(f"| {n} | `{center_formula(factors)}` | **`{expression}`** | **{occupied}** | **{gain_display}** | [{method}]({package}/) |")
    t.append("""
The 29D and 30D lattices are strictly denser than B29 and Q30, respectively,
and match the listed unrestricted B29* and T30 lower bounds. They supersede
the repository’s [shifted Bacher lattices](research/bacher_slip_29_30_20261122/)
by 1.55% and 28.09%, respectively.

The new Hermitian E7 family raises lattice and unrestricted lower bounds
in 34–37D. The 38D coupled lift improves the listed lattice bound and
matches the repository’s earlier periodic bound. Its primitive 37D section
improves the preceding 37D lattice and unrestricted lower bound. In 36D the comparison is the earlier
[metric-deformation lattice](constructions/dim36/) that Cohn lists,
not the weaker historical KP36 packing.

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
        if n in e7_rows:
            r = e7_rows[n]
            sq = F(r["center_density_squared"])
            expression = f"**`{r['center_density']}`**"
            source = f"[{r['method']}]({e7_path[n]}/) (improves Cohn’s listed lattice)"
        elif n in improvements:
            sq = F(improvements[n]["center_density_squared"])
            expression = f"**`{LABELS[n]}`**"
            source = f"[{METHODS[n]}]({PACK}/) (improves Cohn’s listed lattice)"
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
For the new lattices, run `python3 {PACK}/verify.py` and
`python3 {E7}/verify.py` (with `basis.py` for the four Gram determinants).
The 38D package includes its own basis and coset-minimum verification.
The 37D section package replays that parent and checks the section basis.
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
    print("PASS: 48 lattice rows, seven gains against the listed values")


if __name__ == "__main__":
    main()
