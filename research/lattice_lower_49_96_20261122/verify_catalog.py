#!/usr/bin/env python3
"""Build the 49–96 lattice-only table from exact squared center-density data."""
import argparse
from fractions import Fraction as Q
from math import factorial, isqrt, prod
from pathlib import Path
import json
from decimal import Decimal, localcontext

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

def arctan(a, stop):
    return sum(Q((-1)**i, (2*i+1)*a**(2*i+1)) for i in range(stop))

def pi_floor():
    # Even term counts are lower for arctan, odd counts upper.
    p = 16*arctan(5, 60) - 4*arctan(239, 61)
    assert Q(3141592653589793,10**15) < p < Q(22,7)
    return p

def ball_squared(n, p):
    if n % 2 == 0:
        return p**n / factorial(n//2)**2
    return 2**(n+1)*p**(n-1)/prod(range(1,n+1,2))**2

def float_floor(q, digits=12):
    """A rigorously downward-rounded decimal for sqrt(q)."""
    with localcontext() as ctx:
        ctx.prec = 80
        y = (Decimal(q.numerator)/Decimal(q.denominator)).sqrt()
        unit = Decimal(1).scaleb(y.adjusted() - digits + 1)
        v = q / Q(unit)**2
        k = isqrt(v.numerator // v.denominator)
        assert Q(k*unit)**2 <= q < Q((k+1)*unit)**2
        return f"{k*unit:.{digits-1}E}".lower()

def gain(q, p):
    """Downward-rounded 100*(sqrt(q/p)-1) to two decimal places."""
    v = (q/p) * 10**8
    k = isqrt(v.numerator // v.denominator)
    assert Q(k,10**4)**2 <= q/p < Q(k+1,10**4)**2
    assert k>=10000
    return f"{(k-10000)//100}.{(k-10000)%100:02d}%"

def exact(item,n):
    assert item['kind'] in ('published lattice','derived lattice','reported lattice floor',
                            'explicit lattice','lattice existence','audited lattice section')
    assert item['description'] and item['link']
    url=item['link']
    if not url.startswith('https://'):
        assert (ROOT / url.split('#')[0]).is_file(), url
    q=Q(item['center_density_squared'])
    assert q>0
    if 'minimum_squared' in item:
        assert Q(item['minimum_squared'])>0
        assert q==(Q(item['minimum_squared'])/4)**n/Q(item['gram_determinant'])
    return q

def build(data):
    assert [r['dimension'] for r in data['rows']]==list(range(49,97))
    lines=[
        '# Lattice-only sphere-packing lower bounds, 49D–96D','',
        'This catalog separates lattice packings from the [unrestricted 1–96D table](LOWER_BOUNDS_1_96.md).',
        'Each row compares a prior published lattice or an explicitly named derived lattice with',
        'the best verified repository lattice found in this review. A section lattice may be',
        'included even when adding extra translates yields a denser **nonlattice** packing;',
        'only the lattice period counts here. Lattice existence proofs are marked as such.',
        '', 'The center density is δ = (m/4)^(n/2) / sqrt(det G), where m is a proved lower',
        'bound on squared minimum and G is a lattice Gram. Ordinary density is Vₙ δ.',
        'Numerical columns are downward-rounded lower bounds. Gains are rounded down and use',
        'the upper end of a cited rounding interval when supplied; a dash means no',
        'verified increase is claimed. Prior values are the strongest comparisons found',
        'in the cited source audit, not a claim of exhaustive historical priority.','',
        '| Dim. | Prior lattice center density ≥ | Repository lattice center density ≥ | Gain over prior | Current occupied density ≥ |',
        '|---:|---|---|---:|---:|',
    ]
    p=pi_floor()
    improved=[]
    for r in data['rows']:
        n=r['dimension']; before=r['prior']; old=exact(before,n)
        new=r.get('repository'); now=exact(new,n) if new else old
        upper=Q(before.get('center_density_squared_upper',before['center_density_squared']))
        assert upper>=old
        if 'log2_interval' in before:
            u,v=map(Q,before['log2_interval'])
            assert old**u.denominator <= 2**(2*u.numerator)
            assert upper**v.denominator >= 2**(2*v.numerator)
        def cell(item,q):
            tag=' (derived)' if item['kind']=='derived lattice' else ' (reported floor)' if item['kind']=='reported lattice floor' else ' (existence)' if item['kind']=='lattice existence' else ''
            return f"{float_floor(q)} · [{item['description']}]({item['link']}){tag}"
        g='—'
        if new and now>upper:
            g='+'+gain(now,upper) if not r.get('disputed') else '+'+gain(now,upper)+' (repaired comparison only)'
            if not r.get('disputed'): improved.append(n)
        lines.append(f"| {n} | {cell(before,old)} | {cell(new,now) if new else '—'} | {g} | {float_floor(max(old,now)*ball_squared(n,p))} |")
    assert data['strict_comparison_dimensions']==improved
    lines.extend(['','## Scope and replay','',
      '- The 55D thesis target is disputed: the displayed lattice has an explicit shorter vector. The row compares with the [rigorous repair](research/tilted55_20261116/PROOF.md), not the printed target. It is excluded from the improvement count.',
      '- The 89–91D comparison intervals are bounds for a cited derived construction guarantee, not upper bounds on the optimal density. The same distinction applies to historical printed decimal prefixes in other rows.',
      '- The 79D and 86–91D periodic packings and the nine-translate 96D packing are not counted as lattices. In 86–91D the listed bare kernels form single lattices and retain the parent minimum lower bound.',
      '- This catalog makes no lattice-optimality claim. Printed and derived comparisons, exact squared values, and the proof paths are in [catalog.json](research/lattice_lower_49_96_20261122/catalog.json).',
      '', f"Verified improvements over the named non-disputed prior bounds: **{len(improved)} dimensions**.",
      '', 'Run the Python script in [the proof directory](research/lattice_lower_49_96_20261122/verify_catalog.py) to check the table, or add --write to rebuild it.', ''])
    return '\n'.join(lines)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--write',action='store_true')
    args=p.parse_args()
    output=build(json.loads((HERE/'catalog.json').read_text()))
    dest=ROOT/'LATTICE_LOWER_BOUNDS_49_96.md'
    if args.write: dest.write_text(output)
    else: assert dest.read_text()==output
    print('Verified 48 lattice-only rows and exact comparison arithmetic.')
