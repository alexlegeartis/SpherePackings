#!/usr/bin/env python3
"""Build the 49–96 lattice-only table from exact squared center-density data."""
import argparse
from fractions import Fraction as Q
from math import factorial, isqrt, prod
from pathlib import Path
import json
import re
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
    if item['kind']=='derived lattice': derived_tag(item)
    url=item['link']
    if not url.startswith('https://'):
        assert (ROOT / url.split('#')[0]).is_file(), url
    q=Q(item['center_density_squared'])
    assert q>0
    if 'minimum_squared' in item:
        assert Q(item['minimum_squared'])>0
        assert q==(Q(item['minimum_squared'])/4)**n/Q(item['gram_determinant'])
    return q

def derived_tag(item):
    if item['kind'] != 'derived lattice':
        return ''
    years=item['source_years']
    assert years and years==sorted(set(years)) and all(isinstance(y,int) and 1000<=y<=9999 for y in years)
    return ' (derived; ' + '/'.join(map(str,years)) + ' source' + ('s' if len(years)>1 else '') + ')'

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
        'in the cited source audit, not a claim of exhaustive historical priority.',
        'Years in derived rows date the cited inputs, not publication of the derived bound.','',
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
            tag=derived_tag(item) if item['kind']=='derived lattice' else ' (reported floor)' if item['kind']=='reported lattice floor' else ' (existence)' if item['kind']=='lattice existence' else ''
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

BEGIN = '<!-- lattice-lower-1-96:start -->'
END = '<!-- lattice-lower-1-96:end -->'
HEADER = '## Lattice-only lower bounds in 1D–96D'

def readme_high_rows(data):
    p=pi_floor()
    p_hi=16*arctan(5,61)-4*arctan(239,60)
    assert p < p_hi < Q(22,7)
    improved=set(data['strict_comparison_dimensions'])
    lines=[]
    def upper_decimal(value_squared,digits=12):
        x=Decimal(float_floor(value_squared,digits))
        step=Decimal(1).scaleb(x.adjusted()+1-digits)
        high=x+step
        assert Q(high)**2 > value_squared
        return f"{high:.{digits-1}E}".lower()
    for row in data['rows']:
        n=row['dimension']
        if n not in improved: continue
        prev,new=row['prior'],row['repository']
        lo=Q(prev['center_density_squared'])
        hi=Q(prev.get('center_density_squared_upper',prev['center_density_squared']))
        now=Q(new['center_density_squared'])
        assert now>hi and not row.get('disputed')
        old_display=float_floor(lo*ball_squared(n,p))
        if hi>lo:
            old_display=f"[{old_display}, {upper_decimal(hi*ball_squared(n,p_hi))}]"
        elif prev['kind']=='reported lattice floor':
            old_display=f"≥{old_display}"
        else:
            old_display=f"≈{old_display}"
        source=f"[{prev['description']}]({prev['link']})"
        source += derived_tag(prev)
        new_display=float_floor(now*ball_squared(n,p))
        method=f"[{new['description']}]({new['link']})"
        method += derived_tag(new)
        if new['kind']=='lattice existence': method += ' (existence proof)'
        lines.append(f"| {n} | {old_display}<br>{source} | **≥{new_display}**<br>{method} | **+{gain(now,hi)}** |")
    assert len(lines)==len(improved)
    return lines

def write_readme(text, high_rows):
    anchor='<!-- lattice-upper-20261116:start -->'
    if BEGIN in text:
        assert text.count(BEGIN)==text.count(END)==1
        start=text.index(BEGIN)
        end=text.index(END,start)+len(END)
    else:
        # Migrate the former two separate tables without editing their rows.
        legacy='## Lattice-only lower bounds\n'
        assert text.count(legacy)==text.count(anchor)==1
        start=text.index(legacy)
        end=len(text[:text.index(anchor)].rstrip())
    low=[]
    for line in text[start:end].splitlines():
        match=re.match(r'^\| (\d+) \|',line)
        if match and int(match[1])<=48:
            low.append((int(match[1]),line))
    assert low and [n for n,_ in low]==sorted({n for n,_ in low})
    lines=[BEGIN,HEADER,'',
        'Verified improvements over the named prior lattice comparisons. Densities are',
        'occupied-volume fractions. Periodic packings that are not lattices are excluded.',
        '[Full 1D–48D catalog](LATTICE_LOWER_BOUNDS_1_48.md) ·',
        '[Full 49D–96D catalog](LATTICE_LOWER_BOUNDS_49_96.md)',
        '',
        '| n | Prior lattice density / source | Our lattice density / construction | Increase |',
        '|---:|---|---|---:|',
        *[line for _,line in low],*high_rows,'',
        'The 29D and 30D rows improve Cohn’s listed lattice bounds and tie his listed',
        'unrestricted bounds. Cohn already lists the 36D result; the gain is over KP36.',
        'From 49D onward, the gains are rounded down and use the upper end of a prior',
        'comparison interval where applicable. Years in derived rows date the cited',
        'sources, not publication of the derived bound. The catalogs include ties, source',
        'notes and exact center densities. The comparisons do not establish global',
        'optimality or exhaustive historical priority.',END]
    return text[:start]+'\n'.join(lines)+text[end:]

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--write',action='store_true')
    args=p.parse_args()
    output=build(json.loads((HERE/'catalog.json').read_text()))
    dest=ROOT/'LATTICE_LOWER_BOUNDS_49_96.md'
    readme=ROOT/'README.md'
    rows=readme_high_rows(json.loads((HERE/'catalog.json').read_text()))
    original=readme.read_text()
    rendered=write_readme(original,rows)
    if args.write:
        dest.write_text(output)
        readme.write_text(rendered)
    else:
        assert dest.read_text()==output
        assert original==rendered
    print('Verified 48 lattice-only rows, the README summary, and exact comparison arithmetic.')
