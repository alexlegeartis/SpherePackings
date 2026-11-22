#!/usr/bin/env python3
"""Verify the refined 96D MW92 lattice packing existence bound."""
import argparse
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import json, runpy

HERE=Path(__file__).resolve().parent
D=runpy.run_path(str(HERE.parent/'mw92_order3_dichotomy_20261122'/'verify.py'))
PARENT=D['PARENT']; A=D['A']; THRESHOLD=D['PI_HIGH']
B_DIAG=(Q(14,3),Q('5.27266'),Q('5.454682'),Q('5.656144'))
DEPTHS=(Q(34,3),Q('10.741'),Q('10.545318'),Q('10.343856'))
COUNT=(1,3,7,17)

def replay():
    G_A=D['eisenstein_gram']()
    info=D['verify_eisenstein'](G_A,3*A-32)
    G_B=D['gram'](B_DIAG)
    counters=[set() for _ in G_B]; signed_lifts=0
    for z,j,q in D['short'](G_B):
        signed_lifts+=1
        assert z[j]==1 and q+DEPTHS[j]>=16
        prefix=(z[0]%3,)+z[1:j] if j else ()
        counters[j].add(prefix)
    assert tuple(map(len,counters))==COUNT
    for k,d in zip(COUNT[1:],DEPTHS[1:]):
        assert k*THRESHOLD**46*d**46<factorial(46)*PARENT*D['covering']()
    da,db=D['det'](G_A),D['det'](G_B)
    assert da==Q(4129347471480001,16000000000000)
    assert db==Q(19889279550909941977,75000000000000000)
    covol=PARENT**2*max(da,db)
    csq=Q(4**96)/covol
    prev=json.loads((HERE.parent/'mw92_order3_dichotomy_20261122'/'results.json').read_text())['rows'][-1]
    init=json.loads((HERE.parent/'mw92_extensions_20261122'/'results.json').read_text())['rows']
    series_origin=Q(next(x['center_density_squared'] for x in init if x['dimension']==96 and x['kind']=='lattice existence'))
    old=Q(prev['center_density_squared'])
    unres=Q(next(x['center_density_squared'] for x in init if x['kind']=='periodic existence'))
    hist=Q(3**96,2**48)
    assert hist<series_origin<old<csq<unres
    return dict(dimension=96,squared_minimum=16,
                branch_A_gram=[[str(x) for x in r] for r in G_A],
                branch_B_gram=[[str(x) for x in r] for r in G_B],
                branch_A_det=str(da),branch_B_det=str(db),branch_A_checks=info,
                branch_B_depths=list(map(str,DEPTHS)),branch_B_counts=list(COUNT),
                branch_B_short_pairs=signed_lifts,
                max_gram_determinant=str(covol),center_density_squared=str(csq),
                occupied_density_floor=D['floor_sqrt'](csq*D['ball2'](96,D['PI_LOW'])),
                gain_vs_preceding_lattice_percent_floor=D['percent'](csq,old),
                gain_vs_pre_dichotomy_lattice_percent_floor=D['percent'](csq,series_origin),
                gain_vs_BCS96_percent_floor=D['percent'](csq,hist),
                unrestricted_bound_improved=False)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--write',action='store_true')
    args=p.parse_args()
    out=replay()
    f=HERE/'results.json'
    if args.write:f.write_text(json.dumps(out,indent=2)+'\n')
    else:assert out==json.loads(f.read_text())
    print(out['occupied_density_floor'],out['gain_vs_preceding_lattice_percent_floor'],
          out['gain_vs_pre_dichotomy_lattice_percent_floor'])
