"""Replay Poisson population debits against the complete parity census."""
import argparse
from fractions import Fraction
import importlib.util
import json
from math import lcm
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import zipfile

from flint import arb, arb_poly, ctx, fmpq


HERE = Path(__file__).resolve().parent
DISPLAY_SCALE = 10**30


def read_json(path):
    return json.loads(Path(path).read_text())


def rational(value):
    return fmpq(str(Fraction(str(value))))


def real(value):
    return arb(rational(value))


def volume(n):
    return arb.pi()**real(Fraction(n, 2)) / real(Fraction(n + 2, 2)).gamma()


def load_angular(folder):
    loader = importlib.util.spec_from_file_location('parent_repricing', folder/'verify.py')
    module = importlib.util.module_from_spec(loader)
    loader.loader.exec_module(module)
    return module


def poisson_count(record, previous, section):
    from independent_fourier import transform
    from polynomial_tools import sign_half

    ctx.prec = 384
    n = record['dimension']
    scale = rational(record['scale'])
    radius = rational(record['squared_radius'])
    upper = rational(record['density_upper'])
    earlier = rational(record['section_density_upper'])
    target = rational(record['density_threshold'])
    assert n > 1 and scale > 0 and radius >= 1 and upper >= target > 0 and earlier > 0
    assert upper == rational(previous)
    assert target == upper * fmpq(99, 100)
    if n == 25:
        assert section is None and volume(24) <= earlier
    else:
        assert section is not None and rational(section) <= earlier

    determinant = volume(n)/(2**n*upper)
    shortest = 2*volume(n-1)*target/(volume(n)*earlier)
    scaled = 4*arb.pi()**2/scale*shortest**2
    denominator = 10**9
    cut = fmpq(int((scaled*denominator).lower().floor().unique_fmpz()), denominator)
    assert cut > 0

    spatial, dual = transform(n, record['laguerre'])
    if dual(0)<0:
        determinant=volume(n)/(2**n*target)
    assert sign_half(spatial, scale)['pass_sign']
    assert sign_half(-dual, cut)['pass_sign']
    stationary = 2*spatial.derivative()-spatial
    roots = stationary.complex_roots()
    assert sum(int(m) for _,m in roots) == stationary.degree()

    windows = [tuple(map(rational, row)) for row in record['windows']]
    assert windows and all(1 <= a < b <= radius and w>0 for a,b,w in windows)
    assert all(b <= c for (a,b,w),(c,d,z) in zip(windows, windows[1:]))
    values = []
    polynomial = arb_poly([arb(x) for x in spatial.coeffs()])
    for a,b,w in windows:
        left,right = arb(a*scale),arb(b*scale)
        critical = [left,right]+[root.real for root,_ in roots
                                   if root.imag.contains(0) and not(root.real<left or root.real>right)]
        values.extend((polynomial(t)*((arb(scale)-t)/2).exp()/w).lower().fmpq()
                      for t in critical)
    shell = min(values)
    assert shell > 0
    ratio = ((2*arb.pi()/scale)**real(Fraction(n,2))*dual(0)/determinant-spatial(0))
    ratio *= (arb(scale)/2).exp()/(2*shell)
    unit = lcm(*(int(w.q) for a,b,w in windows))
    bound = Fraction(int((ratio*unit).upper().floor().unique_fmpz()), unit)
    assert bound == Fraction(record['opposite_pairs_upper']) and bound >= 0
    return bound


def replay_dimension(angular, archive, n, record):
    parent = load_angular(angular)
    previous_rows = read_json(angular/'certificates.json')
    previous = previous_rows[str(n)]
    section = None
    if n != 25:
        section = previous_rows[str(n-1)]['upper']
        child = subprocess.run([sys.executable, '-B', str(angular/'verify.py'),
                                str(archive), '--worker', str(n-1)],
                               check=True, capture_output=True, text=True)
        if child.stderr:
            print(child.stderr,end='',file=sys.stderr)

    root = parent.benchmark_root(archive,n)
    sys.path[:0]=[str(root),str(archive)]
    import independent_master
    from audit_geometry import check_prefix
    from audit_price import floor

    saved=[]
    original=independent_master.integer_census
    def collect(dim,edges,scores,unit,proof):
        result=original(dim,edges,scores,unit,proof)
        saved.append((dim,edges,scores,unit,proof,result))
        return result
    independent_master.integer_census=collect
    try:
        parent.verify_dimension(archive,n,previous)
    finally:
        independent_master.integer_census=original
    assert len(saved)==2 and all(x[0]==n for x in saved)
    _,edges,scores,unit,proof,profile=saved[-1]
    assert proof==read_json(root/read_json(root/'BOUNDS.json')[str(n)]['certificate'])
    tag=read_json(root/'BOUNDS.json')[str(n)]['certificate'].split('-rigorous',1)[0]
    source=read_json(root/(tag+'.json'))
    prefix=read_json(root/f'prefix-{n}-degree7.json')
    caps={Fraction(x):check_prefix(n,prefix[x]) for x in proof['prefix_ends']}

    theta=record['theta']
    assert theta['dimension']==n
    count=poisson_count(theta, previous['upper'], section)
    windows=[tuple(map(Fraction,row)) for row in theta['windows']]
    assert all(a in edges and b in edges for a,b,w in windows)
    price=Fraction(record['density']['price'])
    assert price>0
    debit=[max([Fraction(0)]+[w for a,b,w in windows if a<=l<r<=b])
           for l,r in zip(edges,edges[1:])]
    charged=[h-floor(unit*price*w) for h,w in zip(scores,debit)]
    modified=original(n,edges,charged,unit,proof)
    prior_census=independent_master.pooled(n,edges,profile,caps)
    new_census=independent_master.pooled(n,edges,modified,caps)
    difference=prior_census-new_census-price*count/((1<<n)-1)
    assert difference>0

    q0=sum(Fraction(c)*(-1)**j for j,c in enumerate(source['d']))
    tau=Fraction(source['tau'])
    assert q0>0 and tau>0
    scale=(real(tau)/8)**real(Fraction(n,2))/real(Fraction(n+2,2)).gamma()/real(q0)
    delta=int((scale*real(difference)*DISPLAY_SCALE).lower().floor().unique_fmpz())
    assert delta>0
    saving=Fraction(delta,DISPLAY_SCALE)
    old=Fraction(previous['upper'])
    upper=max(old-saving,Fraction(theta['density_threshold']))
    pinned=Fraction(read_json(archive/'BOUNDS.json')[str(n)]['upper'])
    expected=record['density']
    assert old==Fraction(expected['previous'])
    assert pinned==Fraction(expected['pinned'])
    assert saving==Fraction(expected['saving'])
    assert upper==Fraction(expected['upper'])<old<pinned
    print(json.dumps(dict(dimension=n,upper=str(upper),previous=str(old),
                          reduction_vs_current_percent=float(100*(1-upper/old)),
                          reduction_vs_pinned_percent=float(100*(1-upper/pinned)))),
          flush=True)


def main():
    table=read_json(HERE/'certificates.json')
    dims=tuple(sorted(map(int,table)))
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive',type=Path)
    parser.add_argument('dimensions',nargs='*',type=int,choices=dims)
    parser.add_argument('--angular-proof',type=Path,default=HERE.parent/'angular_kernel_followup_20261122')
    parser.add_argument('--worker',type=int,choices=dims,help=argparse.SUPPRESS)
    args=parser.parse_args()
    assert __debug__
    angular=args.angular_proof.resolve()
    parent=load_angular(angular)
    if args.worker is not None:
        replay_dimension(angular,args.archive.resolve(),args.worker,table[str(args.worker)])
        return

    def run(directory):
        parent.check_manifest(directory)
        for n in args.dimensions or dims:
            subprocess.run([sys.executable,'-B',str(HERE/'verify.py'),str(directory),
                            '--angular-proof',str(angular),'--worker',str(n)],check=True)

    location=args.archive.resolve()
    if location.is_dir():
        run(location)
    else:
        from hashlib import sha256
        assert sha256(location.read_bytes()).hexdigest()==parent.BENCHMARK_ARCHIVE_SHA256
        with TemporaryDirectory() as temp:
            with zipfile.ZipFile(location) as archive:
                archive.extractall(temp)
            directory=Path(temp)
            if not (directory/'MANIFEST.json').exists():
                [directory]=[p for p in directory.iterdir() if p.is_dir()]
            run(directory)


if __name__=='__main__':
    main()
