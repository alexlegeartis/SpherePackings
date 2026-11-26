"""Exact weight-two input verification for the continuous bank bound."""
from pathlib import Path
import runpy
here = Path(__file__).resolve().parent
base = runpy.run_path(str(here.parent/'qr18_layer_cover_20261126'/'verify.py'))
N,M,Q,four,badQ = (base[k] for k in ('N','M','Q','four','badQ'))
assert len(Q)==512 and len(badQ)==344
for i in range(N):
    for j in range(i):
        y=(1<<i)|(1<<j)
        E=sorted((c^y for c in Q if c not in (0,M)),key=lambda x:(x.bit_count(),x))
        assert four(E) is None
assert badQ=={y for y in range(1<<N) if y.bit_count() in (0,1,2,16,17,18)}
print("All 153 weight-two endpoint failures and entire single-code exception set verified.")
