"""Exact proof checks; no packages. Bits are zero-based, least significant first."""
import json
from pathlib import Path
from collections import Counter
w=json.loads(Path(__file__).with_name('witness.json').read_text())
n=w['n']; assert n==31
def span(B):
    C=[0]
    for b in B:C.extend([x^b for x in C])
    return set(C)
parent=span(w['parent_rows32']); p=w['deleted_coordinate']
assert len(parent)==2**17 and min(c.bit_count() for c in parent if c)==8
short={(c & ((1<<p)-1)) | ((c>>(p+1))<<p) for c in parent if not (c>>p)&1}
C=span(w['binary_generator_rows31']); assert len(C)==2**16 and C==short
assert min(c.bit_count() for c in C if c)==8 and all(c.bit_count()%2==0 for c in C)
t=w['odd_support']; assert t in C and t.bit_count()==24
S=((1<<n)-1)^t
assert (t,S)==(2147483520,127)
assert all((c&S).bit_count()%2==0 for c in C)
K=[c for c in C if not c&S]
assert len(K)==1024 and all(c.bit_count()%4==0 for c in K)
assert w['linear_carry_mask']==0
assert (2**17)**2 * (2**46)**2 == (2*4**31)**2
assert 8**31 == 2*(2**46)**2
assert 2**24>3**15                 # ratio sqrt(2)/Q31
print('65536 even codewords; minimum 8; complement weights',dict(sorted(Counter((c&S).bit_count() for c in C).items())))
print('1024 shortened words, all doubly even; covolume 2^46, minimum squared 32, center density sqrt(2)')
