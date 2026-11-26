"""Count exact contact degrees in the 14D periodic Hamming–Witt packing."""
from collections import Counter
from functools import lru_cache
from itertools import product
from pathlib import Path
from runpy import run_path

d=run_path(str(Path(__file__).with_name("verify.py")))
label, words = d["label"], d["words"]
period=range(6)

@lru_cache(None)
def degree_in_sheet(source_type, center, target_type, target_label, q):
    target=words(target_type,target_label)
    if source_type!=target_type:
        return len(target) if q==1 else 0
    k=4-q
    return (1<<k)*sum((center^w).bit_count()==k for w in target)

# Q(dm,dn)>4 outside [-2,2]^2 because Q≥(dm²+dn²)/2.
offsets=[(dm,dn,dm*dm+dm*dn+dn*dn) for dm,dn in product(range(-2,3),repeat=2)]
hist=Counter()
witt=Counter()
hamming=Counter()
for m,n in product(period,repeat=2):
    src,p=label(m,n)
    neighbors=[(*label(m+dm,n+dn),q) for dm,dn,q in offsets if q<=4]
    for a in words(src,p):
        # Self-translates of one motif word contribute the 24 axial contacts.
        count=24+sum(degree_in_sheet(src,a,t,w,q) for t,w,q in neighbors)
        hist[count]+=1
        (witt if src else hamming)[count]+=1
assert sum(hist.values())==4800 and max(hist)==1620
assert witt[1620]==16 and hamming=={1464:3072}
src,p=label(0,0)
assert 3123 in words(src,p)
parts=Counter({0:24})
for dm,dn,q in offsets:
    if q<=4: parts[q]+=degree_in_sheet(src,3123,*label(dm,dn),q)
assert [parts[i] for i in (0,1,3,4)]==[840,768,12,0]
print(dict(sorted(hist.items())))
print("PASS: max kissing shell 1620; central shell degree split 840+768+12.")
