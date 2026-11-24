"""Two direct polynomial minimal-height certificates, independent of LLL."""
import json
from pathlib import Path
p=47
here=Path(__file__).resolve().parent
def xy(k): return k%p,k//p
def prod(A,B):
    C=[(0,0)]*(len(A)+len(B)-1)
    for j,(a,b) in enumerate(A):
        for k,(c,d) in enumerate(B):
            x,y=C[j+k]
            C[j+k]=((x+a*c-b*d)%p,(y+a*d+b*c)%p)
    return C
for m,pair in json.loads((here/'polynomial_witnesses.json').read_text()).items():
    m=int(m)
    X,Y=[list(map(xy,pair[k])) for k in ('X','Y')]
    assert len(X)-1==16+2*m and len(Y)-1==24+3*m
    assert X[-1]!=(0,0) and Y[-1]!=(0,0)
    XX,YY=prod(prod(X,X),X),prod(Y,Y)
    assert len(XX)==len(YY)
    F=[(0,0)]*len(YY);F[1]=(p-1,0);F[47]=(1,0)
    assert all(((u-a)%p,(v-b)%p)==c for (u,v),(a,b),c in zip(YY,XX,F))
    # Since the affine coordinates are both polynomial there are no finite
    # pole intersections with O. At infinity, w=-z^16 X/(z^24 Y) has order m.
    assert (len(Y)-1)-(len(X)-1)-8==m
    print('Polynomial certificate depth',m,'squared norm',16+2*m)
