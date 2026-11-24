"""Replay rational interval certificates for the paired E8 period lemma."""
from itertools import product
from fractions import Fraction as F
from pathlib import Path
from collections import Counter
import gzip,json
from eight import run

def scalar(parts):
    assert len(parts)==4 and all(isinstance(x,int) and 0<=x<=100 for x in parts)
    u,v,w,t=parts
    q,h=divmod(t,3)
    return u,v,w,q+(h>=1),q+(h>=2),q

def a_of(D,b,c,d):
    return (D-b*c*d)/(b*c+b*d+c*d+3*b*c*d)

def tag(u):
    return tuple(k%2 for k in u)

def per(box,u):
    sq=[n*n for n in u]
    return max(4*(a_of(box[0][1],b,c,d)*sum(u)**2+
                  b*sq[0]+c*sq[1]+d*sq[2]+sum(sq[3:]))
                  for b,c,d in product(*box[1:]))

def thresholds(box,u,v,p,q):
    cc=[x*y for x,y in zip(u,v)]
    aL=a_of(box[0][0],*[x[1] for x in box[1:]])
    aH=a_of(box[0][1],*[x[0] for x in box[1:]])
    m=sum(u)*sum(v)
    C=aL*m+sum(x[0]*c for x,c in zip(box[1:],cc))+sum(cc[3:])
    D=aH*m+sum(x[1]*c for x,c in zip(box[1:],cc))+sum(cc[3:])
    k=max(F(0),C,-D)
    base=max(p/2,q/2,(p+q)/4)
    mid=max(base,p/4+q/2-k,p/2+q/4-k)
    top=max(mid,(p+q)/2-2*k)
    return mid,top

def walk(box,cert,stat):
    op=cert[0];stat[op]+=1
    if op==0:
        assert len(cert)==1 and any(box[i][0]>box[i+1][1] for i in (1,2))
        return
    if op==4:
        _,j,x,y=cert
        assert 0<=j<=3
        l,h=box[j]
        m=(l+h)/2
        assert l<m<h
        walk(box[:j]+[(l,m)]+box[j+1:],x,stat)
        walk(box[:j]+[(m,h)]+box[j+1:],y,stat)
        return
    u=scalar(cert[1]); assert any(tag(u))
    p=per(box,u); L=box[0][0]
    def good(k,g):
        return g>=0 and 2**k*(g/8)**24<=16800*L**4
    # For any r>=8 the product of the carrier's occupancy and min^r
    # is bounded by 2^40 p^r when fixed in one nonzero parity direction.
    for e in (8,24):
        assert max(F(1,4)**e,F(1,16)*F(1,2)**e)<=F(1,256)
    if op==1:
        assert len(cert)==2 and good(40,p)
        return
    assert op==2 and len(cert)==3
    v=scalar(cert[2]);assert any(tag(v)) and tag(v)!=tag(u)
    q=per(box,v)
    s,t=thresholds(box,u,v,p,q)
    p=min(p,q)
    assert s>=p/2  # low-g inherits at most the cut with coefficient 2^40
    assert good(40,min(p,s))
    assert good(36,min(p,t))
    assert good(32,p)
    if not good(40,p):stat['paired essential']+=1
    if not good(36,p):stat['correlation essential']+=1

def main():
    run()
    data=json.loads(gzip.decompress(Path(__file__).with_name('proof.json.gz').read_bytes()))
    for cube,low,high,parts in data:
        cube,low,high=map(F,(cube,low,high))
        stat=Counter()
        cur=low
        for b,cert in parts:
            box=[tuple(map(F,I))for I in b]
            assert box[1:]==[(F(2),cube)]*3 and box[0][0]==cur<box[0][1]
            cur=box[0][1]
            walk(box,cert,stat)
        assert cur==high
        print('PASS',cube,low,high,dict(stat))
if __name__=='__main__':main()
