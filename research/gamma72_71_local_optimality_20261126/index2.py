"""First constructive supercell test, not an all-supercell certificate."""
import re, numpy as np, sys, hashlib

from pathlib import Path
raw=Path(sys.argv[1] if len(sys.argv)>1 else Path(__file__).resolve().parents[1]/"novel_routes_20261113"/"higher_parent"/"neb72.html").read_bytes()
assert hashlib.sha256(raw).hexdigest()=="92daeaf0d93672e82ca5e5bc4479e5775daf59029aa4e983fb464e36ef002044"
html=raw.decode()
def entries(label):
    part=html.split('<a NAME="'+label+'">')[1].split('<p><li>')[0]
    return list(map(int,re.findall(r'-?\d+',re.sub('<[^>]*>',' ',part))))
g=entries('GRAM'); assert g.pop(0)==72
G=np.array(g,dtype=np.int64).reshape(72,72); data=entries('SUBGROUP_GENERATORS');assert data.pop(0)==6
gens=[]
for offset in range(0,len(data),5185):
    assert data[offset]==72
    t=np.array(data[offset+1:offset+5185],dtype=np.int64).reshape(72,72)
    assert np.array_equal(t@G@t.T,G);gens.append(t)
v=np.zeros(72,dtype=np.int64);v[44]=v[51]=1; gv=G@v; assert v@gv==10
rows=np.eye(72,dtype=np.int64)[np.diag(G)==8]
for i in range(3):
    rows=np.unique(np.vstack([rows]+[rows@t for t in gens]),axis=0)
    assert abs(rows).max()<5000
t=rows@gv; S=np.unique(rows[abs(t)==1]*t[abs(t)==1,None],axis=0)
assert ((S@G*S).sum(1)==8).all() and (S@gv==1).all()
Y=10*S-v
def ranks(bits):
    result=[]
    for parity in (0,1):
        selected=Y[(bits&1)==parity][:152]
        mat=selected.copy()%101
        rank=0
        for j in range(72):
            nz=np.flatnonzero(mat[rank:,j])
            if not len(nz):continue
            k=rank+int(nz[0]);mat[[rank,k]]=mat[[k,rank]]
            pivot=mat[rank]*pow(int(mat[rank,j]),-1,101)%101
            mat[rank+1:]=(mat[rank+1:]-mat[rank+1:,j:j+1]*pivot)%101
            rank+=1
        result.append(rank)
    return result
# Indicator x_i of a lift differs from its value on the fixed representative
# only by a constant, irrelevant because both parity subfamilies are checked.
out=[]
for i in range(72):
    r=ranks(S[:,i]);out.append(r)
print("actual distinct oriented contacts sampled",len(S),"of orbit collection",len(rows))
print("coordinate characters full rank both:",sum(a==b==71 for a,b in out),"/72")
if min(min(a) for a in out)<71:print("below full",[(i,x) for i,x in enumerate(out) if min(x)<71])
assert all(a==b==71 for a,b in out)
