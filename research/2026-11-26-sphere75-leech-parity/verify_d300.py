"""Verify fixed 18D discriminant carrier and every index-two thinning. Requires Numpy and Sympy."""
import numpy as np,itertools
from sympy import Matrix
rows=[]
for i,j in itertools.combinations(range(24),2):
 for a,b in itertools.product((-4,4),repeat=2):
  x=[0]*24;x[i]=a;x[j]=b;rows.append(x)
C=[]
for n in range(4096):
 c=0
 for j in range(12):
  if n>>j&1:c^=0xAE3<<j
 c|=(c.bit_count()%2)<<23;C.append(c)
 if c.bit_count()==8:
  S=[i for i in range(24) if c>>i&1]
  for s in range(256):
   if s.bit_count()%2:continue
   x=[0]*24
   for j,i in enumerate(S):x[i]=2-4*(s>>j&1)
   rows.append(x)
 for j in range(24):
  x=[1-2*(c>>i&1) for i in range(24)];x[j]=-3*x[j];rows.append(x)
roots=np.array(rows,dtype=np.int16);del rows
assert len(roots)==196560

import json,tempfile,subprocess,os
from pathlib import Path
from collections import Counter
from fractions import Fraction
workspace=tempfile.TemporaryDirectory(); os.chdir(workspace.name)
dat={"d":5,"D":300,"W_scaled8":[[1,-1,-1,-1,-1,1,1,1,1,1,1,1,-1,-1,1,-1,1,1,3,-1,-1,-1,1,-1],[-1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,-1,1,-1,-1,-1,3,-1,1,-1,-1,1,-1],[-1,1,-3,-1,-1,1,1,-1,-1,1,-1,-1,1,-1,-1,-1,1,-1,-1,1,-1,-1,-1,-1],[0,2,0,2,2,0,0,0,0,0,0,0,0,0,2,0,0,-2,0,-2,0,2,0,2],[1,-3,1,1,-1,-1,-1,1,1,-1,1,-1,-1,1,-1,1,-1,-1,-1,-1,1,1,-1,1]],"char_basis":[[-2,2,2,0,2,-2,0,2,0,0,0,0,2,0,0,0,0,0,0,0,0,0,2,0],[0,-2,0,0,0,0,2,0,0,0,-2,2,0,0,0,-2,0,0,-2,0,0,0,-2,-2],[-1,1,1,-1,-1,-1,-3,-1,1,1,1,-1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,1],[0,0,2,-2,0,-2,0,0,-2,0,0,-2,0,0,0,0,0,-2,-2,0,2,0,0,0],[-1,-1,-1,3,-1,1,-1,-1,-1,1,-1,-1,-1,1,1,1,-1,-1,-1,-1,1,1,1,-1],[-1,-1,-1,-1,1,-1,3,-1,1,-1,-1,1,-1,1,1,-1,1,1,1,1,-1,1,1,1],[-1,-1,1,-1,-1,-1,-1,-1,-1,1,1,1,-1,-3,-1,-1,-1,-1,1,-1,1,-1,1,-1],[0,-2,2,-2,0,0,0,0,-2,0,0,0,0,0,-2,-2,0,0,2,0,0,0,0,-2],[1,3,1,-1,1,-1,-1,1,-1,1,-1,-1,1,-1,-1,1,-1,-1,1,1,-1,1,1,1],[0,-2,0,2,-2,-2,-2,0,0,0,0,0,-2,0,0,-2,0,0,0,0,0,0,0,2],[1,1,-1,1,-1,-1,1,1,1,1,-1,-3,-1,1,1,1,1,1,1,-1,1,-1,1,-1],[1,1,-1,-1,1,1,3,1,-1,1,-1,-1,1,1,-1,1,1,1,1,1,1,1,-1,1],[0,0,2,0,-2,-2,0,0,0,0,0,0,0,0,0,-2,0,0,-2,0,2,0,2,2],[0,0,-2,0,2,0,0,0,-2,0,2,0,2,0,0,0,0,2,0,0,2,0,2,0],[1,-1,1,-1,-1,1,1,1,1,1,1,3,-1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1],[0,-2,0,0,-2,0,-2,0,2,0,0,-2,0,-2,2,0,0,0,0,0,0,0,0,-2],[1,1,-1,-1,1,1,1,3,-1,1,-1,1,-1,-1,-1,-1,1,-1,1,1,1,1,-1,-1],[-1,1,1,1,1,-1,-1,1,-1,-1,-1,1,-1,-3,1,-1,1,-1,1,1,-1,-1,-1,1],[-2,0,2,0,2,0,0,2,0,0,0,0,0,0,2,2,-2,0,2,0,0,0,0,0],[-1,-1,1,-1,-3,1,1,1,1,-1,1,1,1,-1,-1,-1,1,-1,-1,1,-1,-1,1,-1],[-1,-1,1,-1,1,-1,-1,1,1,3,1,1,-1,-1,1,1,1,1,-1,-1,-1,1,1,-1],[0,0,0,0,2,0,0,0,0,-2,2,0,0,0,0,0,0,2,2,2,0,0,2,-2],[0,0,0,0,-2,0,0,2,2,0,0,2,0,-2,0,-2,-2,0,0,2,0,0,0,0],[-1,1,-1,1,-1,-1,-3,1,1,1,1,1,-1,-1,1,-1,1,-1,1,-1,-1,-1,-1,1]]}
Path('base19.data.json').write_text(json.dumps(dat))
import json,numpy as np
from sympy import Matrix
W=np.array(json.load(open('base19.data.json'))['W_scaled8']); G=W@W.T//8; d=len(W); M=Matrix(G.tolist());D=int(M.det());J=np.array(M.adjugate().tolist(),dtype=np.int32)
def inner(b):
 Q=np.indices((2*b+1,)*d,dtype=np.int16).reshape(d,-1).T-b;V=Q.astype(np.int32)@J;ww=(V*Q).sum(1);ix=np.argsort(ww,kind="stable");U,k=np.unique(V[ix]%D,axis=0,return_index=True);return U,Q[ix[k]],ww[ix[k]]
U,Y,sq=inner(3)
if len(U)!=D or max(sq)>=4*D:U,Y,sq=inner(4)

assert D==300
key={a.tobytes():i for i,a in enumerate(roots)}; Vbasis=np.array(dat['char_basis'],dtype=np.int16)
P=roots@W.T//8; tw=16**np.arange(d,dtype=np.int64); val=(P+4)@tw; uy,ic=np.unique(val,return_inverse=True); idu={int(v):i for i,v in enumerate(uy)}; assert all(int(t) in idu for t in (np.eye(d,dtype=int)+4)@tw)
def ec(b,v):
 v=int(v)
 while v:
  k=v.bit_length()-1
  if k in b:v^=b[k]
  else:b[k]=v;return
CM=(roots@Vbasis.T//8%2).astype(np.int64)@(1<<np.arange(24,dtype=np.int64)); B={}
for cm in CM:ec(B,cm)
assert len(B)==24
LL={}; j=idu[int(np.full(d,4)@tw)]
for cm in CM[ic==j]:ec(LL,cm)
assert len(LL)==19
nb=[int(CM[key[z.astype(np.int16).tobytes()]]) for z in W]
with open('base19.tab','w') as out:
 print(d,D,len(nb),*nb,len(LL),*LL.values(),file=out)
 for v,y,w in zip(U,Y,sq):
  t=idu.get(int((y+4)@tw));ms=np.array([]) if t is None else CM[ic==t];B={};ref=int(ms[0]) if len(ms) else 0
  for m in ms:ec(B,int(m)^ref)
  assert len(ms) or w==900
  print(*v,0,w,ref,len(ms),len(B),*B.values(),file=out)
assert max(sq)==900 and len(U)==D and max(sq)<=4*D # Cauchy bounds all profiles of equal or lower energy in [-3,3]^5
tab=open('base19.tab').read().splitlines();m=[int(x) for x in tab[0].split()];nb=np.array(m[3:3+m[2]],dtype=np.uint32);ids={tuple(u):i for i,u in enumerate(U)}
assert all([np.array_equal(list(map(int,row.split()[:d])),U[j]) for j,row in enumerate(tab[1:])])
cb=np.zeros((D,D),dtype=np.uint32)
for i in range(D):
 idx=[ids[tuple(v)] for v in (U[i]-U)%D]; q=(Y[i]-Y-Y[idx])@J; assert not np.any(q%D); cb[i,:]=np.bitwise_xor.reduce(((q//D)&1).astype(np.uint32)*nb,axis=1)
cb.astype('<u4').tofile('carry.bin')

miss=Y[np.where(sq==900)[0][0]]; ex=[(40087,173365),(40087,238),(40087,144483),(40087,160258),(40087,144503),(40087,160376),(40087,181560),(40087,160498),(40087,141204),(40087,174269),(40087,193114),(40087,185940),(40087,185988),(40087,179805),(40087,9514),(40087,146587),(40087,121956),(40087,138215),(126953,13997),(126953,83888)]
bb={};rm=None
for i,j in ex:
 z=roots[i]+roots[j];assert z@z==48 and np.array_equal(z@W.T//8,miss)
 m=int(((z@Vbasis.T//8%2).astype(np.int64))@(1<<np.arange(24,dtype=np.int64)))
 if rm is None:rm=m
 ec(bb,m^rm)
assert len(bb)==19
cpp="#include <bits/stdc++.h>\nusing namespace std;typedef unsigned long long U;const U RHS=1ULL<<63, LOW=RHS-1;\nstruct B{U a[61]={};bool add(U v){while(v&LOW){int k=63-__builtin_clzll(v&LOW);if(a[k])v^=a[k];else{a[k]=v;return true;}}return v==0;}U solve(){U bits=0;for(int k=0;k<61;k++)if(a[k]&&((a[k]>>63)^(__builtin_popcountll(a[k]&bits)&1)))bits|=1ULL<<k;return bits;}};\nstruct K{array<int,6>x;bool operator<(K const&o)const{return x<o.x;}};\nvector<unsigned int>carry;int d,D,nK,nL,F;vector<U>L;B base;map<K,int>ids;vector<K> keys;vector<int>w,num;vector<U>ms;vector<vector<U>>span;vector<int>dif;vector<int>ok;vector<int>s,var;vector<char>twos;double lower_reference;long long nodes;long long bound=1000000;int bestN,target,tmax=0;double bestR=0;mt19937 rng(4321117);\nint pos(int a,int b){return dif[a*D+b];}\nU real(const B& b){for(U v:L){auto B1=b;if(B1.add(v|RHS))return B1.solve();}return 0;}\nbool single(B& b,int v,int y){for(int i=0;i<s.size();i++){int ci=pos(v,s[i]);if(ok[ci]==2)continue;if(twos[i]||ok[ci]==0)return false;for(auto t:span[ci])if(!b.add(t))return false;if(!b.add(RHS|((1ULL<<y)^(1ULL<<var[i])^ms[ci]^carry[v*D+s[i]])))return false;}return true;}\nvoid report(const B&b,int score){if(score<=bestN)return;U h=real(b);if(!h)return;bestN=score;double cd=score*pow(1-F/(4.*D),(24.-d)/2.)/(2*sqrt(D));if(cd>bestR+1e-12){bestR=cd;cout<<setprecision(16)<<\"MIX \"<<d<<\" \"<<D<<\" \"<<F<<\" \"<<score<<\" \"<<cd<<\" \"<<(h&((1<<24)-1));for(int i=0;i<s.size();i++)cout<<\" \"<<s[i]<<\":\"<<(twos[i]?2:int(h>>var[i]&1));cout<<endl;}}\nvoid go(vector<int>const&P,B const&b,int score,int vars,int capDouble){\n int mult=capDouble?2:1;if(++nodes>bound||bestN>=target||score+mult*(int)P.size()<=bestN)return;\n if(score>bestN)report(b,score);\n for(int j=0;j<P.size()&&score+mult*((int)P.size()-j)>bestN;j++){\n  int v=P[j];vector<int> rest;bool full=capDouble;for(int a:s)if(ok[pos(v,a)]!=2)full=false;\n  if(full){for(int k=j+1;k<P.size();k++)if(ok[pos(v,P[k])]==2)rest.push_back(P[k]);\n   s.push_back(v);var.push_back(0);twos.push_back(1);go(rest,b,score+2,vars,capDouble);twos.pop_back();var.pop_back();s.pop_back();}\n  if(bestN>=target||nodes>bound)return;if(vars>=60)continue;\n  B b2=b;if(!single(b2,v,vars)||!real(b2))continue;\n  rest.clear();s.push_back(v);var.push_back(vars);twos.push_back(0);\n  for(int k=j+1;k<P.size();k++)if(ok[pos(v,P[k])]){if(!capDouble&&vars>=27){B tmp=b2;if(!single(tmp,P[k],vars+1)||!real(tmp))continue;}rest.push_back(P[k]);}\n  go(rest,b2,score+1,vars+1,capDouble);twos.pop_back();var.pop_back();s.pop_back();\n }\n}\nint main(int ac,char**av){if(ac>1)bound=atoll(av[1]);cin>>d>>D>>nK;for(int i=0;i<nK;i++){U x;cin>>x;}carry.resize(D*D);{ifstream f(\"carry.bin\",ios::binary);f.read((char*)carry.data(),D*D*4);if(!f)abort();}cin>>nL;L.resize(nL);for(U&x:L)cin>>x;w.resize(D);num.resize(D);ms.resize(D);span.resize(D);\n for(int i=0;i<D;i++){K k;for(int&x:k.x)cin>>x;keys.push_back(k);ids[k]=i;int n;cin>>w[i]>>ms[i]>>num[i]>>n;span[i].resize(n);for(auto&v:span[i])cin>>v;}\n dif.resize(D*D);for(int i=0;i<D;i++)for(int j=0;j<D;j++){K h;for(int k=0;k<6;k++)h.x[k]=(keys[i].x[k]-keys[j].x[k]+D)%D;dif[i*D+j]=ids.at(h);}\n int z=ids.at(K{});lower_reference=d==6?.07508468627929683:.08895810212922774;\n set<int>cuts;for(int i=0;i<D;i++){for(int h:{w[i],w[i]-2*D})if(h>=0&&h*100<=120*D)cuts.insert(h);}vector<int> doms(cuts.begin(),cuts.end());sort(doms.begin(),doms.end(),[&](int a,int b){return abs(double(a)/D-(d==5?8./15:2./3))<abs(double(b)/D-(d==5?8./15:2./3));});\n vector<U>_can(D);for(int i=1;i<D;i++){B b=base;for(U t:span[i])b.add(t);_can[i]=real(b);}\n for(int f:doms){F=f;if(f!=160)continue;ok.assign(D,0);vector<int> P;for(int c=0;c<D;c++){\n  int orig=(num[c]?w[c]:w[c]-2*D);\n  if(c==z||orig<=f)ok[c]=2;\n  else if(w[c]-2*D<=f&&_can[c])ok[c]=1;\n  if(c!=z&&ok[c]==2)P.push_back(c);\n }target=floor(2*sqrt(D)*lower_reference/pow(1-f/(4.*D),(24.-d)/2.)+1e-9)+1;\n target=13;if(P.size()<2)continue;\n vector<int>deg(D);for(int i:P)for(int j:P)deg[i]+=ok[pos(i,j)];\n sort(P.begin(),P.end(),[&](int a,int b){return deg[a]>deg[b];});\n // anchor multiplicity 2; zero translators WLOG for all codes carrying a double class.\n nodes=0;bestN=max(1,(int)floor(2*sqrt(D)*bestR/pow(1-F/(4.*D),(24.-d)/2.)-1e-6));s={z};var={0};twos={1};\n go(P,base,2,24,1);cout<<\"STAT \"<<f<<\" \"<<P.size()<<\" \"<<nodes<<\" \"<<bestN<<\" \"<<target<<endl;\n if(D==300&&f>=145&&f<=225&&d==5){\n vector<int> aa={0,4,288,245,52,287};\n for(int mask=63;mask>=1;mask--){if(mask!=1)continue;\n  s.clear();var.clear();twos.clear();for(int j=0;j<6;j++)if(mask>>j&1){s.push_back(aa[j]);twos.push_back(0);var.push_back(24+var.size());}\n  B bs=base;bool keep=true;for(int j=0;j<s.size();j++)for(int k=0;k<j;k++){int ix=pos(s[j],s[k]);if(ok[ix]==2)continue;if(!ok[ix]){keep=false;break;}for(U a:span[ix])if(!bs.add(a))keep=false;if(!bs.add(RHS|(ms[ix]^carry[s[j]*D+s[k]]^(1ULL<<var[j])^(1ULL<<var[k]))))keep=false;}if(!keep||!real(bs))continue;\n  vector<int> candidates;for(int i=0;i<D;i++){if(find(s.begin(),s.end(),i)!=s.end())continue;\n  B test=bs;if(single(test,i,24+s.size())&&real(test))candidates.push_back(i);}\n  vector<int>scoreC(D);for(int i:candidates){for(int q:s)scoreC[i]+=ok[pos(i,q)]==2;for(int q:candidates)scoreC[i]+=ok[pos(i,q)]==2;}\n  sort(candidates.begin(),candidates.end(),[&](int i,int j){return scoreC[i]>scoreC[j];});\n  nodes=0;go(candidates,bs,s.size(),24+s.size(),0);\n  cout<<\"SINGLE \"<<f<<\" \"<<mask<<\" \"<<candidates.size()<<\" \"<<nodes<<\" \"<<bestN<<\" \"<<target<<endl;\n  if(bestR>lower_reference+1e-10)break;\n }\n }\n\n if(bestR>lower_reference+1e-10)break;\n }\n}"
Path('run.cpp').write_text(cpp)
subprocess.run(['g++','-O3','-std=c++17','run.cpp','-o','proof'],check=True)
with open('base19.tab') as inp:r=subprocess.run(['./proof','1000000'],stdin=inp,capture_output=True,text=True,check=True)
print(r.stdout)
assert 'STAT 160 30 210 12 13' in r.stdout and 'SINGLE 160 1 90 572 12 13' in r.stdout
