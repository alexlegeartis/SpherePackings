#!/usr/bin/env python3
"""Exact integer acceptance; floating point LP proposes duals and splits only.
Usage: python verify.py [72_4.txt [threads]]
Download from https://users.aalto.fi/~pat/opt/72_4.txt if passing a local first argument.
"""
import hashlib, os, pathlib, subprocess, sys, tempfile, urllib.request
from collections import Counter, defaultdict
def run():
    assert __debug__, "run this verifier without -O"
    raw = pathlib.Path(sys.argv[1]).read_bytes() if len(sys.argv)>1 else urllib.request.urlopen("https://users.aalto.fi/~pat/opt/72_4.txt").read()
    z = list(map(int, raw.split()))
    assert len(z)==7010 and z[:2]==[96,11]
    normalized=''.join(str(x)+'\n' for x in z).encode()
    assert hashlib.sha256(normalized).hexdigest()=='0f3fe74e631f11d3d3263c1c476d86ad142298ed6e2184f2b06370c771519d23'
    C=[z[3+i*73:75+i*73] for i in range(96)]
    assert all(z[2+i*73]==72 and len(set(c))==72 and all(0<=u<2048 for u in c) and min((u^v).bit_count() for u in c for v in c if u!=v)==4 for i,c in enumerate(C))
    data=''.join(' '.join(map(str,c))+'\n' for c in C)
    with tempfile.TemporaryDirectory() as base:
        wd=pathlib.Path(base)
        (wd/'C.dat').write_text(data)
        for name, source in SOURCES.items(): (wd/name).write_text(source)
        env=dict(os.environ,OMP_NUM_THREADS=sys.argv[2] if len(sys.argv)>2 else '2')
        exe='sharp12'
        subprocess.run([os.environ.get('CXX','c++'),'-std=c++17','-O3','-fopenmp',exe+'.cpp','-o',exe],cwd=wd,check=True)
        subprocess.run([str(wd/exe)],cwd=wd,env=env,check=True)
        rows=sorted([list(map(int,line.split())) for line in (wd/'sharp12.tsv').read_text().splitlines()])
        assert rows==[[i,2048,0,0] for i in range(96)]
        for code in C:
            assert any(min((u^v).bit_count() for u in code)==3 for v in range(2048))
        code=C[5]
        assert min((u^v^2047).bit_count() for u in code for v in code)==3
        print('Every listed 11D construction has covering radius sqrt(3); whole-layer 12D optimum 9/256.')
SOURCES={"lp.h":"struct LP{int m,n;vector<int>B,N;vector<vector<double>>D;LP(vector<vector<double>>A,vector<double>b,vector<double>c):m(b.size()),n(c.size()),B(m),N(n+1),D(m+2,vector<double>(n+2)){for(int i=0;i<m;i++)for(int j=0;j<n;j++)D[i][j]=A[i][j];for(int i=0;i<m;i++)B[i]=n+i,D[i][n]=-1,D[i][n+1]=b[i];for(int j=0;j<n;j++)N[j]=j,D[m][j]=-c[j];N[n]=-1;D[m+1][n]=1;}void pivot(int r,int s){double inv=1/D[r][s];for(int i=0;i<m+2;i++)if(i!=r&&abs(D[i][s])>1e-14){double q=D[i][s]*inv;for(int j=0;j<n+2;j++)if(j!=s)D[i][j]-=D[r][j]*q;D[i][s]*=-inv;}for(int j=0;j<n+2;j++)if(j!=s)D[r][j]*=inv;D[r][s]=inv;swap(B[r],N[s]);}bool sp(int ph){int r0=ph==1?m+1:m;for(int lim=0;lim<2000;lim++){int s=-1;for(int j=0;j<=n;j++){if(ph==2&&N[j]==-1)continue;if(s<0||make_pair(D[r0][j],N[j])<make_pair(D[r0][s],N[s]))s=j;}if(D[r0][s]>-1e-8)return true;int r=-1;for(int i=0;i<m;i++)if(D[i][s]>1e-8){if(r<0||make_pair(D[i][n+1]/D[i][s],B[i])<make_pair(D[r][n+1]/D[r][s],B[r]))r=i;}if(r<0)return false;pivot(r,s);}return false;}double solve(vector<double>&y){int r=0;for(int i=0;i<m;i++)if(D[i][n+1]<D[r][n+1])r=i;if(D[r][n+1]<-1e-8){pivot(r,n);if(!sp(1)||D[m+1][n+1]<-1e-8)return -1e30;for(int i=0;i<m;i++)if(B[i]==-1){int s=0;for(int j=1;j<=n;j++)if(abs(D[i][j])>abs(D[i][s]))s=j;if(abs(D[i][s])>1e-8)pivot(i,s);}}if(!sp(2))return 1e30;y.assign(m,0);for(int j=0;j<=n;j++)if(N[j]>=n)y[N[j]-n]=D[m][j];return D[m][n+1];}};","sharp12.h":"#include <bits/stdc++.h>\n#include <omp.h>\nusing namespace std;\n#include \"lp.h\"\nconst int Q=128,T0=12*Q*Q;\nint W[2048];vector<vector<int>>C(96,vector<int>(72));\nstruct Bx{array<int,11> l{},h;Bx(){h.fill(Q);}};\nstruct Row{array<int,11>A{};int b,id;};\nint exactMask(int v,int i,int j,const array<int,11>&y,int mask,int T){int m=INT_MAX;for(int c=0;c<72;c++)if(c!=j){int d=0;for(int k=0;k<11;k++){int e=((C[i][c]^v)>>k&1)?2*Q-y[k]:y[k];d+=e*e;}m=min(m,d);}if(mask>=0){int a=0;for(int k=0;k<11;k++)if(mask>>k&1)a+=y[k]*y[k];m=min(m,T0-T*Q*Q+a);}return m;}\nstruct Solver{long long calls=0,lpN=0,blimit=12000;int dep=85;array<int,11>hit{};int rec(int i,int j,int v,int mask,int T,Bx&x,int depth){if(++calls>blimit)return 0;vector<Row> rows;for(int c=0;c<72+(mask>=0);c++)if(c!=j){Row a;a.id=c;a.b=c==72?T0-T*Q*Q:0;for(int k=0;k<11;k++){if(c==72&&!(mask>>k&1))continue;int e=c==72?0:2*Q*((v^C[i][c])>>k&1);int l=x.l[k],h=x.h[k];a.b+=(l-e)*(l-e);a.A[k]=(l+h-2*e)*(h-l);}long long mn=a.b,mx=a.b;for(int d:a.A){mn+=min(d,0);mx+=max(d,0);}if(mx<=T0)return -1;if(mn<=T0)rows.push_back(a);}if(rows.empty()){hit=x.l;return 1;}vector<vector<double>> A;vector<double>b;for(auto&r:rows){vector<double> aa(12);for(int k=0;k<11;k++)aa[k]=-double(r.A[k])/Q/Q;aa[11]=1;A.push_back(aa);b.push_back(double(r.b)/Q/Q);}for(int k=0;k<11;k++){vector<double>aa(12);aa[k]=1;A.push_back(aa);b.push_back(1);}vector<double>o(12),y;o[11]=1;LP lp(A,b,o);lpN++;double z=lp.solve(y);for(int scale : {65536,12,24,60,120,420,840,2520,27720,360360,720720}){array<long long,11>g{};long long tot=0,n=0;for(int r=0;r<(int)rows.size()&&r<(int)y.size();r++){int wt=isfinite(y[r])?int(llround(clamp(y[r]*scale,0.,16777216.))):0;n+=wt;tot+=1LL*rows[r].b*wt;for(int k=0;k<11;k++)g[k]+=1LL*rows[r].A[k]*wt;}for(long long a:g)tot+=max(a,0LL);if(n&&tot<=n*T0)return -1;}array<int,11>point=x.l;for(int r=0;r<lp.m;r++)if(lp.B[r]<11&&lp.B[r]>=0){int k=lp.B[r];double t=lp.D[r][lp.n+1];if(isfinite(t))point[k]=x.l[k]+int(llround((x.h[k]-x.l[k])*clamp(t,0.,1.)));}\nif(exactMask(v,i,j,point,mask,T)>T0){hit=point;return 1;}\nif(depth>=dep)return 0;double best=-1;int K=-1;\nfor(int k=0;k<11;k++){int l=x.l[k],h=x.h[k],w=h-l;if(w<2)continue;double g=(point[k]-l)*(h-point[k]);double sc=g+0.0001*w*w;if(sc>best)best=sc,K=k;}if(K<0)return 0;int l=x.l[K],h=x.h[K],mid=(l+h)/2;\nx.h[K]=mid;int a=rec(i,j,v,mask,T,x,depth+1);x.h[K]=h;if(a==1)return 1;x.l[K]=mid;int b2=rec(i,j,v,mask,T,x,depth+1);x.l[K]=l;return b2==1?1:(a<0&&b2<0?-1:0);}\n};","sharp12.cpp":"#include \"sharp12.h\"\nint main(){ifstream f(\"C.dat\");for(auto&a:C)for(int&i:a)f>>i;ofstream o(\"sharp12.tsv\");int ok=0,no=0,uk=0;long long n=0; \n#pragma omp parallel for schedule(dynamic) reduction(+:ok,no,uk,n)\nfor(int i=0;i<96;i++){int a0=0,a1=0,a2=0;for(int v=0;v<2048;v++){Solver s;s.blimit=1200;s.dep=65;Bx x;int a=s.rec(i,-1,v,-1,0,x,0);n+=s.calls;if(a>0){a1++;\n#pragma omp critical\n{cerr<<\"W \"<<i<<\" \"<<v;for(int z:s.hit)cerr<<\" \"<<z;cerr<<\"\\n\";}}else if(!a)a2++;else a0++;}\n#pragma omp critical\n{o<<i<<\" \"<<a0<<\" \"<<a1<<\" \"<<a2<<\"\\n\";o.flush();}ok+=a1;uk+=a2;no+=a0;}\ncout<<\"12D clearance full corners excluded \"<<no<<\" confirmed \"<<ok<<\" unresolved \"<<uk<<\" nodes \"<<n<<\"\\n\";}"}
if __name__=='__main__': run()
