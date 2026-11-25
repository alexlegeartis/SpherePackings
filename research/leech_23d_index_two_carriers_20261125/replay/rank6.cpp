#include <bits/stdc++.h>
using namespace std; using ll=long long;using wide=__int128;
int mod(int x,int n){return(x%n+n)%n;}
int O[8]={0,1,5,6,7,9,11,23},inc[24],qs[24],ax[64][3],aa[64][64],ss[64][64],wt[2][2][64],gol[4096];
set<pair<int,int>> pts[2][2][64];
int ci(int x,int y,int z){int a[3]={x,y,z},p=mod(x,2),v=32*p,k=0;for(int j=0;j<3;j++){int h=j==2?8:4,u=mod(a[j],h);v|=(u/2)<<j;k+=(a[j]-u)/h;}return v+16*mod(k,2);}
int mul(int x,int n){int y=0;for(n=mod(n,16);n--;y=aa[x][y]);return y;}
void init(){
for(int i:O)inc[i]=1;int best[64];fill(best,best+64,999);for(int x=-7;x<=7;x++)for(int y=-7;y<=7;y++)for(int z=-12;z<=12;z++)if(!mod(x-y,2)&&!mod(y-z,2)){int c=ci(x,y,z),d=2*x*x+2*y*y+z*z;if(d<best[c]){best[c]=d;ax[c][0]=x;ax[c][1]=y;ax[c][2]=z;}}
for(int i=0;i<64;i++)for(int j=0;j<64;j++){auto a=ax[i],b=ax[j];aa[i][j]=ci(a[0]+b[0],a[1]+b[1],a[2]+b[2]);ss[i][j]=ci(a[0]-b[0],a[1]-b[1],a[2]-b[2]);}
for(int i=0;i<4096;i++){int w=0,g=(1<<11)|(1<<9)|(1<<7)|(1<<6)|(1<<5)|(1<<1)|1;for(int j=0;j<12;j++)if(i>>j&1)w^=g<<j;gol[i]=w|(__builtin_parity(w)<<23);}
for(int c=0;c<64;c++){int x=ax[c][0],y=ax[c][1],w=ax[c][2],m=mod(x,2);for(int C:gol){auto R=[&](int i){return m+2*((C>>i)&1);};if(mod(x-R(2),4)||mod(y-R(3),4))continue;int qx=0;for(int j=0;j<24;j++)qx+=qs[j]*R(j);assert(qx%2==0);int r=mod(qx/2,2);for(int t=0;t<2;t++){int z[8],sum=0;for(int j=0;j<8;j++){int sign=j==7?-1:1;z[j]=mod(sign*(R(O[j])+2*t*qs[O[j]])+2,4)-2;sum+=z[j];}
if(mod(2*w-sum,4))continue;int dir=2*w>=sum?4:-4;while(sum!=2*w){int b=0;for(int j=1;j<8;j++)if(dir*z[j]<dir*z[b])b=j;z[b]+=dir;sum+=dir;}
int oc=0,os=0;for(int j=0;j<8;j++){oc+=z[j]*z[j];os+=(j==7?-1:1)*z[j]-2*t*qs[O[j]];}
int ec=0,es=0,flip=99;for(int j=0;j<24;j++)if(!inc[j]&&j!=2&&j!=3){int shift=2*t*qs[j],h=mod(R(j)+shift+2,4)-2;ec+=h*h;es+=h-shift;flip=min(flip,16-8*abs(h));}
int congr=mod(4*m-(es+os+x+y),8);assert(congr==0||congr==4);if(congr==4)ec+=flip;oc=2*oc-w*w;ec*=2;assert(oc>=0);pts[t][r][c].insert({oc,ec});}}
for(int r=0;r<2;r++)for(int t=0;t<2;t++){int mi=99999;for(auto [o,e]:pts[t][r][c])mi=min(mi,o+e);if(mi==99999)throw runtime_error("not index two");wt[t][r][c]=mi;}}
for(int t=0;t<2;t++){map<int,int> h;for(int c=0;c<64;c++)for(int r=0;r<2;r++){h[wt[t][r][c]]++;assert(wt[t][r][c]==wt[t][r][ss[0][c]]);}cerr<<"table t"<<t<<" x16:";for(auto [k,n]:h)cerr<<" "<<k<<":"<<n;cerr<<endl;}
for(int c=1;c<64;c++)assert(min(wt[0][0][c],wt[0][1][c])==64-best[c]);assert(wt[0][0][0]==0&&wt[0][1][0]>=64);
}
struct Site{int u,v,c,r;};struct Aff{ll l,m,b;};struct P{ll x,y,d;};struct Ans{ll n=-1,z=1;P p;int c,r;};
vector<Aff> powers(const vector<Site>&sites,int c,int r,int fa,int fb,int t=1){
vector<Aff>S;auto dom=[](Aff a,Aff b){return a.b-b.b+3*max(0LL,a.l-b.l)+3*max(0LL,a.m-b.m)<=0;};
for(auto h:sites)for(int i=-8;i<=8;i++)for(int j=-8;j<=8;j++){int X=3*i+h.u,Y=3*j+h.v;int cc=ss[c][aa[h.c][aa[mul(4,i)][mul(33,j)]]],rr=r^h.r^(fa&mod(i,2))^(fb&mod(j,2)),w=wt[t][rr][cc];
Aff a={-8LL*X-4LL*Y,-4LL*X-74LL*Y,4LL*X*X+4LL*X*Y+37LL*Y*Y+9LL*w};bool skip=false;for(auto b:S)if(dom(b,a)){skip=true;break;}if(skip)continue;S.erase(remove_if(S.begin(),S.end(),[&](auto b){return dom(a,b);}),S.end());S.push_back(a);}return S;}
P intersect(Aff a,Aff b){ll d=a.l*b.m-a.m*b.l;if(!d)return {0,0,0};ll x=a.m*b.b-a.b*b.m,y=a.b*b.l-a.l*b.b;if(d<0){x=-x;y=-y;d=-d;}if(x<0||y<0||x>3*d||y>3*d)return {0,0,0};ll g=gcd(gcd(x,y),d);return{x/g,y/g,d/g};}
Ans bound(const vector<Site>&sites,int fa,int fb,int t=1){Ans a;for(int c=0;c<64;c++)for(int r=0;r<2;r++){auto S=powers(sites,c,r,fa,fb,t);auto add=[&](P p){if(!p.d)return;ll v=LLONG_MAX;for(auto h:S)v=min(v,h.l*p.x+h.m*p.y+h.b*p.d);ll n=4*p.x*p.x+4*p.x*p.y+37*p.y*p.y+p.d*v,z=144*p.d*p.d;if((wide)n*a.z>(wide)a.n*z){ll g=gcd(n,z);a={n/g,z/g,p,c,r};}};add({0,0,1});add({0,3,1});add({3,0,1});add({3,3,1});auto diff=[](Aff a,Aff b){return Aff{a.l-b.l,a.m-b.m,a.b-b.b};};for(int i=0;i<(int)S.size();i++)for(int j=0;j<i;j++){auto line=diff(S[i],S[j]);for(auto edge:vector<Aff>{{1,0,0},{1,0,-3},{0,1,0},{0,1,-3}})add(intersect(line,edge));for(int k=0;k<j;k++)add(intersect(line,diff(S[i],S[k])));}}return a;}
void report(const string&nm,int f,Ans a){cout<<nm<<" flux"<<f<<" "<<a.n<<"/"<<a.z<<" = "<<setprecision(14)<<1.*a.n/a.z<<" label "<<a.c<<","<<a.r<<" transverse (x3,y3)=("<<a.p.x<<","<<a.p.y<<")/"<<a.p.d<<endl;}
struct Mom{ll n=0;ll a[24]={},m[24][24]={};void put(int *h){n++;for(int i=0;i<24;i++){a[i]+=h[i];for(int j=0;j<24;j++)m[i][j]+=h[i]*h[j];}}};map<tuple<int,int,int>,Mom> mom;
int main(){qs[8]=1;qs[4]=-1;init();struct CG{int i,j,u,v,o,e,c;};vector<CG> cg={{0,0,-1,0,12,48},{0,1,-2,0,23,28},{0,1,-1,0,15,44},{0,1,0,0,15,44},{0,2,-1,-1,15,44},{0,2,-1,-1,31,28},{0,2,1,-1,7,44},{1,2,-1,0,31,28},{1,2,0,0,31,28}};
int sitec[3]={0,32,5};for(auto &v:cg)v.c=aa[ss[sitec[v.j]][sitec[v.i]]][aa[mul(4,v.u)][mul(33,v.v)]];
vector<array<short,28>> rows;auto visit=[&](array<int,24>&z){int sw=0;for(int j=0;j<8;j++)sw+=(j==7?-1:1)*z[O[j]];assert(sw%2==0);int w=sw/2,c=ci(z[2],z[3],w),h[24];for(int i=0;i<24;i++)h[i]=4*z[i];h[2]=h[3]=0;for(int j=0;j<8;j++)h[O[j]]-=(j==7?-1:1)*w;int o=0,e=0;for(int i=0;i<24;i++)if(inc[i])o+=h[i]*h[i];else e+=h[i]*h[i];assert(o%8==0&&e%8==0);o/=8;e/=8;for(auto g:cg)if(g.c==c&&g.o==o&&g.e==e){array<short,28> v;int k=0;for(int a=0;a<24;a++)if(a!=2&&a!=3&&a!=23)v[k++]=h[a];v[21]=3*g.u+g.j-g.i;v[22]=3*g.v+g.j-g.i;v[23]=g.i;v[24]=g.j;v[25]=mod((z[8]-z[4])/2,2);v[26]=g.u;v[27]=g.v;rows.push_back(v);}};
array<int,24>z{};for(int i=0;i<24;i++)for(int j=0;j<i;j++)for(int a:{-4,4})for(int b:{-4,4}){z[i]=a;z[j]=b;visit(z);z[i]=z[j]=0;}
for(int C:gol)if(__builtin_popcount(C)==8){int ix[8],n=0;for(int j=0;j<24;j++)if(C>>j&1)ix[n++]=j;for(int v=0;v<256;v++)if(!__builtin_parity(v)){for(int j=0;j<8;j++)z[ix[j]]=(v>>j&1)?-2:2;visit(z);}for(int j:ix)z[j]=0;}
for(int C:gol){for(int j=0;j<24;j++)z[j]=(C>>j&1)?-1:1;for(int j=0;j<24;j++){z[j]*=-3;visit(z);z[j]/=-3;}}

mt19937 rng(42);shuffle(rows.begin(),rows.end(),rng);for(int f=0;f<4;f++){const int D=69,p=1009;vector<vector<int>>basis(D);int rank=0,tried=0;int inv[p];inv[1]=1;for(int k=2;k<p;k++)inv[k]=p-(p/k)*inv[p%k]%p;
for(auto v:rows){vector<int>a(D);int k=0,sg=(v[25]^((f&1)&mod(v[26],2))^((f>>1)&mod(v[27],2)))?-1:1;for(int b=0;b<3;b++)for(int i=0;i<23;i++)a[k++]=mod(v[i]*(sg*(v[24]==b)-(v[23]==b)),p);assert(k==D);tried++;
for(int i=0;i<D;i++)if(a[i]){if(basis[i].empty()){int x=inv[a[i]];for(int j=i;j<D;j++)a[j]=a[j]*x%p;basis[i]=a;rank++;break;}int q=a[i];for(int j=i;j<D;j++)a[j]=(a[j]+(p-1)*p-q*basis[i][j])%p;}if(rank==D)break;}
cout<<"flux "<<f<<" total "<<rows.size()<<" tried "<<tried<<" rank_mod1009 "<<rank<<"/"<<D<<endl;}
}
