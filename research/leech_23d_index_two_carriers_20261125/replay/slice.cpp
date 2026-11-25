#include <bits/stdc++.h>
using namespace std; using ll=long long;using wide=__int128;
int mod(int x,int n){return(x%n+n)%n;}
int O[8]={0,1,5,6,7,9,11,23},inc[24],qs[24],ax[64][3],aa[64][64],ss[64][64],wt[2][2][64],gol[4096];
set<pair<int,int>> pts[2][2][64];
int ci(int x,int y,int z){int a[3]={x,y,z},p=mod(x,2),v=32*p,k=0;for(int j=0;j<3;j++){int h=j==2?8:4,u=mod(a[j],h);v|=(u/2)<<j;k+=(a[j]-u)/h;}return v+16*mod(k,2);}
int mul(int x,int n){int y=0;for(n=mod(n,16);n--;y=aa[x][y]);return y;}
int wk[64][25];
void init(){
for(auto &v:wk)fill(v,v+25,99999);for(int i:O)inc[i]=1;int best[64];fill(best,best+64,999);for(int x=-7;x<=7;x++)for(int y=-7;y<=7;y++)for(int z=-12;z<=12;z++)if(!mod(x-y,2)&&!mod(y-z,2)){int c=ci(x,y,z),d=2*x*x+2*y*y+z*z;if(d<best[c]){best[c]=d;ax[c][0]=x;ax[c][1]=y;ax[c][2]=z;}}
for(int i=0;i<64;i++)for(int j=0;j<64;j++){auto a=ax[i],b=ax[j];aa[i][j]=ci(a[0]+b[0],a[1]+b[1],a[2]+b[2]);ss[i][j]=ci(a[0]-b[0],a[1]-b[1],a[2]-b[2]);}
for(int i=0;i<4096;i++){int w=0,g=(1<<11)|(1<<9)|(1<<7)|(1<<6)|(1<<5)|(1<<1)|1;for(int j=0;j<12;j++)if(i>>j&1)w^=g<<j;gol[i]=w|(__builtin_parity(w)<<23);}
for(int c=0;c<64;c++){int x=ax[c][0],y=ax[c][1],w=ax[c][2],m=mod(x,2);for(int C:gol){auto R=[&](int i){return m+2*((C>>i)&1);};if(mod(x-R(2),4)||mod(y-R(3),4))continue;int qx=0;for(int j=0;j<24;j++)qx+=qs[j]*R(j);assert(qx%2==0);int r=mod(qx/2,2);for(int t=0;t<2;t++){int z[8],sum=0;for(int j=0;j<8;j++){int sign=j==7?-1:1;z[j]=mod(sign*(R(O[j])+2*t*qs[O[j]])+2,4)-2;sum+=z[j];}
if(mod(2*w-sum,4))continue;int dir=2*w>=sum?4:-4;while(sum!=2*w){int b=0;for(int j=1;j<8;j++)if(dir*z[j]<dir*z[b])b=j;z[b]+=dir;sum+=dir;}
int oc=0,os=0;for(int j=0;j<8;j++){oc+=z[j]*z[j];os+=(j==7?-1:1)*z[j]-2*t*qs[O[j]];}
if(t==0){int ec0=0,ss0=0,fi=99;for(int j=0;j<24;j++)if(!inc[j]&&j!=2&&j!=3&&j!=4&&j!=8){int h=mod(R(j)+2,4)-2;ec0+=h*h;ss0+=h;fi=min(fi,16-8*abs(h));}for(int K=-12;K<=12;K++)if(mod(R(8)-R(4)-2*K,4)==0){int z4=mod(R(4)+K+2,4)-2-K,z8=z4+2*K;int cong=mod(4*m-(os+x+y+ss0+z4+z8),8);assert(cong==0||cong==4);int total=2*(oc+ec0+z4*z4+z8*z8+(cong?fi:0))-w*w;wk[c][K+12]=min(wk[c][K+12],total);}}
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
int main(){qs[8]=1;qs[4]=-1;init();for(int c=0;c<64;c++){for(int r=0;r<2;r++)for(int t=0;t<2;t++){int b=999;for(int k=-12;k<=12;k++)if(mod(k,2)==r)b=min(b,wk[c][k+12]+16*k*t+16*t);assert(b==wt[t][r][c]);}for(auto v:wk[c])cout<<v<<" ";cout<<endl;}}