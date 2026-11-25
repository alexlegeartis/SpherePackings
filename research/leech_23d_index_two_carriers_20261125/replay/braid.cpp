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
using U=unsigned long long;struct V{int x,y,c;};using Bits=array<U,4>;
int ag[256][256],sg[256][256],ww[256],ml[256][17];Bits maskg[70][256];
struct Cn{int i,j,l;};vector<Cn> CN[1600];int g1,g2,K,rad,thr,goal=6;ll cap=500000;
vector<Bits>T;vector<V>B;vector<vector<U>> AD;vector<int> cur,found;ll calls,cases,edgecount;bool capped;
int shift[17][17];
void group(){for(int a=0;a<256;a++){ww[a]=wt[a/128][(a>>6)&1][a&63];for(int b=0;b<256;b++){ag[a][b]=((a^b)&192)|aa[a&63][b&63];sg[a][b]=((a^b)&192)|ss[a&63][b&63];}}
 for(int a=0;a<256;a++){for(int n=-8;n<=8;n++){int h=0;for(int k=0;k<abs(n);k++)h=n<0?sg[h][a]:ag[h][a];ml[a][n+8]=h;}for(int s=0;s<70;s++)for(int b=0;b<256;b++)if(ww[ag[a][b]]>=s)maskg[s][a][b/64]|=1ULL<<(b%64);}
 for(int x=0;x<K;x++)for(int y=0;y<K;y++)for(int i=-7;i<=7;i++)for(int j=-4;j<=4;j++){int X=x+K*i,Y=y+K*j,sup=0;for(int dx:{-rad,rad})for(int dy:{-rad,rad}){int a=X+dx,b=Y+dy;sup=max(sup,4*a*a+4*a*b+37*b*b);}int l=(thr*K*K-sup+K*K-1)/(K*K);if(l>0)CN[x*K+y].push_back({i,j,l});}
 }
bool single(){for(int i=-7;i<=7;i++)for(int j=-4;j<=4;j++)if(i||j){int s=4*i*i+4*i*j+37*j*j;if(s>=thr)continue;if(s+ww[ag[ml[g1][i+8]][ml[g2][j+8]]]<thr)return false;}return true;}
int link(V a,V b){int x=a.x-b.x,y=a.y-b.y,c=sg[a.c][b.c];if(x<0){x+=K;c=ag[c][g1];}if(y<0){y+=K;c=ag[c][g2];}return (T[x*K+y][c/64]>>(c%64))&1;}
void expand(vector<U>P){
 if(cur.size()==goal){found=cur;return;}if(!found.empty()||capped)return;if(++calls>cap){capped=true;return;}int nw=P.size(),color=0;vector<pair<int,int>> ord;vector<U>Z=P;for(;;){bool more=false;for(U v:Z)more|=(v!=0);if(!more)break;color++;vector<U>X=Z;for(int w=0;w<nw;w++)while(X[w]){int j=__builtin_ctzll(X[w]),v=w*64+j;Z[w]&=~(1ULL<<j);X[w]&=~(1ULL<<j);ord.push_back({v,color});for(int h=w;h<nw;h++)X[h]&=~AD[v][h];}}
 for(int h=ord.size()-1;h>=0;h--){auto [v,co]=ord[h];if((int)cur.size()+co<goal)return;cur.push_back(v);vector<U>Q(nw);for(int i=0;i<nw;i++)Q[i]=P[i]&AD[v][i];expand(Q);cur.pop_back();if(!found.empty()||capped)return;P[v/64]&=~(1ULL<<(v%64));}}
void run(){T.assign(K*K,Bits{});B.clear();for(int i=-7;i<=7;i++)for(int j=-4;j<=4;j++)shift[i+8][j+8]=ag[ml[g1][i+8]][ml[g2][j+8]];
 for(int x=0;x<K;x++)for(int y=0;y<K;y++){auto &bits=T[x*K+y];bits={~0ULL,~0ULL,~0ULL,~0ULL};for(auto [i,j,l]:CN[x*K+y]){auto M=maskg[l][shift[i+8][j+8]];for(int h=0;h<4;h++)bits[h]&=M[h];}for(int w=0;w<4;w++){U u=bits[w];while(u){int j=__builtin_ctzll(u);u&=u-1;B.push_back({x,y,64*w+j});}}}
 int n=B.size(),nw=(n+63)/64;assert(!link({0,0,0},{0,0,0}));AD.assign(n,vector<U>(nw));edgecount=0;for(int i=0;i<n;i++)for(int j=0;j<i;j++){int b=link(B[i],B[j]);if(b){AD[i][j/64]|=1ULL<<(j%64);AD[j][i/64]|=1ULL<<(i%64);edgecount++;}}
 vector<U>P(nw,~0ULL);if(n%64)P.back()=(1ULL<<(n%64))-1;found.clear();cur.clear();calls=0;capped=false;expand(P);}
int main(int argc,char**argv){K=atoi(argv[1]);rad=atoi(argv[2]);thr=atoi(argv[3]);string mode=argv[4];if(argc>5)goal=atoi(argv[5])-1;qs[8]=1;qs[4]=-1;if(argc>6){memset(qs,0,sizeof qs);stringstream ins(argv[6]);string k;while(getline(ins,k,',')){int i=stoi(k);qs[abs(i)]=i>=0?1:-1;}}init();group();set<pair<int,int>> filt;if(mode=="filter"){int i,j;while(cin>>i>>j)filt.insert({i,j});}int viable=0,nfound=0,ncap=0;ll ncall=0;for(g1=0;g1<256;g1++)for(g2=0;g2<256;g2++){
 if(mode=="filter"&&!filt.count({g1,g2}))continue;if(mode=="base"&&((g1&63)!=4||(g2&63)!=33))continue;if(!single())continue;viable++;if(mode=="count")continue;
 run();ncall+=calls;ncap+=capped;if(capped||!found.empty())cout<<"OPEN "<<g1<<" "<<g2<<endl;if(!found.empty()){nfound++;if(nfound<=3){cout<<"HIT K="<<K<<" rad="<<rad<<" thr="<<thr<<" flux="<<g1<<","<<g2<<" sites=[0,0;0]";for(int i:found){auto b=B[i];cout<<" ["<<b.x<<","<<b.y<<";"<<b.c<<"]";}cout<<endl;}if(mode=="hit")return 0;}}
cout<<"ALL K "<<K<<" rad "<<rad<<" thr "<<thr<<" goal "<<goal+1<<" viable "<<viable<<" found "<<nfound<<" capped "<<ncap<<" calls "<<ncall<<endl;}
