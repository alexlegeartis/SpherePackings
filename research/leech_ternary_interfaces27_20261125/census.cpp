#include <bits/stdc++.h>
using namespace std;using I=uint32_t;using ll=long long;
int wt(I x){return __builtin_popcount(x);}const I M=(1u<<23)-1;
struct B {I a[23]{};int rank=0;void add(I x){for(int i=22;i>=0;i--)if(x>>i&1){if(a[i])x^=a[i];else{a[i]=x;rank++;break;}}}};
int po(int x,int n){int y=1;while(n--)y=y*x%23;return y;}
vector<I> G,E;set<I> look;
struct Perm{int p[23];I lo[4096],hi[2048];I slow(I x){I y=0;while(x){int i=__builtin_ctz(x);x&=x-1;y|=1u<<p[i];}return y;}void init(){set<int> xx(p,p+23);assert(xx.size()==23);for(int i=0;i<4096;i++)lo[i]=slow(i);for(int i=0;i<2048;i++)hi[i]=slow(i<<12);}I operator()(I x){return lo[x&4095]|hi[x>>12];}} g0,g1;
struct D{int cost,par,alt;};
D tab[2][2][3][2];
void setup(){for(I i=0;i<4096;i++){I x=0;for(int j=0;j<12;j++)if(i>>j&1)x^=2787u<<j;assert(x<=M);if(!(wt(x)&1))E.push_back(x);x|=(wt(x)&1)<<23;assert(!x||(wt(x)>=8&&wt(x)%4==0));G.push_back(x);look.insert(x&M);}assert(E.size()==2048&&look.size()==4096);
for(int i=0;i<23;i++)g0.p[i]=(i+1)%23;g0.init();for(I x:G)assert(look.count(g0(x&M)));
bool ok=false;for(int k:{1,5}){int ki=po(k,21);for(int j=0;j<23;j++){int x=j*k%23;int y=(4*po(x,4)-3*po(x,15)+23*4)%23;g1.p[j]=ki*y%23;}g1.init();ok=true;for(I x:G)if(!look.count(g1(x&M)))ok=false;if(ok)break;}assert(ok);
for(int marked=0;marked<2;marked++)for(int p=0;p<2;p++)for(int a=-1;a<=1;a++)for(int c=0;c<2;c++){
 int best[2]={INT_MAX,INT_MAX};for(int k=-30;k<=30;k++){int y=p+2*c+4*k;if(marked&&y%3)continue;int d=3*y-(marked?24:8)*a;int &b=best[(k+100)&1];b=min(b,d*d);}int e=best[1]<best[0];tab[marked][p][a+1][c]={best[e],e,best[!e]-best[e]};}}
int dist(I S,I neg,int am,int bound=0){I ms[]={neg,M^S,S^neg};int best=100000;for(int p=0;p<2;p++)for(I v:G){
 auto e=tab[1][p][am+1][v>>23];int sum=e.cost,alt=e.alt,par=p^e.par;for(int k=0;k<3;k++){int ones=wt(v&ms[k]),zeros=wt(ms[k])-ones;for(int j=0;j<2;j++){int n=j?ones:zeros;if(!n)continue;auto a=tab[0][p][k][j];sum+=n*a.cost;alt=min(alt,a.alt);par^=(n&1)&a.par;}}best=min(best,sum+(par?alt:0));if(bound&&best<bound)return best;
 }return best;}
int literal(I S,I neg,int am){int best=INT_MAX;for(int p:{0,1})for(I c:G){int dp[]={0,100000};for(int j=0;j<24;j++){int v[2]={100000,100000},a=j==23?3*am:((S>>j)&1)*(1-2*int((neg>>j)&1));for(int k=-6;k<=6;k++){int y=p+2*((c>>j)&1)+4*k;if(j==23&&y%3)continue;int d=3*y-8*a;for(int i=0;i<2;i++)v[(i+k+100)&1]=min(v[(i+k+100)&1],dp[i]+d*d);}copy(v,v+2,dp);}best=min(best,dp[p]);}return best;}

using P=array<int,23>;
P mul(P a,P b){P c;for(int j=0;j<23;j++)c[j]=a[b[j]];return c;}
P inv(P a){P c;for(int j=0;j<23;j++)c[a[j]]=j;return c;}
vector<Perm> stabil(int fixed){P id,p,q;for(int i=0;i<23;i++){id[i]=i;p[i]=g0.p[i];q[i]=g1.p[i];}vector<P> gen{p,q};for(int jj=22;jj>=fixed;jj--){vector<P> transversal;map<int,P> tr;transversal.push_back(id);tr[jj]=id;set<P> gs;for(unsigned i=0;i<transversal.size();i++){P t=transversal[i];for(auto &g:gen){P b=mul(g,t);int x=b[jj];if(!tr.count(x)){tr[x]=b;transversal.push_back(b);}else {P h=mul(inv(tr[x]),b);if(h!=id)gs.insert(h);}}}gen.assign(gs.begin(),gs.end());
 // A few verified elements suffice; a subgroup yields an exact, possibly finer partition.
 if(gen.size()>12){mt19937 r(847+jj);shuffle(gen.begin(),gen.end(),r);gen.resize(12);}}
 vector<Perm> ans;for(auto g:gen){Perm a;copy(g.begin(),g.end(),a.p);a.init();ans.push_back(a);}return ans;}
const int MM=(1<<22)-1;
D tt[2][4][3][3];
void setmark(){for(int p=0;p<2;p++)for(int c=0;c<4;c++)for(int a=-1;a<=1;a++)for(int t=-1;t<=1;t++){
 int best[2]={10000,10000};for(int i=-6;i<=6;i++)for(int j=-6;j<=6;j++){int y=p+2*(c&1)+4*i,z=p+2*(c>>1)+4*j;if((y+z)%3)continue;int h=3*y-8*a,k=3*z+8*a-24*t;best[(i+j+100)&1]=min(best[(i+j+100)&1],h*h+k*k);}int e=best[1]<best[0];tt[p][c][1+a][1+t]={best[e],e,best[!e]-best[e]};}}
int twodist(I s,I neg,int a,int t,int bound=0){I ms[3]={neg,MM^s,s^neg};int nc[3]={wt(ms[0]),wt(ms[1]),wt(ms[2])};int best=100000;
 for(int p=0;p<2;p++)for(I v:G){auto e=tt[p][v>>22][a+1][t+1];int sum=e.cost,alt=e.alt,par=p^e.par;for(int k=0;k<3;k++){int ones=wt(v&ms[k]);for(int j:{0,1}){int n=j?ones:nc[k]-ones;if(!n)continue;auto b=tab[0][p][k][j];sum+=n*b.cost;alt=min(alt,b.alt);par^=(n&1)&b.par;}}best=min(best,sum+(par?alt:0));if(bound&&best<bound)return best;}return best;}
int independent(I s,I neg,int a,int t){int best=100000;for(int p:{0,1})for(I c:G){int dp[2][3]={{0,100000,100000},{100000,100000,100000}};for(int j=0;j<24;j++){int v[2][3];for(auto &r:v)fill(r,r+3,100000);int b=j<22?((s>>j)&1)*(1-2*int(neg>>j&1)):j==22?a:-a+3*t;for(int k=-4;k<=4;k++){int y=p+2*(c>>j&1)+4*k,d=3*y-8*b;for(int l:{0,1})for(int m:{0,1,2})v[(l+k+100)&1][(m+(j<22?0:y)+99)%3]=min(v[(l+k+100)&1][(m+(j<22?0:y)+99)%3],dp[l][m]+d*d);}for(int j=0;j<2;j++)copy(v[j],v[j]+3,dp[j]);}best=min(best,dp[p][0]);}return best;}

// For each off-mark sign orbit, test all three lifts in one fibre.
// Values are nine times squared Euclidean depth. Only fibres with
// some depth strictly exceeding 192/9 need to be counted.
void run(int m){
 int bits=24-m,mask=(1<<bits)-1;
 vector<Perm> gs=m==1?vector<Perm>{g0,g1}:stabil(22);
 vector<char>seen(mask+1);vector<I>q;map<array<int,3>,ll> hist;
 ll total=0,n_orbits=0;
 for(int s=0;s<=mask;s++){
  if(seen[s])continue;seen[s]=1;n_orbits++;q.clear();q.push_back(s);
  for(unsigned j=0;j<q.size();j++)for(auto &g:gs){I y=g(q[j]);if(!seen[y]){seen[y]=1;q.push_back(y);}}
  B b;for(I v:G)if(!(v>>bits))b.add(s&v);
  I free=0;for(int j=0;j<bits;j++)if((s>>j&1)&&!b.a[j])free|=1u<<j;
  for(I neg=free;;neg=(neg-1)&free){
   for(int a=(m==1?0:-1);a<=(m==1?0:1);a++){
    array<int,3>d;auto calc=[&](int t,int bd){return m==1?dist(s,neg,t,bd):twodist(s,neg,a,t,bd);};
    for(int t=-1;t<=1;t++)d[t+1]=calc(t,193);
    ll mult=ll(q.size())<<b.rank;total+=3*mult;
    if(*max_element(d.begin(),d.end())>=193){
     for(int t=-1;t<=1;t++)d[t+1]=calc(t,0);
     sort(d.begin(),d.end());hist[d]+=mult;
    }
   }if(!neg)break;
  }
 }
 assert(total==282429536481LL);assert(hist.size()==2);
 array<int,3> apex{0,288,288},near{32,224,224};
 assert(hist[apex]==1&&hist[near]==(m==1?69876:93150));
 cout<<"mark "<<m<<": "<<n_orbits<<" support orbits; "<<total<<" ternary phases; upper fibres: ";
 for(auto [d,n]:hist)cout<<"["<<d[0]<<","<<d[1]<<","<<d[2]<<"] x "<<n<<"  ";cout<<"\n";
 mt19937 rng(5252+m);for(int i=0;i<12;i++){
  I s=rng()&mask,neg=s&rng();int a=int(rng()%3)-1,t=int(rng()%3)-1;
  if(m==1)assert(dist(s,neg,t)==literal(s,neg,t));
  else assert(twodist(s,neg,a,t)==independent(s,neg,a,t));
 }
}
int main(){setup();setmark();run(1);run(2);cout<<"Independent parity-distance replays passed.\n";}
