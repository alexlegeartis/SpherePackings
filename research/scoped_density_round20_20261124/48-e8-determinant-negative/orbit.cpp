#include<bits/stdc++.h>
using namespace std;using L=int64_t;using U=unsigned;
int wt(U x){return popcount(x);}
int q8(int x){return wt(unsigned(x&(x>>4)))%2;}
int B8(int a,int b){return q8(a^b)^q8(a)^q8(b);}
struct Sub{U mask;int d;vector<int> rep;array<int,32> f;};
vector<Sub>D; vector<pair<int,int>> P[6];
vector<int> qp[6];vector<uint8_t> pol[6];
vector<int16_t> quad[6];
struct Rel{int p,s;vector<int> f;};
vector<Rel> rel;
vector<pair<int,int>> C;
int q(int r,int m,int x){
 int z=0;for(int i=0;i<(int)P[r].size();i++){auto [j,k]=P[r][i];z^=(m>>i&1)&&(x>>j&1)&&(x>>k&1);}return z;
}
int B(int r,int m,int x,int y){return q(r,m,x)^q(r,m,y)^q(r,m,x^y);}
int encode(int a,int b,int m,const vector<int>& v){
 int z=0,i=0;for(auto[j,k]:P[a])z|=(j==k?q(b,m,v[j]):B(b,m,v[j],v[k]))<<i++;
 return z;
}
void init(){
 for(int r=0;r<=5;r++){
  for(int j=0;j<r;j++)for(int k=0;k<=j;k++)P[r].push_back({k,j});
  int n=1<<r,m0=1<<P[r].size();quad[r].resize(m0*n*n);
  for(int m=0;m<m0;m++){
   for(int a=0;a<n;a++){
    int y=0,i=0;
    for(auto[j,k]:P[r]){if(m>>i&1){y^=(a>>k&1)<<j;if(j!=k)y^=(a>>j&1)<<k;}++i;}
    for(int b=0;b<n;b++)quad[r][(m*n+y)*n+b]+=wt(a&b)%2?-n:n;
   }
  }
 }
 set<U> z{1};
 for(int t=0;t<5;t++){auto o=z;for(U m:o)for(int x=1;x<32;x++)if(!(m>>x&1)){
  U k=m;for(int i=0;i<32;i++)if(m>>i&1)k|=1u<<(i^x);z.insert(k);
 }}
 assert(z.size()==374);
 for(U m:z){
  Sub s={m,5-countr_zero(U(wt(m))),{},{}};
  U u=m;for(int x=0;x<5;x++)if(!(u>>(1<<x)&1)){
   s.rep.push_back(1<<x);U w=u;for(int i=0;i<32;i++)if(u>>i&1)w|=1u<<(i^(1<<x));u=w;
  }
  assert(u==4294967295u&&s.rep.size()==s.d);
  for(int a=0;a<(1<<s.d);a++){
   int t=0;for(int i=0;i<s.d;i++)if(a>>i&1)t^=s.rep[i];
   for(int x=0;x<32;x++)if(m>>x&1)s.f[x^t]=a;
  }
  D.push_back(s);
 }
 sort(D.begin(),D.end(),[](auto&a,auto&b){return tie(a.d,a.mask)<tie(b.d,b.mask);});
 for(int i=0;i<374;i++)for(int j=0;j<i;j++)if((D[i].mask&D[j].mask)==D[i].mask){
  vector<int> v;for(int a:D[i].rep)v.push_back(D[j].f[a]);
  Rel u{i,j,{}};
  for(int b=0;b<(1<<P[D[j].d].size());b++)u.f.push_back(encode(D[i].d,D[j].d,b,v));
  rel.push_back(move(u));
 }
}
vector<vector<L>> calc(array<int,5> X){
 vector<vector<L>> A(374);
 for(int i=0;i<374;i++){
  auto &s=D[i];int r=s.d,n=1<<r,T=1<<P[r].size();
  int v[5]={};
  for(int j=0;j<5;j++)for(int k=0;k<r;k++)if(s.f[1<<j]>>k&1)v[k]^=X[j];
  int a[4]={},b[4]={};for(int l=0;l<4;l++)for(int j=0;j<r;j++)a[l]|=(v[j]>>(l+4)&1)<<j,b[l]|=(v[j]>>l&1)<<j;
  auto &av=A[i];av.resize(T);
  for(int m=0;m<T;m++){L z=1;for(int l=0;l<4;l++)z*=quad[r][(m*n+a[l])*n+b[l]];av[m]=z;}
  for(int h=1;h<T;h<<=1)for(int k=0;k<T;k+=2*h)for(int l=0;l<h;l++){L u=av[k+l],v=av[k+l+h];av[k+l]=u+v;av[k+l+h]=u-v;}
  for(auto &w:av){assert(w%T==0);w/=T;}
 }
 for(auto &e:rel)for(int j=0;j<(int)e.f.size();j++)A[e.p][e.f[j]]-=A[e.s][j];
 return A;
}
int gate=8, heavy=5;bool neg=false; int bad=2, one=0;
int R[256][8];
array<int,5> rep(const Sub&,int);

using Vec=array<int,8>;
vector<Vec> local[256];
array<int,8> basis{};
int pairv(const Vec &v,const Vec &w){
 int x=0;for(int j=0;j<8;j++)x+=v[j]*w[j];assert(x%4==0);return (x/4)&1;
}
void localinit(){
 int E[8]={},bo=0;
 for(int u=1;u<256;u++)if(q8(u)){
  int v=u;
  while(v){int j=31-countl_zero(unsigned(v));if(E[j])v^=E[j];else {E[j]=v;basis[bo++]=u;break;}}
 }
 assert(bo==8);
 Vec br[8];
 for(int i=0;i<8;i++)copy(R[basis[i]],R[basis[i]]+8,br[i].begin());
 int bk[256];fill_n(bk,256,-1);
 for(int y=0;y<256;y++){
  int r=0;for(int j=0;j<8;j++)r|=B8(y,basis[j])<<j;assert(bk[r]==-1);bk[r]=y;
 }
 Vec v{};
 auto visit=[&](auto &&f,int i,int rem,int parity,int s)->void{
   if(i==8){
    if((rem!=0&&rem!=8)||s%4)return;
    int y=0;for(int j=0;j<8;j++)y|=pairv(v,br[j])<<j;
    int x=bk[y];
    assert((rem!=0)==q8(x)&&x!=0);
    local[x].push_back(v);
    return;
   }
   for(int z=-4+parity;z<=4;z+=2)if(z*z<=rem){v[i]=z;f(f,i+1,rem-z*z,parity,s+z);}
 };
 visit(visit,0,16,0,0);visit(visit,0,16,1,0);
 for(int i=1;i<256;i++)assert(local[i].size()==(q8(i)?2:16));
}
int absnorm(Vec v){int z=0;for(int i:v)z+=i*i;return z;}
bool shortlow(const Sub&s,int v){
 auto x=rep(s,v);int z=0;vector<int>a;
 for(int u:x){z^=u;if(u)a.push_back(u);}if(z)a.push_back(z);
 for(size_t f=0;f<a.size();f++){
  swap(a[f],a.back());int tail=8*(gate-1);
  for(size_t j=0;j+1<a.size();j++)tail-=q8(a[j])?8:16;
  Vec sum{};
  auto loop=[&](auto &&dfs,size_t p)->bool{
   if(p+1==a.size())return absnorm(sum)<=tail;
   for(const Vec&w:local[a[p]]){
    for(int i=0;i<8;i++)sum[i]+=w[i];
    if(dfs(dfs,p+1))return true;
    for(int i=0;i<8;i++)sum[i]-=w[i];
   }
   return false;
  };
  if(loop(loop,0))return true;
  swap(a[f],a.back());
 }
 return false;
}

bool allowed0(const Sub&s,int v){
 int minraw=0;
 for(int j=0;j<6;j++) {
  int z=s.f[j<5?(1<<j):31];
  if(z)minraw+=q(s.d,v,z)?2:4;
 }
 return minraw&&(minraw>=2*gate||!shortlow(s,v));
}
vector<vector<uint8_t>> acache;
bool allowed(const Sub&s,int v){return acache[&s-D.data()][v];}
array<int,5> rep(const Sub&s,int v){
 int y[5]={};bool seen[256]={};seen[0]=1;
 for(int i=0;i<s.d;i++){
  int cur=0;
  for(int z=1;z<256;z++)if(!seen[z]&&q8(z)==q(s.d,v,1<<i)){
   bool ok=true;for(int j=0;j<i;j++)if(B8(z,y[j])!=B(s.d,v,1<<i,1<<j)){ok=false;break;}
   if(ok){cur=z;break;}
  }
  assert(cur);y[i]=cur;
  bool s2[256];copy(seen,seen+256,s2);for(int j=0;j<256;j++)if(s2[j])seen[j^cur]=1;
 }
 array<int,5> x={};
 for(int i=0;i<5;i++)for(int j=0;j<s.d;j++)if(s.f[1<<i]>>j&1)x[i]^=y[j];
 return x;
}

pair<int,int> label(array<int,5>x){
 int su[32]={};U mask=1;
 for(int i=1;i<32;i++){su[i]=su[i&(i-1)]^x[countr_zero(unsigned(i))];if(!su[i])mask|=1u<<i;}
 int i=0;while(D[i].mask!=mask)++i;
 int k=0,t=0;for(auto[a,b]:P[D[i].d]){
  int u=su[D[i].rep[a]],w=su[D[i].rep[b]];
  t|=(a==b?q8(u):B8(u,w))<<k++;
 }
 return {i,t};
}
vector<vector<pair<int,int>>> orb(const vector<vector<L>>&N,bool dual){
 set<pair<int,int>> seen;vector<vector<pair<int,int>>>out;
 for(int i=0;i<374;i++)for(int j=0;j<(int)N[i].size();j++){
  pair<int,int> p{i,j};
  if(!N[i][j]||(!dual&&!allowed(D[i],j))||seen.count(p))continue;
  vector<pair<int,int>> Q{p}; seen.insert(p);
  for(size_t q=0;q<Q.size();q++){
   auto [k,m]=Q[q];auto wr=rep(D[k],m);
   int a5=dual?0:wr[0]^wr[1]^wr[2]^wr[3]^wr[4];
   for(int t=0;t<heavy-bad;t++){
    array<int,6> x={wr[0],wr[1],wr[2],wr[3],wr[4],0}; x[heavy]=a5;
    swap(x[t],x[t+1]);
    if(dual)for(int j=0;j<heavy;j++)x[j]^=x[heavy]; if(heavy<5)x[heavy]=0;
    auto v=label({x[0],x[1],x[2],x[3],x[4]});
    assert(N[v.first][v.second]>0&&(dual||allowed(D[v.first],v.second)));
    if(seen.insert(v).second)Q.push_back(v);
   }
  }
  out.push_back(move(Q));
 }
 return out;
}

vector<array<int,4>> many(const Sub&s,int v){
 auto x=rep(s,v);
 int n=heavy+1,o=n-bad;
 vector<int> u(x.begin(),x.begin()+heavy);
 u.push_back(accumulate(u.begin(),u.end(),0,bit_xor<int>()));
 int low[6]={};for(int i=0;i<n;i++)if(u[i])low[i]=q8(u[i])?1:2;
 int score[4]={};
 auto kind=[&](int i){return i<o?0:bad==3?(n-i):bad==2?(i==n-2?2:1):one;};
 for(int i=0;i<n;i++)score[kind(i)]+=low[i];
 vector<array<int,4>> out;
 auto insert=[&](array<int,4> p){
  if(p[0]>=8 && (!neg|| 4*p[0]-p[1]>=32))return;
  auto dominates=[&](auto a,auto b){return a[0]<=b[0]&&(neg?a[1]>=b[1]:a[1]<=b[1])&&a[2]<=b[2]&&a[3]<=b[3];};
  for(auto a:out)if(dominates(a,p))return;
  erase_if(out,[&](auto a){return dominates(p,a);});out.push_back(p);
 };
 for(int f=0;f<n;f++){
  int k=kind(f);
  // In the negative chart only A's largest score can help. The
  // least-score cases from other forced slots are still sound.
  int best=(k==1&&neg)?-1:INT_MAX;
  vector<int> p;for(int j=0;j<n;j++)if(j!=f&&u[j])p.push_back(j);
  Vec sum{};
  auto rec=[&](auto&& dfs,int j)->void{
   if(!(k==1&&neg)&&best==low[f])return;
   if(j==(int)p.size()){
    int m=absnorm(sum);assert(m%8==0);m/=8;
    if(k==1&&neg)best=max(best,m);else best=min(best,m);
    return;
   }
   for(const Vec&w:local[u[p[j]]]){
    for(int i=0;i<8;i++)sum[i]+=w[i];
    dfs(dfs,j+1);
    for(int i=0;i<8;i++)sum[i]-=w[i];
   }
  };
  rec(rec,0);
  array<int,4> w={score[0],score[1],score[2],score[3]};
  w[k]+=best-low[f];insert(w);
 }
 return out;
}
template<class T> void put(T c){cout.write(reinterpret_cast<const char*>(&c),sizeof(T));}
int uniform8(){
 heavy=5;bad=0;gate=8;
 for(auto&r:R)for(int&i:r)cin>>i;assert(cin);localinit();
 init();auto N=calc({});
 acache.resize(374);for(int i=0;i<374;i++)acache[i].resize(N[i].size());
 int count=0;L sum=0;
 for(int i=0;i<374;i++)for(int m=0;m<(int)N[i].size();m++){
  assert(N[i][m]>=0);sum+=N[i][m];if(!N[i][m])continue;++count;
  if((acache[i][m]=allowed0(D[i],m)))C.push_back({i,m});
 }
 assert(sum==(1LL<<40));
 auto code=orb(N,false),cha=orb(N,true);
 cerr<<count<<" ordered O8 orbits, "<<C.size()<<" allowed; sym "<<code.size()<<" allowed, "<<cha.size()<<" char orbit rows"<<endl;
 L Z=0;
 for(const auto& v:code){L s=0; for(auto [a,b]:v)s+=N[a][b];cout<<s<<" ";Z+=s;}
 cout<<endl;cerr<<"allowed cardinal "<<Z<<endl;
 for(const auto& V:cha){
  auto [a,b]=V.front();if(N[a][b]==1){assert(a==0&&b==0);continue;}
  auto x=rep(D[a],b);auto M=calc(x);
  L t=0;for(auto&z:M)for(L n:z)t+=n;assert(!t);
  for(const auto& v:code){L s=0;for(auto[i,m]:v)s+=M[i][m];cout<<s<<" ";}
  cout<<"\n";
 }
 return 0;
}

int main(int argc,char **argv){
 if(argc>1&&string(argv[1])=="uniform8")return uniform8();
 if(argc>1)heavy=atoi(argv[1]);assert(2<=heavy&&heavy<=5);
 neg=(argc>3&&string(argv[3])=="negative");
 bad=argc>2?atoi(argv[2]):2;
 one=bad==1&&argc>4?atoi(argv[4]):1;
 assert(bad>=0&&bad<=3&&heavy+1>=bad);
 for(auto&r:R)for(int&i:r)cin>>i;assert(cin);localinit();
 init();auto N=calc({});
 for(int i=0;i<374;i++)for(int j=heavy;j<5;j++)if(D[i].f[1<<j])for(auto &n:N[i])n=0;
 L sum=0;int count=0;
 for(int i=0;i<374;i++)for(L x:N[i]){assert(x>=0);sum+=x;count+=(x>0);}
 assert(sum==(1LL<<(8*heavy)));
 auto orbits=orb(N,true); // dual representatives
 acache.resize(374);for(int i=0;i<374;i++){acache[i].resize(N[i].size(),1);}
 auto code=orb(N,false); // all difference orbits including origin
 assert(code.front().size()==1&&code.front()[0].first==0&&code.front()[0].second==0);
 assert(orbits.front().size()==1&&orbits.front()[0].first==0&&orbits.front()[0].second==0);
 cerr<<count<<" ordered orbits, "<<code.size()-1<<" cols, "<<orbits.size()-1<<" rows"<<endl;
 put(uint32_t(code.size()-1));put(uint32_t(orbits.size()-1));
 for(size_t j=1;j<code.size();j++){
  const auto &v=code[j];L n=0; for(auto[a,b]:v)n+=N[a][b];
  auto [a,b]=v.front();put(n);auto w=many(D[a],b);put(uint16_t(w.size())); for(auto t:w)for(int h=0;h<(bad==3?4:3);h++)put(uint16_t(t[h]));
 }
 for(size_t i=1;i<orbits.size();i++){
  auto [a,b]=orbits[i].front();auto M=calc(rep(D[a],b));
  L tally=0;for(size_t j=1;j<code.size();j++){L n=0;for(auto[c,d]:code[j])n+=M[c][d];put(n);tally+=n;}
  assert(tally==-1);
 }
}
