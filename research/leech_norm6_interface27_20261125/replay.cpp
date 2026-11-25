#include <bits/stdc++.h>
using namespace std;
using V=array<int,24>; using W=unsigned;
int parity(int a){return a&1;} vector<W> words;
void setup(){for(W j=0;j<4096;j++){W x=0;for(int i=0;i<12;i++)if(j>>i&1)x^=2787u<<i;x|=(__builtin_popcount(x)&1)<<23;words.push_back(x);}}
V minusv(V a,const V& b){for(int i=0;i<24;i++)a[i]-=b[i];return a;}
int dot(const V&a,const V&b){int s=0;for(int i=0;i<24;i++)s+=a[i]*b[i];return s;}
int md(int a){return (a%3+3)%3;}
struct Enum{struct O{int y,c,k;};vector<O>opts[2][2][24];vector<pair<V,int>>out;V v;int suff[25],bound;
void walk(int i,int p,W c,int par,int d){if(i==24){if(par==p)out.emplace_back(v,d);return;}
 for(auto a:opts[p][c>>i&1][i]){if(d+a.c+suff[i+1]>bound)break;v[i]=a.y;walk(i+1,p,c,par^parity(a.k),d+a.c);}}
void run(V a,int den,int b){bound=b;out.clear();for(int p=0;p<2;p++)for(int c=0;c<2;c++)for(int i=0;i<24;i++){
 auto &o=opts[p][c][i];o.clear();int kc=(a[i]/den-p-2*c)/4;for(int k=kc-4;k<=kc+4;k++){int y=p+2*c+4*k,z=den*y-a[i];o.push_back({y,z*z,k});}sort(o.begin(),o.end(),[](auto x,auto y){return x.c<y.c;});}
 for(int p=0;p<2;p++)for(auto c:words){suff[24]=0;for(int i=23;i>=0;i--)suff[i]=suff[i+1]+opts[p][c>>i&1][i][0].c;if(suff[0]<=bound)walk(0,p,c,0,0);}}
};
struct Space{V r[24]{};int rank=0;void add(V x){for(auto&v:x)v=md(v);for(int i=0;i<24;i++)if(x[i]){if(r[i][i]){int u=x[i];for(int j=i;j<24;j++)x[j]=md(x[j]-u*r[i][j]);}else{int u=x[i];for(int j=i;j<24;j++)r[i][j]=md(x[j]*u);rank++;return;}}}
 vector<V> ker(){vector<V>b;for(int j=0;j<24;j++)if(!r[j][j]){V a{};a[j]=1;for(int i=23;i>=0;i--)if(r[i][i])a[i]=md(-dot(a,r[i]));b.push_back(a);}return b;}};
void show(V v){for(int i=0;i<24;i++)cout<<v[i]<<(i<23?",":"");}
struct Witness{int n,nn; V v;};
int main(){setup();Enum e;
 vector<Witness> w;auto mk=[&](int n,int nn,unsigned mask,unsigned neg){V v{};for(int i=0;i<24;i++)if(mask>>i&1)v[i]=(neg>>i&1)?-8:8;w.push_back({n,nn,v});};
 mk(0,1,0,0);mk(8,1,1,0);mk(4,1,3,0);mk(12,3,7,0);mk(16,9,31,0);mk(12,1,63,0);mk(14,2,127,0);mk(10,1,383,0);mk(18,36,959,0);
 V s={0,0,0,0,-8,0,0,0,0,0,0,-8,8,0,8,8,8,0,-8,8,0,0,0,8};w.push_back({6,1,s});
 set<pair<int,int>> ids;
 for(auto [n,nmin,v]:w){e.run(v,3,233);int b=INT_MAX,ct=0;Space sp;V base=e.out.at(0).first;for(auto&[y,d]:e.out){sp.add(minusv(y,base));if(d<b)b=d,ct=0;if(d==b)ct++;}
  // A second, coordinate-factorized nearest-distance and multiplicity check.
  int lo=INT_MAX,mult=0;
  for(int p=0;p<2;p++)for(auto c:words){array<pair<int,int>,2> dp{make_pair(0,1),make_pair(100000,0)};
   for(int i=0;i<24;i++){array<pair<int,int>,2> nx{make_pair(100000,0),make_pair(100000,0)};for(int k=-3;k<=3;k++){int a=v[i]-3*(p+2*((c>>i)&1)+4*k);for(int j=0;j<2;j++){auto &out=nx[j^(k&1)];int cand=dp[j].first+a*a;if(cand<out.first)out={cand,dp[j].second};else if(cand==out.first)out.second+=dp[j].second;}}dp=nx;}
   if(dp[p].first<lo)lo=dp[p].first,mult=0;if(dp[p].first==lo)mult+=dp[p].second;
  }
  assert(b==lo&&ct==mult&&b==8*n&&ct==nmin&&ids.insert({n,nmin}).second);
  assert(sp.rank==(n==0?0:n==6?23:24));if(n==6){auto f=sp.ker().at(0);bool col=false;for(int l=1;l<3;l++){bool same=true;for(int j=0;j<24;j++)if(md(l*v[j])!=f[j])same=false;col|=same;}assert(col);}
  cout<<n<<" "<<nmin<<" "<<e.out.size()<<" "<<sp.rank<<"\n";
 }assert(ids.size()==10);}
