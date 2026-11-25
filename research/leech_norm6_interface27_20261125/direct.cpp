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

struct T{unsigned p=0,n=0;T(){}T(const V&a){for(int i=0;i<24;i++){if(md(a[i])==1)p|=1u<<i;if(md(a[i])==2)n|=1u<<i;}}int times(T b)const{return (__builtin_popcount((p&b.p)|(n&b.n))-__builtin_popcount((p&b.n)|(n&b.p))+48)%3;}};
V target[10];vector<T>s[10];T tar[10];int k=0;long long N=0,deep[10]{},orient[10]{},aligned_marks[10]{};
void init(V a){Enum e;e.run(a,3,233);target[k]=a;tar[k]=T(a);V y=e.out.at(0).first;for(auto&[v,d]:e.out)s[k].emplace_back(minusv(v,y));mt19937 gen(34+k);shuffle(s[k].begin(),s[k].end(),gen);k++;}
void ck(V& m){T f(m);N++;for(int j=0;j<k;j++){if(f.times(tar[j]))continue;int color=1;for(auto t:s[j]){color|=1<<f.times(t);if(color==7)break;}if(color!=7){deep[j]++;orient[j]+=3-__builtin_popcount(color);bool c1=true,c2=true;for(int i=0;i<24;i++){if(md(m[i]-target[j][i]))c1=false;if(md(m[i]+target[j][i]))c2=false;}if(c1||c2)aligned_marks[j]++;}}}
int main(){setup();for(unsigned s:{1,3,7,31,63,127,383,959}){V v{};for(int i=0;i<24;i++)if(s>>i&1)v[i]=8;init(v);}V v={0,0,0,0,-8,0,0,0,0,0,0,-8,8,0,8,8,8,0,-8,8,0,0,0,8};init(v);
 for(auto c:words){int w=__builtin_popcount(c);V m{};if(w==12||w==8){vector<int> ix;for(int i=0;i<24;i++)if(c>>i&1)ix.push_back(i);
 for(int sign=0;sign<(1<<w);sign++){for(int j=0;j<w;j++)m[ix[j]]=(sign>>j&1)?-2:2;if(w==12){if(!(__builtin_popcount((unsigned)sign)&1))ck(m);}else if(__builtin_popcount((unsigned)sign)&1){for(int i=0;i<24;i++)if(!(c>>i&1)){m[i]=4;ck(m);m[i]=-4;ck(m);m[i]=0;}}}}
 for(int i=0;i<24;i++)m[i]=(c>>i&1)?-1:1;
 for(int i=0;i<24;i++){m[i]*=5;ck(m);m[i]/=5;
 m[i]*=-3;for(int j=i+1;j<24;j++){m[j]*=-3;for(int l=j+1;l<24;l++){m[l]*=-3;ck(m);m[l]/=-3;}m[j]/=-3;}m[i]/=-3;}}
 cout<<"marks "<<N<<endl;assert(N==16773120);for(int j=0;j<k;j++){cout<<"representative "<<j<<": "<<deep[j]<<" signed marks, "<<orient[j]<<" lifts, "<<aligned_marks[j]<<" aligned marks\n";assert(deep[j]==(j==8?2:0)&&orient[j]==(j==8?4:0)&&aligned_marks[j]==deep[j]);}}
