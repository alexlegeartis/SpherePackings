# Sphere 91: 39D horizontal phase search. Source replay
Source inputs come from alexlegeartis/SpherePackings commit c8bf217549171ce424cb6a03733d18fb7b069e29 via research/e7_hermitian35_20261123/certificate.json and research/e7_lattice38_20261123/{certificate38.json,local4.in,local4.cpp}. Compile with g++ -O3 -mssse3 -std=c++17. 

The common whitespace-separated input “data” is phase_isometry, alpha, horizontal_generators_times_2 (flatten row-major), vertical_generators_times_2 (flatten row-major); the command line then appends local distance output. The scripts use fixed seed 393891, GNU libstdc++ shuffle. Before first and continuous half searches run local4 <local4.in >tab4 and concatenate data + tab4. For quarter search run local8 <local4.in >tab8 and concatenate data + tab8. Calls this round: half 16000, quarter 20000 0 0 0 0 0, cont 12000. Printed progress interval was changed to 500 for half and 200 for quarter; no algorithm change.

## Half-coordinate solver half.cpp
```cpp
#include <bits/stdc++.h>
#include <tmmintrin.h>
using namespace std;
using byte1=uint8_t;using Shape=array<int,5>;using Vec=array<int,3>;
int M[128],invM[128],alpha[5],A[3][5][7],V[3][3],d4[16384],cols[5][128],vp[128][16],ijd[128][128],ijs[128][128],twiceV[128],G[16][5][7],typeg[16],off[16],Synd[5][128],depth[128],QM[128],bs[5],cnt=0, gd1=-1;double gh=100;
alignas(64) byte1 hist[4][8192],base[128][8192],f[8192],z[8192];vector<pair<int,Shape>> elite;mt19937 gen(393891);
int mul(int a,int b){int r=0;while(b){if(b&1)r^=a;b>>=1;a<<=1;if(a&64)a^=67;}return r;}
int sid(int j,int mask){return cols[j][invM[mask]];}
int mod(int x,int n){return (x%n+n)%n;}
int id(Vec x){int a=mod(x[0],4),b=mod(x[1],4),c=mod(x[2],4);a+=4*(1&((x[0]-a+x[1]-b+x[2]-c)/4));return a+8*b+32*c;}
Vec de(int i){return {i%8,i/8%4,i/32};}
int vert(Vec x){int best=100000;for(int a=-4;a<=4;a++)for(int b=-4;b<=4;b++)for(int c=-4;c<=4;c++)if(!((a+b+c)&1)){int s=(x[0]-4*a)*(x[0]-4*a)+(x[1]-4*b)*(x[1]-4*b)+(x[2]-4*c)*(x[2]-4*c);best=min(best,s);}return best;}
inline void xormin(byte1* dst,byte1* src,int o,int w){int hi=o&~15;__m128i mask=_mm_xor_si128(_mm_setr_epi8(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),_mm_set1_epi8(o&15)),W=_mm_set1_epi8(w);for(int a=0;a<8192;a+=16){__m128i v=_mm_loadu_si128((__m128i*)(src+(a^hi)));v=_mm_adds_epu8(_mm_shuffle_epi8(v,mask),W);__m128i d=_mm_loadu_si128((__m128i*)(dst+a));_mm_storeu_si128((__m128i*)(dst+a),_mm_min_epu8(v,d));}}
void dp(Shape R,byte1*out){memset(f,100,sizeof f);f[0]=0;
for(int j=0;j<5;j++){int raw=0;for(int k=0;k<7;k++)raw|=((R[j]>>k)&1)<<(2*k);memset(z,100,sizeof z);
byte1 u[128];for(int b=0;b<128;b++)u[b]=d4[raw^(2*QM[b])];
if(j==0){for(int b=0;b<128;b++)z[cols[0][b]]=u[b];}
else if(j==1){for(int b=0;b<128;b++)for(int a=0;a<128;a++){int t=cols[0][a]^cols[1][b];z[t]=min(z[t],(byte1)(f[cols[0][a]]+u[b]));}}
else for(int b=0;b<128;b++)xormin(z,f,cols[j][b],u[b]);
memcpy(f,z,sizeof f);}memcpy(out,f,sizeof f);}
void prep(Shape R){
bool done[4]={};for(int g=0;g<16;g++){off[g]=0;Shape P{};for(int j=0;j<5;j++){int carry=0;for(int t=0;t<7;t++){int a=((R[j]>>t)&1)-G[g][j][t],p=a&1;P[j]|=p<<t;carry|=(((a-p)/2)&1)<<t;}off[g]^=sid(j,carry);}
int t=typeg[g];if(!done[t]){dp(P,hist[t]);done[t]=true;}}
}
struct Item{int s,v,d;};
double gap(int d){return sqrt(max(0.,8.-d*.25));}
string fmt(Shape R){ostringstream o;for(int x:R)o<<x<<",";return o.str();}
void score(Shape R,int so,int vo,int d,int e,int s2,int v2,int d2,int diag1,int diag2){double half=max({gap(d)+gap(d2),gap(min(diag1,diag2)),sqrt(2.)});double h=half/2;
if(h<gh-1e-9){gh=h;cout<<"BEST "<<cnt<<" "<<(e?"uniform":"parallelogram")<<" R "<<fmt(R)<<" q "<<so<<" "<<vo<<" depths4 "<<d<<" r "<<s2<<" "<<v2<<" "<<d2<<" diagonals "<<diag1<<" "<<diag2<<" avgH "<<h<<" density "<<16/h<<endl;}
}
void setup(){
for(int& x:M)cin>>x;for(int&x:alpha)cin>>x;for(auto&a:A)for(auto&b:a)for(int&x:b)cin>>x;for(auto&a:V)for(int&x:a)cin>>x;for(int &x:d4)cin>>x;
for(int i=0;i<128;i++){invM[M[i]]=i;for(int j=0;j<7;j++)QM[i]|=((M[i]>>j)&1)<<(2*j);for(int j=0;j<5;j++)cols[j][i]=(i&63)|(mul(alpha[j],i&63)<<6)|((i>>6)<<12);}
for(int r=0;r<128;r++){int raw=0;for(int t=0;t<7;t++)raw|=((r>>t)&1)<<(2*t);depth[r]=100;for(int i=0;i<128;i++)depth[r]=min(depth[r],d4[raw^(2*QM[i])]);}
for(int j=0;j<5;j++)for(int k=0;k<7;k++)bs[j]|=(A[1][j][k]&1)<<k;
for(int g=0;g<16;g++){int k=g>>2,l=(g>>1)&1,m=g&1;typeg[g]=(k&1)+2*((l+m)&1);for(int j=0;j<5;j++)for(int t=0;t<7;t++)G[g][j][t]=k*A[0][j][t]+l*A[1][j][t]+m*A[2][j][t];Vec w;for(int a=0;a<3;a++)w[a]=k*V[0][a]+l*V[1][a]+m*V[2][a];for(int v=0;v<128;v++){Vec x=de(v);for(int a=0;a<3;a++)x[a]-=w[a];vp[v][g]=vert(x);}}
for(int i=0;i<128;i++)for(int j=0;j<128;j++){Vec a=de(i),b=de(j),x,y;for(int k=0;k<3;k++)x[k]=a[k]+b[k],y[k]=a[k]-b[k];ijs[i][j]=id(x);ijd[i][j]=id(y);if(i==j)twiceV[i]=id(x);}
prep({0,0,0,0,0});int ch[16]={0,24,24,24,30,26,26,30,24,24,24,16,30,26,26,30};int cv[16]={0,8,8,8,2,6,6,10,8,8,8,16,2,6,6,10};for(int g=0;g<16;g++){assert(hist[typeg[g]][off[g]]==ch[g]);assert(vp[0][g]==cv[g]);}
for(int v=0;v<128;v++){memset(base[v],240,8192);for(int g=0;g<16;g++)xormin(base[v],hist[typeg[g]],off[g],vp[v][g]);}assert(base[0][0]==0);
}
int trial(Shape R){cnt++;prep(R);int sr=0;for(int j=0;j<5;j++)sr^=sid(j,R[j]);int best=-1;Item qi{};static byte1 vals[128][8192];
for(int v=0;v<128;v++){byte1 *dst=vals[v];memset(dst,100,8192);for(int g=0;g<16;g++){xormin(dst,hist[typeg[g]],off[g],vp[v][g]);}
for(int s=0;s<8192;s++)if(dst[s]>best){best=dst[s];qi={s,v,best};}}
if(best>gd1){gd1=best;cout<<"DEPTH "<<cnt<<" depths4 "<<best<<" R "<<fmt(R)<<" phase "<<qi.s<<" "<<qi.v<<endl;}
vector<Item> top;int threshold=max(20,best-2);if(best>=29){double rr=1.7431-gap(best);threshold=max(0,(int)floor(4*(8-rr*rr)));}// include equality / slight rounding conservatively
for(int v=0;v<128;v++)for(int s=0;s<8192;s++)if(vals[v][s]>=threshold)top.push_back({s,v,vals[v][s]});
if(top.empty())return best;
shuffle(top.begin(),top.end(),gen);stable_sort(top.begin(),top.end(),[](const Item&a,const Item&b){return a.d>b.d;});
if(top.size()>1200)top.resize(1200);
for(auto p:top){Vec v3=de(p.v);for(int &x:v3)x*=-3;int d3=vals[id(v3)][p.s];int d2=base[twiceV[p.v]][sr];double h=max({gap(p.d),gap(d2)/2,gap(d3)/3,sqrt(.5)});if(h<gh-1e-9){gh=h;cout<<"BEST "<<cnt<<" uniform R "<<fmt(R)<<" q "<<p.s<<" "<<p.v<<" depths4 "<<p.d<<" "<<d2<<" "<<d3<<" avgH "<<h<<" density "<<16/h<<endl;}}
for(int i=0;i<(int)top.size();i++){auto p=top[i];for(int j=i+1;j<(int)top.size();j++){auto q=top[j];if(gap(p.d)+gap(q.d)>=2*gh-1e-9)break;int s=p.s^q.s;int a=base[ijs[p.v][q.v]][s^sr],b=base[ijd[p.v][q.v]][s];score(R,p.s,p.v,p.d,0,q.s,q.v,q.d,a,b);}}
return best;}
int main(int argc,char**argv){ios::sync_with_stdio(false);cout<<setprecision(13);setup();int lim=argc>1?atoi(argv[1]):100000;set<Shape>seen;
for(int k=0;cnt<lim;k++){Shape R{};for(int j=0;j<5;j++){int typ=k<2?0:(int)(gen()%4);if(typ==0)R[j]=(gen()&1)?82:0;else if(typ==1)R[j]=bs[j]^((gen()&1)?82:0);else{do R[j]=gen()%128;while(depth[R[j]]<4);}} if(k==0)R={0,0,0,0,0};if(k==1)R={82,82,82,82,82};
if(k>=10&&!elite.empty()&&(gen()%100<72)){R=elite[gen()%elite.size()].second;int edits=(gen()%8<6)?1:2;for(int c=0;c<edits;c++){int j=gen()%5,b=gen()%128;switch(gen()%4){case 0:R[j]^=82;break;case 1:R[j]^=bs[j];break;case 2: R[j]=b;break;default:R[j]=gen()%2?82:bs[j];}}}
if(!seen.insert(R).second)continue;int b=trial(R);if(elite.size()<24||b>=elite.back().first){elite.push_back({b,R});shuffle(elite.begin(),elite.end(),gen);stable_sort(elite.begin(),elite.end(),[](auto &a,auto &b){return a.first>b.first;});if(elite.size()>24)elite.pop_back();}if(cnt%25==0)cout<<"PROGRESS "<<cnt<<" max4 "<<gd1<<" bestDensity "<<16/gh<<endl;if(16/gh>18.35927)break;}cout<<"DONE "<<cnt<<endl;}
```

## Quarter-coordinate solver quarter.cpp
```cpp
#include <bits/stdc++.h>
#include <tmmintrin.h>
using namespace std;
using byte1=uint8_t;using Shape=array<int,5>;using Vec=array<int,3>;
int M[128],invM[128],alpha[5],A[3][5][7],V[3][3],d4[2097152],cols[5][128],vp[128][16],ijd[128][128],ijs[128][128],twiceV[128],G[16][5][7],typeg[16],off[16],Synd[5][128],depth[128],QM[128],bs[5],cnt=0, gd1=-1;double gh=100;
alignas(64) byte1 doubled[128][8192];alignas(64) byte1 hist[4][8192],base[128][8192],f[8192],z[8192];vector<pair<int,Shape>> elite;mt19937 gen(393891);
int mul(int a,int b){int r=0;while(b){if(b&1)r^=a;b>>=1;a<<=1;if(a&64)a^=67;}return r;}
int sid(int j,int mask){return cols[j][invM[mask]];}
int mod(int x,int n){return (x%n+n)%n;}
int id(Vec x){int a=mod(x[0],4),b=mod(x[1],4),c=mod(x[2],4);a+=4*(1&((x[0]-a+x[1]-b+x[2]-c)/4));return a+8*b+32*c;}
Vec de(int i){return {i%8,i/8%4,i/32};}
int vert(Vec x){int best=100000;for(int a=-4;a<=4;a++)for(int b=-4;b<=4;b++)for(int c=-4;c<=4;c++)if(!((a+b+c)&1)){int s=(x[0]-4*a)*(x[0]-4*a)+(x[1]-4*b)*(x[1]-4*b)+(x[2]-4*c)*(x[2]-4*c);best=min(best,s);}return best;}
inline void xormin(byte1* dst,byte1* src,int o,int w){int hi=o&~15;__m128i mask=_mm_xor_si128(_mm_setr_epi8(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),_mm_set1_epi8(o&15)),W=_mm_set1_epi8(w);for(int a=0;a<8192;a+=16){__m128i v=_mm_loadu_si128((__m128i*)(src+(a^hi)));v=_mm_adds_epu8(_mm_shuffle_epi8(v,mask),W);__m128i d=_mm_loadu_si128((__m128i*)(dst+a));_mm_storeu_si128((__m128i*)(dst+a),_mm_min_epu8(v,d));}}
void dp(Shape R,byte1*out){memset(f,198,sizeof f);f[0]=0;
for(int j=0;j<5;j++){int raw=0;for(int k=0;k<7;k++)raw|=((R[j]>>(2*k))&3)<<(3*k);memset(z,198,sizeof z);
byte1 u[128];for(int b=0;b<128;b++)u[b]=d4[raw^(4*QM[b])];
if(j==0){for(int b=0;b<128;b++)z[cols[0][b]]=u[b];}
else if(j==1){for(int b=0;b<128;b++)for(int a=0;a<128;a++){int t=cols[0][a]^cols[1][b];z[t]=min(z[t],(byte1)min(250,f[cols[0][a]]+u[b]));}}
else for(int b=0;b<128;b++)xormin(z,f,cols[j][b],u[b]);
memcpy(f,z,sizeof f);}memcpy(out,f,sizeof f);}
void prep(Shape R){
bool done[4]={};for(int g=0;g<16;g++){off[g]=0;Shape P{};for(int j=0;j<5;j++){int carry=0;for(int t=0;t<7;t++){int a=((R[j]>>(2*t))&3)-2*G[g][j][t],p=a&3;P[j]|=p<<(2*t);carry|=(((a-p)/4)&1)<<t;}off[g]^=sid(j,carry);}
int t=typeg[g];if(!done[t]){dp(P,hist[t]);done[t]=true;}}
}
struct Item{int s,v,d;};
double gap(int d){return sqrt(max(0.,8.-d*.0625));}
string fmt(Shape R){ostringstream o;for(int x:R)o<<x<<",";return o.str();}
void score(Shape R,int so,int vo,int d,int e,int s2,int v2,int d2,int diag1,int diag2){double half=max({gap(d)+gap(d2),gap(min(diag1,diag2)),sqrt(2.)});double h=half/2;
if(h<gh-1e-9){gh=h;cout<<"BEST "<<cnt<<" "<<(e?"uniform":"parallelogram")<<" R "<<fmt(R)<<" q "<<so<<" "<<vo<<" depths16 "<<d<<" r "<<s2<<" "<<v2<<" "<<d2<<" diagonals "<<diag1<<" "<<diag2<<" avgH "<<h<<" density "<<16/h<<endl;}
}
void setup(){
for(int& x:M)cin>>x;for(int&x:alpha)cin>>x;for(auto&a:A)for(auto&b:a)for(int&x:b)cin>>x;for(auto&a:V)for(int&x:a)cin>>x;for(int &x:d4)cin>>x;
for(int i=0;i<128;i++){invM[M[i]]=i;for(int j=0;j<7;j++)QM[i]|=((M[i]>>j)&1)<<(3*j);for(int j=0;j<5;j++)cols[j][i]=(i&63)|(mul(alpha[j],i&63)<<6)|((i>>6)<<12);}

for(int j=0;j<5;j++)for(int k=0;k<7;k++)bs[j]|=(A[1][j][k]&1)<<k;
for(int g=0;g<16;g++){int k=g>>2,l=(g>>1)&1,m=g&1;typeg[g]=(k&1)+2*((l+m)&1);for(int j=0;j<5;j++)for(int t=0;t<7;t++)G[g][j][t]=k*A[0][j][t]+l*A[1][j][t]+m*A[2][j][t];Vec w;for(int a=0;a<3;a++)w[a]=k*V[0][a]+l*V[1][a]+m*V[2][a];for(int v=0;v<128;v++){Vec x=de(v);for(int a=0;a<3;a++)x[a]-=w[a];vp[v][g]=4*vert(x);}}
for(int i=0;i<128;i++)for(int j=0;j<128;j++){Vec a=de(i),b=de(j),x,y;for(int k=0;k<3;k++)x[k]=a[k]+b[k],y[k]=a[k]-b[k];ijs[i][j]=id(x);ijd[i][j]=id(y);if(i==j)twiceV[i]=id(x);}
prep({0,0,0,0,0});int ch[16]={0,24,24,24,30,26,26,30,24,24,24,16,30,26,26,30};int cv[16]={0,8,8,8,2,6,6,10,8,8,8,16,2,6,6,10};for(int g=0;g<16;g++){assert(hist[typeg[g]][off[g]]==4*ch[g]);assert(vp[0][g]==4*cv[g]);}
for(int v=0;v<128;v++){memset(base[v],240,8192);for(int g=0;g<16;g++)xormin(base[v],hist[typeg[g]],off[g],vp[v][g]);}assert(base[0][0]==0);
}
int trial(Shape R){cnt++;int sr=0,sy4=0;Shape DB{};for(int j=0;j<5;j++){int a=0,b=0;for(int k=0;k<7;k++){int x=(R[j]>>(2*k))&3;a|=(x>>1)<<k;b|=(x&1)<<k;DB[j]|=(2*(x&1))<<(2*k);}sr^=sid(j,a);sy4^=sid(j,b);}prep(DB);for(int v=0;v<128;v++){memset(doubled[v],240,8192);for(int g=0;g<16;g++)xormin(doubled[v],hist[typeg[g]],off[g],vp[v][g]);}prep(R);int best=-1;Item qi{};static byte1 vals[128][8192];
for(int v=0;v<128;v++){byte1 *dst=vals[v];memset(dst,198,8192);for(int g=0;g<16;g++){xormin(dst,hist[typeg[g]],off[g],vp[v][g]);}
for(int s=0;s<8192;s++)if(dst[s]>best){best=dst[s];qi={s,v,best};}}
if(best>gd1){gd1=best;cout<<"DEPTH "<<cnt<<" depths16 "<<best<<" R "<<fmt(R)<<" phase "<<qi.s<<" "<<qi.v<<endl;}
vector<Item> top;int threshold=max(80,best-5);if(best>=116){double rr=1.7431-gap(best);threshold=max(0,(int)floor(16*(8-rr*rr)));}// include equality / slight rounding conservatively
for(int v=0;v<128;v++)for(int s=0;s<8192;s++)if(vals[v][s]>=threshold)top.push_back({s,v,vals[v][s]});
if(top.empty())return best;
shuffle(top.begin(),top.end(),gen);stable_sort(top.begin(),top.end(),[](const Item&a,const Item&b){return a.d>b.d;});
if(top.size()>1200)top.resize(1200);
for(auto p:top){Vec v3=de(p.v);for(int &x:v3)x*=-3;int d3=vals[id(v3)][p.s^sy4];int d2=doubled[twiceV[p.v]][sr];double h=max({gap(p.d),gap(d2)/2,gap(d3)/3,sqrt(.5)});if(h<gh-1e-9){gh=h;cout<<"BEST "<<cnt<<" uniform R "<<fmt(R)<<" q "<<p.s<<" "<<p.v<<" depths16 "<<p.d<<" "<<d2<<" "<<d3<<" avgH "<<h<<" density "<<16/h<<endl;}}
for(int i=0;i<(int)top.size();i++){auto p=top[i];for(int j=i+1;j<(int)top.size();j++){auto q=top[j];if(gap(p.d)+gap(q.d)>=2*gh-1e-9)break;int s=p.s^q.s;int a=doubled[ijs[p.v][q.v]][s^sr],b=base[ijd[p.v][q.v]][s];score(R,p.s,p.v,p.d,0,q.s,q.v,q.d,a,b);}}
return best;}
int main(int argc,char**argv){ios::sync_with_stdio(false);cout<<setprecision(13);setup();int lim=argc>1?atoi(argv[1]):100000;set<Shape>seen;Shape seed{};for(int j=0;j<5;j++){int m=argc>2+j?atoi(argv[2+j]):82;for(int k=0;k<7;k++)seed[j]|=(2*((m>>k)&1))<<(2*k);}for(int k=0;cnt<lim;k++){Shape R=seed;if(!elite.empty()&&gen()%100<80)R=elite[gen()%elite.size()].second;int edits=1+(gen()%5==0);if(k)for(int c=0;c<edits;c++){int j=gen()%5,i=gen()%7;int frac=(R[j]>>(2*i))&3;int nw=(frac+1+2*(gen()%2))&3;R[j]^=(frac^nw)<<(2*i);}if(!seen.insert(R).second)continue;int b=trial(R);if(elite.size()<24||b>=elite.back().first){elite.push_back({b,R});shuffle(elite.begin(),elite.end(),gen);stable_sort(elite.begin(),elite.end(),[](auto &a,auto &b){return a.first>b.first;});if(elite.size()>24)elite.pop_back();}if(cnt%20==0)cout<<"PROGRESS "<<cnt<<" max16 "<<gd1<<" bestDensity "<<16/gh<<endl;if(16/gh>18.35927)break;}cout<<"DONE "<<cnt<<endl;}
```

## Exact E7 mod-8 nearest-vector table local8.cpp
```cpp
#include<bits/stdc++.h>
using namespace std;int S,T,D[7],L[7][7],x[7],tab[2097152],lim;
int ceildiv(int x,int y){return x>=0?(x+y-1)/y:x/y;}
void f(int j,int rest,int c){
 if(j<0){tab[c]=min(tab[c],(lim-rest)/(S*S*T));assert((lim-rest)%(S*S*T)==0);return;}
 int w=0;for(int i=j+1;i<7;i++)w+=L[i][j]*x[i];int r=sqrt(rest/D[j]);while((r+1)*(r+1)*D[j]<=rest)r++;while(r*r*D[j]>rest)r--;
 for(int z=ceildiv(-r-w,S);z<=-ceildiv(w-r,S);z++){x[j]=z;int y=S*z+w;f(j-1,rest-D[j]*y*y,c|((z&7)<<(3*j)));}
}
int main(){
 cin>>S>>T;for(auto&x:D)cin>>x;for(auto&r:L)for(auto&x:r)cin>>x;
 fill(tab,tab+2097152,10000);lim=S*S*T*100;assert(lim<1000000000);f(6,lim,0);
 map<int,int>C;for(int x:tab)C[x]++;for(auto[k,v]:C)cerr<<k<<" "<<v<<"\n";assert(C.rbegin()->first<10000);
 for(int x:tab)cout<<x<<"\n";
}

```

## Continuous vertical branch cont.cpp
```cpp
#include <bits/stdc++.h>
#include <tmmintrin.h>
using namespace std;
using byte1=uint8_t;using Shape=array<int,5>;using Vec=array<int,3>;
int M[128],invM[128],alpha[5],A[3][5][7],V[3][3],d4[16384],cols[5][128],vp[128][16],ijd[128][128],ijs[128][128],twiceV[128],G[16][5][7],typeg[16],off[16],Synd[5][128],depth[128],QM[128],bs[5],cnt=0, gd1=-1;double gh=100;
alignas(64) byte1 hist[4][8192],base[128][8192],f[8192],z[8192];vector<pair<int,Shape>> elite;mt19937 gen(393891);
int mul(int a,int b){int r=0;while(b){if(b&1)r^=a;b>>=1;a<<=1;if(a&64)a^=67;}return r;}
int sid(int j,int mask){return cols[j][invM[mask]];}
int mod(int x,int n){return (x%n+n)%n;}
int id(Vec x){int a=mod(x[0],4),b=mod(x[1],4),c=mod(x[2],4);a+=4*(1&((x[0]-a+x[1]-b+x[2]-c)/4));return a+8*b+32*c;}
Vec de(int i){return {i%8,i/8%4,i/32};}
int vert(Vec x){int best=100000;for(int a=-4;a<=4;a++)for(int b=-4;b<=4;b++)for(int c=-4;c<=4;c++)if(!((a+b+c)&1)){int s=(x[0]-4*a)*(x[0]-4*a)+(x[1]-4*b)*(x[1]-4*b)+(x[2]-4*c)*(x[2]-4*c);best=min(best,s);}return best;}
inline void xormin(byte1* dst,byte1* src,int o,int w){int hi=o&~15;__m128i mask=_mm_xor_si128(_mm_setr_epi8(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),_mm_set1_epi8(o&15)),W=_mm_set1_epi8(w);for(int a=0;a<8192;a+=16){__m128i v=_mm_loadu_si128((__m128i*)(src+(a^hi)));v=_mm_adds_epu8(_mm_shuffle_epi8(v,mask),W);__m128i d=_mm_loadu_si128((__m128i*)(dst+a));_mm_storeu_si128((__m128i*)(dst+a),_mm_min_epu8(v,d));}}
void dp(Shape R,byte1*out){memset(f,100,sizeof f);f[0]=0;
for(int j=0;j<5;j++){int raw=0;for(int k=0;k<7;k++)raw|=((R[j]>>k)&1)<<(2*k);memset(z,100,sizeof z);
byte1 u[128];for(int b=0;b<128;b++)u[b]=d4[raw^(2*QM[b])];
if(j==0){for(int b=0;b<128;b++)z[cols[0][b]]=u[b];}
else if(j==1){for(int b=0;b<128;b++)for(int a=0;a<128;a++){int t=cols[0][a]^cols[1][b];z[t]=min(z[t],(byte1)(f[cols[0][a]]+u[b]));}}
else for(int b=0;b<128;b++)xormin(z,f,cols[j][b],u[b]);
memcpy(f,z,sizeof f);}memcpy(out,f,sizeof f);}
void prep(Shape R){
bool done[4]={};for(int g=0;g<16;g++){off[g]=0;Shape P{};for(int j=0;j<5;j++){int carry=0;for(int t=0;t<7;t++){int a=((R[j]>>t)&1)-G[g][j][t],p=a&1;P[j]|=p<<t;carry|=(((a-p)/2)&1)<<t;}off[g]^=sid(j,carry);}
int t=typeg[g];if(!done[t]){dp(P,hist[t]);done[t]=true;}}
}
struct Item{int s,v,d;};
double gap(int d){return sqrt(max(0.,8.-d*.25));}
string fmt(Shape R){ostringstream o;for(int x:R)o<<x<<",";return o.str();}
void score(Shape R,int so,int vo,int d,int e,int s2,int v2,int d2,int diag1,int diag2){double half=max({gap(d)+gap(d2),gap(min(diag1,diag2)),sqrt(2.)});double h=half/2;
if(h<gh-1e-9){gh=h;cout<<"BEST "<<cnt<<" "<<(e?"uniform":"parallelogram")<<" R "<<fmt(R)<<" q "<<so<<" "<<vo<<" depths4 "<<d<<" r "<<s2<<" "<<v2<<" "<<d2<<" diagonals "<<diag1<<" "<<diag2<<" avgH "<<h<<" density "<<16/h<<endl;}
}
void setup(){
for(int& x:M)cin>>x;for(int&x:alpha)cin>>x;for(auto&a:A)for(auto&b:a)for(int&x:b)cin>>x;for(auto&a:V)for(int&x:a)cin>>x;for(int &x:d4)cin>>x;
for(int i=0;i<128;i++){invM[M[i]]=i;for(int j=0;j<7;j++)QM[i]|=((M[i]>>j)&1)<<(2*j);for(int j=0;j<5;j++)cols[j][i]=(i&63)|(mul(alpha[j],i&63)<<6)|((i>>6)<<12);}
for(int r=0;r<128;r++){int raw=0;for(int t=0;t<7;t++)raw|=((r>>t)&1)<<(2*t);depth[r]=100;for(int i=0;i<128;i++)depth[r]=min(depth[r],d4[raw^(2*QM[i])]);}
for(int j=0;j<5;j++)for(int k=0;k<7;k++)bs[j]|=(A[1][j][k]&1)<<k;
for(int g=0;g<16;g++){int k=g>>2,l=(g>>1)&1,m=g&1;typeg[g]=(k&1)+2*((l+m)&1);for(int j=0;j<5;j++)for(int t=0;t<7;t++)G[g][j][t]=k*A[0][j][t]+l*A[1][j][t]+m*A[2][j][t];Vec w;for(int a=0;a<3;a++)w[a]=k*V[0][a]+l*V[1][a]+m*V[2][a];for(int v=0;v<128;v++){Vec x=de(v);for(int a=0;a<3;a++)x[a]-=w[a];vp[v][g]=vert(x);}}
for(int i=0;i<128;i++)for(int j=0;j<128;j++){Vec a=de(i),b=de(j),x,y;for(int k=0;k<3;k++)x[k]=a[k]+b[k],y[k]=a[k]-b[k];ijs[i][j]=id(x);ijd[i][j]=id(y);if(i==j)twiceV[i]=id(x);}
prep({0,0,0,0,0});int ch[16]={0,24,24,24,30,26,26,30,24,24,24,16,30,26,26,30};int cv[16]={0,8,8,8,2,6,6,10,8,8,8,16,2,6,6,10};for(int g=0;g<16;g++){assert(hist[typeg[g]][off[g]]==ch[g]);assert(vp[0][g]==cv[g]);}
for(int v=0;v<128;v++){memset(base[v],240,8192);for(int g=0;g<16;g++)xormin(base[v],hist[typeg[g]],off[g],vp[v][g]);}assert(base[0][0]==0);
}
int trial(Shape R){cnt++;prep(R);int sr=0;for(int j=0;j<5;j++)sr^=sid(j,R[j]);int best=-1;Item qi{};static byte1 vals[128][8192];
for(int v=0;v<128;v++){byte1 *dst=vals[v];memset(dst,100,8192);for(int g=0;g<16;g++){xormin(dst,hist[typeg[g]],off[g],vp[v][g]);}
for(int s=0;s<8192;s++)if(dst[s]>best){best=dst[s];qi={s,v,best};}}
if(best>gd1){gd1=best;cout<<"DEPTH "<<cnt<<" depths4 "<<best<<" R "<<fmt(R)<<" phase "<<qi.s<<" "<<qi.v<<endl;}
vector<Item> top;int threshold=max(20,best-2);if(best>=29){double rr=1.7431-gap(best);threshold=max(0,(int)floor(4*(8-rr*rr)));}// include equality / slight rounding conservatively
for(int v=0;v<128;v++)for(int s=0;s<8192;s++)if(vals[v][s]>=threshold)top.push_back({s,v,vals[v][s]});
if(top.empty())return best;
shuffle(top.begin(),top.end(),gen);stable_sort(top.begin(),top.end(),[](const Item&a,const Item&b){return a.d>b.d;});
if(top.size()>1200)top.resize(1200);
for(auto p:top){Vec v3=de(p.v);for(int &x:v3)x*=-3;int d3=vals[id(v3)][p.s];int d2=base[twiceV[p.v]][sr];double h=max({gap(p.d),gap(d2)/2,gap(d3)/3,sqrt(.5)});if(h<gh-1e-9){gh=h;cout<<"BEST "<<cnt<<" uniform R "<<fmt(R)<<" q "<<p.s<<" "<<p.v<<" depths4 "<<p.d<<" "<<d2<<" "<<d3<<" avgH "<<h<<" density "<<16/h<<endl;}}
for(int i=0;i<(int)top.size();i++){auto p=top[i];for(int j=i+1;j<(int)top.size();j++){auto q=top[j];if(gap(p.d)+gap(q.d)>=2*gh-1e-9)break;int s=p.s^q.s;int a=base[ijs[p.v][q.v]][s^sr],b=base[ijd[p.v][q.v]][s];score(R,p.s,p.v,p.d,0,q.s,q.v,q.d,a,b);}}
return best;}

constexpr int NN=256,SS=NN*NN,TARG=474514;using Wt=array<byte1,16>;
struct Cand{Shape R;int s;Vec x;int d;}; map<Shape,vector<Cand>> winners;map<Shape,int> srp;
vector<Wt> boxweights;vector<vector<int>> eqs;int nshape=0,global1=0,here=0,remain=0,boxes=0;Shape current;map<array<int,4>,Cand> spool;
Vec cen[16];int fracSyn(Shape R){int s=0;for(int j=0;j<5;j++)s^=sid(j,R[j]);return s;}
Wt weight(int s){Wt w;for(int g=0;g<16;g++)w[g]=hist[typeg[g]][s^off[g]];return w;}
pair<int,int> boxcost(Vec x,Vec c,int r){int d[3][2];for(int j=0;j<3;j++)for(int p=0;p<2;p++){int w=abs(x[j]-c[j]-p*2*NN)%(4*NN);d[j][p]=min(w,4*NN-w);}int lo=INT_MAX,hi=INT_MAX;
for(int a=0;a<2;a++)for(int b=0;b<2;b++){int x=d[0][a],y=d[1][b],z=d[2][a^b];lo=min(lo,x*x+y*y+z*z);hi=min(hi,(x+r)*(x+r)+(y+r)*(y+r)+(z+r)*(z+r));}return {lo,hi};}
int exact(const Wt&w,Vec x){int d=INT_MAX;for(int g=0;g<16;g++)d=min(d,(int)w[g]*(SS/4)+boxcost(x,cen[g],0).first);return d;}
void process(Vec x,int rad,const vector<int>&ix){boxes++;int lo[16],hi[16];for(int g=0;g<16;g++){auto p=boxcost(x,cen[g],rad);lo[g]=p.first;hi[g]=p.second;}vector<int>next;for(int k:ix){int d=INT_MAX,u=INT_MAX;for(int g=0;g<16;g++){int w=boxweights[k][g]*(SS/4);d=min(d,w+lo[g]);u=min(u,w+hi[g]);}here=max(here,d);if(d>global1){global1=d;cout<<"CLEARANCE "<<nshape<<" "<<fmt(current)<<" s "<<eqs[k][0]<<" v "<<x[0]<<","<<x[1]<<","<<x[2]<<" d65536 "<<d<<" approx "<<double(d)/SS<<endl;}
if(d>=int(7.12*SS)&&spool.size()<50000){for(int i=0;i<min(12,(int)eqs[k].size());i++){int s=eqs[k][i];array<int,4>key={s,x[0]/64,x[1]/64,x[2]/64};auto jt=spool.find(key);if(jt==spool.end()||jt->second.d<d)spool[key]={current,s,x,d};}}
if(u>=TARG){if(rad==1||((d>=TARG)&&rad<=4)){remain++;continue;}next.push_back(k);}}
if(next.empty())return;int r=rad/2;for(int a:{-r,r})for(int b:{-r,r})for(int c:{-r,r})process({x[0]+a,x[1]+b,x[2]+c},r,next);}
double fg(int d){return sqrt(max(0.,8.-double(d)/SS));}
vector<Cand> go(Shape R){current=R;here=remain=boxes=0;spool.clear();prep(R);boxweights.clear();eqs.clear();map<Wt,int> ids;
for(int s=0;s<8192;s++){Wt w=weight(s);int k=ids.size();auto e=ids.emplace(w,k);if(e.second){boxweights.push_back(w);eqs.emplace_back();}eqs[e.first->second].push_back(s);}vector<int>I(boxweights.size());iota(I.begin(),I.end(),0);
for(int i=64;i<4*NN;i+=128)for(int j=64;j<2*NN;j+=128)for(int k=64;k<2*NN;k+=128)process({i,j,k},64,I);
vector<Cand>v;for(auto&p:spool)v.push_back(p.second);shuffle(v.begin(),v.end(),gen);stable_sort(v.begin(),v.end(),[](auto&a,auto&b){return a.d>b.d;});if(v.size()>480)v.resize(480);if(nshape<20||!v.empty()||nshape%200==0)cout<<"GEOM "<<nshape<<" "<<fmt(R)<<" weightTypes "<<boxweights.size()<<" boxes "<<boxes<<" unclosed "<<remain<<" maxCentre65536 "<<here<<" pool "<<v.size()<<endl;return v;}
string xy(Vec x){return to_string(x[0])+","+to_string(x[1])+","+to_string(x[2]);}
void testpair(vector<Cand> const &U,vector<Cand>const &V,Shape R,Shape T){Shape d{},cp{},cm{};for(int j=0;j<5;j++)d[j]=R[j]^T[j],cp[j]=R[j]&T[j],cm[j]=T[j]&~R[j];int sp=fracSyn(cp),sm=fracSyn(cm);prep(d);for(auto&p:U)for(auto&q:V){double gaps=fg(p.d)+fg(q.d);if(gaps>=2*gh-1e-12)break;Vec a,b;for(int j=0;j<3;j++)a[j]=p.x[j]+q.x[j],b[j]=p.x[j]-q.x[j];int a2=exact(weight(p.s^q.s^sp),a),b2=exact(weight(p.s^q.s^sm),b);
double h=.5*max({gaps,fg(min(a2,b2)),sqrt(2.)});if(h<gh-1e-10){gh=h;cout<<"WITNESS cycle "<<fmt(R)<<' '<<p.s<<' '<<xy(p.x)<<" | "<<fmt(T)<<' '<<q.s<<' '<<xy(q.x)<<" clearances65536 "<<p.d<<' '<<q.d<<' '<<a2<<' '<<b2<<" avgH "<<h<<" density "<<16/h<<endl;}}}
void testuniform(Shape R,vector<Cand>const &U){prep({0,0,0,0,0});Wt w=weight(fracSyn(R));prep(R);for(auto&p:U){Vec v=p.x,u=p.x;for(int j=0;j<3;j++){u[j]*=-3;v[j]*=2;}int a=exact(w,v),b=exact(weight(p.s),u);double h=max({fg(p.d),fg(a)/2,fg(b)/3,sqrt(.5)});if(h<gh-1e-10){gh=h;cout<<"WITNESS uniform "<<fmt(R)<<" s "<<p.s<<" v "<<xy(p.x)<<" clearances65536 "<<p.d<<' '<<a<<' '<<b<<" avgH "<<h<<" density "<<16/h<<endl;}}}
int main(int argc,char**argv){ios::sync_with_stdio(false);cout<<setprecision(14);setup();for(int g=0;g<16;g++){int k=g>>2,l=(g>>1)&1,m=g&1;for(int j=0;j<3;j++)cen[g][j]=(k*V[0][j]+l*V[1][j]+m*V[2][j])*(NN/2);}int lim=argc>1?atoi(argv[1]):5000;set<Shape>seen;vector<pair<int,Shape>> e;gh=1.;for(int k=0;nshape<lim;k++){Shape R{};if(k==0){}else if(k==1)R.fill(82);else if(k==2)copy(bs,bs+5,R.begin());else if(k<36){int selector=k-4;for(int j=0;j<5;j++)R[j]=(selector&(1<<j))?82:0;}else {for(int j=0;j<5;j++)R[j]=gen()%3==0?82:gen()%3==0?bs[j]:gen()%128;if(!e.empty()&&gen()%3!=0){R=e[gen()%e.size()].second;int j=gen()%5;R[j]=gen()%3==0?82:gen()%128;}}
if(!seen.insert(R).second)continue;nshape++;auto pool=go(R);if(!pool.empty()){testuniform(R,pool);testpair(pool,pool,R,R);for(auto &a:winners)testpair(pool,a.second,R,a.first);if(winners.size()<60)winners[R]=pool;}if(e.size()<25||here>=e.back().first){e.emplace_back(here,R);shuffle(e.begin(),e.end(),gen);stable_sort(e.begin(),e.end(),[](auto&a,auto&b){return a.first>b.first;});if(e.size()>25)e.pop_back();}if(16/gh>18.35927){cout<<"CANDIDATE STOP"<<endl;break;}}cout<<"FINISHED "<<nshape<<" density "<<16/gh<<" maxCentre65536 "<<global1<<endl;}
```
