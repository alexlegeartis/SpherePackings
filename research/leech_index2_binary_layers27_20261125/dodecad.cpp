#include <bits/stdc++.h>
using namespace std;using I=uint32_t;const I A=(1u<<24)-1;
int wt(I x){return __builtin_popcount(x);}int pa(I x){return __builtin_parity(x);}
struct Basis{I r[24]{};bool add(I x){for(int i=23;i>=0;i--)if(x>>i&1){if(r[i])x^=r[i];else{r[i]=x;return 1;}}return 0;} I red(I x){for(int i=23;i>=0;i--)if(x>>i&1)x^=r[i];return x;}}Gbase;
vector<I>G,P,B;I D,C;int ei,ej;vector<vector<I>> errors(4096);int col[24];int syn(I x){int v=0;while(x){int k=__builtin_ctz(x);x&=x-1;v^=col[k];}return v;}I p0[2048],b0[2048],bv[1024];int e0[2048],be[1024],bP[2048],bS[1024];vector<unsigned char> F(1<<24);
int q(I v){return ((wt(v&D)/2)&1)+2*((wt(v&C)/2)&1);}
int par(I v){return pa(v&D)+2*pa(v&C);}
struct Pt{I p,b;int e,h;};
Pt pt(int a){int h=a>>23,k=(a>>12)&2047,s=(a>>2)&1023;I p=p0[k],b=b0[k]^bv[s];int e=e0[k]^be[s]^par(bv[s]&b0[k])^(a&3);if(h){e^=par(b&1);b^=1;}return{p,b,e,h};}
int direct(int a){auto[p,b,e,h]=pt(a);int best=200;for(I r:G){int target=q(r);for(int od=0;od<2;od++){int sum=0,z=e^od,alt[2]={999,999};for(int k=0;k<24;k++){int j=(D>>k&1)?0:1;int v=2*int(p>>k&1)+4*int(b>>k&1)+h-4*int(r>>k&1)-2*od;int u=((v+4)&7)-4;z^=(((v-u)/8)&1)<<j;alt[j]=min(alt[j],64-16*abs(u));sum+=u*u;}z^=target;sum+=(z&1?alt[0]:0)+(z&2?alt[1]:0);best=min(best,sum/4);}}return best;}
int main(){auto start=chrono::steady_clock::now();for(int k=0;k<4096;k++){I v=0;for(int j=0;j<12;j++)if(k>>j&1)v^=2787u<<j;v|=pa(v)<<23;G.push_back(v);if(wt(v)==12&&(v&1))D=v;if(!(v&1)&&Gbase.add(v))P.push_back(v);}C=A^D;ei=__builtin_ctz(D);ej=__builtin_ctz(C);assert(ei==0&&P.size()==11&&Gbase.add(A));Basis E=Gbase;
for(int j=0;j<24;j++){I h=(1u<<j)^(1u<<((D>>j&1)?ei:ej));if(E.add(h))B.push_back(h);I x=Gbase.red(1u<<j);int v=0,i=0;for(int z=0;z<24;z++)if(!Gbase.r[z])v|=int(x>>z&1)<<i++;col[j]=v;}assert(B.size()==10);
for(int s=0;s<1024;s++){for(int j=0;j<10;j++)if(s>>j&1)bv[s]^=B[j];be[s]=q(bv[s]);bS[s]=syn(bv[s]);}for(int k=0;k<2048;k++){for(int i=0;i<24;i++){int t=0;for(int j=0;j<11;j++)if(k>>j&1)t+=int(P[j]>>i&1)+2*int((i==ei||i==ej)&&((wt(P[j]&D)/2)&1));p0[k]|=(t&1)<<i;b0[k]|=((t>>1)&1)<<i;e0[k]^=((t>>2)&1)<<((D>>i&1)?0:1);}assert(par(b0[k])==q(p0[k]));bP[k]=syn(b0[k]);}
vector<int> shortSyn(1<<24);for(int i=1;i<(1<<24);i++) shortSyn[i]=shortSyn[i&(i-1)]^col[__builtin_ctz(i)];for(int w=0;w<=8;w++)for(int i=0;i<(1<<24);i++)if(wt(i)==w)errors[shortSyn[i]].push_back(i);shortSyn.clear();shortSyn.shrink_to_fit();
auto cut=[&](int &mask,int e,int raw,int alt0,int alt1){if(raw>=24)return;for(int j=0;j<4;j++)if((mask>>j&1)&&raw+(((j^e)&1)?alt0:0)+(((j^e)&2)?alt1:0)<24)mask&=~(1<<j);};
vector<int>rq;for(I r:G)rq.push_back(q(r));int n[2]={};for(int k=0;k<2048;k++){I p=p0[k];int w=wt(p);for(int s=0;s<1024;s++){int a=(k<<12)|(s<<2);auto[p,b,e,h]=pt(a);int mask=15;
for(int od=0;od<2&&mask;od++){I u=od?(A^p):p,t=od?(b^u):b;int ee=od?(e^1^par(u&~b)):e;int W=wt(u);if(W>=24)continue;int a0=(u&D)?8:16,a1=(u&C)?8:16;for(int ri=0;ri<4096&&mask;ri++){I r=G[ri],bad=(t^r)&~u;int raw=W+4*wt(bad);if(raw>=24)continue;int c=ee^par(u&t&~r)^rq[ri];cut(mask,c,raw,(bad&D)?0:a0,(bad&C)?0:a1);}}
for(int j=0;j<4;j++)if(mask>>j&1)F[a+j]=1,n[0]++;a|=1<<23;auto u=pt(a);I d=u.b;mask=15;for(I er:errors[bP[k]^bS[s]^col[0]]){int ww=wt(er),raw=6+2*ww;if(raw>=24)break;I r=p^d^er;int e1=u.e^par(d&~r)^q(r),alt1=(er&D)?4:12,alt2=(er&C)?4:12;cut(mask,e1,raw,alt1,alt2);if(!mask)break;r=d^er;int e2=u.e^1^par((A^p)&~d&r)^par(p&d&~r)^q(r);cut(mask,e2,raw,alt1,alt2);if(!mask)break;}
for(int j=0;j<4;j++)if(mask>>j&1)F[a+j]=1,n[1]++;}
if(k%256==255)cerr<<" k "<<k<<" "<<n[0]<<" "<<n[1]<<" sec "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<endl;}
mt19937 rng(765);for(int i=0;i<500;i++){int a=rng()&0xffffff;if(bool(F[a])!=(direct(a)>=24)){cerr<<"FAIL "<<a<<" "<<direct(a)<<" "<<int(F[a])<<endl;return 2;}}
assert(n[0]==24553 && n[1]==24576 && D==6678945);
ofstream file("fardod.bin",ios::binary);for(int i=0;i<(1<<24);i++)if(F[i])file.write((char*)&i,4);return !file.good();}
