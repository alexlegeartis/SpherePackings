#include<bits/stdc++.h>
using namespace std;using V=vector<int>;int pc(int x){return __builtin_popcount(x);}int F[16],re[65536],ds[65536];V qual,rs;long long parallelograms,fulltests,geometry,rects,filtered;struct Z{int i,j,k,n;};vector<Z> outer;
int main(int ac,char**av){int ty=atoi(av[1]),gg=atoi(av[2]);const int cols[]={1,2,3,3,5,6,7,7,9,10,11,11,13,14,15,15};for(int a=0;a<16;a++)for(int i=0;i<16;i++)F[a]|=(pc(a&(ty==2?cols[i]:ty==1?i:(8+(i&7))))&1)<<i;for(int b=0;b<65536;b++){int d=20,r=65536;for(int u:F){r=min(r,b^u);d=min(d,pc(b^u));}ds[b]=d;re[b]=r;if(b==r){rs.push_back(b);if(d>=7&&(pc(b)&1))qual.push_back(b);}}
for(int i=-5;i<=5;i++)for(int j=-5;j<=5;j++)for(int k=-5;k<=5;k++){int n=(i+j)*(i+j)+(i+k)*(i+k)+(j+k)*(j+k);if(n>0&&n<=12&&((i>0)||(i==0&&j>0)||(i==0&&j==0&&k>0)))outer.push_back({i,j,k,n});}
sort(outer.begin(),outer.end(),[](Z a,Z b){return a.n<b.n;});V good(65536);
for(int A=0;A<65536;A++)if(pc(A)==gg){V bank;for(int b:rs){int d=20;for(int u:F)d=min(d,pc((b^u)&~A));good[b]=4*d+gg>=27;if(good[b])bank.push_back(b);}if(bank.empty())continue;geometry++;
for(int u:qual)for(int v:qual)if(v>u&&ds[u^v]>=6){rects++;for(int b:bank)if(good[re[b^u]]&&good[re[b^v]]&&good[re[b^v^u]]){parallelograms++;
// Costs where the binary part is not U are sheet-independent. Only eight choices of central-bit generators remain.
for(int sht=0;sht<8;sht++){bool ok=true;fulltests++;for(auto [i,j,k,n]:outer){int vec[16];for(int l=0;l<16;l++)vec[l]=2*i*((u>>l)&1)+(j+k)*(((A>>l)&1)+2*((b>>l)&1))+2*k*((v>>l)&1);vec[0]+=4*(i*(sht&1)+(j+k)*((sht>>1)&1)+k*(sht>>2));int best=999;
for(int w:F){int rad=0,kp=0,fix=16;for(int l=0;l<16;l++){int x=vec[l]-2*((w>>l)&1),r=((x+2)&3)-2;rad+=r*r;kp+=(x-r)/4;fix=min(fix,16-8*abs(r));}best=min(best,rad+(kp&1)*fix);}if(2*best+5*n<64){ok=false;break;}}
if(ok){ofstream f("fccWIN.txt");f<<ty<<" "<<A<<" "<<u<<" "<<v<<" "<<b<<" "<<sht<<"\n";cout<<"WIN "<<ty<<" "<<A<<" "<<u<<" "<<v<<" "<<b<<" "<<sht<<endl;return 0;}}}}}
cout<<"FCC2x2 type="<<ty<<" m="<<gg<<" A="<<geometry<<" rect="<<rects<<" shells="<<parallelograms<<" full="<<fulltests<<endl;}
