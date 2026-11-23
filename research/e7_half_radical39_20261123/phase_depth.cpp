#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
using namespace std;
const int bound[7]={12,7,9,5,10,7,4};
const int pairs[6][2]={{0,1},{0,2},{2,3},{0,4},{4,5},{5,6}};
int frac[128],v[7],M[128],alpha[5],checked=0;
int mul(int x,int y){
 int z=0;while(y){if(y&1)z^=x;x<<=1;if(x&64)x^=67;y>>=1;}return z;
}
int pw(int x,int j){
 int k=1;while(j){if(j&1)k=mul(k,x);x=mul(x,x);j>>=1;}return k;
}
int q(int x){int n=pw(x,9);return n^pw(n,2)^pw(n,4);}
void gen(int i){
 if(i==7){
   int w=0,m=0;
   for(int j=0;j<7;j++){w+=2*v[j]*v[j];m|=(((v[j]-((82>>j)&1))/2)&1)<<j;}
   for(auto& z:pairs)w-=2*v[z[0]]*v[z[1]];
   if(w<=14){assert(w==6||w==14);++checked;frac[m]=min(frac[m],w);}
   return;
 }
 for(int x=-bound[i];x<=bound[i];x++)if((x&1)==(82>>i&1)){v[i]=x;gen(i+1);}
}
int main(){
 array<bool,128> seen{};
 for(int &x:M){cin>>x;assert(cin&&x>=0&&x<128&&!seen[x]);seen[x]=true;}
 for(int &x:alpha){cin>>x;assert(cin&&x>=0&&x<64&&pw(x,8)==x);}
 assert(M[0]==0&&M[64]==82);
 fill(frac,frac+128,30);gen(0);
 assert(checked==632);
 array<int,2> cc{};
 for(int j=0;j<128;j++){
   int a=j&63, q0=q(a^60);assert(q0==0||q0==1);
   assert(frac[M[j]]==6+8*q0);cc[q0]++;
 }
 assert(cc[0]==56&&cc[1]==72);
 int K[5][128],C[2][128];
 for(int j=0;j<128;j++){
   int x=j&63,e=j>>6,q0=q(x);
   C[1][j]=frac[M[j]];
   C[0][j]=x?16-8*(q0^e):24*e;
   for(int i=0;i<5;i++)K[i][j]=x|(mul(alpha[i],x)<<6)|(e<<12);
 }
 array<int,6> maxima{};
 for(int m=0;m<32;m++){
  array<int,8192>D,G;
  D.fill(1000);D[0]=0;
  for(int i=0;i<5;i++){
   G.fill(1000);for(int s=0;s<8192;s++)if(D[s]<1000)for(int j=0;j<128;j++){
    int k=s^K[i][j];G[k]=min(G[k],D[s]+C[m>>i&1][j]);}
   D=G;
  }
  auto w=__builtin_popcount(unsigned(m));
  auto h=*max_element(D.begin(),D.end());maxima[w]=max(maxima[w],h);
  if(w==3)assert(count(D.begin(),D.end(),26)==72&&count(D.begin(),D.end(),18)==8120);
  cout<<"pattern "<<m<<" maximum quarter-units "<<h<<"\n";
 }
 assert((maxima==array<int,6>{24,22,20,26,24,30}));
 cout<<"PASS all 262144 half-radical/integral classes\n";
}
