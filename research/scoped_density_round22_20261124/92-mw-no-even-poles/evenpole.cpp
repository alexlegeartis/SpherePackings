#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <iostream>
using namespace std;
constexpr int p=47;
int mod(int x){x%=p;return(x>=0?x:x+p);}
int powm(int a,int n){int y=1;while(n){if(n&1)y=y*a%p;a=a*a%p;n>>=1;}return y;}
int main(){
 array<int,23> roots{};int nr=0;for(int x=1;x<p;++x)if(powm(x,23)==1)roots[nr++]=x;
 assert(nr==23);
 array<int,25> inv{};for(int n=1;n<=24;++n){inv[n]=powm(3*n,p-2);assert(mod(inv[n]*3*n)==1);}
 array<uint64_t,26> stop{};uint64_t all=0;
 for(uint32_t mask=0;mask<(1u<<23);++mask){
  // Seven roots in A. The complement makes B after removing s^3.
  if(popcount(mask)!=7)continue;
  ++all;array<int,17>B{};B[0]=1;int d=0;
  for(int j=0;j<23;++j)if(!((mask>>j)&1u)){
   int r=roots[j];++d;for(int i=d;i;--i)B[i]=mod(B[i]-r*B[i-1]);
  }
  assert(d==16);
  array<int,25>C{};C[0]=1;
  int n=1;for(;n<=24;++n){
   int z=0;for(int i=1;i<=d&&i<=n;++i)z+=(5*i-3*n)*B[i]*C[n-i];
   C[n]=mod(inv[n]*mod(z));if(n>=14&&C[n])break;
  }
  ++stop[n];
 }
 assert(all==245157);
 assert(stop[14]==239476 && stop[15]==5612 && stop[16]==69);
 for(int n=14;n<=25;++n)cout<<n<<" "<<stop[n]<<'\n';
 assert(stop[25]==0);
}
