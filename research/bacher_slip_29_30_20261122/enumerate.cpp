#include <bits/stdc++.h>
using namespace std;
using Z=__int128_t;
constexpr int d=28;constexpr long scale=1u<<20;
long den,r[d],v[d],L[d][d],G[d][d],P[d],Q[d];
long long nodes;long leaves,mn,bd;
bool zero;
long fl(long a,long b){return a>=0?a/b:-(b-1-a)/b;}
long ce(long a,long b){return -fl(-a,b);}
void go(int i,long B){
 ++nodes;
 if(i<0){
  if(zero&&all_of(v,v+d,[](long x){return !x;}))return;
  long a=0;for(int j=0;j<d;j++)for(int k=0;k<d;k++)a+=v[j]*G[j][k]*v[k];
  if(a<=bd){mn=min(mn,a);++leaves;}
  return;
 }
 long t=0,k=L[i][i];
 for(int j=i+1;j<d;j++)t+=L[i][j]*v[j];
 Z z=(Z(B+1)*Q[i]-1)/(Z(scale)*P[i]);
 long b=sqrtl((long double)z);while(Z(b+1)*(b+1)<=z)++b;while(Z(b)*b>z)--b;
 long lo=ce(-b-t+k*r[i],den*k),hi=fl(b-t+k*r[i],den*k);
 for(long j=lo;j<=hi;j++){
  v[i]=den*j-r[i];assert(abs(v[i])<100000);
  long a=k*v[i]+t;long nb=B-Z(scale)*P[i]*Z(a)*a/Q[i];assert(nb>=0);
  go(i-1,nb);
 }
}
int main(int c,char**argv){
 ifstream f("ldl_fast.txt");for(int i=0;i<d;i++)f>>P[i]>>Q[i];for(auto&l:L)for(auto&i:l)f>>i;
 for(auto&l:G)for(auto&i:l)f>>i;
 int id;
 while(cin>>id>>den>>bd){
  assert(den>=1&&den<=24&&bd<4608&&bd>=0);
  for(int i=0;i<d;i++){cin>>r[i];assert(0<=r[i]&&r[i]<den);}
  zero=id==0;nodes=leaves=0;mn=bd+1;
  go(d-1,bd*scale);
  cout<<id<<" nodes "<<nodes<<" leaves "<<leaves<<" min_numerator "<<mn<<endl;
 }
}
