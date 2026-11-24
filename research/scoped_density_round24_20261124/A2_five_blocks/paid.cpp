#include <bitset>
#include <iostream>
#include <vector>
#include <fstream>
#include <array>
#include <algorithm>
#include <cstdint>
using namespace std;
constexpr int N=4096;
using B=bitset<N>;using W=array<int,5>;using V=vector<int>;
array<B,N>E;
vector<W>A;
V D; vector<W>C0;
bool ok(const W&a,const W&b){
 int z=0;for(int i=0;i<5;i++){bool p=a[i]/3!=b[i]/3;bool q=a[i]%3!=b[i]%3;z+=p?(q?6:18):(q?24:0);}return z>=72;
}
int removed,n,m;long long ct,limit=500000;int best=18;int tbest=-1;
int nice(const V&V0){
 int n=0;for(int a:V0)for(int z:A[D[a]])if(z/3){n++;break;}return n;
}
void save(const V&V0,const char* f){
 ofstream o(f);
 for(int i=0;i<18;i++)if(!(m>>i&1)){for(int c:C0[i])o<<c<<" ";o<<"\n";}
 for(int a:V0){for(int c:A[D[a]])o<<c<<" ";o<<"\n";}
}
void dfs(B P,V& C){
 if(++ct>limit)return;
 if(P.none()){
  int t=C.size()+18-removed;
  if(t>best){best=t;cout<<"SIZE "<<t<<" drop "<<removed<<" in mask "<<m<<" nodes "<<ct<<endl;save(C,"best-paid.txt");}
  if(t>=18 && nice(C)>tbest){tbest=nice(C);cout<<"TIE "<<t<<" mixed "<<tbest<<" drop "<<removed<<" mask "<<m<<endl;save(C,"mixed-paid.txt");}
  return;
 }
 V u,k;B S=P;
 for(int b=1;S.any();b++){
  B T=S;while(T.any()){int v=T._Find_first();T.reset(v);S.reset(v);T&=~E[v];u.push_back(v);k.push_back(b);}
 }
 for(int i=u.size()-1;i>=0;i--){
  if(int(C.size())+k[i]+18-removed<best||ct>limit)return;
  int v=u[i];C.push_back(v);dfs(P&E[v],C);C.pop_back();P.reset(v);
  if(best>=20)return;
 }
}
int main(){
 ifstream h("ter5.txt");W w;
 while(h>>w[0]>>w[1]>>w[2]>>w[3]>>w[4])C0.push_back(w);
 for(int i=0;i<248832;i++){int j=i;W X;for(int t=4;t>=0;t--)X[t]=j%12,j/=12;A.push_back(X);}
 ifstream f("donors.txt");int rounds=0;
 while(f>>removed>>n>>m){
  D.resize(n);for(int &i:D)f>>i;
  if(n>N)continue;
  if(removed>12)limit=80000;else if(removed>10)limit=200000;
  vector<pair<int,int>> K;
  for(int i:D){int z=0;for(int j:D)z+=ok(A[i],A[j]);K.push_back({-z,i});}
  sort(K.begin(),K.end());for(int i=0;i<n;i++)D[i]=K[i].second;
  for(int i=0;i<n;i++){E[i].reset();for(int j=0;j<n;j++)if(ok(A[D[i]],A[D[j]]))E[i].set(j);}
  B P;for(int i=0;i<n;i++)P[i]=1;
  ct=0; V Q;dfs(P,Q);
  rounds++;
  cout<<"chart "<<rounds<<" "<<removed<<" size "<<n<<" nodes "<<ct<<" best "<<best<<" mixed "<<tbest<<endl;
  if(best>=20)break;
 }
}
