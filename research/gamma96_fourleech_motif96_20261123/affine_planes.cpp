#include<bits/stdc++.h>
using namespace std;using U=uint32_t;
int main(){
 U F[24];for(auto &v:F)cin>>v;
 int n;cin>>n;assert(n==4095);
 unordered_map<U,vector<U>> S,E;vector<U>A;
 unordered_map<U,U> D;
 map<int,int> rc;
 for(int i=0;i<n;i++){
  U a;cin>>a;A.push_back(a);
  auto &s=S[a];for(int j=0;j<24;j++){U v;cin>>v;s.push_back(v);assert(v);}
  sort(s.begin(),s.end());assert(adjacent_find(s.begin(),s.end())==s.end());
  U lhs=0;for(U k=a;k;k&=k-1)lhs^=F[countr_zero(k)];D[a]=lhs;
  U P[24]={};for(U w:s)for(U v=w^s[0];v;v^=P[bit_width(v)-1]){
    U &u=P[bit_width(v)-1];if(!u){u=v;break;}
  }
  for(U v:P)if(v)E[a].push_back(v);
  rc[E[a].size()]++;
 }
 for(auto [a,b]:rc)cerr<<"rank "<<a<<" "<<b<<endl;
 sort(A.begin(),A.end());
 auto orth=[&](U a,U b){U z=D.at(a);for(U v:E.at(b))if(popcount(z&v)%2)return false;return true;};
 auto chk=[&](U a,U b,U c){
  auto &u=S.at(a),&v=S.at(b),&w=S.at(c);
  for(U x:u)for(U y:v)if(binary_search(u.begin(),u.end(),x^y)||binary_search(v.begin(),v.end(),x^y)||binary_search(w.begin(),w.end(),x^y))return false;
  return true;
 };
 int yes=0,all=0;
 ofstream out("unused_good_affine.txt"), structout("admissible_affine.txt");
 for(U a:A){
  for(U b:A)if(b>a&&(a^b)>b&&orth(a,b)&&orth(b,a)&&orth(a,a^b)&&orth(a^b,a)&&orth(b,a^b)&&orth(a^b,b)){
   U c=a^b;
   auto parity=[&](U x,U y){return popcount(D.at(x)&S.at(y)[0])%2;};
   int da=parity(a,b),db=parity(b,a),dc=parity(c,a);
   if(da!=parity(a,c)||db!=parity(b,c)||dc!=parity(c,b)||(da^db)!=dc)continue;
   all++;structout<<a<<" "<<b<<" "<<c<<"\n";
   yes+=(da==0&&db==0&&dc==0);
  }
 }
 cerr<<"admissible planes "<<all<<" homogeneous "<<yes<<endl;
}
