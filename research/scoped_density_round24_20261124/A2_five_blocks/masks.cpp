#include <bitset>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
using namespace std;using W=array<int,5>;using P=array<int,18>;
int main(){
 ifstream f("ter5.txt");vector<W>C; W row;
 while(f>>row[0]>>row[1]>>row[2]>>row[3]>>row[4])C.push_back(row);
 int ids[243];fill(ids,ids+243,-1);
 auto idx=[](W a){int i=0;for(auto x:a)i=3*i+x;return i;};
 for(int j=0;j<18;j++)ids[idx(C[j])]=j;
 vector<array<int,3>>S;array<int,3>ps={0,1,2};
 do{S.push_back(ps);}while(next_permutation(ps.begin(),ps.end()));
 vector<P>A;W w={0,1,2,3,4};
 do{
  for(int h=0;h<7776;h++){
   int k=h;array<int,5>t;for(int&i:t)i=k%6,k/=6;
   P g;bool ok=true;
   for(int j=0;j<18;j++){
    int i=0;for(int e=0;e<5;e++)i=3*i+S[t[e]][C[j][w[e]]];
    g[j]=ids[i];if(g[j]<0){ok=false;break;}
   }
   if(ok)A.push_back(g);
  }
 }while(next_permutation(w.begin(),w.end()));
 cerr<<"Aut "<<A.size()<<endl;
 vector<bool> seen(1<<18);
 for(int m=0;m<(1<<18);m++)if(!seen[m]&&__builtin_popcount(m)>=6&&__builtin_popcount(m)<=14){
   for(auto&g:A){int z=0;for(int j=0;j<18;j++)if(m>>j&1)z|=1<<g[j];seen[z]=true;}
   cout<<m<<"\n";
 }
}
