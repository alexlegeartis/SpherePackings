#include <algorithm>
#include <array>
#include <vector>
#include <iostream>
#include <cassert>
using namespace std;
int pc(int x){return __builtin_popcount((unsigned)x);}
int main(){
 int g[]={0,1,3,2};vector<int>C;array<bool,1024>cs{},B{};
 for(string v:{"01112","03110","21310","21132"})for(int k=0;k<5;k++)for(int s:{1,3}){int w=0;for(int i=0;i<5;i++)w|=g[(v[(k+i)%5]-48)*s%4]<<(2*i);if(!cs[w])C.push_back(w);cs[w]=1;}assert(C.size()==40);
 for(int a=0;a<1024;a++){int d=10;for(int c:C)d=min(d,pc(c^a));B[a]=d==3;}assert(count(B.begin(),B.end(),true)==112);
 for(int v=1;v<1024;v++)if(pc(v)==4){int n=0;for(int a=0;a<1024;a++)n+=B[a]&&B[a^v];bool pairs=true;for(int i=0;i<5;i++)pairs&=(((v>>(2*i))&3)%3==0);assert((n>=40)==pairs);}
 // The affine pentacode group includes all dihedral permutations of the five pairs. Verify against C, using a translating offset if necessary.
 for(int s:{1,-1})for(int k=0;k<5;k++){vector<int>q;for(int c:C){int z=0;for(int i=0;i<5;i++)z|=((c>>(2*i))&3)<<(2*((10+k+s*i)%5));q.push_back(z);}bool ok=false;for(int d:C){bool yes=true;for(int x:q)yes&=cs[x^q[0]^d];ok|=yes;}assert(ok);}
 long tests=0,solutions=0;
 for(int src:{1,2})for(int dst:{1,2}){int ii[]={0,1,2*src,2*src+1},oo[]={0,1,2*dst,2*dst+1};vector<int>I,J;for(int i=0;i<10;i++){if(find(ii,ii+4,i)==ii+4)I.push_back(i);if(find(oo,oo+4,i)==oo+4)J.push_back(i);}
 const int v=(3|3<<(2*dst));vector<int>good;for(int a=0;a<1024;a++)if(B[a]&&B[a^v])good.push_back(a);
 int rows[]={0,1,2,3};
 do{do{for(int eps=0;eps<16;eps++){
  vector<int>z;for(int c:C){int f=0;for(int i=0;i<4;i++){int u=(eps>>i)&1;for(int j=0;j<4;j++)u+=((pc(rows[i]&j)%2)?-1:1)*(1&(c>>ii[j]));int fl=(u>=0?u/2:-((-u+1)/2));f|=(fl&1)<<oo[i];}for(int i=0;i<6;i++)f|=((c>>I[i])&1)<<J[i];z.push_back(f);}
  for(int t:good){bool ok=true;int b=t^z[0];for(int f:z)if(!(B[f^b]&&B[f^b^v])){ok=false;break;}if(ok){solutions++;if(solutions==1){cout<<"possible src="<<src<<" dst="<<dst<<" eps="<<eps<<" beta="<<b<<"\n";}}}
  tests++;
 }}while(next_permutation(J.begin(),J.end()));}while(next_permutation(rows+1,rows+4));
 }
 cout<<"candidate tests "<<tests<<" feasible roundings "<<solutions<<endl;
}
