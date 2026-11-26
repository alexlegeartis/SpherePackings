// Additional higher-block half-integral isometries. Single 4x4 triality component handled separately.
#include <algorithm>
#include <array>
#include <vector>
#include <iostream>
#include <numeric>
#include <cassert>
#include <functional>
#include <set>
using namespace std;
int pc(int x){return __builtin_popcount((unsigned)x);}
int main(){
 int g[]={0,1,3,2};vector<int>C;array<bool,1024>cs{},B{};
 for(string v:{"01112","03110","21310","21132"})for(int k=0;k<5;k++)for(int s:{1,3}){int w=0;for(int i=0;i<5;i++)w|=g[(v[(k+i)%5]-48)*s%4]<<(2*i);if(!cs[w])C.push_back(w);cs[w]=1;}assert(C.size()==40);
 for(int a=0;a<1024;a++){int d=10;for(int c:C)d=min(d,pc(c^a));B[a]=d==3;}
 long long mats=0,tests=0,sol=0;int graphs=0;
 for(int k:{3,4,5}){
 vector<int> choices=k==3?vector<int>{7,11}:(k==4?vector<int>{15}:vector<int>{31});
 for(int si:choices)for(int so:choices){
 vector<int>I,O,I0,O0;for(int i=0;i<10;i++){((si>>(i/2))&1?I:I0).push_back(i);((so>>(i/2))&1?O:O0).push_back(i);}
 vector<int>A(k),deg(k);
 function<void(int)> inc=[&](int c){
 if(c<k){for(int u=0;u<k;u++)for(int v=u+1;v<k;v++)if(deg[u]<2&&deg[v]<2){A[c]=(1<<u)|(1<<v);deg[u]++;deg[v]++;inc(c+1);deg[u]--;deg[v]--;}return;}
 graphs++;
 vector<int>H={0};for(int h:A){int v=0;for(int r=0;r<k;r++)if((h>>r)&1)v|=3<<(O[2*r]);if(find(H.begin(),H.end(),v)==H.end()){int L=H.size();for(int i=0;i<L;i++)H.push_back(H[i]^v);}}
 vector<int>good;array<bool,1024> goodset{};for(int a=0;a<1024;a++){bool ok=true;for(int v:H)if(!B[a^v]){ok=false;break;}if(ok){good.push_back(a);goodset[a]=true;}}
 if(good.size()<40)return;
 vector<vector<int>>W(2*k,vector<int>(2*k));vector<int>rsum(2*k);
 function<void(int,int)> rec=[&](int r,int used){
 if(r<2*k){int pos[4],p=0;for(int cc=0;cc<2*k;cc++)if(A[cc/2]>>(r/2)&1)pos[p++]=cc;assert(p==4);int mask=0;for(int i:pos)mask|=1<<i;
 int free=mask&used;if(free)free&=free-1; // gauge: first previously fixed column positive; all previously unseen columns positive
 int u=free;while(true){for(int i:pos)W[r][i]=((u>>i)&1)?-1:1;
 bool ok=true;for(int q=0;q<r;q++){int v=0;for(int i:pos)v+=W[r][i]*W[q][i];if(v){ok=false;break;}}
 if(ok)rec(r+1,used|mask);
 if(!u)break;u=(u-1)&free;}return;}
 mats++;
 int low[40],par[40];for(int j=0;j<40;j++){low[j]=par[j]=0;for(int ri=0;ri<2*k;ri++){int z=0;for(int cc=0;cc<2*k;cc++)z+=W[ri][cc]*((C[j]>>I[cc])&1);int fl=z>=0?z/2:-((-z+1)/2);low[j]|=(fl&1)<<O[ri];par[j]|=(z&1)<<O[ri];}}
 do{int f0[40];for(int j=0;j<40;j++){f0[j]=low[j];for(unsigned i=0;i<I0.size();i++)f0[j]|=((C[j]>>I0[i])&1)<<O0[i];}
 for(int e=0;e<(1<<(2*k));e++){int E=0;for(int i=0;i<2*k;i++)E|=((e>>i)&1)<<O[i];int v0=f0[0]^(par[0]&E);
 for(int a0:good){int b=a0^v0;bool ok=true;for(int j=1;j<40;j++)if(!goodset[b^f0[j]^(par[j]&E)]){ok=false;break;}if(ok){sol++;if(sol==1){cout<<"possible k="<<k<<" src="<<si<<" dst="<<so<<" eps="<<e<<" beta="<<b<<endl;}}}
 tests++;
 }}while(next_permutation(O0.begin(),O0.end()));};
 rec(0,0);
 };inc(0);}
 cout<<"up to k="<<k<<" matrices="<<mats<<" tests="<<tests<<" graphs="<<graphs<<" survivors="<<sol<<endl;}
 return sol?1:0;
}
