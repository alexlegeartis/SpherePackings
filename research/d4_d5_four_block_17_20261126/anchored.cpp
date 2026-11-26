#include <algorithm>
#include <array>
#include <vector>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
using namespace std;
int bas[8]={31181,58116,71012,35581,40360,35279,18,258}, G[17],gr[131072],rc[17],prows[8]; uint64_t ok[512][4];bool good[4096];int how[20],maxsize=8,goals=12,done=0;long long nodes=0;vector<int> cur,solution;
int dep4(int x){return x<4?(x?2:0):1;}
int dep(int x){int p=x>>12,w=__builtin_popcount(p&15);w+=w%2;return dep4(x&15)+dep4((x>>4)&15)+dep4((x>>8)&15)+(w?w/2:!!p*2);}
bool addbasis(int x){for(int j=16;j>=0;j--)if((x>>j)&1){if(rc[j])x^=rc[j];else{rc[j]=x;return true;}}return false;}
bool dfs(const vector<int>&a,int l){
 if(l==0){solution=cur;return true;}if(a.size()<l)return false; nodes++;
 // strict proper greedy coloring in the compatibility graph, last color upper
 vector<int>v=a,order,colors,t;int color=0;
 while(!v.empty()){color++;vector<int> col; t.clear();for(int x:v){bool can=true;for(int y:col)if(good[x^y]){can=false;break;}if(can){col.push_back(x);order.push_back(x);colors.push_back(color);}else t.push_back(x);}v.swap(t);}
 for(int j=order.size()-1;j>=0;j--){if(colors[j]<l)return false;int x=order[j];t.clear();for(int k=0;k<j;k++)if(good[x^order[k]])t.push_back(order[k]);cur.push_back(x);if(dfs(t,l-1))return true;cur.pop_back();}
 return false;
}
void scan(int mask,const int *kernel){ // rows in RREF over binary coefficient space of B
 done++;uint64_t cos[8][4]={};int rep[8];int n=0;for(int x=0;x<256;x++)if(!(x&mask))rep[n++]=x;
 int elem[32]={};for(int j=0;j<5;j++)for(int i=0;i<(1<<j);i++)elem[i+(1<<j)]=elem[i]^kernel[j];
 for(int j=0;j<8;j++)for(int r:elem){int x=rep[j]^r;cos[j][x>>6]|=1ULL<<(x%64);}
 vector<int>A;memset(good,0,sizeof(good));
 for(int u=0;u<512;u++)for(int j=0;j<8;j++){bool pass=true;for(int k=0;k<4;k++)if(cos[j][k]&~ok[u][k]){pass=false;break;}if(pass){A.push_back(u*8+j);good[u*8+j]=true;}}
 // contains all nonzero inside B, so at least 7 options. Solve needed sizes only.
 for(int sz=goals;sz>maxsize;sz--){cur.clear();if(dfs(A,sz-1)){maxsize=sz;printf("FOUND after %d subspaces: %d kernels",done,sz);for(int k=0;k<5;k++)printf(" %d",gr[kernel[k]]);printf(" representatives 0");for(int x:solution){int u=x>>3,j=x&7;printf(" %d",gr[(u<<8)^rep[j]]);}puts("");fflush(stdout);break;}}
 if(done%5000==0){fprintf(stderr,"%d spaces best %d nodes %lld\n",done,maxsize,nodes);fflush(stderr);}
}
int main(int argc,char**argv){if(argc>1)goals=atoi(argv[1]);int j=0;for(int x:bas){if(!addbasis(x))return 1;G[j++]=x;}for(int k=0;k<17;k++)if(addbasis(1<<k))G[j++]=1<<k;
 for(int i=1;i<131072;i++)gr[i]=gr[i&(i-1)]^G[__builtin_ctz(i)];
 for(int i=1;i<131072;i++)if(dep(gr[i])>=4)ok[i>>8][(i>>6)&3]|=1ULL<<(i&63);
 for(int mask=0;mask<256;mask++)if(__builtin_popcount(mask)==5){
 int piv[5],np=0;for(int c=0;c<8;c++)if((mask>>c)&1)piv[np++]=c;
 vector<pair<int,int>> free;for(int r=0;r<5;r++)for(int c=piv[r]+1;c<8;c++)if(!((mask>>c)&1))free.emplace_back(r,c);
 for(int a=0;a<(1<<(int)free.size());a++){int kernel[5];for(int r=0;r<5;r++)kernel[r]=1<<piv[r];for(size_t k=0;k<free.size();k++)if((a>>k)&1)kernel[free[k].first]^=1<<free[k].second;scan(mask,kernel);if(maxsize>=goals)return 0;}}
 printf("ALL %d kernels checked; max %d, nodes %lld\n",done,maxsize,nodes);}

