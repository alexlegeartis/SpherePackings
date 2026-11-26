// Enumerate the double cosets of the coordinate permutations modulo the affine Best code symmetries.
// The proof itself verifies the generators against the exact 40 words.
#include <algorithm>
#include <array>
#include <vector>
#include <map>
#include <iostream>
#include <numeric>
#include <cassert>
#include <cstdint>
using namespace std; using P=array<int,10>;using H=array<bool,1024>;
int pc(int m){return __builtin_popcount((unsigned)m);}
int perm(int m,const P&p){int z=0;for(int i=0;i<10;i++)if((m>>i)&1)z|=1<<p[i];return z;}
P comp(P p,P q){P z;for(int i=0;i<10;i++)z[i]=p[q[i]];return z;}
int rankp(P p){int r=0,used=0;for(int i=0;i<10;i++){r=r*(10-i)+p[i]-pc(used&((1<<p[i])-1));used|=1<<p[i];}return r;}
bool dis(const vector<int>& s,int i,int k,int used){
 if(k==0)return true;for(;i+k<=(int)s.size();i++)if(!(s[i]&used)&&dis(s,i+1,k-1,used|s[i]))return true;return false;
}
int main(){int G[]={0,1,3,2};vector<int>C;H Cs{};for(string v:{"01112","03110","21310","21132"})for(int k=0;k<5;k++)for(int s:{1,3}){int w=0;for(int i=0;i<5;i++)w|=G[((v[(k+i)%5]-48)*s)%4]<<(2*i);if(!Cs[w])C.push_back(w);Cs[w]=1;}assert(C.size()==40);
 vector<P> aut;
 // Supplied candidate group: permutations of five pairs in D5, and arbitrary within-pair swaps; retain only checked affine automorphisms.
 for(int s:{1,-1})for(int k=0;k<5;k++)for(int flips=0;flips<32;flips++){P p;for(int i=0;i<5;i++)for(int j=0;j<2;j++)p[2*i+j]=2*((10+k+s*i)%5)+(j^((flips>>i)&1));int a=perm(C[0],p);
 for(int y:C){int x=a^y;bool ok=true;for(int v:C)if(!Cs[perm(v,p)^x]){ok=false;break;}if(ok){aut.push_back(p);break;}}}
 P p;iota(p.begin(),p.end(),0);vector<bool> checked(3628800,false);int classes=0,nperm=0,nontriv=0,mind=1024;map<int,int> sizes;vector<int> supports(1024);iota(supports.begin(),supports.end(),0);sort(supports.begin(),supports.end(),[](int x,int y){return pc(x)==pc(y)?x<y:pc(x)<pc(y);});
 do{int ri=rankp(p);if(checked[ri])continue;classes++;int ne=0;for(auto&a:aut){auto t=comp(a,p);for(auto&b:aut){int z=rankp(comp(t,b));if(!checked[z]){checked[z]=1;ne++;}}}
 nperm+=ne;H d{};for(int x:C)for(int y:C)d[x^perm(y,p)]=true;int ds=count(d.begin(),d.end(),true);mind=min(mind,ds);sizes[ds]++;
 for(int a=0;a<1024;a++)if(!d[a]){vector<int>S;for(int x:supports)if(pc(x)<=4&&d[a^x])S.push_back(x);if(!dis(S,0,4,0)){cerr<<"missing "<<ri<<","<<a<<endl;return 1;}nontriv++;}
 }while(next_permutation(p.begin(),p.end()));assert(nperm==3628800);cout<<"aut="<<aut.size()<<" double-coset reps="<<classes<<" perms="<<nperm<<" nontrivial corner proofs="<<nontriv<<" min difference population="<<mind<<endl;}

