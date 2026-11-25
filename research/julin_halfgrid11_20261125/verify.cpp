// C++17, exact verification for 11D Julin/Steiner half-grid extensions.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>
using namespace std;
using Word=uint32_t; using Bits=vector<uint64_t>;
static int pc(Word t){return __builtin_popcount(t);}
static Word rot(Word t){return (t>>1)|((t&1)<<10);}
static vector<Word> cycle(Word t){vector<Word> v={t};for(Word u=rot(t);u!=t;u=rot(u))v.push_back(u);return v;}
static int distance4(Word a,Word b){Word x=a^b,h=x&2047;return pc(h)+4*pc((x>>11)&~h);}
static int first(const Bits& b){for(unsigned j=0;j<b.size();++j)if(b[j])return int(64*j+__builtin_ctzll(b[j]));return -1;}
struct Checker{
 vector<Bits> E;int words;
 bool has(Bits S,int goal){
  int n=0;for(auto w:S)n+=__builtin_popcountll(w);if(n<goal)return false;if(goal==1)return true;
  Bits rem=S;vector<int> seq,col;int c=0;
  while(first(rem)>=0){Bits p=rem;++c;for(int v;(v=first(p))>=0;){
   rem[v/64]&=~(1ull<<(v%64));p[v/64]&=~(1ull<<(v%64));for(int j=0;j<words;++j)p[j]&=~E[v][j];
   seq.push_back(v);col.push_back(c);}}
  for(int k=int(seq.size())-1;k>=0;--k){if(col[k]<goal)return false;int v=seq[k];Bits T=S;
   for(int j=0;j<words;++j)T[j]&=E[v][j];if(has(T,goal-1))return true;S[v/64]&=~(1ull<<(v%64));}
  return false;
 }
 bool test(vector<Word> v,int g){
  int n=v.size();vector<int> deg(n),idx(n);iota(idx.begin(),idx.end(),0);
  for(int i=0;i<n;i++)for(int j=0;j<i;j++)if(distance4(v[i],v[j])>=16){deg[i]++;deg[j]++;}
  stable_sort(idx.begin(),idx.end(),[&](int i,int j){return deg[i]>deg[j];});
  words=(n+63)/64;E.assign(n,Bits(words));
  for(int i=0;i<n;++i)for(int j=0;j<i;++j)if(distance4(v[idx[i]],v[idx[j]])>=16){E[i][j/64]|=1ull<<(j%64);E[j][i/64]|=1ull<<(i%64);}
  Bits all(words,~0ull);if(n%64)all.back()=(1ull<<(n%64))-1; return has(all,g);
 }
};
int main(){
 Word seed=0;for(int j:{0,1,3,4,5,9})seed|=1u<<j;auto rows=cycle(seed);rows.push_back(0);set<Word> reps,base;
 for(int i=0;i<12;i++)for(int j=0;j<i;j++){Word x=rows[i]^rows[j];auto o=cycle(x);reps.insert(*min_element(o.begin(),o.end()));base.insert(x);}
 assert(base.size()==66&&reps.size()==6);vector<vector<Word>> orbs;for(auto r:reps)orbs.push_back(cycle(r));
 vector<Word> old(base.begin(),base.end());for(int i=0;i<5;i++)old.push_back(3u<<(2*i));old.push_back(1023);
 for(int i=0;i<72;i++)for(int j=0;j<i;j++)assert(pc(old[i]^old[j])>=4);
 vector<vector<Word>> by(64);
 for(Word h=0;h<2048;++h){Word rest=2047^h;
  for(Word b=rest;;b=(b-1)&rest){int m=0;
   for(int o=0;o<6;++o){for(auto x:orbs[o])if(pc(h)+4*pc((b^x)&rest)<16){m|=1<<o;break;}if(pc(m)>3)break;}
   if(pc(m)<=3)for(Word t=h;;t=(t-1)&h){by[m].push_back(h|((b|t)<<11));if(t==0)break;}
   if(b==0)break;
  }
 }
 Checker ck;int tests=0;
 for(int m=0;m<64;++m)if(pc(m)==3){vector<Word> v;for(int s=0;s<64;s++)if(!(s&~m))v.insert(v.end(),by[s].begin(),by[s].end());
   sort(v.begin(),v.end());assert(adjacent_find(v.begin(),v.end())==v.end());bool exists=ck.test(v,40);
   cout<<"removed "<<m<<" candidates "<<v.size()<<" forty_additions "<<exists<<"\n";if(exists)return 1;tests++;}
 assert(tests==20);cout<<"All 20 three-orbit retentions certified; attainable size is 72.\n";
}
