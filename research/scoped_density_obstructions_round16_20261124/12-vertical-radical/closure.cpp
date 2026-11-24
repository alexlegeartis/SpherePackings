#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cassert>
using namespace std;
using W = uint64_t;
int n,h,f,goal,passall;
W S, vacancies, full, fixedholes;
unsigned long long total=0;
vector<vector<W>> tests;
vector<int> ext;
inline bool forced(int v,W known) {
 for(W p:tests[v]) if((p&known)==p) return true;
 return false;
}
void check(W cut) {
 W alive=full&~(cut|fixedholes), known=S;
 if(passall)known=alive;
 else {
  bool flag;
  do {
   W old=known;
   W candidates=alive&~known;
   while(candidates) {
    W b=candidates&-candidates; candidates^=b;
    if(forced(__builtin_ctzll(b),known)) known|=b;
   }
   flag=(old!=known);
  } while(flag);
 }
 int count=__builtin_popcountll(known);
 for(int v:ext) {
  if(forced(v,known)) ++count;
  if(count>=goal)return;
 }
 cerr<<"counterexample mask "<<cut<<", in-shell closure "<<__builtin_popcountll(known)<<", total "<<count<<endl;
 exit(1);
}
void rec(W avail,int left,W cut) {
 if(!left) {++total; check(cut); return;}
 while(__builtin_popcountll(avail)>=left) {
  W b=avail&-avail;avail^=b;
  rec(avail,left-1,cut|b);
 }
}
int main() {
 cin>>n>>h>>f>>goal>>passall>>S>>vacancies>>fixedholes;
 assert(n==243&&h<=11&&f==36&&__builtin_popcountll(vacancies)>=h);
 full=(1ull<<f)-1;
 tests.resize(n);
 for(int i=0;i<n;i++) {
  int k;cin>>k;tests[i].resize(k);
  for(auto& a:tests[i]) {cin>>a; assert(__builtin_popcountll(a)==4);}
  if(i>=f&&k)ext.push_back(i);
 }
 sort(ext.begin(),ext.end(),[](int a,int b){return tests[a].size()>tests[b].size();});
 if(passall) {
  assert(!fixedholes);
  for(int i=0;i<f;i++)if(!(S>>i&1))assert(forced(i,S));
 }
 rec(vacancies,h,0);
 cout<<"PASS "<<total<<" deletion patterns, guaranteed "<<goal<<" potential positions"<<endl;
}
