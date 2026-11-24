#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
using bits=vector<uint64_t>;
const int N=531441; using pt=array<int,12>;
pt V[N]; bool good[N]; int S[12];
int diff(int x,int y) {
   int r=0; for(int j=0;j<12;j++) { int a=V[x][j]-V[y][j]; if(a<0) a+=3; r+=a*S[j]; }
   return r;
}
vector<bits> adj;
int n,m; unsigned long long nodes;
bool remove(bits& P,int x) { bool h=P[x>>6]>>(x&63)&1; P[x>>6]&=~(1ull<<(x&63)); return h; }
bool any(const bits &P) { for(auto x:P) if(x) return 1; return 0; }
bool search(bits P,int k) {
    nodes++;
    if(k==0) return true;
    vector<int> order,colors;
    int color=0;
    bits W=P;
    while(any(W)) {
        color++;
        bits A=W;
        for(int i=0;i<m;i++) while(A[i]) {
            int w=__builtin_ctzll(A[i])+64*i;
            order.push_back(w); colors.push_back(color);
            remove(W,w);remove(A,w);
            for(int j=i;j<m;j++) A[j]&=~adj[w][j];
        }
    }
    for(int p=(int)order.size()-1;p>=0;p--) {
        if(colors[p]<k) return false;
        int v=order[p];bits R=P;
        for(int i=0;i<m;i++) R[i]&=adj[v][i];
        if(search(R,k-1)) return true;
        remove(P,v);
    }
    return false;
}
int main() {
   int t,a,b;cin>>t>>a>>b;
   vector<int> L(a),E(b);for(auto& x:L) {cin>>x;good[x]=1;} for(auto& x:E) cin>>x;
   S[0]=1;for(int j=1;j<12;j++)S[j]=3*S[j-1];
   for(int i=0;i<N;i++)for(int j=0;j<12;j++) V[i][j]=i/S[j]%3;
   for(int e:E) {
       vector<int> Q; for(int x:L) if(good[diff(x,e)]) Q.push_back(x);
       vector<pair<int,int>> rank;
       for(int x:Q) {
           int deg=0;for(int y:Q) deg+=good[diff(x,y)];
           rank.push_back({-deg,x});
       }
       sort(rank.begin(),rank.end());Q.clear();for(auto [r,x]:rank) Q.push_back(x);
       n=Q.size();m=(n+63)/64;adj.assign(n,bits(m));
       for(int i=0;i<n;i++)for(int j=0;j<i;j++) if(good[diff(Q[i],Q[j])]) {
           adj[i][j>>6]|=1ull<<(j&63);adj[j][i>>6]|=1ull<<(i&63);
       }
       bits P(m,~0ull);if(n%64)P.back()=(1ull<<(n%64))-1;
       nodes=0;
       auto found=search(P,t-2);
       cout<<e<<" "<<n<<" "<<nodes<<" "<<found<<endl;
       if(found) return 1;
   }
}
