#include<bits/stdc++.h>
using namespace std;using I=int64_t;using J=__int128_t;using U=uint64_t;using W=__uint128_t;constexpr int N=40,M=80;constexpr I S=1LL<<20;
struct Record{W p;int q;array<short,N> z;};
const U allParity=(1ULL<<N)-1;
I floorq(I a,I b){return a>=0?a/b:-((-a+b-1)/b);}
int last(W x){U a=x>>64;return a?127-__builtin_clzll(a):63-__builtin_clzll(U(x));}
struct Lat{
 int h,G[N][N]{},z[N]{}; I mu[N][N]{},d[N]{},E[N]{},cen[N]{},wide[N]{}; W basis[N]{},prot[M]{};U dec[M]{},nodes=0,maxn=0;
 vector<Record> pool;
 function<void(W,int,const array<short,N>&)> visit;
 int bound;bool hit=false,cap=false,collect=false;array<short,N>wit;int wq;
 Lat(int a,bool small=false):h(a){
  string nm=(small?"jointK":"fixed")+to_string(h)+".in";ifstream r(nm);
  for(auto&row:G)for(int&x:row)r>>x;
  for(int i=0;i<N;i++){
   W q=0;for(int j=0;j<M;j++){int t;r>>t;if(t&1)q|=W(1)<<j;}basis[i]=q;
  }
  for(I &x:d)r>>x;for(auto&row:mu)for(I&x:row)r>>x;for(I&x:wide)r>>x;assert(r);
  for(int i=0;i<N;i++){assert(d[i]>0&&wide[i]>=0&&wide[i]<1000);for(int j=0;j<N;j++)assert(abs(mu[i][j])<(1LL<<33));}
  for(int i=N-2;i>=0;i--)E[i]=E[i+1]+wide[i+1];
  for(int i=0;i<N;i++){
   W k=basis[i];U w=1ULL<<i;
   while(k){int j=last(k);if(prot[j]){k^=prot[j];w^=dec[j];}else{prot[j]=k;dec[j]=w;break;}}
   assert(k||small);
  }
 }
 U label(W r)const{
  U s=0;while(r){int j=last(r);if(!prot[j])return 1ULL<<63;r^=prot[j];s^=dec[j];}
  return s;
 }
 void go(int i,J b,bool started,U tar){
  if(hit||cap)return;
  if(++nodes>maxn){cap=true;return;}
  if(i<0){
   if(!started)return;
   array<short,N>r;long long q=0;
   for(int j=0;j<N;j++){assert(abs(z[j])<1000);r[j]=z[j];if(z[j])for(int k=0;k<N;k++)q+=1LL*z[j]*G[j][k]*z[k];}
   if(q>bound)return;
   assert(q>=14);
   if(collect){
    W p=0;for(int j=0;j<N;j++)if(z[j]&1)p^=basis[j];if(visit)visit(p,int(q),r);else pool.push_back({p,int(q),r});
   }else{hit=true;wq=q;wit=r;}
   return;
  }
  int step=1,want=0;if(tar<=allParity){step=2;want=(tar>>i)&1;}
  I lo=started?-wide[i]:0,hi=wide[i],cc=cen[i];
  lo+=((want-lo)%step+step)%step;hi-=((hi-want)%step+step)%step;
  if(lo>hi)return;
  I a=floorq(cc-S*want+S*(step==1?0:1),S*step)*step+want;
  if(step==1)a=floorq(cc+S/2,S);
  a=clamp(a,lo,hi);I dir=cc>S*a?1:-1,m=max(a-lo,hi-a)/step;
  for(I t=0;t<=2*m&&!hit&&!cap;t++){
   I u=a+(t&1?dir:-dir)*step*((t+1)/2);if(u<lo||u>hi)continue;
   I dis=max<I>(0,abs(cc-S*u)-E[i]);J cost=J(d[i])*dis*dis;
   if(cost>b)break;
   z[i]=u;
   if(u)for(int j=0;j<i;j++)cen[j]-=u*mu[i][j];
   go(i-1,b-cost,started||u,tar);
   if(u)for(int j=0;j<i;j++)cen[j]+=u*mu[i][j];
  }
 }
 struct TrieNode{int a[2]{-1,-1};int bound=0,src=-1;};
 vector<TrieNode> trie;
 const vector<Record>* other;
 int seen_id=-1;
 int congruent_cap(U v,int b)const{
  // For b in shared parity, its norm is constant mod8 in the entire b+2L class.
  long long q=0;
  for(int i=0;i<N;i++){
   long long t=0;for(int j=0;j<N;j++)if(v>>j&1)t+=G[i][j];
   assert(t%2==0);
   if(v>>i&1)q+=t;
  }
  assert(q%4==0);
  return b-int((b-q)%8+8)%8;
 }
 void trie_add(U a,int q,int src){
  int k=0;trie[k].bound=max(trie[k].bound,q);
  for(int i=N-1;i>=0;i--){
   bool b=(a>>i)&1;int t=trie[k].a[b];
   if(t<0){trie[k].a[b]=t=trie.size();trie.emplace_back();}
   k=t;trie[k].bound=max(trie[k].bound,q);
  }
  assert(trie[k].src<0);trie[k].src=src;
 }
 void walk(int i,J used,bool started,int t){
  if(hit||cap)return;if(++nodes>maxn){cap=true;return;}
  const J all=J(trie[t].bound)*S*S*S;
  if(used>all)return;
  if(i<0){
   if(!started)return;
   int q=0;for(int j=0;j<N;j++)if(z[j])for(int k=0;k<N;k++)q+=z[j]*G[j][k]*z[k];
   if(q>trie[t].bound)return;
   assert(q>=14&&trie[t].src>=0);
   for(int j=0;j<N;j++)wit[j]=z[j];
   seen_id=trie[t].src;hit=true;wq=q;return;
  }
  I lo=started?-wide[i]:0,hi=wide[i],cc=cen[i],a=clamp(floorq(cc+S/2,S),lo,hi),dir=cc>S*a?1:-1,m=max(a-lo,hi-a);
  for(I step=0;step<=2*m&&!hit&&!cap;step++){
   I u=a+(step&1?dir:-dir)*((step+1)/2);if(u<lo||u>hi)continue;
   I di=max<I>(0,abs(cc-S*u)-E[i]);J nw=used+J(d[i])*di*di;
   if(nw>all)break;
   int dst=trie[t].a[u&1];if(dst<0||nw>J(trie[dst].bound)*S*S*S)continue;
   z[i]=u;
   if(u)for(int j=0;j<i;j++)cen[j]-=u*mu[i][j];
   walk(i-1,nw,started||u,dst);
   if(u)for(int j=0;j<i;j++)cen[j]+=u*mu[i][j];
  }
 }
 void batch(W p,const vector<Record>& R,int q,U limit){
  trie.clear();trie.emplace_back();seen_id=-1;
  const U rr=label(p);assert(rr<=allParity);
  long long qr=0;U tags=0;
  for(int i=0;i<N;i++){
   long long z=0;for(int j=0;j<N;j++)if(rr>>j&1)z+=G[i][j];
   assert(z%2==0);
   if((z/2)&1)tags|=U(1)<<i;
   if(rr>>i&1)qr+=z;
  }
  assert(qr%4==0);
  // A shortcut cap reduction: store exact residue form per occurrence (norms are integer).
  for(int i=0;i<int(R.size());i++){
   const U a=label(p^R[i].p);
   if(a==0||a>allParity)continue;
   // No short candidate of parity 0 is needed: such a pair reduces to an already tested pure point.
   // The simpler common bound still suffices for verification; no exact congruence theory assumed.
   const int b=q-R[i].q;
   const long long r0=qr-R[i].q+4*__builtin_parityll((a^rr)&tags);
   const int bound=b-int((b-r0)%8+8)%8;
   if((i%997)==0)assert(bound==congruent_cap(a,b));
   trie_add(a,bound,i);
  }
  nodes=0;maxn=limit;hit=cap=false;
  if(trie[0].bound>0)walk(N-1,0,false,0);
 }
 void run(int q,U tar,U limit,bool cl=false){
  collect=cl;bound=q;nodes=0;maxn=limit;hit=cap=false;
  go(N-1,J(q)*S*S*S,false,tar);
 }
};

struct Codec{
 W q[M]{};U code[M]{};W basis[38]{};
 Codec(const Lat& D){
  int n=0;
  for(W z:D.basis){
   W x=z;U u=0;
   while(x){
    int k=last(x);
    if(q[k]){x^=q[k];u^=code[k];}
    else{
     assert(n<38);basis[n]=z;code[k]=u^(1ULL<<n);q[k]=x;n++;break;
    }
   }
  }
  assert(n==38);
 }
 U to(W v)const{
  U z=0;while(v){int k=last(v);assert(q[k]!=0);v^=q[k];z^=code[k];}
  return z;
 }
 W back(U a)const{
  W z=0;while(a){int k=__builtin_ctzll(a);assert(k<38);z^=basis[k];a&=a-1;}
  return z;
 }
};
struct Hash{size_t operator()(U u)const{u*=0x9e3779b9a1245791ULL;return u^(u>>32);}};
struct Automorphism{
 U table[5][256]{};
 Automorphism(const Codec&C,const char* file){
  ifstream f(file);W p[M]{};
  for(W &v:p){U u,w;f>>u>>w;v=(W(w)<<64)|u;}assert(f);
  U d[40]{};
  for(int i=0;i<38;i++){
   W b=C.basis[i],z=0;
   while(b){int k=last(b);b^=W(1)<<k;z^=p[k];}
   d[i]=C.to(z);
  }
  for(int j=0;j<5;j++)for(int u=0;u<256;u++)for(int i=0;i<8;i++)if(u>>i&1)table[j][u]^=d[8*j+i];
 }
 U operator()(U x)const{
  U v=0;for(int i=0;i<5;i++){v^=table[i][x&255];x>>=8;}
  return v;
 }
};
int swapK[40][40]{};
int main(int argc,char**argv){
 const U samplecap=argc>1?stoull(argv[1]):3500000000ULL;
 const U cp=argc>2?stoull(argv[2]):60000000ULL;
 const int focus=argc>3?atoi(argv[3]):-1;
 const int focus_count=argc>4?atoi(argv[4]):20;
 Lat L[2]={Lat(0),Lat(1)}; Lat K[2]={Lat(0,true),Lat(1,true)};
 Codec C(K[0]);
 Automorphism A(C,"aut0.txt"),J(C,"aut1.txt"),U0(C,"aut2.txt");
 ifstream sw("swap41.txt");for(auto&r:swapK)for(int&i:r)sw>>i;assert(sw);
 // Exhaust norm16 on K+, and generate K− under the exact U swap, verified in setup_swap.py.
 K[0].run(16,1ULL<<63,2000000000ULL,true);
 assert(!K[0].cap && K[0].pool.size()==77028);
 for(const Record &t:K[0].pool){
  array<short,40>x{};
  for(int j=0;j<40;j++){int v=0;for(int i=0;i<40;i++)v+=t.z[i]*swapK[i][j];assert(abs(v)<1000);x[j]=v;}
  W r=0;long long w=0;
  for(int i=0;i<40;i++){
   if(x[i]&1)r^=K[1].basis[i];
   for(int j=0;j<40;j++)w+=1LL*x[i]*K[1].G[i][j]*x[j];
  }
  assert(w==16&&r==C.back(U0(C.to(t.p))));K[1].pool.push_back({r,16,x});
 }
 cerr<<"norm16 complete both "<<K[0].nodes<<" "<<K[0].pool.size()<<endl;
 unordered_set<U,Hash> s16[2],s20[2];
 vector<U> r16[2];
 for(int i=0;i<2;i++)for(auto a:K[i].pool){U z=C.to(a.p);assert(s16[i].insert(z).second);r16[i].push_back(z);}
 if(samplecap){
  long long nraw=0;U qtest=C.to(K[0].pool[0].p);
  for(int i=0;i<38;i++){
   U c=1ULL<<i,t=c;
   for(int j=0;j<21;j++)t=A(t);assert(t==c);
   assert(J(J(J(c)))==c && A(J(c))==J(A(c)));
  }
  K[0].visit=[&](W p,int q,const array<short,40>&z){
   if(q!=20)return;
   nraw++;
   U c=C.to(p);
   if(s20[0].count(c))return;
   U v=c;
   for(int j=0;j<3;j++){
    U w=v;
    for(int i=0;i<21;i++){s20[0].insert(w);s20[1].insert(U0(w));w=A(w);}
    v=J(v);
   }
   if((nraw&((1<<18)-1))==0)cerr<<"sample "<<nraw<<" orbit20 "<<s20[0].size()<<" nodes "<<K[0].nodes<<endl;
  };
  K[0].run(20,1ULL<<63,samplecap,true);
  cerr<<"sample cap "<<K[0].cap<<" norms20 "<<nraw<<" orbit "<<s20[0].size()<<" nodes "<<K[0].nodes<<endl;
 }
 ifstream f("newlabels.txt");U lo,hi;int k,n=0,nskip=0,nchecked=0,nshortpure=0;
 while(f>>k>>lo>>hi){
  if(focus>=0&&k!=focus)continue;
  W p=(W(hi)<<64)|lo;bool hit=false,lost=false;U nodes=0;
  U r=C.to(p);
  if(samplecap){
   bool rejected=false;
   for(int h=0;h<2&&!rejected;h++)for(U a:r16[h]){
    const U b=r^a; if(s16[1-h].count(b)||s20[1-h].count(b)){rejected=true;break;}
   }
   if(rejected){nskip++;continue;}
  }
  auto out=[&](int h,const Record* x){
   cout<<"BLOCK "<<k<<" side "<<h<<" q "<<L[h].wq+(x?x->q:0)<<" nodes "<<nodes<<" ";
   for(short v:L[h].wit)cout<<v<<' ';
   for(int j=0;j<40;j++)cout<<(x?x->z[j]:0)<<' ';
   cout<<endl;
  };
  for(int h=0;h<2&&!hit;h++){
   U t=L[h].label(p);assert(t!=0&&t<=allParity);
   L[h].run(L[h].congruent_cap(t,36),t,cp);
   nodes+=L[h].nodes;lost|=L[h].cap;
   if((hit=L[h].hit)){out(h,nullptr);nshortpure++;}
  }
  if(hit)continue;
  if(nchecked++>=focus_count)break;
  cout<<"CHECK "<<k<<" "<<nskip<<" "<<nshortpure<<" pure nodes "<<nodes<<endl;
  for(int h=0;h<2&&!hit;h++){
   Lat&D=L[1-h];D.batch(p,K[h].pool,36,cp*40);
   nodes+=D.nodes;lost|=D.cap;
   if((hit=D.hit))out(1-h,&K[h].pool.at(D.seen_id));
   else cout<<"HALF "<<k<<" "<<h<<" nodes "<<D.nodes<<" cap "<<D.cap<<endl;
  }
  if(!hit){
   cout<<"RESULT "<<k<<" nodes "<<nodes<<" "<<(!lost?"EXACT_NO_VECTOR":"UNRESOLVED")<<endl;
   if(!lost)break;
  }
 }
 cerr<<"FIN screened away "<<nskip<<" pure "<<nshortpure<<" checked "<<nchecked<<endl;
}
