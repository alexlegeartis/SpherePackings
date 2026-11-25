#include <bits/stdc++.h>
using namespace std; using I=uint32_t;
int wt(I x){return __builtin_popcount(x);} int pa(I x){return __builtin_parity(x);} const I ALL=(1u<<24)-1;
struct Basis{I r[24]{};bool add(I x){for(int i=23;i>=0;i--)if(x>>i&1){if(r[i])x^=r[i];else {r[i]=x;return 1;}}return 0;} I red(I x){for(int i=23;i>=0;i--)if((x>>i&1)&&r[i])x^=r[i];return x;}}; Basis br;
int pack(I x){x=br.red(x);int k=0,s=0;for(int j=0;j<24;j++)if(!br.r[j])s|=int(x>>j&1)<<k++;return s;}
vector<I> R,P,B; vector<int> col; vector<I> rb(1024),bcarry(1024),bval(4096);vector<int>ecarry(1024),bsyn(4096);vector<unsigned char> F(1<<24),BF(1<<22);
int S(I x){int s=0;while(x){int j=__builtin_ctz(x);x&=x-1;s^=col[j];}return s;}
void dist(I kept,vector<unsigned char>& d){fill(d.begin(),d.end(),99);d[0]=0;for(int j=0;j<24;j++){int k=col[j],w=(kept>>j)&1,bit=k&-k;for(int b=0;b<8192;b+=2*bit)for(int i=b;i<b+bit;i++){int h=i^k,u=d[i],v=d[h];d[i]=min(u,v+w);d[h]=min(v,u+w);}}}
struct Point{I p,b;int e,half;};
Point getp(int addr){int h=addr>>23,k=(addr>>13)&1023,n=(addr>>1)&4095;I p=rb[k],d=bval[n],b=bcarry[k]^d;int e=(addr&1)^ecarry[k]^((wt(d)/2)&1)^pa(d&bcarry[k]);if(h){e^=b&1;b^=1;}return {p,b,e,h};}
int direct(int a){auto [p,b,e,h]=getp(a);int res=500;for(I r:R)for(int odd=0;odd<2;odd++){int sum=0,alt=1000,z=e^odd;for(int i=0;i<24;i++){int v=2*int((p>>i)&1)+4*int((b>>i)&1)+h-4*int((r>>i)&1)-2*odd;int u=((v+4)&7)-4;z^=((v-u)/8)&1;sum+=u*u;alt=min(alt,64-16*abs(u));}res=min(res,(sum+(z?alt:0))/4);}return res;}
int main(int ac,char**av){int mark=ac>1?atoi(av[1]):3; auto start=chrono::steady_clock::now();vector<I> G;for(I k=0;k<4096;k++){I w=0;for(int j=0;j<12;j++)if((k>>j)&1)w^=2787u<<j;w|=pa(w)<<23;G.push_back(w);if(!w)continue;assert(wt(w)>=8&&wt(w)%4==0);if(!(w&1)&&!pa(w&mark)&&br.add(w))P.push_back(w);}assert(P.size()==10&&br.add(ALL));for(I w:G)if(!pa(w&mark))R.push_back(w);assert(R.size()==2048);Basis whole=br;for(int j=1;j<24;j++){I w=1u^(1u<<j);if(whole.add(w))B.push_back(w);}assert(B.size()==12);for(int j=0;j<24;j++)col.push_back(pack(1u<<j));for(int i=0;i<4096;i++){for(int j=0;j<12;j++)if(i>>j&1)bval[i]^=B[j];bsyn[i]=S(bval[i]);} for(int k=0;k<1024;k++){for(int i=0;i<24;i++){int v=0;for(int j=0;j<10;j++)if(k>>j&1)v+=int((P[j]>>i)&1);rb[k]|=(v&1)<<i;bcarry[k]|=((v>>1)&1)<<i;ecarry[k]^=(v>>2)&1;}assert(!pa(bcarry[k]));}
vector<unsigned char> fd(8192),da(8192),db(8192);dist(ALL,fd);vector<vector<I>> emin(8192);vector<int>count(25);for(auto v:fd)count[v]++;for(int i=0;i<25;i++)if(count[i])cerr<<" d"<<i<<"="<<count[i];cerr<<endl;
for(int s=0;s<8192;s++)if(fd[s]>=7&&(fd[s]&1)){I v=0;int m=0;for(int j=0;j<24;j++)if(!br.r[j]){if(s>>m&1)v|=1u<<j;m++;}for(I r:R)if(wt(v^r)==fd[s])emin[s].push_back(v^r);}
int ni=0,nh=0,nro=0;for(int k=0;k<1024;k++){I p=rb[k];int w=wt(p),sc=S(bcarry[k]);dist(ALL^p,da);dist(p,db);for(int i=0;i<4096;i++){
int base=(k<<13)|(i<<1),s=sc^bsyn[i];auto [q,b,e,h]=getp(base);int D1=da[s],D2=db[s],lo=min(w+4*D1,24-w+4*D2);unsigned mask=0; if(lo>=24)mask=3;else if((w+4*D1>=24||(w==16&&D1==0))&&(24-w+4*D2>=24||(w==8&&D2==0))){mask=3;for(I r:R){if(w==16&&D1==0&&!((r^b)&(ALL^p))){int c=e^pa(p&b&~r);mask&=~(1u<<c);}
if(w==8&&D2==0&&!((r^b)&p)){// x-z = complement(p)+2(b xor complement(p))+4*(e+1+popcount(complement(p)&~b))
I t=ALL^p,bb=b^t;int ee=e^1^pa(t&~b);int c=ee^pa(t&bb&~r);mask&=~(1u<<c);}
if(!mask)break;}}
if(mask&1)F[base]=1,ni++;if(mask&2)F[base+1]=1,ni++;if(mask==3)BF[base>>1]=1,nro++;
base|=(1<<23);Point u=getp(base);I d=u.b;int ss=s^col[0],dd=fd[ss];if(6+2*dd>=24)mask=3;else if(dd<7)mask=0;else{mask=3;for(I err:emin[ss]){I r1=p^d^err, r2=d^err;int c1=u.e^pa(d&~r1);int c2=u.e^1^pa((ALL^p)&(~d)&r2)^pa(p&d&~r2);mask&=~((1u<<c1)|(1u<<c2));if(!mask)break;}}
if(mask&1)F[base]=1,nh++;if(mask&2)F[base+1]=1,nh++;}}
cerr<<"mark "<<mark<<" far integral "<<ni<<" far half "<<nh<<" robust "<<nro<<" sec "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<endl;
mt19937 rng(123);for(int t=0;t<1000;t++){int a=rng()&0xffffff;bool pred=F[a];int v=direct(a);if((v>=24)!=pred){cerr<<"FAIL random "<<a<<" "<<v<<" "<<pred<<endl;return 2;}}
assert((mark==3 && ni==50648 && nh==47104)||(mark==15 && ni==28096 && nh==24576));
ofstream file(string("far")+to_string(mark)+".bin",ios::binary);for(int i=0;i<(1<<24);i++)if(F[i]) file.write((char*)&i,4);return !file.good();}
