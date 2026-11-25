#include <bits/stdc++.h>
using namespace std; using U=uint64_t;
struct G{int D=64,add[64][64],sub[64][64],wt[64];array<int,3>P[64];U T[65][64]; int mod(int x,int y){return(x%y+y)%y;}int id(array<int,3>x){int p=mod(x[0],2),v=p*32,k=0;for(int j=0;j<3;j++){int d=j==2?8:4;assert(mod(x[j],2)==p);int s=mod(x[j],d);v|=(s/2)<<j;k+=(x[j]-s)/d;}return v+16*mod(k,2);}
G(){for(int i=0;i<64;i++){int p=i/32;P[i]={p+2*(i&1),p+2*((i>>1)&1),p+2*((i>>2)&3)};P[i][0]+=4*((i>>4)&1);assert(id(P[i])==i);}for(int i=0;i<64;i++)for(int j=0;j<64;j++){auto a=P[i],b=P[i];for(int k=0;k<3;k++){a[k]+=P[j][k];b[k]-=P[j][k];}add[i][j]=id(a);sub[i][j]=id(b);}vector<int>gol,oct;for(int i=0;i<4096;i++){int word=0,gen=(1<<11)|(1<<9)|(1<<7)|(1<<6)|(1<<5)|(1<<1)|1;for(int j=0;j<12;j++)if(i>>j&1)word^=gen<<j;word|=__builtin_parity(word)<<23;gol.push_back(word);if(__builtin_popcount(word)==8)oct.push_back(word);}assert(oct.size()==759);
vector<int>o,ext;for(int k=0;k<24;k++)((oct[0]>>k&1)?o:ext).push_back(k);int u=ext[0],v=ext[1],nsh=0;int old[64];fill(old,old+64,100000);for(int u=-10;u<=10;u++)for(int v=-10;v<=10;v++)for(int w=-16;w<=16;w++)if(mod(u-v,2)==0&&mod(v-w,2)==0){int i=id({u,v,w});old[i]=min(old[i],2*u*u+2*v*v+w*w);}int smallest[64];fill(smallest,smallest+64,9999);
auto visit=[&](array<int,24>&z){int so=0;for(int j=0;j<8;j++)so+=(j==7?-1:1)*z[o[j]];assert(so%2==0);int x=z[u],y=z[v],w=so/2;int c=id({x,y,w}),norm=64-(2*x*x+2*y*y+w*w);assert(norm>=0);if(norm>0)smallest[c]=min(smallest[c],norm);nsh++;};
array<int,24> z{};for(int i=0;i<24;i++)for(int j=0;j<i;j++)for(int a:{-4,4})for(int b:{-4,4}){z[i]=a;z[j]=b;visit(z);z[i]=z[j]=0;}
for(int w:oct){int ix[8],q=0;for(int j=0;j<24;j++)if(w>>j&1)ix[q++]=j;for(int m=0;m<256;m++)if(!__builtin_parity(m)){for(int j=0;j<8;j++)z[ix[j]]=(m>>j&1)?-2:2;visit(z);}for(int j:ix)z[j]=0;}
for(int c:gol){for(int j=0;j<24;j++)z[j]=(c>>j&1)?-1:1;for(int j=0;j<24;j++){z[j]*=-3;visit(z);z[j]/=-3;}}
assert(nsh==196560);for(int j=1;j<64;j++){assert(smallest[j]<9999);if(smallest[j]<9999){assert(smallest[j]==64-old[j]);wt[j]=smallest[j];}else wt[j]=min(64,96-old[j]);}wt[0]=0;
for(int s=0;s<=64;s++)for(int t=0;t<64;t++){U z=0;for(int r=0;r<64;r++)if(wt[add[t][r]]>=s)z|=1ULL<<r;T[s][t]=z;}map<int,int> cnt;for(int j:wt)cnt[j]++;cerr<<"weight profile ";for(auto[k,v]:cnt)cerr<<k<<":"<<v<<" ";cerr<<endl;}
int mul(int i,int n){n=mod(n,16);int r=0;while(n--){r=add[r][i];}return r;}}g;


int main(){int A=4,B=33,u[3]={0,1,2},c[3]={0,32,5};int best=999999,edges[3][3]{};set<array<int,4>> tight;for(int a=0;a<3;a++)for(int b=0;b<3;b++){int bs=999999;for(int i=-15;i<=15;i++)for(int j=-15;j<=15;j++){if(i==0&&j==0&&a==b)continue;int l=g.add[g.add[g.mul(A,i)][g.mul(B,j)]][g.sub[c[b]][c[a]]];int x=3*i+u[b]-u[a],y=3*j+u[b]-u[a];int n=4*x*x+4*x*y+37*y*y+9*g.wt[l];bs=min(bs,n);best=min(best,n);if(n==576){edges[a][b]++;tight.insert({u[b]-u[a],x,y,g.wt[l]});}}cout<<a<<","<<b<<"="<<bs<<"/144 "<<"at "<<edges[a][b]<<" pairs"<<endl;} cout<<"global "<<best<<"/144; colors: twice site32="<<g.add[32][32]<<" != 5; square plane det=9/16; center=1/2"<<endl;assert(best==576&&g.add[32][32]!=5);for(auto z:tight){for(auto i:z)cout<<i<<" ";cout<<endl;}}

