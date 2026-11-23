#include<bits/stdc++.h>
using namespace std;
int main(){
 int lg[12],g0[]={7,1,4,4,2,7,3,3,4,7,2,1};
 for(int i=0;i<12;i++){lg[i]=1<<23;for(int j=0;j<12;j++)lg[i]|=(g0[j]&1)<<(i+j);}
 vector<int> G{0};
 for(int a:lg){int n=G.size();for(int i=0;i<n;i++)G.push_back(G[i]^a);}
 assert(set<int>(G.begin(),G.end()).size()==4096);
 assert(count(G.begin(),G.end(),0xffffff)==1);
 for(int a:G){assert(__builtin_popcount(a)%4==0);assert(a==0||__builtin_popcount(a)>=8);}
 for(int a:lg)for(int b:lg)assert(__builtin_popcount(a&b)%2==0);
 int sy[24]={};for(int j=0;j<24;j++)for(int i=0;i<12;i++)if(lg[i]>>j&1)sy[j]|=1<<i;
 std::array<short,16384>a,b; std::array<long long,1024> hist={};

 int pr[2][2][4];
 for(int g=0;g<2;g++)for(int c=0;c<2;c++)for(int i=0;i<4;i++){int q=g+2*c+4*i;q=min(q,16-q);pr[g][c][i]=q*q;}
 int N0=0; vector<unsigned char> tab(1<<25,255);vector<int> idx(1<<24,-1); int ci=0;
 for(int g=0;g<2;g++)for(int c:G)if(c<(c^0xffffff)){
  idx[c]=ci++; vector<int> out(8192,999);
  for(int cc:{c,c^0xffffff}){
   a.fill(12345);a[0]=0;
   for(int j=0;j<24;j++){
    int s=sy[j],*m=pr[g][cc>>j&1];
    for(int p=0;p<4096;p++)for(int d=0;d<4;d++){
     int v0=a[4*p+d]+m[0];
     int v1=a[4*(p^s)+((d+3)&3)]+m[1];
     int v2=a[4*p+((d+2)&3)]+m[2];
     int v3=a[4*(p^s)+((d+1)&3)]+m[3];
     b[4*p+d]=min(min(v0,v1),min(v2,v3));
    }
    a.swap(b);
   }
   for(int p=0;p<4096;p++)for(int d=0;d<2;d++)
    out[2*p+d]=min(out[2*p+d],int(a[4*p+2*d]));
  }
  for(int p=0;p<4096;p++)for(int d=0;d<2;d++){
   int q=out[2*p+d];if(q>=999)continue;assert(q<1024);hist[q]++;N0++;tab[((ci-1)*4096+p)*2+d]=q;
   
  }
 }
 cerr<<"count "<<N0<<endl;
 assert(N0==(1<<24));
 for(int j=0;j<1024;j++)if(hist[j])cout<<j<<" "<<hist[j]<<endl;
 int basis[24][24];ifstream s0("odd_basis.txt");for(auto&x:basis)for(int&y:x)s0>>y;assert(s0);
 array<int,24>x={};vector<unsigned char> O(1<<24,255);
 for(unsigned i=0;i<(1u<<24);i++){
  int g=x[0]&1,c=0,hh=0,v=0;
  for(int j=0;j<24;j++){
   assert((x[j]&1)==g);
   int z=(x[j]-g)&15;c|=((z>>1)&1)<<j;
   int a=z>>2;v+=a;if(a&1)hh^=sy[j];
  }
  if(c>>23)c^=0xffffff;
  assert(v%2==0&&idx[c]>=0);
  int j=(i^(i>>1));int val=tab[((idx[c]-2048+g*2048)*4096+hh)*2+(v%4)/2];
  assert(val<90&&O[j]==255);O[j]=val;
  if(i+1<(1u<<24)){
    int b=__builtin_ctz(i+1),t=(j>>b&1)?-1:1;
    for(int j=0;j<24;j++)x[j]=(x[j]+t*basis[b][j])&15;
  }
 }
 array<long long,1024> vh={};for(int z:O)vh[z]++;assert(vh==hist);



 const int W[2][9]={{72,0,1,0,0,0,0,0,0},{13210188,205176,86860,16512,2448,7032,6732,0,0}};
 const int keys[9]={0,24,32,40,48,56,64,72,80};
 const int N=1<<24;
 vector<int64_t> P(N);
 for(auto &C:W){
  int w[90];fill(w,w+90,-1);for(int j=0;j<9;j++)w[keys[j]]=C[j];
  for(int j=0;j<N;j++){assert(w[O[j]]>=0);P[j]=w[O[j]];}
  for(int d=1;d<N;d*=2)for(int j=0;j<N;j+=2*d)for(int k=0;k<d;k++){
   int64_t a=P[j+k],b=P[j+d+k];P[j+k]=a+b;P[j+d+k]=a-b;
  }
  assert(*min_element(P.begin(),P.end())>=0);
  cout<<"POSITIVE DEFINITE "<<C[0]<<" "<<P[0]<<" "<<*min_element(P.begin(),P.end())<<endl;
 }
}
