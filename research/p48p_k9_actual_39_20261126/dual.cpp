#include <cstdio>
#include <vector>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <algorithm>
#include <unordered_map>
#include <string>
using namespace std;const int n=39,S=16384;int M[n][n],V[n][n],L[n][n],D[n],H[n],E[n],v[n],t[n],cut;long long nodes=0,leaves=0,hist[200],hyp=0;FILE *fo;unordered_map<string,int> seen;
void leaf(){int q=0;for(int i=0;i<n;i++)if(v[i]){int s=0;for(int j=0;j<n;j++)s+=v[j]*M[i][j];q+=v[i]*s;}if(q<=0||q>cut)return;leaves++;hist[q]+=2;string k(n,0);for(int i=0;i<n;i++){int s=0;for(int z=0;z<n;z++)s+=v[z]*V[z][i];k[i]=(s%30+30)%30;}auto it=seen.find(k);if(it==seen.end()||q<it->second){seen[k]=q;fwrite(k.data(),1,39,fo);fwrite(&q,4,1,fo);for(char &z:k)z=(30-int(z))%30;seen[k]=q;}}
void rec(int i,long long av,bool zero){nodes++;if(i<0){leaf();return;}int y=t[i],m=y>=0?-(y+S/2)/S:(-y+S/2)/S,dir=(-y>=m*S?1:-1);
for(int z=0;;z++){int a=m+(z+1)/2*(z%2?dir:-dir),u=max(0,abs(a*S+y)-E[i]); long long rem=av-1LL*u*u*D[i];if(rem<0)break;if((zero&&a<0)||abs(a)>H[i])continue;v[i]=a;int* row=L[i];if(a)for(int j=0;j<i;j++)t[j]+=row[j]*a;rec(i-1,rem,zero&&a==0);if(a)for(int j=0;j<i;j++)t[j]-=row[j]*a;}}
int main(int argc,char**argv){cut=argc>1?atoi(argv[1]):100;freopen("exactdual.dat","r",stdin);for(int i=0;i<n;i++)for(int j=0;j<n;j++)scanf("%d",&M[i][j]);for(int i=0;i<n;i++)for(int j=0;j<n;j++)scanf("%d",&V[i][j]);for(int i=0;i<n;i++)for(int j=0;j<n;j++)scanf("%d",&L[i][j]);for(int i=0;i<n;i++)scanf("%d",&D[i]);for(int i=0;i<n;i++)scanf("%d",&H[i]);for(int i=0;i<n;i++)scanf("%d",&E[i]);fo=fopen("low.dat","wb");rec(n-1,1LL*cut*S*S*S,1);fclose(fo);printf("cut %d nodes %lld nonzero vectors %lld classes %d\n",cut,nodes,2*leaves,int(seen.size())); for(int i=1;i<=cut;i++)if(hist[i])printf("norm30 %d count %lld\n",i,hist[i]);}
