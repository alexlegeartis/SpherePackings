#include<bits/stdc++.h>
using namespace std;using I=int64_t;using J=__int128_t;
constexpr int U=84;constexpr I S=1LL<<20;
I G[U][U],D[U],M[U][U],W[U],E[U],cen[U];
int z[U],n,goal,printed=0;uint64_t nodes=0,limit, cnt[128];
bool cut=false;
I fd(I a,I b){return a>=0?a/b:-((-a+b-1)/b);}
void visit(int i,J bud,bool live){
 if(cut)return;if(++nodes>limit){cut=true;return;}
 if(i<0){
   if(!live)return;
   I q=0;for(int j=0;j<n;j++)if(z[j])for(int k=j;k<n;k++)if(z[k])q+=(j==k?1:2)*z[j]*G[j][k]*z[k];
   if(q>goal)return;
   assert(q>0 && (q%2==0));cnt[q]++;
   if(q<=10 && !printed){cout<<"witness "<<q;for(int i=0;i<n;i++)cout<<" "<<z[i];cout<<endl;printed=1;}
   return;
 }
 const I lo=live?-W[i]:0,hi=W[i], a=clamp(fd(cen[i]+S/2,S),lo,hi),p=cen[i]>=S*a?1:-1;
 for(I j=0;j<=2*max(a-lo,hi-a);j++){
   I t=a+(j&1?p:-p)*((j+1)/2);if(t<lo||t>hi)continue;
   I d=max<I>(0,abs(cen[i]-S*t)-E[i]);J q=J(D[i])*d*d;if(q>bud)break;
   z[i]=t;if(t)for(int c=0;c<i;c++)cen[c]-=t*M[i][c];
   visit(i-1,bud-q,live||t);
   if(t)for(int c=0;c<i;c++)cen[c]+=t*M[i][c];
   if(cut)return;
 }
}
int main(int argc,char**argv){
 cin>>n>>goal;assert(n<=84);
 for(int i=0;i<n;i++)for(int j=0;j<n;j++)cin>>G[i][j];
 for(int i=0;i<n;i++)cin>>D[i];
 for(int i=0;i<n;i++)for(int j=0;j<n;j++)cin>>M[i][j];
 for(int i=0;i<n;i++)cin>>W[i];assert(cin);
 for(int i=n-2;i>=0;i--)E[i]=E[i+1]+W[i+1];
 limit=argc>1?stoull(argv[1]):50000000ULL;
 visit(n-1,J(goal)*S*S*S,false);
 cout<<"cap "<<cut<<" nodes "<<nodes;
 for(int i=2;i<=goal;i+=2)cout<<" "<<2*cnt[i];
 cout<<endl;
}
