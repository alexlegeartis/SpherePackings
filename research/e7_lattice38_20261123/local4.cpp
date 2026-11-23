#include<bits/stdc++.h>
using namespace std;int S,T,D[7],L[7][7],x[7],tab[16384],lim;
int ceildiv(int x,int y){return x>=0?(x+y-1)/y:x/y;}
void f(int j,int rest,int c){
 if(j<0){tab[c]=min(tab[c],(lim-rest)/(S*S*T));assert((lim-rest)%(S*S*T)==0);return;}
 int w=0;for(int i=j+1;i<7;i++)w+=L[i][j]*x[i];int r=sqrt(rest/D[j]);while((r+1)*(r+1)*D[j]<=rest)r++;while(r*r*D[j]>rest)r--;
 for(int z=ceildiv(-r-w,S);z<=-ceildiv(w-r,S);z++){x[j]=z;int y=S*z+w;f(j-1,rest-D[j]*y*y,c|((z&3)<<(2*j)));}
}
int main(){
 cin>>S>>T;for(auto&x:D)cin>>x;for(auto&r:L)for(auto&x:r)cin>>x;
 fill(tab,tab+16384,10000);lim=S*S*T*40;assert(lim<1000000000);f(6,lim,0);
 map<int,int>C;for(int x:tab)C[x]++;for(auto[k,v]:C)cerr<<k<<" "<<v<<"\n";assert(C.rbegin()->first<10000);
 for(int x:tab)cout<<x<<"\n";
}
