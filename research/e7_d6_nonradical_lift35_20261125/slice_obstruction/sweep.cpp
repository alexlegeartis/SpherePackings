#include <bits/stdc++.h>
using namespace std;
using U=unsigned char;
const int M[128]={0,1,7,6,8,9,15,14,27,26,28,29,19,18,20,21,38,39,33,32,46,47,41,40,61,60,58,59,53,52,50,51,44,45,43,42,36,37,35,34,55,54,48,49,63,62,56,57,10,11,13,12,2,3,5,4,17,16,22,23,25,24,30,31,82,83,85,84,90,91,93,92,73,72,78,79,65,64,70,71,116,117,115,114,124,125,123,122,111,110,104,105,103,102,96,97,126,127,121,120,118,119,113,112,101,100,98,99,109,108,106,107,88,89,95,94,80,81,87,86,67,66,68,69,75,74,76,77},B[6]={88,4613,5130,6396,7485,5014},A[4]={1,14,15,22};
const int F[7][7]={{2,-1,-1,0,-1,0,0},{0,3,-1,0,-1,0,0},{0,0,4,-3,-2,0,0},{0,0,0,5,-2,0,0},{0,0,0,0,4,-5,0},{0,0,0,0,0,3,-4},{0,0,0,0,0,0,1}},D[7]={30,10,5,3,3,5,40};
U r4[16384],costs[128][128];int mv[128],S[4][128],xv[7],tag[64];using E=array<U,8192>;
int wt(int a){return __builtin_popcount((unsigned)a);}
int gf(int x,int y){int r=0;while(y){if(y&1)r^=x;x*=2;if(x&64)x^=67;y>>=1;}return r;}
void en(int k,int b){if(k<0){int ix=0;for(int i=0;i<7;i++)ix|=(xv[i]&3)<<(2*i);r4[ix]=min<int>(r4[ix],(1440-b)/60);return;}int t=0;for(int j=k+1;j<7;j++)t+=F[k][j]*xv[j];int h=(int)sqrt(b/D[k]);while((h+1)*(h+1)*D[k]<=b)h++;while(h*h*D[k]>b)h--;for(int v=(int)ceil((-h-t)/(double)F[k][k]);v<=(int)floor((h-t)/(double)F[k][k]);v++){int q=D[k]*(t+v*F[k][k])*(t+v*F[k][k]);xv[k]=v;en(k-1,b-q);}}
void init(){memset(r4,99,sizeof r4);en(6,1440);assert(*max_element(r4,r4+16384)<=24);for(int i=0;i<128;i++){mv[M[i]]=i;for(int j=0;j<4;j++)S[j][i]=(i&63)|(gf(A[j],i&63)<<6)|((i>>6)<<12);}for(int b=0;b<128;b++)for(int t=0;t<128;t++){int ix=0;for(int i=0;i<7;i++)ix|=(((b>>i&1)+2*(M[t]>>i&1))&3)<<(2*i);costs[b][t]=r4[ix];}for(int t=0;t<64;t++)for(int j=0;j<6;j++)if(t>>j&1)tag[t]^=B[j];}
E add(const E&a,int j,int b){E q;q.fill(96);for(int t=0;t<128;t++){int k=S[j][t];U w=costs[b][t];for(int s=0;s<8192;s++)q[s]=min<U>(q[s],a[s^k]+w);}return q;}
E empty_(){E e;e.fill(90);e[0]=0;return e;}

constexpr int QD=128; using P=array<int,6>;extern int u0;
int ex[2][2][64],pflag[64],par[64];long long hist[32]={},un=0,good=0,exam=0,pairkill=0,scount=0;bool two=false;int T1=13*QD*QD/2,T2=2*QD*QD,Rlim=2;long long best=0; vector<array<int,7>> suspect; time_t deadline;
struct Site{int term;P r;};
int val(const P& x,int h,int w,int cutoff,bool single,std::vector<Site>* ss=nullptr){int best=9999999;
for(int m=0;m<2;m++){int cost[6][2],inc[6][2],rr[6][2];int hh=h*w;for(int i=0;i<6;i++)for(int l=0;l<2;l++){int q=w*x[i]-m*QD/2;int n=(int)lround((-q-QD*l)/(2.*QD));int r=q+QD*(l+2*n);int a=abs(r)+hh;rr[i][l]=r;cost[i][l]=a*a;inc[i][l]=(2*QD-abs(r)+hh)*(2*QD-abs(r)+hh)-a*a;}
for(int k=0;k<64;k++){int p=par[k],carry=pflag[k],v=ex[w-1][m][k],md=9999999,im=-1;P r;for(int i=0;i<6;i++){int bit=p>>i&1;r[i]=rr[i][bit];int n=(r[i]-(w*x[i]-m*QD/2)-QD*bit)/(2*QD);carry^=(n&1);v+=cost[i][bit];if(inc[i][bit]<md){md=inc[i][bit];im=i;}}
if(carry){v+=md;r[im]+=r[im]>=0?-2*QD:2*QD;}best=min(best,v);if(ss){int d=ex[w-1][m][k];for(int u:r)d+=u*u; if(d<=T1+12*QD*h&&d>T1-18*QD*h)ss->push_back({d,r});}if(single&&best<=cutoff)return best;
}}
return best;}
bool convexpair(const P&x,int h){vector<Site> a;val(x,h,1,T1,false,&a);for(int i=0;i<(int)a.size();i++)for(int j=0;j<i;j++){
 // Sum of two costs at their midpoint mixture, max over box; search rational lambda candidates at interval endpoints and each of six L1 kinks.
 for(int k=-1;k<6;k++){int n=1,d=2;if(k>=0){int u=a[i].r[k],v=a[j].r[k];if(u==v)continue;d=v-u;n=v;if(d<0)d=-d,n=-n;if(n<0||n>d)continue;}
 long long t=(long long)n*(a[i].term-T1)+(d-n)*(a[j].term-T1)+6LL*h*h*d;
 for(int l=0;l<6;l++)t+=2LL*h*abs(n*a[i].r[l]+(d-n)*a[j].r[l]);if(t<=0)return true;
 }}return false;}
void work(P &x,int h,int level){if(time(0)>=deadline){cerr<<"TIME reached"<<endl;exit(2);}++exam;if(val(x,h,1,T1,true)<=T1 ||(two&&val(x,h,2,T2,true)<=T2)){hist[level]++;return;}
if(h<=QD/4&&convexpair(x,h)){pairkill++;hist[level]++;return;}
int d=val(x,0,1,T1,false);if(d>best)best=d;if(d>T1){int e=val(x,0,2,T2,false);if(++good<=20){cout<<"FIRST y";for(int z:x)cout<<" "<<z;cout<<" d1 "<<d<<"/"<<QD*QD<<" d2 "<<e<<"/"<<QD*QD<<endl;}if(e>T2){cout<<"WIN u="<<u0<<endl;exit(0);}}
if(h<=Rlim){un++;array<int,7>null{};copy(x.begin(),x.end(),null.begin());null[6]=h;suspect.push_back(null);return;}
for(int k=0;k<64;k++){P y;for(int i=0;i<6;i++)y[i]=x[i]+((k>>i&1)?h/2:-h/2);work(y,h/2,level+1);}
}

E EX[2][2];int off[2][2],u0;
void gen(const array<int,4>&f0){for(int w=1;w<=2;w++)for(int m=0;m<2;m++){int of=m?2212:0;auto e=empty_();for(int j=0;j<4;j++){int f=0,c=0;for(int i=0;i<7;i++){int z=w*((f0[j]>>i)&1)+(m&&j?(82>>i)&1:0);f|=(z&1)<<i;c|=(z>>1&1)<<i;}of^=S[j][mv[c]];e=add(e,j,f);}EX[w-1][m]=e;off[w-1][m]=of;}}
int main(int argc,char**argv){int secs=argc>1?atoi(argv[1]):115;Rlim=argc>2?atoi(argv[2]):4;two=argc>3?atoi(argv[3]):false;deadline=time(0)+secs;int nprof=argc>4?atoi(argv[4]):1;int setting=argc>5?atoi(argv[5]):0;int start=argc>6?atoi(argv[6]):0;
init();for(int k=0;k<64;k++){int p=par[k]=(k&31)|(__builtin_parity(k&31)<<5);pflag[k]=(k>>5)^((__builtin_popcount(p)/2)&1);}
vector<int> rep;set<int>a;for(int i=0;i<8192;i++)if(!a.count(i)){rep.push_back(i);for(int k=0;k<64;k++)a.insert(i^tag[k]);}cerr<<"coset reps "<<rep.size()<<endl;
mt19937 rng(6060);for(int f=0;f<start+nprof;f++){array<int,4> b={21,60,41,40};if(f&&setting==0)b[(f-1)/127]^=(f-1)%127+1;if(setting==1)for(int j=0;j<4;j++)b[j]=rng()%128;if(f<start)continue;
gen(b);int ngoodprof=0;long long e0=exam,u1=un,p1=pairkill;
for(int u:rep){for(int w=0;w<2;w++)for(int m=0;m<2;m++)for(int k=0;k<64;k++)ex[w][m][k]=(QD*QD/4)*EX[w][m][off[w][m]^(w?0:u)^tag[k]];u0=u;for(int bank=0;bank<2;bank++){P x;x.fill(QD);x[0]+=2*QD*bank;work(x,QD,0);}if(un>u1&&f==0){cerr<<"unresolved at profile "<<f<<" u "<<u<<" count "<<un<<endl;u1=un;}}
cerr<<"PROFILE "<<f<<" "<<b[0]<<" "<<b[1]<<" "<<b[2]<<" "<<b[3]<<" exactExam "<<exam-e0<<" un "<<un<<" pairClosed "<<pairkill-p1<<" first "<<good<<" bestMid "<<best<<endl;}
cerr<<"FINAL exam "<<exam<<" unresolved "<<un<<" good "<<good<<" pairkill "<<pairkill<<endl;}
