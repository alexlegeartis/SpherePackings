#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <vector>
using Z=int64_t; using W=__int128_t;
constexpr Z D=1000000000000LL;
Z fl(W a,W b){assert(b>0);return a>=0?a/b:-((-a+b-1)/b);}
Z ce(W a,W b){return -fl(-a,b);}
struct I {
    Z l,h;
    I(Z n=0):l(n*D),h(n*D){}
    I(Z l_,Z h_):l(l_),h(h_){assert(l<=h);}
    static I rat(Z a,Z b){return {fl(W(a)*D,b),ce(W(a)*D,b)};}
    I operator-()const{return {-h,-l};}
    friend I operator+(I a,I b){return {a.l+b.l,a.h+b.h};}
    friend I operator-(I a,I b){return a+-b;}
    friend I operator*(I a,I b){
        auto [u,v]=std::minmax({W(a.l)*b.l,W(a.l)*b.h,W(a.h)*b.l,W(a.h)*b.h});
        return {fl(u,D),ce(v,D)};
    }
    friend I operator/(I a,I b){assert(b.l>0);return a*I(fl(W(D)*D,b.h),ce(W(D)*D,b.l));}
    I sq()const {
        W a=W(l)*l,b=W(h)*h;
        return {l<=0&&h>=0?0:fl(std::min(a,b),D),ce(std::max(a,b),D)};
    }
};
struct Box{I p,u;int depth;};
void verify(int n,I k,std::array<int,6> cfs){
    const int den=(n+4)*(n+6);
    const I one(1),third=I::rat(1,3),half=I::rat(1,2),fourth=I::rat(1,4),fifth=I::rat(1,5);
    const I c3=I::rat(4,n+4),c4a=I::rat(24,den),c4b=I::rat(6*n,den),c5a=I::rat(192,den*(n+8)),c5b=I::rat(80*n,den*(n+8));
    const I m2coeff=I::rat(2,n+2),N=I::rat(n+2,2);
    const I c0=I::rat(cfs[0],10000),c1=I::rat(cfs[1],10000),c2=I::rat(cfs[2],10000),
      d0=I::rat(cfs[3],10000),d1=I::rat(cfs[4],10000),d2=I::rat(cfs[5],10000);
    std::vector<Box> todo{{I(0,D),I(0,D),0}};
    Z checks=0,v3=0,v5=0; int maxdepth=0;
    while(!todo.empty()){
        const auto [p,u,depth]=todo.back();todo.pop_back();++checks;
        const I b=one-p,w=u*p*b,L=k/(one+k*p);
        const I slow(p.l,p.l),shigh(p.h,p.h);
        const I A((k*slow/(one+k*slow)).l,(k*shigh/(one+k*shigh)).h);
        const I R=L*(half+third*A+fourth*A.sq()+fifth*A.sq()*A/(one-A));
        const I r=R*w,lo=-p+u*b,hi=b-u*p;
        if((R*(r+c3*hi)).h<m2coeff.l){++v3;continue;}
        const I c=c0+c1*p+c2*w,d=-w*(d0+d1*p+d2*w);
        const I Z4=c5a*b+(r-2*c)*c4a;
        const I Z3=c5b*w+(c.sq()-2*d-2*c*r)*c3;
        const I K=(r-2*c)*c4b*w+2*c*d+r*(c.sq()-2*d)+R*d.sq()*N;
        bool ok=true;
        if(Z4.h<0){
            const Z a=Z4.h,base=(K+Z4*w).h;
            for(Z t:{Z3.l,Z3.h}){
                W zden=-2*a,znum=W(t)*D;
                Z z;
                if(znum<lo.l*zden)z=lo.l;
                else if(znum>hi.h*zden)z=hi.h;
                else {ok &= W(-4*a)*base+W(t)*t<0;continue;}
                ok &= W(base)*D*D+W(a)*z*z+W(t)*z*D<0;
            }
        } else {
            Z z=std::max(p.h,b.h);
            W v=std::max({W(Z3.l)*lo.l,W(Z3.l)*hi.h,W(Z3.h)*lo.l,W(Z3.h)*hi.h});
            ok=W(K.h)*D*D+W(Z4.h)*z*z+v*D<0;
        }
        if(ok){++v5;continue;}
        if(depth>=40){
            printf("FAIL n=%d checks=%lld depth=%d p=%lld:%lld u=%lld:%lld\n",n,(long long)checks,depth,(long long)p.l,(long long)p.h,(long long)u.l,(long long)u.h);fflush(stdout);assert(0);
        }
        maxdepth=std::max(maxdepth,depth+1);
        if(depth%2==0){
            Z m=(p.l+p.h)/2;assert(p.l<m&&m<p.h);
            todo.push_back({I(p.l,m),u,depth+1});todo.push_back({I(m,p.h),u,depth+1});
        }else{
            Z m=(u.l+u.h)/2;assert(u.l<m&&m<u.h);
            todo.push_back({p,I(u.l,m),depth+1});todo.push_back({p,I(m,u.h),depth+1});
        }
    }
    assert(checks==2*(v3+v5)-1);
    printf("PASS dimension=%d checks=%lld degree3=%lld degree5=%lld depth=%d\n",n,(long long)checks,(long long)v3,(long long)v5,maxdepth);fflush(stdout);
}
int main(){
    verify(48,I(1),{1575,5517,19024,5375,10270,-28895});
    verify(72,I::rat(9,10),{1929,-1924,11708,2775,402,-10479});
}
