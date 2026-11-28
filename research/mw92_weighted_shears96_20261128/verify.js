// Exact replay for the 96D coupled Eisenstein tails. Node.js, no dependencies.
// Parent MW92 and the equal-ball threshold A are linked in README.md.
'use strict';
const assert=(x,m='assertion')=>{if(!x)throw Error(m);};
function gcd(a,b){while(b)[a,b]=[b,a%b];return a<0n?-a:a;}
function F(n,d=1n){n=BigInt(n);d=BigInt(d);assert(d!==0n);const g=gcd(n,d);return [d<0n?-n/g:n/g,d<0n?-d/g:d/g];}
const add=([p,q],[r,s])=>F(p*s+q*r,q*s), sub=([p,q],[r,s])=>F(p*s-q*r,q*s);
const mul=([p,q],[r,s])=>F(p*r,q*s), div=([p,q],[r,s])=>F(p*s,q*r);
const pow=([p,q],k)=>[p**k,q**k], lt=([p,q],[r,s])=>p*s<r*q;
const eq=([p,q],[r,s])=>p===r&&q===s;
const pstr=([p,q])=>p.toString()+'/'+q.toString();
const A=F(1100098526n,100000000n),s=sub(F(16),A),t=sub(F(16),mul(F(3),s));
assert(lt(F(0),t)&&lt(F(16),mul(F(4),s)));
assert(eq(sub(A,div(t,F(3))),F(32,3)));
const norm=(i,j)=>i*i-i*j+j*j;
const B=F(537,50); // 10.74
assert(lt(mul(F(3),pow(B,46n)),pow(A,46n)));
const B2=F(1074136,100000),dd=F(3583,1000000),p=F(61433,100000),p2=sub(F(1),p),Sf=F(16,3),ef=F(3869092,1000000);
assert(lt(mul(F(3),pow(B2,46n)),pow(A,46n)));
const VV=add(add(mul(F(34,9),pow(p,2n)),mul(F(34,12),pow(p2,2n))),mul(B2,mul(p,p2)));
const Wf=sub(sub(F(16),mul(F(3),ef)),mul(F(48),pow(dd,2n)));
assert(lt(Wf,VV));assert(lt(F(0),Wf));
const ni=(x,y)=>sub(add(pow(x,2n),pow(y,2n)),mul(x,y));const rho=[add(F(2,3),mul(F(2),dd)),add(F(1,3),dd)];
function em([a,b],[c,d]){return [sub(mul(a,c),mul(b,d)),sub(add(mul(a,d),mul(b,c)),mul(b,d))];}
const fixed=new Set(), fixedBounds={};
let Df=F(0);
for(const [i,j,n] of [[1,0,1],[1,-1,3],[2,0,4]]){
 const [x,y]=em(rho,[F(i),F(j)]);
 assert(lt(F(-2),x)&&lt(x,F(2))&&lt(F(-2),y)&&lt(y,F(2)));
 assert(lt(F(16),mul(F(9,2),Sf)));
 let count=0;
 for(let a=-4;a<=4;a++)for(let b=-4;b<=4;b++){
 const Q=add(mul(F(n),ef),mul(Sf,ni(add(F(a),x),add(F(b),y))));
 if(!lt(Q,F(16)))continue;
 count++;
 if(n===1){fixed.add(((a+b)%3+3)%3);let d=sub(F(16),Q);if(lt(Df,d))Df=d;}
 else if(n===3){assert(a===-1&&b===0);assert(eq(sub(F(16),Q),Wf));}
 else assert(false);
 }
 fixedBounds[n]=count;
}
assert(fixed.size===3&&fixedBounds[3]===1&&fixedBounds[4]===0);
assert(lt(F(16),mul(F(5),ef)));
assert(lt(mul(F(15),pow(Df,46n)),pow(A,46n)));
const df=pow(mul(F(4),ef),2n);
const dg=mul(F(9,16),mul(pow(s,2n),pow(F(514,125),2n)));
const old=F(149769,625);
const periodic=F(300153807922100058481n,1440000000000000000n);
assert(lt(periodic,dg)&&lt(dg,df)&&lt(df,old));
function atan(z,k){let v=F(0);for(let i=0n;i<k;i++)v=add(v,F((-1n)**i,(2n*i+1n)*z**(2n*i+1n)));return v;}
const pi=sub(mul(F(16),atan(5n,40n)),mul(F(4),atan(239n,41n)));
const fact48=Array.from({length:48},(_,i)=>BigInt(i+1)).reduce((a,b)=>a*b,1n);
const center2=div(F(2n**192n),mul(F(47n**14n),df));
const occ2=mul(center2,pow(div(pow(pi,48n),F(fact48)),2n));
function isqrt(n){let a=0n,b=1n<<BigInt(Math.ceil(n.toString(2).length/2)+1);while(b-a>1n){let c=(a+b)/2n;if(c*c<=n)a=c;else b=c;}return a;}
const floor=(v,k)=>isqrt(v[0]*10n**(2n*k)/v[1]);
const frac=(z,d)=>[z/10n**d,(z%10n**d).toString().padStart(Number(d),'0')].join('.');
const occ=floor(occ2,37n), gain=floor(div(old,df),8n)-10n**8n;
const bcsRatio2=mul(center2,F(2n**48n,3n**96n));
const bcsGain=floor(bcsRatio2,8n)-10n**8n;
assert(occ>5938937064720705n);
assert(lt(F((10n**8n+gain)**2n,10n**16n),div(old,df)));
console.log(JSON.stringify({Df:pstr(Df),Wf:pstr(Wf),VV:pstr(VV),branchCounts:fixedBounds,old:pstr(old),fixedDet:pstr(df),genericDet:pstr(dg),occupied:frac(occ,15n)+'e-22',gainPercent:frac(gain,6n)},null,2));
