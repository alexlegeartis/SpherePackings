#!/usr/bin/env node
// Sphere188 94D tube cover. No dependencies. Larger 8D covering uses SV Theorem 1.2.
'use strict';
const assert=(x,m='assertion')=>{if(!x)throw Error(m);};
const gcd=(x,y)=>{while(y)[x,y]=[y,x%y];return x<0n?-x:x;};
const rat=(x,y)=>{let d=gcd(x,y);return `${x/d}/${y/d}`;};
const det=G=>G[0][0]*(G[1][1]*G[2][2]-G[1][2]**2n)-G[0][1]*(G[0][1]*G[2][2]-G[0][2]*G[1][2])+G[0][2]*(G[0][1]*G[1][2]-G[0][2]*G[1][1]);
const sections=[[[2,1029,864],[44,245,1920,2161]],[[864,1995,94],[2161,1880,277,3]],[[1363,1995,2160],[2161,336,245,3]],[[2160,1029,29],[44,277,329,2161]],[[29,1995,2],[2161,40,1781,3]],[[94,1029,1363],[44,1781,7,2161]],[[0,987],[1,0,0,1]],[[0,21],[46,0,0,46]]];
// F_47[omega], omega^2 + omega + 1 = 0; encode a+47*b.
const pack=(a,b)=>(a%47+47)%47+47*((b%47+47)%47);
const decode=x=>[x%47,Math.floor(x/47)];
function add(x,y){let [a,b]=decode(x),[c,d]=decode(y);return pack(a+c,b+d);}
function mul(x,y){let [a,b]=decode(x),[c,d]=decode(y);return pack(a*c-b*d,a*d+b*c-b*d);}
const neg=x=>mul(46,x);
function inv(x){assert(x!==0);let y=1;for(let n=2207;n;n>>=1,x=mul(x,x))if(n&1)y=mul(y,x);return y;}
function trim(x){while(x.length&&!x.at(-1))x.pop();return x;}
function padd(x,y){return trim(Array.from({length:Math.max(x.length,y.length)},(_,i)=>add(x[i]??0,y[i]??0)));} 
function pmul(x,y){let z=Array(x.length+y.length-1).fill(0);for(let i=0;i<x.length;i++)for(let j=0;j<y.length;j++)z[i+j]=add(z[i+j],mul(x[i],y[j]));return trim(z);}
function psub(x,y){return padd(x,y.map(neg));}
function gcddeg(x,y){x=x.slice();y=y.slice();while(y.length){const a=inv(y.at(-1));while(x.length>=y.length){let d=x.length-y.length,u=mul(x.at(-1),a);for(let j=0;j<y.length;j++)x[d+j]=add(x[d+j],neg(mul(u,y[j])));trim(x);}[x,y]=[y,x];}return x.length-1;}
function valid(Q){const [x,y]=Q;assert(x.length<=17&&y.length<=25);let a=Array(49).fill(0);a[0]=a[48]=1;assert(!psub(pmul(y,y),padd(pmul(pmul(x,x),x),a)).length);}
function dot(P,Q){const x=psub(P[0],Q[0]),y=psub(P[1],Q[1]);return !x.length&&!y.length?16:8-gcddeg(x,y)-Math.min(x.length?17-x.length:Infinity,y.length?25-y.length:Infinity);}
// Eight nonzero vertices of a regular simplex, written in u=t^8.
function expand(x){const z=Array(8*(x.length-1)+1).fill(0);x.forEach((c,i)=>z[8*i]=c);return trim(z);}
const pts=sections.map(a=>a.map(expand));pts.forEach(valid);
const H=pts.map(a=>pts.map(b=>dot(a,b)));
for(let i=0;i<8;i++)for(let j=0;j<8;j++)assert(H[i][j]===(i===j?16:8),'simplex Gram');

// A literal extra A8 glue root in the same F47[omega](t) field.
const extra=[[864,1029,2],[2161,1920,245,44]].map(expand);
valid(extra); assert(dot(extra,extra)===16);
const expected=[0,0,0,0,8,0,8,8];
assert(pts.every((v,i)=>dot(v,extra)===expected[i]),'A8 glue pairings');
// 3*extra = 3*(v4+v6+v7) - sum v_i. Equality follows by zero squared norm.
const c=expected.map(v=>v===8?2:-1);
const pairing=(a,b)=>a.reduce((s,x,i)=>s+x*b.reduce((t,y,j)=>t+y*H[i][j],0),0);
assert(pairing(c,c)===9*16 && c.reduce((z,v,i)=>z+v*expected[i],0)===3*16);
assert(9*dot(extra,extra)-6*c.reduce((z,v,i)=>z+v*expected[i],0)+pairing(c,c)===0);
assert(c.every(v=>((v%3)+3)%3===2)); // nontrivial index-three glue, covolume 8^4
const f=(n,d=1n)=>{n=BigInt(n);d=BigInt(d);assert(d>0n);let g=gcd(n,d);return [n/g,d/g];};
const A=(a,b)=>f(a[0]*b[1]+a[1]*b[0],a[1]*b[1]), M=(a,b)=>f(a[0]*b[0],a[1]*b[1]), S=(a,b)=>A(a,f(-b[0],b[1])),V=(a,b)=>f(a[0]*b[1],a[1]*b[0]),P=(a,j)=>f(a[0]**BigInt(j),a[1]**BigInt(j)), L=(a,b)=>a[0]*b[1]<a[1]*b[0], E=(a,b)=>a[0]*b[1]===a[1]*b[0];
const s0=f(1100098526n,100000000n),p=f(14,3),s=S(f(16),s0),D0=f(107840207,10000000),D1=f(1079109,100000),Db=f(1078516,100000);
const Q=d=>S(f(16),d);
const Df=d=>M(p,S(Q(d),V(p,f(4))));
const Dg=M(f(3,4),P(s,2)), old=Df(D0), newD=Df(D1), fallback=Df(Db);
assert(E(old,f(850355653n,45000000n)));
assert(E(Dg,f(187426112780529507n,10000000000000000n)));
assert(L(Dg,newD)&&L(newD,fallback)&&L(fallback,old),'strict improvements, fixed dominates');
for(let D of [D1,Db]){
 const q=Q(D);
 assert(L(f(4),S(q,V(p,f(4)))),'j>=2 safe');
 assert(L(f(0),S(D,M(f(2),p))) && L(f(16),A(q,M(f(6),p))),'only j=±1 phase types');
}
assert(L(f(16),M(f(4),s))&&L(f(0),S(f(16),M(f(3),s))),'generic stability');
assert(E(S(s0,V(S(f(16),M(f(3),s)),f(3))),f(32,3)));
// Explicit rational backup: A8* at scale 64/51 relative to the sqrt8E8 coordinates.
const scale=f(64,51),radiusA8=M(M(f(8),P(scale,2)),f(20,27)),rateA8=M(f(3),P(f(51,64),8));
assert(L(radiusA8,M(f(2),p))&&L(rateA8,f(1,2)));
assert(L(M(A(f(2),rateA8),P(Db,46)),P(s0,46)),'backup coverage');
// Schuermann--Vallentin lattice plus rational approximation. Prove strict rate <.42575
// using pi>3.14159 (itself checked by the exact alternating Machin formula).
function atan(q,h){let a=f(0);for(let k=0n;k<h;k++)a=A(a,f((-1n)**k,(2n*k+1n)*q**(2n*k+1n)));return a;}
const piLower=S(M(f(16),atan(5n,40n)),M(f(4),atan(239n,41n)));
assert(L(f(314159,100000),piLower),'pi bound');
const V8=V(P(f(314159,100000),4),f(24));
const rateSV=f(1703,4000);
assert(L(M(f(32013,10000),P(f(8),4)),M(M(rateSV,V8),P(M(f(2),p),4))),'strict density margin for commensurable cover');
assert(L(M(A(f(2),rateSV),P(D1,46)),P(s0,46)),'main covering contradiction');
assert(L(M(P(f(500437,500000),2),newD),old),'gain >.0874%');
const oldFloor=f(8078244769088839n,10n**36n),newFloor=f(808530552528770n,10n**35n);
assert(L(M(P(newFloor,2),newD),M(P(oldFloor,2),old)),'rounded-down occupied floor');
console.log(JSON.stringify({genericDet:rat(...Dg),fixedDet:rat(...newD),fallbackDet:rat(...fallback),priorDet:rat(...old),occupiedDensityFloor:'8.08530552528770e-21',gainPercent:'>0.0874'},null,2));
