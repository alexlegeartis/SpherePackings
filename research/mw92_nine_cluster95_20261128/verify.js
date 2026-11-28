// Sphere168 95D exact arithmetic, ECMAScript/Node.js, no dependencies.
// Parent geometric and CFR theorems are stated in README.md.
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
function squareNorm(a){return a.reduce((s,x,i)=>s+x*a.reduce((t,y,j)=>t+y*H[i][j],0),0);}
// No nonzero F2 relation on the first three outer vertices. For a
// relation choose signs [+], [+,-], [+,+,-] according to support size.
for(let mask=1;mask<8;mask++){
 const a=Array.from({length:8},(_,j)=>j<3&&((mask>>j)&1)?1:0),k=a.reduce((s,x)=>s+x,0);
 if(k>=2)a[a.lastIndexOf(1)]=-1;
 assert(squareNorm(a)===(k===3?32:16)&&squareNorm(a)<4*16);
}
// The F3 Gram kernel consists of the constants. The chosen images
// are all eight nonzero points of F3^2, whose sum is 0. Therefore
// they kill every possible relation in L/3L, even if not all vanish.
const img=Array.from({length:8},(_,j)=>[(j+1)%3,Math.floor((j+1)/3)]);
for(let c=0;c<3**8;c++){
 const v=Array.from({length:8},(_,i)=>Math.floor(c/3**i)%3);
 const kernel=H.every(r=>r.reduce((s,x,j)=>s+x*v[j],0)%3===0);
 assert(kernel===v.every(x=>x===v[0]),'mod3 kernel');
 if(kernel)for(let t=0;t<2;t++)assert(img.reduce((s,r,j)=>s+r[t]*v[j],0)%3===0);
}
for(let k of [3,8]){
 const sum=H.slice(0,k).reduce((s,row)=>s+row.slice(0,k).reduce((s,x)=>s+x,0),0);
 assert(sum===8*k*(k+1));
 for(let i=0;i<k;i++)assert(16*(k+1)**2+sum-2*(k+1)*H[i].slice(0,k).reduce((s,x)=>s+x,0)===sum,'equal centroid radii');
}
// Tail matrices, scaled by 600000000. All inequalities are exact.
const M=600000000n, A=6600591156n, D2=6470412420n;
const s=16n*M-A,p=14n*M/3n,q=16n*M-D2,t=16n*M-3n*s;
const Dg=6317700000n,Df=6328284000n;
const Gg=[[s,s/2n,-s/2n],[s/2n,s,-s/2n],[-s/2n,-s/2n,16n*M-Dg]];
const Gf=[[p,p/2n,p/2n],[p/2n,q,q/2n],[p/2n,q/2n,16n*M-Df]];
assert(s%2n===0n&&p%2n===0n&&q%2n===0n);
const w=D2-2n*p;
const leRootSum=(a,b,d)=>{
 const z=d-a-b;
 assert(a>0n&&b>0n&&z>0n&&z*z>4n*a*b,'strict surd gate');
};
function demands(G,fixed){
 const delta=det(G);
 assert(G[0][0]>0n&&G[0][0]*G[1][1]-G[0][1]**2n>0n&&delta>0n,'positive G');
 for(let i=0;i<3;i++){const j=(i+1)%3,k=(i+2)%3;assert(16n*M*(G[j][j]*G[k][k]-G[j][k]**2n)<64n*delta,'inverse bound');}
 let out=new Map();
 for(let i=-7;i<=7;i++)for(let j=-7;j<=7;j++)for(let k=0;k<=7;k++){
  const v=[i,j,k].map(BigInt);
  if(!k&&j<0||!k&&!j&&i<=0)continue;
  let Q=0n;for(let u=0;u<3;u++)for(let h=0;h<3;h++)Q+=v[u]*v[h]*G[u][h];
  if(Q>=16n*M)continue;
  if(!k){
   if(fixed){if(!j)assert(Math.abs(i)===1&&Q+34n*M/3n>=16n*M);else assert(j===1&&Q+((i%3+3)%3===1?w:D2)>=16n*M);}
   else {let n=i*i+i*j+j*j;assert((n===1||n===3)&&Q+(n===1?A:t)>=16n*M);}
   continue;
  }
  assert(k===1,'only k=1');
  let ph=[fixed?(i%3+3)%3:i,j].join(','),f=16n*M-Q;
  if(!out.has(ph)||out.get(ph)<f)out.set(ph,f);
 }
 return out;
}
const keys=(...xs)=>new Set(xs.map(x=>x.join(',')));
const near=keys([0,0],[1,0],[0,1],[-1,1],[1,-1],[1,1]);
const far=keys([0,-1],[-1,0],[-1,2],[0,2],[2,-1],[2,0]);
const g=demands(Gg,false),f=demands(Gf,true);
assert(g.size===near.size+far.size);
for(let ph of near)assert(g.get(ph)>0n&&g.get(ph)<=Dg,'generic large demand');
for(let ph of far)assert(g.get(ph)===Dg-2n*s,'generic weak demand');
leRootSum(6n*M,Dg-2n*s,Dg);
assert(15n*Dg**46n<2n*A**46n);
const tiny=keys([0,1],[1,1],[0,-2],[1,-2]),mid='1,0';
assert(f.size===11&&f.get(mid)===Df-2n*p);
for(let ph of tiny)assert(f.get(ph)===Df-2n*q,'fixed small demand');
for(let [ph,b] of f){assert(b>0n&&b<=Df);}
leRootSum(4n*M,Df-2n*p,Df);
leRootSum(64n*M,9n*(Df-2n*q),9n*Df);
assert(125n*Df**46n<18n*A**46n);
const dg=det(Gg),df=det(Gf),oldn=3578813470795833160453541n,oldd=50000000000000000000000n;
assert(dg<df&&df*oldd<oldn*M**3n);
// Machin with exact rational arithmetic, for an outward occupied-density bound.
const zfrac=(n,d)=>{let g=gcd(n,d);return [n/g,d/g];};
const plus=([n,d],[p,q])=>zfrac(n*q+p*d,d*q);
function atan(q,h){let a=[0n,1n];for(let k=0n;k<h;k++)a=plus(a,[(-1n)**k,(2n*k+1n)*q**(2n*k+1n)]);return a;}
let [x,y]=atan(5n,40n),[u,v]=atan(239n,41n);
const pil=plus([16n*x,y],[-4n*u,v]),[pi,piD]=pil;
const dfact=Array.from({length:48},(_,i)=>BigInt(2*i+1)).reduce((a,b)=>a*b,1n);
const densN=(2n**143n*pi**47n)**2n*M**3n;
const densD=(dfact*47n**7n*piD**47n)**2n*df;
function isqrt(n){let a=0n,b=1n<<BigInt(Math.ceil(n.toString(2).length/2)+1);while(b-a>1n){let c=(a+b)/2n;if(c*c<=n)a=c;else b=c;}return a;}
const floor=isqrt(densN*10n**72n/densD);
const gf=isqrt(oldn*M**3n*10n**16n/(oldd*df))-10n**8n;
const result={previous_det:rat(oldn,oldd),generic_det:rat(dg,M**3n),fixed_det:rat(df,M**3n),worst_det:rat(df,M**3n),
 generic_coverage_centers_per_L_cell:'15/2',fixed_coverage_centers_per_L_cell:'125/18',
 occupied_density_floor:`${floor/10n**15n}.${(floor%10n**15n).toString().padStart(15,'0')}e-21`,
 gain_percent_floor:`${gf/10n**6n}.${(gf%10n**6n).toString().padStart(6,'0')}`};
assert(floor>2129324206844210n);
assert(oldn*M**3n*10n**16n>(10n**8n+gf)**2n*oldd*df);
console.log(JSON.stringify(result,null,2));
