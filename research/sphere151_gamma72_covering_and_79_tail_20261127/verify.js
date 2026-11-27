'use strict';
const assert=require('node:assert/strict');

const md=(a,m)=>((a%m)+m)%m,sc=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0),sq=a=>sc(a,a);
function zShort(limit){
 const out=[];
 for(let a=-5;a<=5;a++)for(let b=-5;b<=5;b++)for(let c=-5;c<=5;c++){
  const y=[a,b,c],norm=2*sq(y);if(norm>=limit||(a+b+c)%2)continue;
  const u=(a+b-c)/2,v=(a+c-b)/2,w=(b+c-a)/2;if(w%2)continue;
  const H=[2*u+v,v,v,v],h=[];
  function f(i,S,par){
   if(i===4){if(par%2||S===0)return;
    let x=h.map((e,j)=>(e-H[j])/2),p=(x[0]+x[1]+x[2]+x[3])/2,d=[x[0],x[0]+x[1],p-x[3],p,u,v,w/2];out.push([S,d]);return;
   }
   for(let e=-4+md(H[i],2);e<=4;e+=2){if(S+3*e*e>=limit)continue;h.push(e);f(i+1,S+3*e*e,par+(e-H[i])/2);h.pop();}
  }
  f(0,norm,0);
 }
 return out;
}
const short=zShort(64),two=short.filter(q=>q[0]<24),three=short.filter(q=>q[0]<32);

const mask=z=>z.reduce((s,x,i)=>s+((x&1)<<i),0),pc=x=>{let n=0;while(x){n++;x&=x-1;}return n&1;};
const read=(L,p)=>L.reduce((s,c,i)=>s+(pc(c&p)<<i),0);
const P2=new Set(two.map(([n,z])=>mask(z))),P3=new Set(three.map(([n,z])=>mask(z)));

function Kern(){
 const V=new Map(),K=[];
 function go(K,start){
   const s=K.slice().sort((a,b)=>a-b),k=s.join();if(V.has(k))return;V.set(k,s);
   for(let x=start;x<128;x++)if(!K.includes(x)){
    const a=K.map(y=>x^y);if(a.some(y=>P3.has(y)))continue;go([...K,...a],x+1);
   }
 }
 go([0],1);return [...V.values()];
}
function inspect(K){
 const f=p=>Math.min(...K.map(x=>x^p)),R=[...new Set([...P2].map(f))],S=new Set(R),A=R.map(a=>R.filter(b=>a!==b&&S.has(f(a^b))));
 for(let i=0;i<R.length;i++){
   const N=A[i];
   if(N.length>14)return ['degree',K,R[i],N.length];
   for(let a of N){
    const b=N.filter(c=>c!==a && S.has(f(a^c)));
    if(b.length>1)return ['link',K,R[i],a,b];
   }
   for(let j=0;j<i;j++)if(A[j].filter(a=>N.includes(a)).length>4)return ['codeg',K,R[i],R[j]];
 }
 return ['survive',K,R.length,A.map(k=>k.length)];
}
const KK=Kern(),I=KK.map(inspect);

assert.equal(short.length,2154);assert.equal(two.length,60);assert.equal(three.length,180);
assert.equal(Math.min(...short.map(x=>x[0])),16);
assert.equal(P2.size,30);assert.equal(P3.size,90);
const expected={link:37,codeg:1},counts={};
assert.equal(KK.length,38);
for(let r of I)counts[r[0]]=(counts[r[0]]||0)+1;
assert.deepEqual(counts,expected);
let dead=0;
for(let i=0;i<128;i++)if(!P3.has(i)&&i)for(let j=1;j<i;j++)if(!P3.has(j)&&!P3.has(i^j))dead++;
assert.equal(dead,0);
assert(!short.some(([q,z])=>mask(z)===0));
const N=2n**79n*41n**11n*234n*(936n**79n),D=4n*3235n**79n; // D_* = N/D
assert(81n*D<16n*N && 2n*N<11n*D);
assert(N*16n*1000n**2n>81n*D*1041n**2n && N*16n*1000n**2n<81n*D*1043n**2n);
const example=[
  [[0,0,0,0],[0,0,0]],
  [[0,0,0,0],[-2,-2,0]],
  [[0,0,0,0],[-2,0,-2]],
  [[-1,-1,-1,1],[-2,0,0]],
  [[-1,-1,1,-1],[-2,0,0]]
];
const std=[[0,0,0,0],[2,0,0,0],[1,1,1,1],[3,1,1,1]],U=[0,1,6,7];
for(let [h,[a,b,c]] of example){
 assert((a+b+c)%2===0);
 const k=U.indexOf(md((a+b-c)/2,2)+2*md((a+c-b)/2,2)+4*md((b+c-a)/2,2));
 assert(k>=0);let d=h.map((v,i)=>(v-std[k][i])/2);
 assert(d.every(Number.isInteger)&&d.reduce((s,x)=>s+x,0)%2===0);
}
for(let i=0;i<5;i++)for(let j=0;j<i;j++){
 const [h,y]=example[i],[g,x]=example[j],r=3*sq(h.map((v,i)=>v-g[i]))+2*sq(y.map((v,i)=>v-x[i]));
 assert(r>0 && (i===4&&j===3?r<32:r<24));
}
// Exact normalized center-density bounds for k=1,2,3.
const p=[[0n,1n],[1n,1n],[4n,3n],[2n,1n]],V=[0n,2n,3n,4n];
for(let k=1;k<=3;k++){
 const n=BigInt(k),[u,v]=p[k],lead=2n**n;
 function compare(a,b,strict=true){assert(strict?a*v*lead<b*u:a*v*lead<=b*u);}
 function A(x,y){return [u*y**(72n+n),v*x**72n*(8n*y-6n*x)**n];}
 function B(x,y,c,cmax){return [BigInt(c)**2n*4n**n*y**(72n+n),V[k]**2n*x**72n*(8n*y-BigInt(cmax)*x)**n];}
 compare(...A(1n,1n),false);compare(...A(5n,4n));
 compare(...B(5n,4n,4,5));compare(...B(3n,2n,4,5));
 compare(...B(3n,2n,4096,4));compare(...B(79n,40n,4096,4));
 compare(2n**72n*u*40n**72n,v*79n**72n*8n**n);
}
console.log('PASS: all scalar bounds 73–75; exact 79 gate and two scopes (all additive type 0; all labels on nominal type 1).');

