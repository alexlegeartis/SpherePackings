'use strict';
const assert=require('node:assert/strict');
const mod=(x,n)=>(x%n+n)%n;
const pc=x=>{let n=0;while(x){n++;x&=x-1;}return n;};
const cols=Array.from({length:12},(_,i)=>(1<<(2+(i>>2)))+(i%4));
const syn=new Uint8Array(4096);
for(let x=1;x<4096;x++)syn[x]=syn[x&(x-1)]^cols[31-Math.clz32(x&-x)];
const depth=new Uint8Array(131072);
const C=[];
for(let m=0;m<4096;m++){
 if(syn[m]===0)C.push(m);
 let d=new Int8Array(32);d.fill(99);d[0]=0;
 for(let i=0;i<12;i++){
  let a=cols[i],cost=(m>>i&1)?0:1;
  for(let j=0;j<32;j++){
   if(j>(j^a))continue;
   let l=d[j],r=d[j^a];d[j]=Math.min(l,r+cost);d[j^a]=Math.min(r,l+cost);
  }
 }
 for(let s=0;s<32;s++){let z=pc(m)+4*d[s];depth[32*m+s]=z;assert(z<=12);}
}
assert.equal(C.length,128);
assert.equal(Math.min(...C.slice(1).map(pc)),4);
const diff=(a,b)=>(((a^b)>>5)<<5)|((a^b^syn[(a>>5)&(b>>5)]^syn[b>>5])&31);
const D=(a,b)=>depth[diff(a,b)];
const P=x=>{
 let c=x.map(v=>mod(v,2));
 let m=c[0]*15+c[1]*240+c[2]*3840;
 let s=4*c[1]+8*c[2]+16*c[0]+mod(x.reduce((s,v)=>s+Math.floor(v/2),0),2);
 return m*32+s;
};
for(let i=0;i<64;i++){
 let x=[i%4,(i>>2)%4,i>>4];if(mod(x[0]+x[1]+x[2],2))continue;
 for(let a=-2;a<=2;a++)for(let b=-2;b<=2;b++)for(let c=-2;c<=2;c++){
  if(mod(a+b+c,2))continue;
  let r=a*a+b*b+c*c;if(r===0||r>=8)continue;
  const y=x.map((v,i)=>v+[a,b,c][i]);
  assert(D(P(x),P(y))+2*r>=16);
 }
}
let best=0;
for(let x=0;x<=16;x++)for(let y=0;y<=16;y++){
 const a=[x,16-x,y,16-y];
 for(let u=0;u<=16;u++)for(let v=0;v<=16;v++){
  const b=[u,v,16-u,16-v];
  for(let p=0;p<=16;p++)for(let q=0;q<=16;q++){
   const c=[p,q,16-q,16-p];let cap=0;
   for(let i=0;i<4;i++)cap+=Math.min(144,2*Math.min(a[i]*b[i],a[i]*c[i],b[i]*c[i]));
   best=Math.max(best,cap);
  }
 }
}
assert.equal(best,512);
const S=[],T=[];
for(let x=0;x<depth.length;x++){
 if(depth[x]>=10)T.push(x);
 if(depth[x]>=11)S.push(x);
 if(depth[x]===12)assert.equal(syn[x>>5],0);
}
assert.equal(S.length,972);
assert.equal(T.length,4956);
// A,B,C are the opposite pairs of an octahedron.
// Translate the first vertex to 0 and test each candidate for its opposite.
let nPairs=0,largest=0;
for(let a of T){
 const L=S.filter(x=>D(a,x)>=11);largest=Math.max(largest,L.length);
 if(L.length<4)continue;
 const n=L.length,A=new Uint8Array(n*n);
 for(let i=0;i<n;i++)for(let j=0;j<i;j++)
  A[i*n+j]=A[j*n+i]=D(L[i],L[j]);
 for(let i=0;i<n;i++)for(let j=0;j<i;j++)if(A[i*n+j]>=10){
  nPairs++;
  const Q=[];for(let k=0;k<n;k++)if(A[i*n+k]>=11&&A[j*n+k]>=11)Q.push(k);
  for(let u=0;u<Q.length;u++)for(let v=0;v<u;v++)assert(A[Q[u]*n+Q[v]]<10);
 }
}
assert.equal(nPairs,827136);
assert.equal(largest,416);
const W=[0,123390,130572,8177,104841,26525,46264,30818,46245,84831,99956,84803];
assert.equal(W.length,new Set(W).size);
for(let i=0;i<W.length;i++)for(let j=0;j<i;j++)assert(D(W[i],W[j])>=10);
// At t=17/16 first-neighbour requirement >= 3; axial requirement >9/4.
const a=17n,b=16n;
assert(4n*a**4n*b-b**5n>=3n*a**5n);
assert(4n*(4n*a**4n*b-2n*b**5n)>9n*a**5n);
// At t=3/2 first requirement >5/2, axial requirement >9/4.
assert(2n*(4n*3n**4n*2n-2n**5n)>5n*3n**5n);
assert(4n*(4n*3n**4n*2n-2n*2n**5n)>9n*3n**5n);
// (t,c)=(8/5,1/8), density ratio >1 iff c t^4 <1.
assert(8n**4n<8n*5n**4n);
console.log('PASS: 128 words and FCC equality; four-colour row cap; half-grid depths; octahedron obstruction; 12-point near-deep witness and exact scale gates');
