// Exhaust every generator subset directly in the literal Euclidean quotient Z^31/(4D31).
// Dynamic program for distance to the even-carry period, independent of code splits.
'use strict';
const W=require('./witness.json'), n=W.n;
const assert=x=>{if(!x)throw Error('verification failed')};
const wt=x=>{let r=0;while(x){x&=x-1;r++;}return r};
const B=W.binary_generator_rows31.map(b=>Array.from({length:n},(_,i)=>2*((b>>>i)&1)+(i===0?4*((wt(b)>>1)&1):0)));
B.push(Array.from({length:n},(_,i)=>((W.odd_support>>>i)&1)+(i===0?4:0)));
const v=Array(n).fill(0),N=2**B.length,seen=new Set(),hist={};
const sq8=x=>{x=(x%8+8)%8;return Math.min(x,8-x)**2;};
let prev=0,min=1e9;
for(let j=0;j<N;j++){
 const u=j^(j>>>1);
 if(j){const b=31-Math.clz32(u^prev),sgn=(u>>>b)&1?1:-1;for(let i=0;i<n;i++)v[i]+=sgn*B[b][i];}prev=u;
 let ev=0,od=1e9,k=0n,p=0;
 for(const x of v){
  let a=sq8(x),b=sq8(x+4); // add even / odd multiples of four locally
  [ev,od]=[Math.min(ev+a,od+b),Math.min(od+a,ev+b)];
  k=4n*k+BigInt((x%4+4)%4);p^=(Math.floor(x/4)&1);
 }
 k=2n*k+BigInt(p);assert(!seen.has(k));seen.add(k);
 if(j){min=Math.min(min,ev);hist[ev]=(hist[ev]||0)+1;}else assert(ev===0);
}
assert(seen.size===131072&&min===32);
// Doubling each even generator lies in 4D31. Find the subset equal to twice the odd one.
// This proves the enumerated residue classes form a group, so the quotient index is exact.
const key=v=>{let k=0n,p=0;for(const x of v){k=4n*k+BigInt((x%4+4)%4);p^=Math.floor(x/4)&1;}return 2n*k+BigInt(p)};
const even=new Set();v.fill(0);prev=0;
for(let j=0;j<N/2;j++){
 let u=j^(j>>>1);if(j){let b=31-Math.clz32(u^prev),sgn=(u>>>b)&1?1:-1;for(let i=0;i<n;i++)v[i]+=sgn*B[b][i];}prev=u;even.add(key(v));
}
assert(B.slice(0,-1).every(g=>key(g.map(x=>2*x))===0n));
assert(even.has(key(B.at(-1).map(x=>2*x)))&&!even.has(key(B.at(-1))));
assert(2n*4n**31n/BigInt(N)===2n**46n);
console.log('131072 distinct quotient classes, closed. Minimum nonzero coset squared norm',min,'period minimum 32; min-class counts',JSON.stringify(hist));
