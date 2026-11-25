'use strict';
const assert=require('node:assert/strict'), w=require('./witness.json');
const wt=x=>{let n=0;for(;x;x&=x-1)n++;return n};
const C=[0];for(const g of w.binary_generator_rows31)C.push(...C.map(x=>x^g));
const S=(2**31-1)^w.odd_support;
const count=C.filter(c=>wt(c)===8).length, hist=Array(8).fill(0);
for(const c of C)hist[wt(c&S)]++;
assert.equal(count,681);assert.deepEqual(hist,[1024,0,21504,0,35840,0,7168,0]);
const H=Array.from({length:31},(_,i)=>w.binary_generator_rows31.reduce((s,g,k)=>s|((g>>>i&1)<<k),0));
function zeroSupports(k){let n=0;function go(i,d,s){if(!d){n+=s===0;return}
 for(let q=i;q<=31-d;q++)go(q+1,d-1,s^H[q]);}go(0,k,0);return n;}
const dual=Array.from({length:7},(_,i)=>zeroSupports(i+1));
assert.deepEqual(dual,[0,0,0,0,0,0,31]);
const neighbors=4*(31*30/2)+count*128+hist[0]*24+hist[2]*2;
assert.equal(neighbors,156612);console.log({neighbors,dual,complementWeightCounts:hist});
