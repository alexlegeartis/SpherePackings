'use strict';
const assert=require('node:assert/strict');
const T=[[7,8,[[1,1],[12,11]]],[6,16,[[12,11],[120,91]]],[5,64,[[120,91],[3,2]]],[4,128,[[3,2],[39,20]]],[3,1024,[[39,20],[5,2]]],[2,2048,[[5,2],[3,1]]],[0,32768,[[3,1]]]];
for(const [k,M,I] of T)for(const [p,q] of I){const z=8n*BigInt(q)-BigInt(k*p);assert(z>0n&&z*BigInt(p)**14n>=BigInt(M/8)**2n*BigInt(q)**15n);}
console.log('PASS: 13 exact slab endpoints');
