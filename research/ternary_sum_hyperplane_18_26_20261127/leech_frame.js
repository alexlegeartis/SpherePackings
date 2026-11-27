'use strict';
const assert=require('node:assert/strict');
const W=[[4,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-4,0,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[3,-1,1,1,1,-1,-1,-1,1,-1,1,-1,1,1,1,1,1,1,1,1,1,1,1,-1],[1,-3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,0,2,0,0,2,0,0,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0],[2,2,0,2,0,0,0,0,2,0,2,0,2,0,0,2,0,0,0,0,0,0,2,0]];
const U=[[2,0,-2,0,0,-2,0,0,-2,2,-2,2,2,0,0,0,0,0,0,0,0,0,0,0],[0,0,-2,-2,0,0,0,-2,2,-2,0,2,0,2,0,0,0,0,0,0,0,0,0,2],[0,0,0,0,4,0,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,4,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]];
const T=[[0,0,0,0,0,-1],[1,0,-1,0,0,1],[-1,-1,1,0,0,0],[-1,-1,0,1,0,0],[-1,0,0,1,-1,0],[-1,0,0,0,0,0]];
const B=[[4,3,5,6,4,2],[3,6,6,9,6,3],[5,6,10,12,8,4],[6,9,12,18,12,6],[4,6,8,12,10,5],[2,3,4,6,5,4]];
const J=[[2,1,0,1,0,-1],[1,2,1,0,1,-1],[0,1,2,-1,1,-1],[1,0,-1,2,-1,0],[0,1,1,-1,2,-1],[-1,-1,-1,0,-1,2]];
const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
const G=W.map(w=>W.map(v=>dot(w,v)/8));
assert.deepEqual(G,T.map(r=>T.map(s=>dot(r,B.map(b=>dot(b,s))))));
const det=(M)=>M.length===0?1:M[0].reduce((s,x,i)=>s+(-1)**i*x*det(M.slice(1).map(r=>r.filter((_,j)=>j!==i))),0);
assert.equal(Math.abs(det(T)),1);
assert.equal(det(G),243);
for(let i=0;i<6;i++)for(let j=0;j<6;j++) assert.equal(dot(G[i],J.map(r=>r[j])),3*Number(i===j));
const q=[0,-2,2,2,1,0], p=J.map(r=>dot(r,q));
assert.equal(dot(p,q),6);
const pc=n=>{n-=(n>>>1)&1431655765;n=(n&858993459)+((n>>>2)&858993459);return ((n+(n>>>4)&252645135)*16843009)>>>24;};
const roots=[];
for(let i=0;i<24;i++)for(let j=0;j<i;j++)for(let a of [-4,4])for(let b of [-4,4]){
 let x=Array(24).fill(0);x[i]=a;x[j]=b;roots.push(x);
}
for(let n=0;n<4096;n++){
 let b=0;for(let j=0;j<12;j++)if(n>>j&1)b^=0xAE3<<j;
 b|=(pc(b)&1)<<23;
 if(pc(b)===8){
  let S=[];for(let j=0;j<24;j++)if(b>>j&1)S.push(j);
  for(let s=0;s<256;s++)if(!(pc(s)&1)){
   let x=Array(24).fill(0);for(let j=0;j<8;j++)x[S[j]]=2-4*(s>>j&1);roots.push(x);
  }
 }
 for(let k=0;k<24;k++)roots.push(Array.from({length:24},(_,j)=>(j===k?-3:1)*(1-2*(b>>j&1))));
}
assert.equal(roots.length,196560);
for(let r of W.concat(U))assert(roots.some(x=>x.every((v,i)=>v===r[i])));
for(let i=0;i<6;i++)for(let j=0;j<6;j++)assert.equal(dot(U[i],W[j]),8*Number(i===j));
const shell=roots.filter(r=>W.every((w,i)=>dot(r,w)===8*q[i]));
assert.equal(shell.length,27);
const z=shell.map(r=>r.map((x,j)=>3*x-W.reduce((s,w,i)=>s+p[i]*w[j],0)));
const visited=new Set();
let nt=0;
for(let i=0;i<27;i++){
 assert.equal(dot(z[i],z[i]),144);
 if(visited.has(i))continue;
 let Q=z.map((v,j)=>dot(v,z[i])===-72?j:-1).filter(x=>x>=0);
 assert.equal(Q.length,2);
 const tri=[i,...Q];
 assert.equal(dot(z[Q[0]],z[Q[1]]),-72);
 for(let j=0;j<27;j++)if(!tri.includes(j)) for(let k of tri)assert.equal(dot(z[k],z[j]),0);
 for(let j of tri)visited.add(j);nt++;
}
assert.equal(nt,9);
assert.equal(3**11,243*27**2);
console.log('PASS: explicit primitive Leech rank-six E6* section; 27 shortest contacts form nine orthogonal A2 triangles; period index27');
