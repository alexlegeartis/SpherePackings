"use strict";
const assert=require("node:assert/strict");
const result=(()=>{
const n=8,unit=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>i===j?2:0));
const a=[[1,-1,-1,-1,-1,-1,-1,1],unit[0].map((x,i)=>x+unit[1][i]),unit[1].map((x,i)=>x-unit[0][i])];
for(let j=1;j<=5;j++)a.push(unit[j+1].map((x,i)=>x-unit[j][i]));
const dot=(v,w)=>v.reduce((s,x,i)=>s+x*w[i],0);
const G=a.map(r=>a.map(t=>dot(r,t)/4));
function det(M) {
 if(!M.length)return 1;
 let d=0;
 for(let j=0;j<M.length;j++)if(M[0][j])
  d+=(j%2?-1:1)*M[0][j]*det(M.slice(1).map(r=>r.filter((_,i)=>i!==j)));
 return d;
}
if(det(G)!==1)throw new Error("not a unimodular E8 basis");
const mod=(x)=>((x%3)+3)%3,P=[1,3,9,27,81,243,729,2187],N=6561;
const q=new Int16Array(N);q.fill(32767);
const x=Array(n).fill(0);
for(let p=0;p<2;p++){
 const rec=(i,Q,S)=>{
  if(i<n){for(let z=-5;z<=5;z++)if((z+10)%2===p&&Q+z*z<=32){x[i]=z;rec(i+1,Q+z*z,S+z);}return;}
  if((S%4)!==0)return;
  const m=a.map(r=>dot(r,x)/4);
  if(!m.every(Number.isInteger))throw new Error("nonintegral");
  const k=m.reduce((s,y,j)=>s+mod(y)*P[j],0);q[k]=Math.min(q[k],Q);
 };rec(0,0,0);
}
const h={};for(const y of q)h[y]=(h[y]??0)+1;
const first=q.findIndex(q=>q===32), seen=new Uint8Array(N), Q=[first];seen[first]=1;
for(let p=0;p<Q.length;p++){
 const code=Q[p],m=P.map(t=>Math.floor(code/t)%3);
 for(let i=0;i<n;i++){
  const w=m.reduce((s,v,j)=>s+P[j]*mod(v-m[i]*G[i][j]),0);
  if(!seen[w]){seen[w]=1;Q.push(w);}
 }
}
return {G,h,first,orbit:Q.length,only:Q.every(t=>q[t]===32)};
})();
assert.deepEqual(result.h,{0:1,8:240,16:2160,24:2240,32:1920});
assert(result.only);
assert.equal(result.orbit,1920);
console.log("PASS: the 1920 deepest sqrt(3)E8 discriminant classes are one Weyl orbit.");
