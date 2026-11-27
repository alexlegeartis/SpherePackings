'use strict';
const assert=require('node:assert/strict');
const obs=[];
for(let n=2;n<=7;n++)for(let p=0;p<=1;p++){
 const V=p?[1,3,5,7]:[0,1,2,3,4],goal=(p?8:2)*n;
 function go(A,m,r){
  if(A.length===V.length){
   if(m||r)return;
   const S=new Set();
   for(let i=0;i<A.length;i++)for(let j=i;j<A.length;j++)
    if(A[i]&&A[j]>Number(i===j)){
     S.add((V[i]+V[j])/(p?2:1));S.add((V[j]-V[i])/(p?2:1));
    }
   const D=[];for(let j=1;j*j<=2*n;j++)if(!S.has(j))D.push(j);
   if(D.length)obs.push([n,p,A,D]);
   return;
  }
  const v=V[A.length];
  for(let i=0;i<=m&&i*v*v<=r;i++)go(A.concat(i),m-i,r-i*v*v);
 }
 go([],16,goal);
}
assert.deepEqual(obs,[[2,0,[15,0,1,0,0],[1]],[2,1,[16,0,0,0],[2]],
 [4,0,[14,0,2,0,0],[1]],[5,0,[6,10,0,0,0],[3]],
 [6,0,[4,12,0,0,0],[3]],[6,0,[13,0,3,0,0],[1,3]],
 [7,0,[2,14,0,0,0],[3]]]);
assert(2n*15n**15n < 16n**15n);
console.log('PASS: complete small-normal magnitude patterns and exact density comparison');
