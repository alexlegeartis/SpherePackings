'use strict';
// Conway–Sloane, Quaternary constructions (1994), pentacode cyclic-negation form.
const Z=['01112','03110','21310','21132'],G=[0,1,3,2],C=new Set();
for(const p of Z)for(let k=0;k<5;k++)for(const s of [1,3]){
 let u=0;for(let i=0;i<5;i++)u|=G[(+p[(i+k)%5]*s)%4]<<(2*i);C.add(u);
}
const pc=x=>{let n=0;for(;x;x&=x-1)n++;return n;};
if(C.size!==40)throw Error('motif size');
for(let x of C)for(let y of C)if(x!==y&&pc(x^y)<4)throw Error('code distance');
const D=new Set([...C].flatMap(x=>[...C].map(y=>x^y)));
if(D.size!==261)throw Error('difference count');
let nontrivial=0;
for(let a=0;a<1024;a++)if(!D.has(a)){
 let S=[...D].map(d=>a^d).filter(s=>pc(s)<=3).sort((u,v)=>pc(u)-pc(v));
 function f(i,k,u){
  if(k===4)return true;
  for(let j=i;j<S.length;j++)if(!(S[j]&u)&&f(j+1,k+1,u|S[j]))return true;
  return false;
 }
 if(!f(0,0,0))throw Error('missing disjoint supports '+a);
 nontrivial++;
}
console.log({vertices:1024,differenceWords:D.size,fourDisjointWitnessCases:nontrivial});
// If w_i is distance from nearest corner a_i, then w_i in [0,1/2].
// Cost to d equals sum w_i^2 plus sum_{i in a xor d} (1-2w_i).
// If a in D the result is <=2.5 directly. Otherwise the checked four
// disjoint-support witnesses give min <= mean.
// sum (w_i^2+(1-2w_i)/4) <= 10/4. The cube center attains it.
