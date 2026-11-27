"use strict";
const A=(t)=>{if(!t)throw Error("assertion failed");};
const range=(n)=>Array.from({length:n},(_,i)=>i);
const pc=(x)=>{let a=0;for(;x;x&=x-1)a++;return a;};
const syn=(w,n)=>range(n).reduce((a,i)=>a^((w>>i&1)?16+i:0),0);
const checks=(n)=>range(1<<n).map(w=>syn(w,n));
const hr12=checks(12),hr16=checks(16);
const H=range(4096).filter(w=>!hr12[w]),C=range(65536).filter(w=>!hr16[w]);
A(H.length===128&&C.length===2048);
A(Math.min(...H.filter(x=>x).map(pc))===4);
A(Math.min(...C.filter(x=>x).map(pc))===4);
const residues=new Set([1,3,4,5,9]);
const rows=[0,...range(11).map(i=>range(11).reduce((a,j)=>a|(residues.has((j-i+11)%11)?0:2<<j),0))];
const W=new Set;
for(let i=0;i<12;i++)for(let j=0;j<i;j++){const w=rows[i]^rows[j];A(pc(w)===6);W.add(w);W.add(w^4095);}
for(let i=0;i<6;i++){W.add(3<<(2*i));W.add((3<<(2*i))^4095);}
A(W.size===144);
const sum=new Set;
for(const x of W)for(const y of W)sum.add(x^y);
A([...sum].every(x=>x===0||pc(x)>=4));
for(let x=0;x<4096;x++)if([0,4,8,12].includes(pc(x)))A(sum.has(x));
const even=new Set;
for(const w of W)for(const x of H)even.add(w^x);
A(even.size===2048&&[...even].every(x=>pc(x)%2===0));
const dist=Array(32).fill(99);
for(let w=0;w<4096;w++)dist[hr12[w]]=Math.min(dist[hr12[w]],pc(w));
A(range(16).every(s=>dist[s]===(s===0?0:2)));

const mod=(x,m)=>((x%m)+m)%m;
function label(z){
 const e=mod(z[0]+z[1],2);
 const u=[z[0]-e,z[1],z[2]-e,z[3]];
 const k=[(u[0]+u[1])/2,(u[0]-u[1])/2,(u[2]+u[3])/2,(u[2]-u[3])/2];
 A(k.every(Number.isInteger));
 const parity=mod(k.reduce((a,b)=>a+b,0),2);
 let low=k.reduce((a,b,i)=>a^(mod(b,2)?i:0),0);
 if(parity)low=[0,2,1,3][low];
 return [e,low^(parity?28:0)];
}
const box=(b)=>{
 const v=[];for(const a of b)for(const c of b)for(const d of b)for(const e of b)v.push([a,c,d,e]);return v;
};
const sq=z=>z.reduce((a,b)=>a+b*b,0), evenZ=z=>mod(z.reduce((a,b)=>a+b,0),2)===0;
const sites=box(range(4)).filter(evenZ);
const short=box([-2,-1,0,1,2]).filter(z=>evenZ(z)&&0<sq(z)&&sq(z)<8);
A(sites.length===128);
for(const z of sites){
 const s=label(z);
 for(let j=0;j<4;j++){const t=[...z];t[j]+=4;A(String(label(t))===String(s));}
 for(const d of short){
  const t=label(z.map((v,i)=>v+d[i]));
  const h=s[0]===t[0]?dist[s[1]^t[1]]:3;
  A(sq(d)+2*h>=8);
 }
}
A(label([0,0,0,0])[1]===0);
A((label([1,1,0,0])[1]^label([1,-1,0,0])[1])!==label([2,0,0,0])[1]);
const x2=range(16).map(i=>(i<8?1:0)+((i===0||i===8)?2:0));
function bankDist4(offset){
 let m=Infinity;
 for(const w of C){
  let q=0;
  for(let i=0;i<16;i++){
   const a=mod(x2[i]-offset-2*(w>>i&1),4),d=Math.min(a,4-a);
   q+=d*d;
  }
  m=Math.min(m,q);
 }
 return m;
}
A(bankDist4(0)===12&&bankDist4(1)===12);
console.log("PASS: Witt/Hamming difference sets, nonlinear 16D equality packing, and 17D equality shift; "+sites.length*short.length+" short normal checks");

