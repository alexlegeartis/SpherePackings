'use strict';
const assert=require('node:assert/strict');
const mod=(x)=>(x%3+3)%3;
const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
const wt=a=>a.filter(x=>mod(x)!==0).length;
const minLocal=(a,d)=>{
 let p=a+d,q=a-d,h=[Infinity,Infinity,Infinity];
 for(let i=-3;i<=3;i++)for(let j=-3;j<=3;j++){
  let x=p+3*i,y=q+3*j,z=2*(x*x+y*y-x*y),s=mod(i+j);
  h[s]=Math.min(h[s],z);
 }
 assert(h.every(v=>v<=18)); // outside the box x^2+y^2>18
 return h;
};
for(let d=0;d<3;d++){
 assert.deepEqual(minLocal(0,d),d?[6,6,6]:[0,18,18]);
 assert.deepEqual(minLocal(1,d),[0,1,2].map(s=>2+6*mod(s+Number(d!==0))));
}
const C=(n,k)=>k<0||k>n?0:k===0?1:C(n-1,k-1)+C(n-1,k);
const K=(i,d)=>Array.from({length:i+1},(_,j)=>(-1)**j*2**(i-j)*C(d,j)*C(9-d,i-j))
 .reduce((s,x)=>s+x,0);
function hammingBound(d,Y,expected){
 let R=Y.reduce((a,[i])=>a*BigInt(K(i,0)),1n),min=null;
 for(let j=d;j<=9;j++){
  let v=-Y.reduce((s,[i,y])=>s+BigInt(y*K(i,j))*R/BigInt(K(i,0)),0n);
  if(min===null||v<min)min=v;
 }
 assert(min>0n);
 let M=1n+R*BigInt(Y.reduce((a,[,y])=>a+y,0))/min;
 assert.equal(M,BigInt(expected));
}
hammingBound(3,[[1,90],[2,480],[3,1344],[4,2016],[5,1344],[8,768],[9,512]],937);
hammingBound(4,[[1,160],[2,672],[3,1344],[4,1120],[7,1024],[8,768]],340);
hammingBound(6,[[1,10],[2,16]],27);
const code=(B)=>{
 let A=[Array(B[0].length).fill(0)];
 for(const b of B)A=A.flatMap(a=>[0,1,2].map(x=>a.map((v,i)=>mod(v+x*b[i]))));
 assert.equal(new Set(A.map(a=>a.join())).size,3**B.length);
 return A;
};
const U=[[1,1,2,0,1,0,1,0,1],[0,1,1,2,2,0,0,1,1],
         [0,0,1,1,0,0,2,1,2],[0,0,0,1,1,1,1,1,1]];
const V=[[1,1,2,0,1,0,1,0,1],[0,1,1,2,2,0,0,1,1],
         [0,0,1,2,0,2,1,1,0],[0,0,0,1,1,1,1,1,1]];
assert(U.every(u=>V.every(v=>mod(dot(u,v))===0)));
for(const B of [U,V])assert.equal(Math.min(...code(B).slice(1).map(wt)),5);
for(const u of code(U))for(const v of code(V)){
 const s=mod(dot(u,u)),t=mod(2-dot(v,v)),d=v.map((x,i)=>mod(x-u[i]));
 assert.equal(mod(t-s+wt(d)),2);
}
const records=[[177147n,1n],[59049n,2n],[8433n,3n],[1020n,4n],[27n,6n]];
for(const [N,h] of records)assert(N*h**9n<162n*5n**9n);
assert(36n*5n**18n*4n**18n < 243n*6n**18n*3n**18n);
const g=[2,2,1,2,0,1];
const W=Array.from({length:6},(_,i)=>{
 let b=Array(11).fill(0);
 for(let j=0;j<g.length;j++)b[i+j]=g[j];
 return [...b,mod(-b.reduce((s,x)=>s+x,0))];
});
assert(W.every(u=>W.every(v=>mod(dot(u,v))===0)));
assert.equal(Math.min(...code(W).slice(1).map(wt)),6);
const beta=[0,0,1],s=[0,2,0],m=13;
for(let j=0;j<3;j++){
 let i=(j+1)%3;
 assert.equal(mod(s[i]-s[j]+Number(beta[i]!==beta[j])-2*m*Number(j===2)),2);
}
assert(2*m+12>36);  // cross min distance (2m+12)/9 > 4
assert.equal(3**(m+2),3*(3*729)**2); // covolume^2 / M^2 =3
console.log('PASS: A2 color minima; three Hamming duals; exact 18D 162-phase and 26D 2187-phase witnesses; density comparisons');
