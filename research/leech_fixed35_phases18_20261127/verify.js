'use strict';
const assert=require('node:assert/strict');
const A=[[2,0,-1,0,0,0],[0,2,0,-1,0,0],[-1,0,2,-1,0,0],
         [0,-1,-1,2,-1,0],[0,0,0,-1,2,-1],[0,0,0,0,-1,2]];
const B=[[4,3,5,6,4,2],[3,6,6,9,6,3],[5,6,10,12,8,4],
         [6,9,12,18,12,6],[4,6,8,12,10,5],[2,3,4,6,5,4]];
const det=(M)=>M.length===0?1:M[0].reduce((acc,x,i)=>acc+((-1)**i)*x*det(
 M.slice(1).map(r=>r.filter((_,j)=>j!==i))),0);
const dot=(x,y)=> x.reduce((a,c,i)=>a+c*y[i],0);
const prod=(M,x)=> M.map(r=>dot(r,x));
const mod3=(x)=>(x%3+3)%3;
const norm=(M,x)=>dot(x,prod(M,x));
assert.equal(det(A),3);
assert.equal(det(B),243);
assert.equal(det(A.slice(1).map(r=>r.slice(1))),4);
for(let j=0;j<6;j++) assert.deepEqual(prod(A,B.map(r=>r[j])),
 A.map((_,i)=>3*Number(i===j)));
assert(B.every(r=>r.every((v,j)=>mod3(v- r[0]*B[0][j])===0))); // B00=1 (mod 3)
function enumBox(width,f,x=[]){
 if(x.length===width.length){f(x);return;}
 for(let j=-width[x.length];j<=width[x.length];j++)enumBox(width,f,x.concat(j));
}
let k4=0;
enumBox([1,1,1,1,1,1], x => {
 const n=norm(B,x);
 if(n===4)k4++;
 if(x.some(v=>v!==0))assert(n>=4);
});
assert.equal(k4,54);
// Projected minimal-shell Lagrange moments: 7 distinct nodes given last count.
const shell=[6480,104976,65610,19440,0,0,k4];
let num=196560n,den=1n;
for(let j=0;j<=5;j++){
 const total=shell.reduce((s,n,i)=>s+BigInt(n)*BigInt(i)**BigInt(j),0n);
 assert.equal(total*den,num);
 num*=6n*BigInt(3+j); den*=BigInt(12+j);
}
const mins={};
enumBox(B.map((r,i)=>Math.floor(Math.sqrt(2*r[i]))), x=>{
 const h=prod(A,x), n=dot(h,x);
 if(n>6)return;
 const label=h.map(mod3).join('');
 mins[label]=Math.min(mins[label]??Infinity,n);
});
const M=Object.values(mins);
assert.equal(M.length,243);
assert.deepEqual([0,2,4,6].map(n=>M.filter(x=>x===n).length),[1,72,90,80]);
const V=Array.from({length:243},(_,i)=>{
 const r=[0];for(let j=0;j<5;j++){r.push(i%3);i=Math.floor(i/3);}return r;
});
const iso=V.filter(v=>v.some(x=>x>0)&&mod3(norm(A,v)/2)===0);
assert.equal(iso.length,80);
const spec={};
for(const h of V){
 const n0=iso.filter(v=>mod3(dot(v,h))===0).length;
 const ev=(3*n0-80)/2;
 spec[ev]=(spec[ev]||0)+1;
}
assert.deepEqual(spec,{80:1,8:90,'-1':80,'-10':72});
const clique=[[0,0,0,0,0,0]];
for(let i of [0,2,3,4,5]){
 let v=clique.at(-1).slice();v[i]++;clique.push(v);
}
for(let i=0;i<6;i++)for(let j=0;j<i;j++)assert.equal(
 norm(A,clique[i].map((x,k)=>x-clique[j][k])),2);
for(const [m,k] of [[243n,6n],[27n,8n],[1n,12n]])
 assert(m*k**9n<6n*10n**9n);
// D_opt^2 = 6^2/243 * (5/6)^18; B*_18^2 = (3/4)^18.
assert(36n*5n**18n*4n**18n < 243n*6n**18n*3n**18n);
console.log('PASS: E6 matrices; 243 minimum norms; exact 11-design jets; isotropic spectrum; six-clique; strict density comparisons');
