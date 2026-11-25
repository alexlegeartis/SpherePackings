'use strict';
// Exact integer mesh replay. Run with node verify.js (single process).
const assert=require('node:assert/strict'),q=9,K=81,tau=540;
const word=Array.from({length:K},(_,a)=>[0,...Array.from({length:4},(_,j)=>Math.floor(a/3**j)%3)]);
const add=(a,b,sgn=1)=>{let w=word[a].map((x,j)=>(3+x+sgn*word[b][j])%3);return w.slice(1).reduce((s,x,i)=>s+3**i*((3+x-w[0])%3),0)};
const minus=word.map((_,a)=>word.map((_,b)=>add(a,b,-1))),plus=word.map((_,a)=>word.map((_,b)=>add(a,b)));
function ok(a,b){const f=[0,0,0];for(let i=0;i<5;i++)f[(3+word[a][i]-word[b][i])%3]++;return Math.max(...f)<=2;}
let cl=[];
function grow(tail,cand){if(tail.length===6){cl.push(tail);return;}
 for(let j=0;j<cand.length;j++){if(tail.length+cand.length-j<6)break;
 const v=cand[j];grow([...tail,v],cand.slice(j+1).filter(w=>ok(v,w)));}}
grow([0],word.map((_,i)=>i).filter(i=>ok(0,i)));assert.equal(cl.length,12);
const mk=x=>{let a=[0,0,0];for(const z of x)a[z>>5]|=1<<(z%32);return a},lst=new Map;
for(const a of cl)for(const b of cl){const dd=new Set(a.flatMap(i=>b.map(j=>minus[i][j])));
 // Enumerate global translations explicitly, then deduplicate by masks.
 for(let t=0;t<K;t++){let m=mk([...dd].map(u=>plus[u][t]));lst.set(m.join(),m);}}
const dds=[...lst.values()];
let local=[];for(let a=q;a>=0;a--)for(let b=q-a;b>=0;b--){let c=q-a-b;if(a<b||a<c)continue;
 let p=[a,b,c];local.push(p.map(x=>2*(q*q+a*a+b*b+c*c-2*q*x)));}
let profiles=0;function verify(typ){profiles++;let good=[0,0,0];
 for(let u=0;u<K;u++){let best=1e9;for(let s=0;s<3;s++){let val=0;for(let j=0;j<5;j++)val+=local[typ[j]][(word[u][j]+s)%3];best=Math.min(best,val);}
 if(best<=tau)good[u>>5]|=1<<(u%32);}
 for(const m of dds)assert.ok((good[0]&m[0])||(good[1]&m[1])||(good[2]&m[2]));}
function run(xs,lo){if(xs.length===5){verify(xs);return;}for(let i=lo;i<local.length;i++)run([...xs,i],i);}
run([],0);assert.equal(profiles,42504);
// The triangle centroid has squared distance 4/9 to each of the three local labels.
assert.equal(5*4*3*q*q/9,tau);
console.log({normalized18SiteCodes:cl.length,translatedDifferenceMasks:dds.length,profiles,maxSquaredMeshDistance:'20/9',continuousUpper:'2000/729'});