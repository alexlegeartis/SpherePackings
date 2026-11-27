"use strict";
const assert=require("node:assert/strict");
const log=(...v)=>console.log(...v);
const mod=(x,p)=>((x%p)+p)%p;
function det(G) {
 const n=G.length;
 if(n===0)return 1;
 if(n===1)return G[0][0];
 let s=0;
 for(let j=0;j<n;j++) if(G[0][j])
   s+=((j&1)?-1:1)*G[0][j]*det(G.slice(1).map(r=>r.filter((_,i)=>i!==j)));
 return s;
}
function adj(G) {
 const n=G.length;
 return G.map((r,i)=>r.map((_,j)=>(((i+j)&1)?-1:1)*det(G.filter((_,a)=>a!==j).map(t=>t.filter((_,b)=>b!==i)))));
}
function rankD(d){
 return Array.from({length:d},(_,i)=>Array.from({length:d},(_,j)=>{
  if(i===j)return 6;
  const u=Math.min(i,j),v=Math.max(i,j);
  return (v===u+1&&v!==d-1)||(u===d-3&&v===d-1)?-3:0;
 }));
}
const extend=(G,m)=>G.map((r,i)=>[...r,m[i]]).concat([[...m,6]]);
const d4=rankD(4), d5=rankD(5), d6=rankD(6);
const g6=extend(d5,[-1,0,0,0,3]);
const g7=extend(g6,[-3,0,0,0,0,-2]);
const g8=extend(g7,[-3,0,0,0,0,-2,3]);

function atlas(G, cap=6) {
 const n=G.length,B=[],D=[],m=Array(n).fill(0),result=new Map();
 for(let l=0;l<n;l++){
  const H=G.slice(0,l+1).map(r=>r.slice(0,l+1));
  B.push(adj(H));D.push(det(H));
 }
 let count=0;
 function visit(r){
  const A=B[r],mm=D[r], t=0;
  let L=0,C=0;for(let i=0;i<r;i++){
   L+=A[i][r]*m[i];
   for(let j=0;j<r;j++)C+=m[i]*m[j]*A[i][j];
  }
  for(let k=-6;k<=6;k++){
   const q=A[r][r]*k*k+2*L*k+C;
   if(q>cap*mm)continue;
   m[r]=k;
   if(r<n-1){visit(r+1);continue;}
   count++;
   const h=A.map(row=>mod(row.reduce((s,v,j)=>s+v*m[j],0),mm)).join(",");
   if(!result.has(h)||result.get(h).q>q)result.set(h,{q,m:[...m]});
  }
 }
 visit(0);
 if(result.size!==D[n-1])throw new Error("missing "+JSON.stringify([G,D[n-1],result.size]));
 const values=[...result.values()].sort((a,b)=>b.q-a.q);
 return {d:D[n-1],count,values,g:G};
}

function roots(G,det,adj){
 const n=G.length,D=det(G), A=adj(G),cap=6;
 const B=G.map((_,i)=>{let k=0;while((k+1)**2*D<=cap*A[i][i])k++;return k;});
 const tail=[], ds=[];
 for(let i=0;i<n;i++){
  const P=G.slice(0,i).map(r=>r.slice(0,i)),d=det(P),J=adj(P);
  ds.push(d);
  const M=[];
  for(let j=i;j<n;j++){
   let row=[];for(let k=i;k<n;k++){
    let x=d*G[j][k];
    for(let a=0;a<i;a++)for(let b=0;b<i;b++)x-=G[j][a]*J[a][b]*G[b][k];
    row.push(x);
   }M.push(row);
  }tail.push(M);
 }
 const result=[], a=Array(n).fill(0);
 const dot=(v,w)=>v.reduce((s,x,i)=>s+x*w.reduce((t,y,j)=>t+y*G[i][j],0),0);
 function go(i){
  if(i<0){if(a.find(x=>x!==0)>0&&dot(a,a)===6)result.push([...a]);return;}
  const M=tail[i];
  let L=0,C=0;
  for(let j=1;j<M.length;j++){
   L+=a[i+j]*M[0][j];
   for(let k=1;k<M.length;k++)C+=a[i+j]*a[i+k]*M[j][k];
  }
  for(let t=-B[i];t<=B[i];t++){
   if(M[0][0]*t*t+2*L*t+C>cap*ds[i])continue;
   a[i]=t;go(i-1);
  }
 }
 go(n-1);return {R:result,dot};
}
function hasA3(R,dot){
 for(let i=0;i<R.length;i++)for(let j=0;j<i;j++)if(Math.abs(dot(R[i],R[j]))===3){
  for(let k=0;k<R.length;k++)if(k!==i&&k!==j&&Math.abs(dot(R[i],R[k]))===3&&dot(R[j],R[k])===0)return true;
  for(let k=0;k<R.length;k++)if(k!==i&&k!==j&&Math.abs(dot(R[j],R[k]))===3&&dot(R[i],R[k])===0)return true;
 }
 return false;
}
function hasD4(R,dot){
 for(let i=0;i<R.length;i++){
  const S=R.filter((r,j)=>j!==i&&Math.abs(dot(r,R[i]))===3);
  for(let a=0;a<S.length;a++)for(let b=0;b<a;b++)if(dot(S[a],S[b])===0)
   for(let c=0;c<b;c++)if(dot(S[a],S[c])===0&&dot(S[b],S[c])===0)return true;
 }return false;
}

function eq(G,H){
 const {R,dot}=roots(G,det,adj), R2=R.concat(R.map(v=>v.map(x=>-x))), n=G.length;
 const out=[];
 function go(i){
  if(i===n)return Math.abs(det(out))===1;
  for(const v of R2)if(out.every((w,j)=>dot(w,v)===H[i][j])){out.push(v);if(go(i+1))return true;out.pop();}
  return false;
 }
 return go(0)?out:undefined;
}
function gate(d,m,r){
 d=BigInt(d);m=BigInt(m);r=BigInt(r);
 return BigInt(({1:2,2:3,3:4,4:4,5:4,6:3,7:2,8:1})[Number(r)])*m**r<=21870n*d**(r-1n)*2n**r;
}

const keyOf=(A,D,m)=>A.map(row=>mod(row.reduce((s,v,j)=>s+v*m[j],0),D)).join(",");
const mincheck=(A,D,P,m,q)=>{
 for(let k=2;k*k*q<6*D;k++){
  const t=P.get(keyOf(A,D,m.map(x=>x*k)));
  if(t+k*k*q<6*D)return false;
 }
 return true;
};

const eqTo=(G,H)=>Boolean(eq(G,H));
const A3=rankD(3);
const C2=[[6,2],[2,6]], A2=[[6,3],[3,6]];
const C3=[[6,-2,-2],[-2,6,-2],[-2,-2,6]];
const T3=[[6,-3,-3],[-3,6,2],[-3,2,6]];
const T4a=extend(T3,[-1,-2,-2]), T4b=extend(T3,[-1,-1,3]);
const T5=extend(T4a,[-3,2,2,-2]), T6=extend(T5,[-3,1,1,-1,1]);
const rootsD4=(G)=>{const t=roots(G,det,adj);return hasD4(t.R,t.dot);};
const rootsA3=(G)=>{const t=roots(G,det,adj);return hasA3(t.R,t.dot);};
const certify=(G)=>assert(G.every((r,i)=>r[i]===6&&r.every((x,j)=>Math.abs(x)<=6&&x===G[j][i])));
function eligible(G) {
 certify(G);
 const o=atlas(G),D=o.d,A=adj(G), rank=9-G.length;
 assert(o.values[0].q<6*D);
 assert(!gate(D,8*D-o.values[0].q,rank),"higher lift may be below the Hermite cutoff");
 const P=new Map(o.values.map(t=>[keyOf(A,D,t.m),t.q]));
 return o.values.filter(t=>gate(D,6*D-t.q,rank)).map(t=>{
  const q=6*D-t.q,H=extend(G,t.m);
  assert.equal(det(H),q);
  assert(mincheck(A,D,P,t.m,q));
  return {...t,H,q,D};
 });
}
// Mordell's bound gamma_9^7 <= gamma_8^8 = 2^8 suffices to force min=6.
assert(8n**63n>2n**72n*21870n**7n);
// Rootless-in-the-scaled-sense branch.
const L1=eligible([[6]]);
assert.deepEqual(L1.map(t=>Math.abs(t.m[0])).sort(),[2,2,3]);
for(const x of L1)if(Math.abs(x.m[0])!==3)assert(eqTo(C2,x.H));
const badC2=eligible(C2).filter(t=>t.m.every(x=>Math.abs(x)!==3));
assert.equal(badC2.length,1);
for(const x of badC2)assert(eqTo(C3,x.H));
const L3=eligible(C3);
assert.equal(L3.length,18);
for(const x of L3)
 assert([...x.m,x.m.reduce((s,v)=>s+v,0)].some(v=>Math.abs(v)===3));
// √3A2 but no √3A3.
const L2=eligible(A2);
assert.equal(L2.length,8);
for(const x of L2)assert(eqTo(A3,x.H)||eqTo(T3,x.H));
assert.equal(det(T3),120);
const aD5=atlas(rankD(5));
assert.deepEqual([...new Set(aD5.values.map(x=>x.q))].slice(0,4), [45,37,36,32].map(x=>81*x));
for(const x of eligible(T3))
 assert(rootsA3(x.H)||eqTo(T4a,x.H)||eqTo(T4b,x.H));
assert.deepEqual([det(T4a),det(T4b),det(T5),det(T6)],[400,448,1280,3072]);
for(const G of [T4a,T4b])for(const x of eligible(G))assert(rootsA3(x.H)||eqTo(T5,x.H));
for(const x of eligible(T5))assert(rootsA3(x.H)||eqTo(T6,x.H));
assert.equal(atlas(T6).values[0].q,9728);
assert.equal(eligible(T6).length,0);
// √3A3 but no √3D4, up to isometry.
const st=[[0,0,-3],[0,-3,2],[0,-3,1]].map(m=>extend(A3,m));
const starts=eligible(A3);
assert.equal(starts.length,19);
for(const x of starts)assert(rootsD4(x.H)||st.some(g=>eqTo(g,x.H)));
function advance(deck){
 const next=[];let total=0, withD4=0;
 for(const G of deck)for(const x of eligible(G)){
  total++;
  if(rootsD4(x.H)){withD4++;continue;}
  const {R}=roots(x.H,det,adj);
  if(!next.some(t=>t.D===x.q&&t.roots===R.length&&eqTo(t.G,x.H)))
   next.push({D:x.q,roots:R.length,G:x.H});
 }
 return {total,withD4,next};
}
let deck=st;
for(const [t,w,rank,D] of [[60,22,5,[1320,1350,1350,1428]],[85,56,6,[3375,3648,3795,3920]],[0,0,7,[]]]){
 if(rank===7) assert.deepEqual(deck.map(G=>{const v=atlas(G);return [v.d,v.values[0].q];}),[[3375,12150],[3648,13184],[3795,13398],[3920,14112]]);
 const o=advance(deck);
 assert.equal(o.total,t);assert.equal(o.withD4,w);
 assert.deepEqual(o.next.map(x=>x.D).sort((x,y)=>x-y),D);
 log("rank "+rank+": "+t+" eligible, "+w+" contain scaled D4; residual classes:",o.next);
 deck=o.next.map(x=>x.G);
}
log("PASS: every even integral rank-nine lattice with det <=21870 contains scaled D4; apply the root-seed theorem.");
