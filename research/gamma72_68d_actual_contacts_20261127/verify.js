"use strict";
const fs=require("fs");
const D=JSON.parse(fs.readFileSync(__dirname+"/data.json","utf8"));
const ok=(p,msg="assertion")=>{if(!p)throw new Error(msg);};
const N=24;
const range=n=>Array.from({length:n},(_,i)=>i);
const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
const tr=A=>A[0].map((_,j)=>A.map(r=>r[j]));
const prod=(A,B)=>A.map(r=>tr(B).map(c=>dot(r,c)));
const eq=(A,B)=>JSON.stringify(A)===JSON.stringify(B);
const add=(a,b)=>a.map((x,i)=>x+b[i]);
const sub=(a,b)=>a.map((x,i)=>x-b[i]);
const scale=(a,b)=>a.map(x=>x*b);
const md=(a,b)=>((a%b)+b)%b;
const {F,A,B,catalogue_gram:G,normals68:R,witnesses:W,colour_groups:C}=D;
const I=range(N).map(i=>range(N).map(j=>Number(i===j)));
const Z=I.map(r=>scale(r,0));
const stack=(X,Y,Z)=>X.map((r,i)=>r.concat(Y[i],Z[i]));
const T=stack(I,I,A).concat(stack(Z,B,B),stack(Z,Z,I.map(r=>scale(r,2))));
const GF=stack(F,Z,Z).concat(stack(Z,F,Z),stack(Z,Z,F));
ok(eq(prod(prod(T,GF),tr(T)),G.map(r=>scale(r,2))),"catalogue factorization");
ok(eq(prod(A,A),A.map((r,i)=>sub(r,scale(I[i],2)))),"A polynomial");
ok(eq(A.map((r,i)=>add(r,B[i])),I),"B=I-A");
ok(eq(prod(prod(A,F),tr(A)),F.map(r=>scale(r,2))),"A norm doubling");
ok(eq(prod(prod(B,F),tr(B)),F.map(r=>scale(r,2))),"B norm doubling");
if(process.argv[2]){
  const s=fs.readFileSync(process.argv[2],"utf8");
  const frag=s.split('<a NAME="GRAM"><STRONG>GRAM</STRONG></a><br>')[1].split("<p><li>")[0];
  const num=frag.replace(/<[^>]+>/g," ").match(/-?\d+/g).map(Number);
  ok(num[0]===72 && eq(num.slice(1),G.flat()),"live/archive catalogue agreement");
}
const PR=prod(R,F);
const gram=prod(PR,tr(R)).map(r=>scale(r,2));
ok(eq(gram,range(4).map(i=>range(4).map(j=>i===j?8:-2))),"normal Gram 10I-2J");
ok(eq(R.map(r=>[23,11,4,3].map(i=>r[i])),range(4).map(i=>range(4).map(j=>Number(i===j)))),"primitive minor");
const key=v=>{const s=v.reduce((x,y)=>x+y,0);return v.map(x=>md(x+s,10)).join("");};
const phases=new Map();
for(const a of range(10))for(const b of range(10))for(const c of range(10))for(const d of range(10)){
  const v=[a,b,c,d];phases.set(key(v),v);
}
ok(phases.size===2000,"whole discriminant group");
const upper=new Map();
for(const [b,w] of W){
  ok(eq(PR.map(r=>dot(r,w)),b),"witness normal inner products");
  const q=dot(w,F.map(r=>dot(r,w)));
  let ambient;
  if(q===4 || q===6)ambient=8;
  else if(q===8 || q===10)ambient=10;
  else if(q===16 && w.every(z=>z%2===0))ambient=8;
  else throw Error("invalid Leech witness");
  const h=dot(b,b)+b.reduce((x,y)=>x+y,0)**2;
  ok(h%2===0 && 5*ambient-h/2>0,"horizontal upper credit");
  for(const a of [b,scale(b,-1)]){
    const k=key(a), old=upper.get(k)??40;
    upper.set(k,Math.min(old,5*ambient-h/2));
  }
}
const old=[[0,0,0,0],...range(4).map(i=>range(4).map(j=>Number(i===j)))];
const oldKeys=new Set(old.map(key));ok(oldKeys.size===5,"old motif");
const U=(v,w)=>upper.get(key(sub(v,w)))??40;
ok(U(old[1],old[0])===39,"inherited upper contact at mu=39/5");
const verify=(d,n)=>ok(BigInt(n)*BigInt(d)**34n<=5n*39n**34n,"exact relative density");
let report=[];
for(let d=39;d>=33;d--){
  const remain=new Set([...phases.keys()].filter(k=>!oldKeys.has(k)&&old.every(a=>U(phases.get(k),a)>=d)));
  const size=remain.size;
  let groups=C[d];
  if(groups){
    for(const row of groups)for(let i=0;i<row.length;i++){
      ok(remain.delete(row[i]),"proper colour input partition");
      for(let j=0;j<i;j++)ok(U(phases.get(row[i]),phases.get(row[j]))<d,"colour class incompatibility");
    }
    ok(remain.size===0,"colouring covers every candidate");
  }
  let cap=5+(groups?groups.length:size);
  verify(d,cap);report.push([d,size,cap]);
}
ok(JSON.stringify(report)==="[[39,0,5],[38,0,5],[37,42,18],[36,95,28],[35,213,53],[34,417,422],[33,594,599]]","frozen cap table");
verify(32,2000);
console.log("PASS: exact 68D Gamma72 graph on all 2000 dual phases; [5 mu ceiling, addable candidates, size bound]");
console.log(report);
