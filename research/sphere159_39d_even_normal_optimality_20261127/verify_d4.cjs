
"use strict";
const assert=require("node:assert/strict");
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
function DAtlas(G,d){
 const B=adj(G),D=det(G),st=G.length-d;
 const M=[12,420,10920,185640][st];
 const res=new Map(), a=Array(d).fill(0), j=[];
 function emit(Q){
  if(Q>6*M)return;
  const m=a.slice(0,-1).map((x,l)=>(x-a[l+1])/2);
  m.push((a[d-2]+a[d-1])/2,...j);
  const x=B.map(row=>row.reduce((s,v,l)=>s+v*m[l],0));
  const f=x.reduce((s,y,l)=>s+y*m[l],0);
  assert.equal(f*M,Q*D);
  const k=x.map(v=>mod(v,D)).join(",");
  if(!res.has(k)||res.get(k)[0]>Q)res.set(k,[Q,[...a],[...j],m]);
 }
 function rec(i,q,dot){
  if(i<d){
   for(let t=-8;t<=8;t++)if(mod(t,2)===p){
    const u=q+t*t;
    if(u>72)continue;
    a[i]=t;rec(i+1,u,dot+t*(i?3:1));
   }
   return;
  }
  if(st===0){emit(q);return;}
  for(let j1=-6;j1<=6;j1++){
   j[0]=j1;
   const w=12*j1-dot,Q1=35*q+w*w;
   if(Q1>2520)continue;
   if(st===1){emit(Q1);continue;}
   for(let j2=-6;j2<=6;j2++){
    j[1]=j2;
    const z=70*j2+35*a[0]+3*w,Q2=26*Q1+z*z;
    if(Q2>65520)continue;
    if(st===2){emit(Q2);continue;}
    for(let j3=-6;j3<=6;j3++){
     j[2]=j3;
     const y=52*(j3-j2)+z;
     emit(17*Q2+35*y*y);
    }
   }
  }
 }
 for(var p=0;p<2;p++)rec(0,0,0);
 const vals=[...res.values()], n=vals.length;
 assert.equal(n,D);
 const levels=[...new Set(vals.map(r=>r[0]))].sort((a,b)=>b-a);
 const top=levels.slice(0,3);
 return {D,n,M,top,pools:top.map(t=>vals.filter(x=>x[0]===t))};
}
const A4=DAtlas(d4,4),A5=DAtlas(d5,5),
      D6=DAtlas(d6,6),A6=DAtlas(g6,5),A7=DAtlas(g7,5),A8=DAtlas(g8,5);
assert.deepEqual([A4.D,A5.D,A6.D,A7.D,A8.D],[324,972,2835,6318,12393]);
assert.equal(A4.top[0],36); assert.equal(A4.top[1],28);
assert.equal(A4.pools[0].length,3);
for(const r of A4.pools[0]){
 assert(r[1].map(Math.abs).sort((a,b)=>b-a).join()==="6,0,0,0" ||
        r[1].every(x=>Math.abs(x)===3));
}
assert.deepEqual(A5.top,[45,37,36]);
assert.deepEqual(D6.top,[54,46,38]);
for(const r of A5.pools[2])assert(r[1].map(Math.abs).sort((a,b)=>b-a).join()==='6,0,0,0,0');
for(const r of D6.pools[0]) assert(r[1].every(x=>Math.abs(x)===3));
for(const r of D6.pools[1]) assert.deepEqual(r[1].map(Math.abs).sort((a,b)=>a-b),[1,3,3,3,3,3]);
for(const r of A5.pools[0]) assert(r[1].every(x=>Math.abs(x)===3));
for(const r of A5.pools[1]) assert.deepEqual(r[1].map(Math.abs).sort((a,b)=>a-b),[1,3,3,3,3]);
assert.deepEqual(A6.top.slice(0,2),[1584,1376]);
assert.deepEqual(A7.top.slice(0,2),[44100,41860]);
assert.equal(A8.top[0],633360);
const sig = rows => rows.map(r=>JSON.stringify(r.slice(1,3))).sort();
assert.deepEqual(sig(A6.pools[0]),[
 [[-6,0,0,0,0],[-2]], [[-6,0,0,0,0],[1]], [[-4,0,0,0,0],[-3]], [[2,0,0,0,0],[-3]]
].map(r=>JSON.stringify(r)).sort());
assert.deepEqual(sig(A7.pools[0]),[
 [[-6,0,0,0,0],[-2,3]], [[-6,0,0,0,0],[1,3]]
].map(r=>JSON.stringify(r)).sort());
assert.equal(A5.pools[0].length,2);
function tell(name,o){
 console.log(name+" det="+o.D+" top="+o.top.join(",")+"/"+o.M+" deepest="+o.pools[0].length);
}
for(const [i,o] of [A4,A5,D6,A6,A7,A8].entries())tell(["D4","D5","D6","T6","T7","T8"][i],o);
const r6=(()=>{
const M=2187;
const pow=[1,3,9,27,81,243];
const mod=(x,p)=>((x%p)+p)%p;
const rep=new Int16Array(M);rep.fill(32767);
const tie=new Int16Array(M);
for(let p=0;p<2;p++) {
 const rec=(i,q,s,j)=>{
  if(i===5){for(let w=p-12;w<=12;w+=2){
   if(mod(s+(w-p)/2,2)!==0)continue;
   const Q=q+w*w;if(Q>144)continue;
   const k=j+243*mod(w,9);
   if(Q<rep[k]){rep[k]=Q;tie[k]=1;}else if(Q===rep[k])tie[k]++;
  }return;}
  for(let x=p-8;x<=8;x+=2){const r=q+3*x*x;if(r>144)continue;
   rec(i+1,r,s+(x-p)/2,j+pow[i]*mod(x,3));}
 };rec(0,0,0,0);
}
const hist={};for(const q of rep)hist[q]=(hist[q]??0)+1;
return{rep:Array.from(rep),tie:Array.from(tie),hist};})();
const r7=(()=>{
const M=4374,N=2187;
const pow=[1,3,9,27,81,243,729];
const mod=(x,p)=>((x%p)+p)%p;
const rep=new Int16Array(M); rep.fill(32767);
const tie=new Int16Array(M);
const a=Array(6).fill(0);
for(let p=0;p<2;p++) {
  const rec=(i,q,s,j)=>{
    if(i===6){for(let w=-5;w<=5;w++) {
      if(mod(s+w,2)!==0)continue;
      const Q=q+2*w*w;if(Q>54)continue;
      const k=j+729*mod(w,3)+N*mod(w-p,2);
      if(Q<rep[k]){rep[k]=Q;tie[k]=1;}else if(Q===rep[k]) tie[k]++;
    }return;}
    for(let x=p-8;x<=8;x+=2){const r=q+x*x;if(r>54)continue;
      a[i]=x;rec(i+1,r,s+(x-p)/2,j+pow[i]*mod(x,3));}
  };
  rec(0,0,0,0);
}
const hist={};for(const q of rep) hist[q]=(hist[q]??0)+1;
return {rep:Array.from(rep),tie:Array.from(tie),hist};})();
assert.deepEqual(r6.hist,{0:1,16:54,24:72,40:432,48:270,64:486,72:240,88:432,96:144,112:54,144:2});
assert.deepEqual(r7.hist,{0:1,6:56,8:126,14:576,16:756,22:756,24:728,30:672,32:576,38:126,54:1});
assert.deepEqual(r6.rep.flatMap((q,i)=>q===144?[i]:[]),[729,1458]);
assert.deepEqual(r7.rep.flatMap((q,i)=>q===54?[i]:[]),[2187]);
console.log("E6",r6.hist);
console.log("E7",r7.hist);
// Each call proves det(seed) times its rank-r Hermite determinant > 21870.
function check(D,m,n,r,k){
  D=BigInt(D);m=BigInt(m);n=BigInt(n);r=BigInt(r);k=BigInt(k);
  assert(D*k*m**r>21870n*(2n*n)**r);
}
check(324,11,3,5,4);
check(972,10,3,4,4);
check(2916,17,6,3,4);
check(2187,26,9,3,4);
check(4374,17,6,2,3);
check(2835,286,105,3,4);
check(6318,13,6,2,3);
assert.equal(633360*17,58*185640);
assert.equal(12393*44,32076*17);
assert(32076>21870);
assert.equal(3**8*10,21870*3);
console.log("PASS: all discriminant maxima, exceptional extensions and Hermite gates.");
