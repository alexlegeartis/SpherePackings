'use strict';
const assert=require('node:assert/strict');
const fs=require('node:fs'),path=require('node:path'),zlib=require('node:zlib');
const {atlas39,P,source,imposs,graph,mod,histA,normalSmall5}=require('./verify.js');

function DN(n){
 const L=3**n,C=[];
 for(let p=0;p<2;p++)for(let i=0;i<L;i++){let j=i,a=[];for(let k=0;k<n;k++){a.push(2*(j%3)+p);j=Math.floor(j/3)}for(let b=0;b<2;b++)C.push(a.map((v,k)=>v+(k===0?6*b:0)))}
 return C;
}
function keyD(d){
 let p=mod(d[0],2),j=0,m=1,t=0;
 for(let x of d){let r=mod(x,6);if((r&1)!==p)throw 'mixed parity';j+=m*((r-p)>>1);m*=3;t+=(x-r)/6;}
 return p*2*m+2*j+mod(t,2);
}
function kpD(n,C,H,P){
 const out=[],A=P.map(r=>r.slice(0,n)),T=P.map(r=>r.slice(n)),tail=new Int8Array(8-n);
 const dot=(a,b)=>a.reduce((s,t,i)=>s+t*b[i],0);
 for(const u of C){
  let best=Infinity,ex,YY=A.map(a=>dot(a,u)),p=u[0]&1,ss=u.reduce((s,x)=>s+x,0);
  const rec=(i,k,s)=>{
   if(i<tail.length){for(let a=-8+p;a<=8;a+=2){const z=k+5*a*a;if(z>360||z>=best)continue;tail[i]=a;rec(i+1,z,s+a);}return;}
   if(s%4)return;
   const X=T.map((a,i)=>YY[i]+dot(a,tail));if(X.some(x=>x%4))throw 'map';
   const jj=mod(-X[0]/4,3);let idx=0;
   for(let i=7;i>0;i--){const x=X[i]-X[0];if(x%12)throw 'class';idx=idx*3+mod(x/12,3)}
   idx*=10;
   for(let j=-13+mod(jj+13,3);j<=13;j+=3){let z=k+2*j*j;if(z>360||z>=best)continue;
    const id=idx+mod(j,10),h=2*H[id]+z;if(h<best){best=h;ex=[[...tail],j,id];}
   }
  };
  rec(0,0,ss);out.push({h:best,w:ex});
 }
 return out;
}

const gram=JSON.parse(fs.readFileSync(path.join(source,'gram.json')));
const gens=JSON.parse(fs.readFileSync(path.join(source,'gens.json')));
const c=zlib.gunzipSync(fs.readFileSync(path.join(source,'witness.gz')));
const H=atlas39(new Int8Array(c.buffer,c.byteOffset,c.byteLength),gram,gens);
const C=DN(4), A=kpD(4,C,H,P), witness=[0,319,322,109];
assert(C.length===324 && C.every((u,i)=>keyD(u)===i));
assert.deepEqual(histA(A),{"0":1,"180":3,"220":48,"240":32,"260":72,"280":24,"300":96,"320":24,"340":24});
const d=(u,v)=>A[keyD(u.map((x,i)=>x-v[i]))].h;
assert(C.every((u,i)=>d(C[0],u)===A[i].h));
for(let i=0;i<witness.length;i++)for(let j=0;j<i;j++){
  const a=C[witness[i]],b=C[witness[j]];
  assert.equal(d(a,b),340);
  assert.equal(normalSmall5(a.map((x,k)=>x-b[k])),4);
}
const cap=16n*340n**44n;
assert(360n**44n<cap);  // singleton strictly worse
for(let h of [...new Set(A.map(a=>a.h))].filter(h=>h>0).sort((a,b)=>b-a)){
 const a=BigInt(h)**44n;
 let lo=1,hi=325;
 while(lo<hi){const m=(lo+hi)>>1;if(BigInt(m)**2n*a>cap)hi=m;else lo=m+1;}
 const X=C.filter(u=>d(u,C[0])>=h);
 if(lo>X.length+1)continue;
 const r=imposs(graph(X,(u,v)=>-d(u,v),-h),lo-1);
 assert(!r.found);
 console.log([h,lo,X.length,r.visits]);
}
// squared center-density comparison against Chen-Hu-Li-Wang-Wu
assert(cap*2n**44n*5n**43n*11n**46n < 157n**44n*240n**44n*324n);
console.log('44D D4 exact canonical optimum',witness,'squared min 17/3; equals old T44. PASS');

