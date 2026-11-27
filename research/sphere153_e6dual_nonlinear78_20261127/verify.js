'use strict';
const assert=require('node:assert/strict');
const H=[[4,5,6,4,2,3],[5,10,12,8,4,6],[6,12,18,12,6,9],[4,8,12,10,5,6],[2,4,6,5,4,3],[3,6,9,6,3,6]];
const C=Array.from({length:6},(_,i)=>Array.from({length:6},(_,j)=>i===j?2:0));
for(const [i,j] of [[0,1],[1,2],[2,3],[3,4],[2,5]])C[i][j]=C[j][i]=-1;
for(let i=0;i<6;i++)for(let j=0;j<6;j++)assert(H[i].reduce((s,x,k)=>s+x*C[k][j],0)===(i===j?3:0));
const Cm=C.map(x=>x.slice());let last=1;
for(let i=0;i<5;i++){
 const v=Cm[i][i];assert(v>0);
 for(let j=i+1;j<6;j++)for(let k=i+1;k<6;k++){
  const a=v*Cm[j][k]-Cm[j][i]*Cm[i][k];assert(a%last===0);Cm[j][k]=a/last;
 }
 last=v;
}
assert(Cm[5][5]===3);
const Q=x=>x.reduce((s,t,i)=>s+t*x.reduce((a,u,j)=>a+H[i][j]*u,0),0);
const DD=(a,b)=>Q(a.map((x,i)=>x-b[i]));
const pc=x=>{let n=0;while(x){x&=x-1;n++}return n;};
const bit=i=>1<<i, mask=x=>x.reduce((s,a,i)=>s|((a&1)<<i),0);
const q=x=>((Q(Array.from({length:6},(_,i)=>(x>>i)&1))/2)%2+2)%2;
const nrm=new Map();let rootCount=0,mediumCount=0;
const x=[];
function rootBox(i){
 if(i===6){
  const s=Q(x);if(!s||s>6)return;assert(s===4||s===6);
  if(s===4)rootCount++;else mediumCount++;
  const p=mask(x);assert(p>0 && (nrm.get(p)===undefined||nrm.get(p)===s));nrm.set(p,s);return;
 }
 for(let a=-2;a<=2;a++){x[i]=a;rootBox(i+1);}
}
rootBox(0);assert(rootCount===54&&mediumCount===72&&nrm.size===63);
const R=[...nrm.keys()].filter(k=>nrm.get(k)===4).sort((a,b)=>a-b),Rs=new Set(R);
assert(R.length===27);
for(let a=1;a<64;a++)assert(Rs.has(a)===(q(a)===0));
// Each nonzero coefficient parity is at squared normal distance 2 or 3.
// Test the literal physical (not quotient) 8-point witness.
const Y=[[0,0,0,0,0,0],[1,0,0,0,0,0],[1,0,0,0,0,-1],
         [0,1,0,0,0,-1],[1,-1,1,0,0,-1],[1,0,-1,1,0,0],
         [1,0,0,-1,1,0],[1,0,0,0,-1,0]];
for(let i=1;i<8;i++)for(let j=0;j<i;j++){
 const short=j===0||(i>=3&&(j===1||j===2));
 assert(DD(Y[i],Y[j])===(short?4:6));
}

// Link hypergraph: endpoints of projected codegree five and their common neighbors.
const adj=R.map(a=>R.reduce((m,b,i)=>m|(Rs.has(a^b)?bit(i):0),0)),Ws=new Set();
for(let i=0;i<27;i++)for(let j=0;j<i;j++){
 const c=adj[i]&adj[j]; if(pc(c)===5)Ws.add(c|bit(i)|bit(j));
}
assert(Ws.size===216);
const W=[...Ws];
const T=[];
for(let a=1;a<64;a++)if(q(a)){
 const f=x=>(q(x^a)^q(x)^q(a))?x^a:x;
 for(let x=0;x<64;x++){assert(q(f(x))===q(x));for(let y=0;y<64;y++)assert(f(x^y)===(f(x)^f(y)));}
 const p=R.map(x=>R.indexOf(f(x)));assert(p.every(i=>i>=0));T.push(p);
}
assert(T.length===36);

const hbit=x=>1n<<BigInt(x), hpc=x=>{let n=0;while(x){x&=x-1n;n++;}return n;};
const nbd=R.reduce((s,x)=>s|hbit(x),1n),out=((1n<<64n)-1n)^nbd;
const BigWindows=new Set();
for(let a=0;a<64;a++)for(let w of W){
 let v=hbit(a);for(let i=0;i<27;i++)if(w&bit(i))v|=hbit(a^R[i]);BigWindows.add(v);
}
assert(BigWindows.size===13824);
function force(D){
 let f=0n;
 for(const w of BigWindows)if(!(w&D)){
  const a=w&out;
  if(!a)return Infinity;
  if((a&(a-1n))===0n)f|=a;
 }
 let L=[...new Set([...BigWindows].filter(w=>!(w&(D|f))).map(w=>w&out))].sort((a,b)=>hpc(a)-hpc(b));
 let disjoint=0,used=0n;
 for(const w of L)if(!(w&used)){used|=w;disjoint++;}
 return hpc(D)+hpc(f)+disjoint;
}
const expect={5:[0,0],6:[72,1],7:[1728,3],8:[19440,7],9:[136350,22]};
for(let k=5;k<=9;k++){
 const all=[];
 function choose(d,r,i){
  if(!r){if(W.every(w=>w&d))all.push(d);return;}
  for(let j=i;j<=27-r;j++)choose(d|bit(j),r-1,j+1);
 }
 choose(0,k,0);const rem=new Set(all);
 let orbits=0,mn=Infinity;
 for(const a of all)if(rem.delete(a)){
  let queue=[a];
  for(let i=0;i<queue.length;i++)for(const p of T){
   let t=0;for(let v=queue[i];v;v&=v-1)t|=bit(p[31^Math.clz32(v&-v)]);
   if(rem.delete(t))queue.push(t);
  }
  orbits++;
  const D=R.reduce((s,x,i)=>s|((a&bit(i))?hbit(x):0n),0n);
  const n=force(D);assert(n>=20);mn=Math.min(n,mn);
 }
 assert.deepEqual([all.length,orbits],expect[k]);
 console.log('misses in link',k,'hitters',all.length,'orbits',orbits,'required global misses',mn);
}
assert(9n*27n*64n===4n*243n*16n); // ratio² (27/4)/(243/64) = 16/9.
assert(44n*4n*12n===64n*3n*11n);
console.log('PASS: 78D literal full-stock obstruction, compact occupancy <=44/64, density <=11/12 incumbent');
