const pc=x=>{x=x-((x>>>1)&0x55555555);x=(x&0x33333333)+((x>>>2)&0x33333333);return Math.imul((x+(x>>>4))&0xf0f0f0f,0x01010101)>>>24};

const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
const tr=A=>A[0].map((_,i)=>A.map(r=>r[i]));
const mm=(A,B)=>{const T=tr(B);return A.map(r=>T.map(t=>dot(r,t)));};
const I=Array.from({length:24},(_,i)=>Array.from({length:24},(_,j)=>+(i===j)));
const p24=(v)=>v.reduce((s,x,i)=>s+((x&1)<<i),0);
const mp=(x,M)=>{let y=0;for(let a=x;a;a&=a-1)y^=M[31-Math.clz32(a&-a)];return y;};
function invert(A){
 const V=A.map(p24), W=I.map(p24);
 for(let i=0;i<24;i++){
  const j=V.findIndex((v,k)=>k>=i&&(v>>i&1));if(j<0)throw Error('singular');
  [V[i],V[j]]=[V[j],V[i]];[W[i],W[j]]=[W[j],W[i]];
  for(let k=0;k<24;k++)if(k!==i&&(V[k]>>i&1)){V[k]^=V[i];W[k]^=W[i];}
 }
 return W;
}
function gf(A,B){const R=B.map(p24);return A.map(a=>Array.from({length:24},(_,j)=>mp(p24(a),R)>>j&1));}
function gamma(c,vals){
 const G=Array.from({length:72},(_,i)=>vals.slice(72*i,72*i+72));
 const U=c.leech_basis_change,J=c.leech_frame, F0=G.slice(0,24).map(r=>r.slice(0,24).map(x=>x/2));
 const X=G.slice(0,24).map(r=>r.slice(48,72));
 const if0=invert(F0),iu=invert(U);
 const A=U.map(r=>mp(mp(mp(p24(r),X.map(p24)),if0),iu));
 const B=A.map((x,i)=>x^(1<<i));
 const F=mm(mm(U,F0),tr(U)),FM=F.map(p24), fwd=mm(F,tr(J));
 if(mm(fwd,J).some((row,i)=>row.some((x,j)=>x!==8*+(i===j))))throw Error('frame');
 const add4=(a,b)=>[(a[0]^b[0]),(a[1]^b[1]^(a[0]&b[0]))];
 let R=[[0,0]];
 for(const r of fwd){
  const a=[p24(r),p24(r.map(x=>x>>1))],a2=add4(a,a);
  const T=new Map();for(const t of R)for(const b of [[0,0],a,a2,add4(a2,a)]){const p=add4(t,b);T.set(p.join(','),p);}
  R=[...T.values()];
 }
 if(R.length!==8192)throw Error('mod 4');
 const M=new Map(),Jcols=tr(J);
 let n=0;
 function reg(w){
  const u=Jcols.map(y=>dot(y,w));
  if(u.some(x=>x%8))return;
  n++;const x=p24(u.map(x=>x/8)), a=mp(x,A), b=mp(x,B);
  if(!M.has(a))M.set(a,new Set());M.get(a).add(b);
 }
 for(let i=0;i<24;i++)for(let j=i+1;j<24;j++)for(const a of [-4,4])for(const b of [-4,4]){
  const w=Array(24).fill(0);w[i]=a;w[j]=b;reg(w);
 }
 for(const [lo,hi] of R){
  if(!lo&&pc(hi)===8){
   const X=[];for(let m=hi;m;m&=m-1)X.push(31-Math.clz32(m&-m));
   for(let z=0;z<256;z++){const w=Array(24).fill(0);X.forEach((j,i)=>w[j]=z>>i&1?-2:2);reg(w);}
  }else if(lo===0xffffff){
   const v=Array.from({length:24},(_,j)=>hi>>j&1?-1:1);
   for(let j=0;j<24;j++){const w=v.slice();w[j]*=-3;reg(w);}
  }
 }
 if(n!==196560||M.size!==4095||M.has(0)||[...M.values()].some(s=>s.size!==24))throw Error('M4');
 return {A,B,FM,M4:[...M].sort((a,b)=>a[0]-b[0]).map(([a,V])=>[a,[...V]])};
}

const assert=require('node:assert/strict'),fs=require('node:fs'),crypto=require('node:crypto');
const cbuf=fs.readFileSync(process.argv[2]),gbuf=fs.readFileSync(process.argv[3]);
const blob=x=>crypto.createHash('sha1').update(Buffer.from('blob '+x.length+'\0')).update(x).digest('hex');
assert.equal(blob(cbuf),'5cc5386748aa477917fa7fd113a9c54a343800f5');
assert.equal(crypto.createHash('sha256').update(gbuf).digest('hex'),'1791fab8c152302db2652ec66f227eb473ce4007f934431ef6bbb2ab2f2d87ab');
const d=gamma(JSON.parse(cbuf),JSON.parse(gbuf)),R=new Map(d.M4);
const u=720,v=8359937,plane=[u,v,u^v];
assert(plane.every(a=>R.get(a).length===24 && mp(a,d.A)===a));
for(const a of plane)for(const b of plane)if(a!==b){
 const sa=new Set(R.get(a)), vb=mp(b,d.FM);
 for(const x of R.get(b))assert(!sa.has(x),'4+4 with B-total zero');
 for(const x of R.get(a))assert.equal(pc(vb&x)&1,0,'4+6 with B-total zero');
}
const parity = (p) => p.reduce((s,x,i)=>s+((x&1)<<i),0)+(1<<p.length)*((p.reduce((s,x)=>s+x,0)/2)&1);
const ep=(n)=>Array.from({length:1<<n},(_,i)=>i).filter(i=>!(pc(i)&1));
const unit=(n,i,s=1)=>Array.from({length:n},(_,j)=>i===j?s:0);
function roots(n) {
 const res=new Set;
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)for(const a of [-1,1])for(const b of [-1,1]){
  const x=unit(n,i,a);x[j]=b;res.add(parity(x));
 }
 return res;
}
for(const n of [4,5,6]){
 const r=roots(n);
 for(const c of ep(n))for(const high of [0,1])assert.equal(r.has(c+(high<<n)),pc(c)===2);
}
const qn=(n,c)=>{const w=pc(c&((1<<n)-1));return w || (c?4:0);};
const Q=ep(4).flatMap(x=>[x,x+16]);
assert.deepEqual(Q.filter(x=>qn(4,x)===4).sort((a,b)=>a-b),[15,16,31]);
const F4mul=[0,2,3,1];
const f=c=>{
 const b=c&15, v=(b^((b>>3)*15))&3, u=(b>>3)+2*(c>>4), w=F4mul[v];
 return [u,v,w];
};
function rank(p){const a=new Set(p.filter(Boolean));return a.size===0?0:a.size===1?1:2;}
assert.equal(new Set(Q.map(c=>f(c).join(','))).size,16);
const depth=[0,4,6];
for(const x of Q)for(const y of Q){
 const a=f(x),b=f(y),c=a.map((z,i)=>z^b[i]);
 if(x!==y)assert(qn(4,x^y)+depth[rank(c)]>=8);
}
const d5=roots(5),cross=new Set;
for(let j=0;j<4;j++)for(const s of [-1,1]){
 const x=unit(4,j,s);x[0]--;cross.add(parity(x));
}
assert.equal(cross.size,8);
assert([...cross].every(x=>cross.has(x^16)));
const embed=(c,k)=>{
 const b=Array.from({length:4},(_,i)=>c>>i&1);
 if(parity(b)!==c)b[0]+=2;
 b[0]+=k;b.push(k);
 return parity(b);
};
const P=Q.flatMap(c=>[embed(c,0),embed(c,1)]);
assert.equal(new Set(P).size,32);
for(const c of Q)for(const e of Q)for(const k of [0,1])for(const l of [0,1]){
 assert.equal(d5.has(embed(c,k)^embed(e,l)),k===l?roots(4).has(c^e):cross.has(c^e));
}
console.log('PASS: pinned good 2-plane, all rank-2 base contacts, both doubled root graphs, and exact D4-pair snap graph.');
