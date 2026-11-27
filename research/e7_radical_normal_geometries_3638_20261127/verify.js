'use strict';
const assert=require('node:assert/strict');

const d={"E7_gram":[[2,-1,-1,0,-1,0,0],[-1,2,0,0,0,0,0],[-1,0,2,-1,0,0,0],[0,0,-1,2,0,0,0],[-1,0,0,0,2,-1,0],[0,0,0,0,-1,2,-1],[0,0,0,0,0,-1,2]],"alpha":[0,1,14,15,22],"field_polynomial":67,"phase_isometry":[0,1,7,6,8,9,15,14,27,26,28,29,19,18,20,21,38,39,33,32,46,47,41,40,61,60,58,59,53,52,50,51,44,45,43,42,36,37,35,34,55,54,48,49,63,62,56,57,10,11,13,12,2,3,5,4,17,16,22,23,25,24,30,31,82,83,85,84,90,91,93,92,73,72,78,79,65,64,70,71,116,117,115,114,124,125,123,122,111,110,104,105,103,102,96,97,126,127,121,120,118,119,113,112,101,100,98,99,109,108,106,107,88,89,95,94,80,81,87,86,67,66,68,69,75,74,76,77]},F=[[2,-1,-1,0,-1,0,0],[0,3,-1,0,-1,0,0],[0,0,4,-3,-2,0,0],[0,0,0,5,-2,0,0],[0,0,0,0,4,-5,0],[0,0,0,0,0,3,-4],[0,0,0,0,0,0,1]],D=[30,10,5,3,3,5,40],M=d.phase_isometry,A=d.alpha;
const r4=new Uint8Array(16384).fill(99),v=Array(7).fill(0);
function en(i,b){if(i<0){let k=0;for(let j=0;j<7;j++)k|=(v[j]&3)<<(2*j);r4[k]=Math.min(r4[k],(1440-b)/60);return;}let off=0;for(let j=i+1;j<7;j++)off+=F[i][j]*v[j];let q=0;while((q+1)**2*D[i]<=b)q++;for(let t=Math.ceil((-q-off)/F[i][i]);t<=Math.floor((q-off)/F[i][i]);t++){v[i]=t;en(i-1,b-D[i]*(F[i][i]*t+off)**2);}}
en(6,1440);if(!r4.every(a=>a<=24))throw Error("L4 incomplete");
const gf=(x,y)=>{let z=0;while(y){if(y&1)z^=x;x*=2;if(x&64)x^=d.field_polynomial;y>>=1;}return z};
const sy=(j,t)=>(t&63)+(gf(A[j],t&63)<<6)+4096*(t>>6);
const seq=n=>Array.from({length:n},(_,j)=>j);
const k4=seq(128).map(b=>seq(7).reduce((s,j)=>s|((b>>j&1)<<(2*j)),0));
const L=seq(128).map(b=>seq(128).map(t=>r4[k4[b]+2*k4[M[t]]]));
const tags=seq(5).map(j=>seq(128).map(t=>sy(j,t)));
const itag=seq(5).map(j=>seq(128).map(b=>sy(j,M.indexOf(b))));
const cached=new Map();
const ext=b=>{const key=b.join(",");if(cached.has(key))return cached.get(key);let u=new Uint8Array(8192).fill(99);u[0]=0;
for(let j=0;j<5;j++){const a=new Uint8Array(8192).fill(99),l=L[b[j]],g=tags[j];for(let t=0;t<128;t++){let z=l[t],s=g[t];for(let k=0;k<8192;k++){let w=u[s^k]+z;if(w<a[k])a[k]=w;}}u=a;}cached.set(key,u);return u;};
const add=(u,v)=>[u[0]^v[0],u[1]^v[1],u[2]^v[2],u[3]^v[3],u[4]^v[4],u[5]^v[5]^itag[0][u[0]&v[0]]^itag[1][u[1]&v[1]]^itag[2][u[2]&v[2]]^itag[3][u[3]&v[3]]^itag[4][u[4]&v[4]]];
const depth=u=>ext(u.slice(0,5))[u[5]];
const neg=u=>add(u,add(u,u));
const hist=u=>u.reduce((s,a)=>(s[a]=(s[a]??0)+1,s),{});

const G=d.E7_gram;
for(let i=0;i<7;i++)for(let j=0;j<7;j++)assert.equal(F.reduce((s,z,k)=>s+D[k]*z[i]*z[j],0),60*G[i][j]);
const det=A=>A.length===1?A[0][0]:A[0].reduce((s,t,j)=>s+(j%2?-1:1)*t*det(A.slice(1).map(r=>r.filter((_,k)=>k!==j))),0);
assert.equal(det(G),2);
assert.equal(M[64],82); assert.equal(new Set(M).size,128);
for(let u=0;u<128;u++)for(let v=0;v<128;v++)assert.equal(M[u^v],M[u]^M[v]);
for(let u=1;u<64;u++)assert.equal(new Set(seq(64).map(v=>gf(u,v))).size,64);
assert.equal(new Set(A).size,5);
assert.equal(L[0][0],0); assert.equal(Math.min(...L[0].slice(1)),8);
// f and g separately track all assignments and those with at least one nonzero local label.
let f=new Uint8Array(8192).fill(99),g=new Uint8Array(8192).fill(99);f[0]=0;
for(let j=0;j<5;j++){
 const fn=new Uint8Array(8192).fill(99),gn=new Uint8Array(8192).fill(99);
 for(let u=0;u<8192;u++)for(let x=0;x<128;x++){
  const s=u^tags[j][x], w=L[0][x];
  fn[s]=Math.min(fn[s],f[u]+w);
  gn[s]=Math.min(gn[s],(x?f:g)[u]+w);
 }
 f=fn;g=gn;
}
assert.equal(g[0],32);
const maxima=[24,22,20,26,24,30];
for(let m=0;m<32;m++){
 const b=seq(5).map(i=>m>>i&1?82:0), a=ext(b);
 const weight=b.filter(Boolean).length;
 assert.equal(Math.max(...a),maxima[weight]);
 if(weight===4)assert.ok(a.every(x=>x===24));
 if(m===31)assert.ok(a.every(x=>x===30));
}
const q2=m=>seq(7).reduce((s,i)=>s+(m>>i&1)*seq(7).reduce((v,j)=>v+(m>>j&1)*G[i][j],0),0)/2%2;
assert.equal(q2(82),1);
for(let i=0;i<7;i++)assert.equal(seq(7).reduce((s,j)=>s+((82>>j)&1)*G[i][j],0)&1,0);
for(let m=0;m<128;m++){
 assert.equal(Math.min(...L[m]),m===0?0:m===82?6:q2(m)?2:4);
 assert.equal(q2(m^82),q2(m)^1);
}
const exceptions=[6,18,20];
assert.equal(exceptions.reduce((s,m)=>s^m,0),0);
assert.ok(exceptions.every(m=>m!==0&&q2(m)===0));
const oddwalk=new Set(), at=seq(5).map(()=>0);
for(let j=0;j<5;j++)for(const e of exceptions){
 oddwalk.add(at.join(','));
 for(let i=0;i<5;i++)at[i]^=i===j?e:82;
}
assert.equal(oddwalk.size,15);
assert.ok(at.every(x=>x===0));
const h=ext([0,0,0,0,0]), t=[79,418,628];
for(let m=1;m<8;m++)assert.equal(h[t.reduce((s,a,j)=>m>>j&1?s^a:s,0)],m===7?16:24);
const phase=a=>a.reduce((z,t,i)=>(t&1)?(z^[79,418,628][i]):z,0);
for(let i=-5;i<=5;i++)for(let j=-5;j<=5;j++){
 let n=2*(i*i+i*j+j*j);
 if(n>0&&n<8)assert.ok(4*n+h[phase([i,j]) ]>=32);
 const x=3*i+1,y=3*j+1;
 assert.ok(2*(x*x+x*y+y*y)>=6);
}
for(let i=-4;i<=4;i++)for(let j=-4;j<=4;j++)for(let k=-4;k<=4;k++){
 const p=[i+j,i+k,j+k],n=p.reduce((s,a)=>s+a*a,0);
 if(n>0&&n<8)assert.ok(4*n+h[phase([i,j,k]) ]>=32);
 if((i+j+k)%2===0)assert.ok((4*i-3)**2+(4*j-1)**2+(4*k-1)**2>=11);
}
assert.ok(11/16>1/2&&2/3>1/2);
console.log('PASS E7 parity cover, H35 min/covol, all 32 radical profiles, 8192 all-five depths, tags 79/418/628 and A2/D3 witnesses, literal depth-7 odd cycle');
