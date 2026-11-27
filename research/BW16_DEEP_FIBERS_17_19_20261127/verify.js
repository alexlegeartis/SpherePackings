'use strict';
const assert=require('node:assert/strict');

const pc=x=>{let c=0;while(x){c++;x&=x-1;}return c;};
const C=Array.from({length:32},(_,a)=>{let w=0;for(let i=0;i<16;i++)if((pc(i&(a&15))+(a>>4))&1)w|=1<<i;return w;});
const U=Array.from({length:1<<16},(_,w)=>w).filter(w=>C.every(c=>pc(c&w)%2===0));
const id=new Int16Array(1<<16).fill(-1),W=[];
for(let w of U)if(id[w]<0){for(let c of C)id[w^c]=W.length;W.push(w);}
const N=256,decode=i=>[i>>7,(i>>6)&1,W[i&63]];
function add(i,j){
 const [p,b,w]=decode(i),[q,d,v]=decode(j);
 return (p^q)*128+(b^d^(pc(w&v)&1))*64+id[w^v];
}
function depth(i){
 if(!i)return 0;const [p,b,w]=decode(i);
 if(p)return 2+b;if(w===0)return 2;
 return Math.min(...C.map(c=>pc(w^c)))/2;
}
const D=Array.from({length:N},(_,i)=>depth(i));
function find9(){
 const A=Array.from({length:N},(_,i)=>i).filter(i=>D[i]===3);
 function rec(L,T){
   if(T.length===9)return T;
   if(L.length+T.length<9)return null;
   for(let i=0;i<L.length;i++) {const a=L[i],r=rec(L.slice(i+1).filter(x=>D[add(x,a)]===3),[...T,a]);if(r)return r;}return null;
 }
 return rec(A,[0]);
}

assert.equal(U.length,2048);
assert.equal(W.length,64);
assert.deepEqual([0,2,3].map(d=>D.filter(x=>x===d).length),[1,135,120]);
for(let i=0;i<N;i++){assert.equal(add(i,i),0);assert.equal(add(i,0),i);}
const deep=[0,12,22,26,192,206,210,221,252];
for(let i=0;i<deep.length;i++)for(let j=0;j<i;j++)assert.equal(D[add(deep[i],deep[j])],3);
const vec=(p,w,k,b)=>Array.from({length:16},(_,i)=>{
 if(p){const u=(w>>i)&1? -1:1;return i===k?-3*u:u;}
 return ((w>>i)&1)?(((b>>i)&1)?-2:2):0;
}); // vector y in units 1/sqrt(8)
for(let a=0;a<N;a++)if(D[a]===3){
  const [p,b,w]=decode(a),S=[];
  for(let c of C){
   const u=c^w;
   if(p){for(let i=0;i<16;i++)S.push(vec(1,u,i,0));continue;}
   if(pc(u)!==6)continue;
   for(let v=u;;v=(v-1)&u){if(pc(v)%2===b)S.push(vec(0,u,0,v));if(!v)break;}
  }
  assert.equal(S.length,512);
  const seen=new Set(S.map(v=>v.join()));
  assert.equal(seen.size,512);
  const sum=Array.from({length:16},()=>new Int32Array(16));
  for(let v of S){
   assert(v.reduce((s,x)=>s+x*x,0)===24);
   assert(seen.has(v.map(x=>-x).join()));
   const u=v.map(x=>(x-p)/2);
   let k=0;for(let i=0;i<16;i++)if(u[i]&1)k|=1<<i;
   assert.equal(id[k],a&63);
   assert.equal((u.reduce((s,x)=>s+x,0)-pc(k))/2&1,b);
   for(let i=0;i<16;i++)for(let j=0;j<=i;j++)sum[i][j]+=v[i]*v[j];
  }
  for(let i=0;i<16;i++)for(let j=0;j<=i;j++)assert.equal(sum[i][j],i===j?768:0);
}
assert(9n*3n**8n<4n**8n);
assert(256n*2n**8n===4n**8n);
const mod=(x,m)=>((x%m)+m)%m;
const R17=x=>[0,192,0,206][mod(x[0],4)];
const R18=x=>deep[mod(x[0],2)+2*mod(x[1],2)];
const R19=x=>{
 const p=mod(x[0],2)+2*mod(x[1],2),c=mod(x.reduce((s,v)=>s+Math.floor(v/2),0),2);
 return deep[2*p+c];
};
for(let i=0;i<4;i++)for(let d=-1;d<=1;d+=2)assert(1+D[add(R17([i]),R17([i+d]))]>=4);
for(let i=0;i<4;i++){let a=[i%2,i>>1];for(let x=-2;x<=2;x++)for(let y=-2;y<=2;y++){
 const q=x*x+x*y+y*y;if(q>0&&q<4)assert(q+D[add(R18(a),R18([a[0]+x,a[1]+y]))]>=4);
}}
for(let i=0;i<64;i++){
 const a=[i%4,(i>>2)%4,i>>4];if(a.reduce((s,x)=>s+x,0)&1)continue;
 for(let x=-2;x<=2;x++)for(let y=-2;y<=2;y++)for(let z=-2;z<=2;z++){
  const q=x*x+y*y+z*z;if(q===0 || q>=8 || (x+y+z)%2)continue;
  assert(q+2*D[add(R19(a),R19([a[0]+x,a[1]+y,a[2]+z]))]>=8);
 }
}
const target=[[0n,1n],[1n,256n],[3n**18n,2n**36n],[13n**19n,3n**18n*5n**21n]];
const pk=[[0n,1n],[1n,1n],[4n,3n],[2n,1n]]; // squares of unit-minimum maximum point densities
const vol=[[0n,1n],[2n,1n],[157n,50n],[314n,75n]]; // lower bounds on unit ball volumes, pi > 3.14
const A=(k,t,u)=>[pk[k][0]*u**(16n+BigInt(k)),256n*pk[k][1]*t**16n*(4n*u-3n*t)**BigInt(k)];
const B=(k,t,u)=>[81n*(4n**BigInt(k))*vol[k][1]**2n*u**(16n+BigInt(k)),256n*vol[k][0]**2n*t**16n*(4n*u-2n*t)**BigInt(k)];
const RC=(k,t,u)=>[256n*pk[k][0]*u**16n,4n**BigInt(k)*pk[k][1]*t**16n];
function bound(k,x,strict){const y=target[k],a=x[0]*y[1],b=x[1]*y[0];assert(strict?a<b:a<=b);}
for(let [k,a,b,c,d] of [[1,13n,10n,15n,8n],[2,32n,25n,7n,4n]]){
 bound(k,A(k,1n,1n),k!==1);bound(k,A(k,a,b),true);
 bound(k,B(k,a,b),true);bound(k,RC(k,c,d),true);
 assert(c*(16n+BigInt(k))<=32n*d); // B bound decreasing on its interval
}
bound(3,A(3,1n,1n),true);bound(3,A(3,73n,60n),true);
bound(3,B(3,51n,40n),true);bound(3,RC(3,31n,20n),true);
assert(31n*19n<32n*20n);
// T(t)=16 delta_19 t^8(4−2t)^(3/2), increasing on candidate interval.
// At its first endpoint T² exceeds all three normal stock cutoffs.
const num=256n*target[3][0]*73n**16n*(240n-146n)**3n, den=target[3][1]*60n**19n;
assert(num>128n*den && num>125n*den && 16n*num>1331n*den);
const cube=[];
for(let p of [0,1])for(let [y,z] of [[0,0],[1,0],[-1,0],[0,1],[0,-1]])cube.push([p,y,z]);
for(let a of cube)for(let b of cube)assert(a.reduce((s,x,i)=>s+(x-b[i])**2,0)<=5);
console.log('PASS: complete BW16 palette, 120 deep shells, nine-clique, three explicit stacks, exact aperture and density comparisons');

