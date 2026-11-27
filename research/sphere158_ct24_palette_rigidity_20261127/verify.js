"use strict";
// Sphere158 exact finite check. Run: node verify.js
const assert = require("node:assert/strict");
const r = (() => {
'use strict';
const N=729,pow=[1,3,9,27,81,243],V=Array.from({length:N},(_,v)=>pow.map(p=>Math.floor(v/p)%3));
const Q=V.map(v=>v.reduce((s,x)=>s+x*x,0)%3),tag=Q.map(q=>q===2?-1:q);
const neg=V.map(v=>v.reduce((s,x,i)=>s+(x?3-x:0)*pow[i],0));
const add=new Uint16Array(N*N);
for(let a=0;a<N;a++)for(let b=0;b<N;b++){let s=0;for(let j=0;j<6;j++)s+=((V[a][j]+V[b][j])%3)*pow[j];add[a*N+b]=s;}
const E=Array.from({length:N},(_,t)=>{
 if(!t)return [];
 const z=[];for(let x=0;x<N;x++)if(!x||x===t||tag[add[t*N+neg[x]]]>tag[x])z.push(x);
 return z;
});
const H={},opts=[];
for(let t=1;t<N;t++)H[E[t].length]=(H[E[t].length]||0)+1;
for(const s of [1,4,13]){
 const em=new Uint8Array(N);for(const i of E[s])em[i]=1;
 let max=0,arg=null,cnt=0;
 for(let t=1;t<N;t++){
  const w=add[t*N+neg[s]],L=E[t];
  for(let a=1;a<N;a++){
   const B=a*N,v=add[B+w];if(v===0||tag[v]>tag[a])continue;cnt++;
   let c=0;for(let l=0;l<L.length;l++)c+=em[add[B+L[l]]];
   if(c>max){max=c;arg=[a,t];}
  }
 }
 opts.push([s,max,arg,cnt]);
}
return {H,opts};

})();
assert.deepEqual(r, {
 H: {244: 504, 245: 224},
 opts: [[1, 135, [39, 1], 353099], [4, 135, [14, 4], 353099], [13, 135, [5, 13], 353100]]
});
for (let n=440; n<=729; n++) for (let a=n-Math.floor(14*(729-n)/27); a<=376; a++) {
 const k=Math.max(2,n-a);
 assert((135*a-38879)*k>110*(729-a));
}
const cap = j => {
 let ans=j*439;
 for (let ell=2;ell<=j;ell++) ans=Math.max(ans, Math.floor(911*ell/2)+ell+(j-ell)*439);
 for (let n=440;n<=729;n++) ans=Math.max(ans, n+(j-1)*(n<=486?439:Math.min(439,2*(729-n))));
 return ans;
};
assert.equal(cap(3),1387); assert.equal(cap(4),1826);
assert(2n*1387n**2n<3n*1458n**2n);
assert(1387<1458);
assert(1387*1000<952*1458);
assert.equal(Math.floor(27*1826/4),12325);
assert(12325n**2n<2n*8748n**2n);
console.log(JSON.stringify(r));
console.log("PASS: outsider overlaps, rigidity arithmetic, colour caps and density comparisons.");
