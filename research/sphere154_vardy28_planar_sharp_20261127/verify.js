'use strict';
const assert=require('node:assert/strict');
const n=14,m=(1<<n)-1;
const pc=x=>{let s=0;while(x){s++;x&=x-1;}return s;};
const col=Array.from({length:n},(_,j)=>((471<<j)|(471>>(n-j)))&m);
const tr=Array.from({length:n},(_,j)=>col.reduce((s,c,k)=>s|(((c>>j)&1)<<k),0));
const Id=Array.from({length:n},(_,i)=>1<<i),synC=[...col,...Id],synD=[...Id,...tr];
const lin=(L,x)=>L.reduce((s,a,i)=>s^((x>>i)&1?a:0),0);
function dd(g){
 const out=new Uint8Array(1<<n).fill(255);out[0]=0;
 for(let c of g){
  assert(pc(c)%2===1); // all syndrome columns change parity
  const t=out.slice();for(let i=0;i<out.length;i++)out[i]=Math.min(out[i],t[i^c]+1);
 }
 return out;
}
const dC=dd(synC),dD=dd(synD),rad=[];
let minC=28,minD=28;
for(let x=1;x<1<<n;x++){
 const a=lin(col,x),b=lin(tr,x),u=pc(a)+pc(x),v=pc(b)+pc(x);
 minC=Math.min(minC,u);minD=Math.min(minD,v);
 assert(u%2===0&&v%2===0);
 if(lin(tr,a)===x)rad.push([x,a,u]);
}
assert.deepEqual([minC,minD],[8,8]);
assert.deepEqual(rad,[[5461,5461,14],[10922,10922,14],[16383,16383,28]]);
const hist=d=>Array.from({length:7},(_,i)=>d.filter(x=>x===i).length);
for(const d of [dC,dD]){
 assert.deepEqual(hist(d),[1,28,378,3276,7735,4888,78]);
 const deep=[...d.keys()].filter(i=>d[i]===6);
 assert(deep.every(x=>pc(x)%2===0));
 for(const a of deep)for(const b of deep)assert(d[a^b]!==6);
}
const deep=[...dC.keys()].filter(x=>dC[x]===6);
let tested=0;
for(let i=0;i<28;i++)for(let j=0;j<i;j++)for(let k=0;k<j;k++){
 const a=synC[i],b=synC[j],c=synC[k],U=[0,a,b,c,a^b,a^c,b^c,a^b^c];
 assert(new Set(U).size===8);
 // Any candidate affine cube with all entries at least 5 contains a depth-six even syndrome.
 for(const s of deep){tested++;assert(U.some(g=>dC[s^g]<=4));}
}
assert(tested===3276*78);

// Literal square attainment. Syndrome columns store x followed by y.
const A=[715,15538],B=[12300,6534];
const cs=([x,y])=>[dC[lin(col,x)^y],dD[lin(tr,y)^x]];
assert.deepEqual(cs(A),[6,6]);assert.deepEqual(cs(B),[6,6]);
assert.deepEqual(cs([A[0]^B[0],A[1]^B[1]]),[4,4]);
assert((pc(A[0])+pc(A[1]))%2===0&&(pc(B[0])+pc(B[1]))%2===0);
// Opposite-parity cap angle coefficient: 2*pi+4 < 6*sqrt3 using pi<22/7.
assert(72**2<3*(6*7)**2);
// Low-radius area bounds with pi replaced by 3; 8(mu/32)^15 times normal density.
// Below 9: at most 2*2^14 normal sites in a disk of diameter sqrt(mu).
assert(1048576n*9n**15n<3n*9n*32n**15n);
// 9..18: at most 2 floor(2^14/29)=1128 sites in a disk of diameter sqrt(mu-8).
assert(2*Math.floor(2**14/29)===1128);
for(const mu of [9n,18n])assert(36096n*mu**15n<3n*(mu-8n)*32n**15n);
// 18..26: at most 4 sites in a disk of diameter sqrt(mu-16).
for(const mu of [18n,26n])assert(128n*mu**15n<3n*(mu-16n)*32n**15n);
console.log('PASS exact Vardy syndrome radii, all 3276 three-punctures, four square phases and density bounds');
