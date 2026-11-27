'use strict';
const assert = require('node:assert/strict');
const dual = require('./dual.json');
const Z = [0,3,4,7,8,12], K = [1,56,63,72,63,1];
const P = [
 [1,56,63,72,63,1], [1,-56,63,-72,63,1], [1,0,-7,0,7,-1],
 [1,-8,-1,8,-1,1], [1,8,-1,-8,-1,1], [1,0,9,0,-9,-1]
];
const Y0 = [3,3,-1,-1,-1,-1,-1,-1], C = [1,1,1,1,-1,-1,-1,-1];
const residue = (r) => {
  const s = r[0];
  const a = r.map((x,i) => (x-s*C[i])/2);
  assert(a.every(Number.isInteger));
  let q = s & 1, p = 0;
  for (let i=1;i<=6;i++) {p += a[i]; q |= (p & 1) << i;}
  return q;
};
assert.equal(residue(Y0),81); // 2w in E7, supplied as doubled vector.
const minimum=Array(256).fill(Infinity), cmin=Array(256).fill(null);
for(let e=0;e<=1;e++) for(let b=0;b<2;b++){
  const mod=(3*e+2*b)%4;
  const L=Array.from({length:19},(_,i)=>i-9).filter(i=>((i-mod)%4)===0);
  function go(y,sum,n){
    if(y.length===7){
      let z=-sum;
      if(!L.includes(z)||n+z*z>96)return;
      const y8=y.concat(z),r=y8.map((x,i)=>(x-e*Y0[i])/2);
      assert(r.every(Number.isInteger));
      let code=128*e+residue(r);
      if(n+z*z<minimum[code]){minimum[code]=n+z*z;cmin[code]=y8;}
      return;
    }
    for(const v of L)if(n+v*v<=96)go(y.concat(v),sum+v,n+v*v);
  }
  go([],0,0);
}
assert(minimum.every(n=>Z.includes(n/8)));
const d=minimum.map(n=>Z.indexOf(n/8));
assert.deepEqual(K,K.map((_,j)=>d.filter(x=>x===j).length));
assert.equal(minimum[0],0);
const parity=(n)=>{let p=0;while(n){p^=1;n&=n-1;}return p;};
const groups=P.map(()=>[]);
for(let k=0;k<128;k++)for(let neg=0;neg<2;neg++){
  const row=Z.map(()=>0), im=Z.map(()=>0), t=parity(k&81);
  const chi=(x)=> [((parity(k&(x&127))^(neg*(x>>7)))? -1:1),t*(x>>7)];
  for(let x=0;x<256;x++){
    const [v,i]=chi(x);
    (i?im:row)[d[x]]+=v;
  }
  assert(im.every(x=>x===0));
  const h=P.findIndex(r=>r.every((x,i)=>x===row[i]));
  assert(h>=0);
  groups[h].push(chi);
}
assert.deepEqual(groups.map(g=>g.length),[1,1,72,63,63,56]);
for(let h=0;h<6;h++) for(let x=0;x<256;x++){
  let a=0,b=0;
  for(const chi of groups[h]){const [v,i]=chi(x); if(i)b+=v; else a+=v;}
  assert.equal(b,0);
  assert.equal(K[d[x]]*a,groups[h].length*P[h][d[x]]);
}
const O=[[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]];
const permAverageTimesR=(h,s)=>O.reduce((q,o)=>q+o.reduce(
 (p,j,i)=>p*BigInt((504/K[s[j]])*P[h[i]][s[j]]),1n),0n);
const R=6n*504n**3n;
const cutoff=(M,t)=>{
  const lhs=M*M*BigInt(t)**21n*2500n*2n**82n*3n**23n;
  const rhs=39n**2n*43n**21n*4096n**2n*16n**21n;
  assert(lhs<rhs);
  return Number(M)/4096*(t/16)**10.5/(43**10.5/(2**41*3**11.5));
};
const gate={8:9,10:11,12:12,13:13,14:14,15:15,16:16};
const expected={8:274627,10:29906,12:9696,13:4096,14:2645,15:1192,16:663};
assert.deepEqual(Object.keys(dual),Object.keys(gate));
const trivial=cutoff(256n**3n,6);
const clique=cutoff(256n**3n/8n,7);
console.log('t≤6 elementary ratio',trivial,'; t=7 clique ratio',clique);
let max=Math.max(trivial,clique);
for(const [t,ys] of Object.entries(dual)){
  assert(ys.length>0&&ys.every(([h,y])=>/^[0-5]{3}$/.test(h)&&h!=='000'&&Number.isSafeInteger(y)&&y>0));
  let hmin=null;
  for(let i=0;i<6;i++) for(let j=i;j<6;j++) for(let k=j;k<6;k++){
    if(Z[i]+Z[j]+Z[k]<Number(t))continue;
    const v=-ys.reduce((q,[h,y])=>q+BigInt(y)*permAverageTimesR([...h].map(Number),[i,j,k]),0n);
    if(hmin===null||v<hmin)hmin=v;
  }
  assert(hmin>0n);
  const M=1n+R*ys.reduce((q,[_,y])=>q+BigInt(y),0n)/hmin;
  assert.equal(M,BigInt(expected[t]));
  const relative=cutoff(M,gate[t]);max=Math.max(max,relative);
  console.log('t gate',t,'; min dual integer',hmin.toString(),'; M≤',String(M),
      '; comparison up through t=',gate[t],'; ratio',relative);
}
console.log('PASS: E7 finite geometry, 6 Fourier sums and inverse sums, 7 exact duals, ratios < 39/50; numerical maximum display=',max);
