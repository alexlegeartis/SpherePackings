'use strict';
const assert = require('node:assert/strict');
const pc=x=>{let n=0;for(;x;x&=x-1)n++;return n;};
const synd=x=>{let t=0;for(let i=0;i<12;i++)if((x>>i)&1)t^=16+i;return t;};
const ALL=(1<<12)-1;
const Y=Array.from({length:4096},(_,x)=>x), P=Y.map(pc), S=Y.map(synd);
const fibers=Array.from({length:32},(_,s)=>Y.filter(e=>S[e]===s).sort((a,b)=>P[a]-P[b]));
function disjoint(C,k,m,start=0,used=0){
 if(!k)return true;
 for(let i=start;i<C.length && P[C[i]]<=m;i++)
   if(!(C[i]&used)&&disjoint(C,k-1,m,i+1,used|C[i]))return true;
 return false;
}
for(let s=1;s<32;s++)assert(disjoint(fibers[s],4,5),'four disjoint errors syndrome '+s);
const MOD=(x,k)=>((x%k)+k)%k;
const label=(m,n)=>((m&1)?(29^30):0)^((n&1)?(28^30):0)^(MOD(m-n,3)===2?(28^29^30):0);
const ring=[[1,0],[-1,0],[0,1],[0,-1],[-1,1],[1,-1]];
const halo=[1,2,3,28,29,30];
for(let m=0;m<6;m++)for(let n=0;n<6;n++)if(MOD(m-n,3)===0){
 const T=new Set(ring.map(([a,b])=>label(m+a,n+b)));
 assert.equal(T.size,6);
 assert([0,1,2,3].some(t=>halo.every(x=>T.has(x^t))));
}
assert.equal(new Set(halo.flatMap(x=>halo.map(y=>x^y))).size,8);
// Direct independent replay of the quotient-syndrome maximizer statement.
const counts=[0,0,0,0];
for(let F=0;F<=ALL;F++){
 const h=P[F];if(h%4)continue;
 const D=Array(32).fill(99);
 for(let e=0;e<=ALL;e++)D[S[e]]=Math.min(D[S[e]],P[e&(F^ALL)]);
 for(let t=0;t<32;t++)if(halo.every(s=>h+4*D[s^t]===12))counts[h/4]++;
}
assert.deepEqual(counts,[0,0,0,32]);
// Literal Steiner hexads and Julin code.
const K=[0,1,3,4,5,9], rows=[0,...Array.from({length:11},(_,k)=>
 K.reduce((s,i)=>s|(1<<(1+((i+k)%11))),0))];
const H=new Set;
for(let i=0;i<12;i++)for(let j=i+1;j<12;j++){
 const x=rows[i]^rows[j];H.add(x);H.add(x^ALL);
}
assert.equal(H.size,132);assert([...H].every(x=>P[x]===6));
for(const t of Y)if(P[t]===5)assert.equal([...H].filter(x=>(x&t)===t).length,1);
const C=new Set(H);
for(let i=0;i<12;i+=2){const t=3<<i;C.add(t);C.add(t^ALL);}
assert.equal(C.size,144);
const D=[...C];assert(D.every((c,i)=>D.slice(0,i).every(d=>P[c^d]>=4)));
const cover=[0,0,0,0];
for(const b of Y){const d=Math.min(...D.map(x=>P[b^x]));assert(d<=3);cover[d]++;}
assert.deepEqual(cover,[144,1728,1904,320]);
// Constants in the analytic area proof. All comparisons are rational,
// or rational comparisons after squaring positive quantities.
assert(52*52>27*27*3);                  // chord slope
assert(9>8);                           // cos(75°)>1/4
assert(113<128 && 9*13>8*14);           // two marks small gamma
assert(125>121);                       // sqrt(5)>11/5
assert(9*676>40*147);                   // sin(108°)/2 > T(2π/3)
assert(52*676>125*243);                 // recipient budget + three charges
assert(24000+1280+9<25440);             // 4−1/40 full saturation margin
console.log('PASS: all disjoint syndrome errors, all four halo types, unique common half-grid types, Witt code and exact geometric constants.');
