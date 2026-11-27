'use strict';
const ok=require('node:assert/strict');
const w=[1,3,9,27], mod=x=>(x%3+3)%3;
const digits=x=>w.map(k=>Math.floor(x/k)%3);
const key=a=>a.reduce((s,x,i)=>s+w[i]*mod(x),0);
const plus=(a,b,t=1)=>key(digits(a).map((x,i)=>x+t*digits(b)[i]));
const neg=a=>plus(0,a,-1);
function tag(v){const [a,b,c,d]=v;return key([a,a+b,(a+b+c-d)/2,(a+b+c+d)/2]);}
const norm=v=>v.reduce((s,x)=>s+x*x,0);
const ds=Array.from({length:81},()=>[]),roots=[];
for(let a=-2;a<=2;a++)for(let b=-2;b<=2;b++)for(let c=-2;c<=2;c++)for(let d=-2;d<=2;d++){
    const v=[a,b,c,d],n=norm(v);
    if((a+b+c+d)%2!==0||n>6)continue;
    ds[tag(v)].push(v);
    if(n===2)roots.push(v);
}
ok.equal(roots.length,24);
ok.ok(ds.every(a=>a.length)); // all 81 residues realized with norm <=6
const nd=ds.map(r=>Math.min(...r.map(norm)));
const counts={}; for(const n of nd)counts[n]=(counts[n]??0)+1;
ok.deepEqual(counts,{0:1,2:24,4:24,6:32});
const q=(a,b)=>new Set([plus(a,a),a,b,neg(b),plus(a,b),plus(a,b,2),plus(neg(a),b),plus(neg(a),b,2)]);
const planes=new Map();
for(let a=1;a<81;a++)for(let b=1;b<81;b++)if(b!==a&&b!==neg(a)){
    const s=[...q(a,b)].sort((x,y)=>x-y);ok.equal(s.length,8);
    planes.set(s.join(','),s);
}
ok.equal(planes.size,130);
ok.equal([...planes.values()].filter(s=>s.every(a=>nd[a]===6)).length,8);
const t=[...q(5,46)];
ok.ok(t.every(a=>nd[a]===6));
ok.equal(tag([2,-1,-1,0]),5);
ok.equal(tag([1,-1,0,2]),46);
const hs=new Set([5,neg(5)]);
const up=v=>{const [a,b,c,d]=v; return 3*(a*a+b*b+c*c)-(a+b+c)**2;};
// In units 1/192: squared norm under E=(33/32)P_U+(63/64)P_V.
// 192 v^T E v = 189 ||v||² + 3 (3||P_U v||²).
const en=v=>189*norm(v)+3*up(v);
for(const v of roots)ok.ok(en(v)>=384);
for(const c of t)for(const v of ds[c]){
    ok.equal(norm(v),6);
    ok.equal(en(v),hs.has(c)?1188:1152);
}
ok.equal(16*128*33*63,231*9*32*64); // 128*(33/32)*(63/64)/9=231/16
ok.ok(16*128/9<231);
const D=300153807922100058481n,W=1440000000000000000n;
ok.ok(81n*D*10000n>128n*128n*W*10303n); // ratio² exceeds 1.0303
ok.ok(81n*D>128n*128n*W);
ok.ok(1074136270n*3n>32n*100000000n);
ok.ok(106743n*3n>32n*10000n);

// Hermitian backup: p=5, r=191/50, theta=(2+omega)/3.
// Integer q below is 9*N(i+theta) for unit j=1.
const E=(a,b)=>a*a-a*b+b*b;
const sh={};
for(let a=-4;a<=4;a++)for(let b=-4;b<=4;b++){
  const m=E(3*a+2,3*b+1);
  if(191*9+250*m<16*450)sh[m]=(sh[m]??0)+1;
}
ok.deepEqual(sh,{3:3,12:3,21:6});
ok.equal(16*150-191*3-250,1577); // D=1577/150
ok.equal(16*50-3*191,227);
ok.ok(16*50<3*191+250);
ok.ok(16*150<4*191*3+250);
ok.ok(D*100000000n>14325n*14325n*W*100n); // sqrt(D)>14.325
ok.ok(1051475582n*150n>1577n*100000000n);
// sqrt(Dstrong) - 4/sqrt(3) > sqrt(16-r-35/3).
// A simple rational overbound on RHS terms: 4/sqrt(3)<231/100, sqrt(77/150)<72/100.
ok.ok(3*231*231>16*10000);
ok.ok(150*72*72>77*10000);
ok.ok(1577*10000>150*(231+72)**2);
ok.equal(21**3%47,2); // exact constant in cubic-descent R_+; -21 gives -2
console.log('PASS 81 D4 residues, eight 9-site planes, exact 9-colour tie, conditional D96 +1.51%, Hermitian finite branch gates');

