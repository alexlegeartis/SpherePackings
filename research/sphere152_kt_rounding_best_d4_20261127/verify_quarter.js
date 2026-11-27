'use strict';
const assert=require('node:assert/strict');
const M=(x,n)=>(x%n+n)%n;
const V=[[0,0],[1,0],[0,1],[1,1]];
function prod(i,j){const [a,b]=V[i],[c,d]=V[j];return M(a*c-b*d,2)+2*M(a*d+b*c-b*d,2);}
const H=[];
for(let a=0;a<4;a++)for(let b=0;b<4;b++)for(let c=0;c<4;c++){
    H.push([a,b,c,a^b^c,a^prod(2,b)^prod(3,c),a^prod(3,b)^prod(2,c)]);
}
const hSet=new Set(H.map(x=>x.join(''))),
    hBase=H.filter(x=>x.slice(0,3).reduce((s,t)=>s+t,0)===1).slice(0,3);
assert.equal(hBase.length,3);
const weight=h=>h.filter(x=>x!==0).length;
const perms=[];
function perm(p,rem){if(!rem.length){perms.push(p);return;}
  for(const i of rem)perm([...p,i],rem.filter(j=>j!==i));
}
perm([],[0,1,2,3,4,5]);
const reps=[H[0],H.find(h=>weight(h)===4),H.find(h=>weight(h)===6)];
// Exhibit monomial maps from the chosen parity word to each other of its weight.
for(const w of H){
  const b=reps[[0,4,6].indexOf(weight(w))];
  assert(b);
  let found=false;
  for(const p of perms){
    if(!w.every((x,j)=>(x===0)===(b[p[j]]===0)))continue;
    const maps=w.map((x,j)=>x?[1,2,3].filter(t=>prod(t,b[p[j]])===x):[1,2,3]);
    const c=[];
    function search(j){
      if(j===6)return hBase.every(g=>hSet.has(p.map((l,i)=>prod(c[i],g[l])).join('')));
      for(const z of maps[j]){c[j]=z;if(search(j+1))return true;}
      return false;
    }
    if(search(0)){found=true;break;}
  }
  assert(found,'monomial orbit '+w);
}
const q=(a,b)=>a*a-a*b+b*b, costs=[];
for(let a=0;a<8;a++)for(let b=0;b<8;b++){
  costs[a+8*b]=V.map(([u,v])=>{
    let n=10000;
    for(let i=-2;i<=2;i++)for(let j=-2;j<=2;j++)n=Math.min(n,q(a-4*u+8*i,b-4*v+8*j));
    return n;
  });
}
const raw=[...Array(64).keys()];
const parity=k=>M(k,2)+2*M(Math.floor(k/8),2);
const opp=k=>M(-M(k,8),8)+8*M(-Math.floor(k/8),8);
let max=0,total=0, witness=null;
const allMax=[];
for(const w of reps){
  const options=w.map((s,j)=>raw.filter(k=>parity(k)===s&&(j<3?(M(k,8)<4&&Math.floor(k/8)<4):(k<=opp(k)))));
  const x=[]; let current=0, count=0; const rows=[];
  function go(i){
    if(i!==6){for(const c of options[i]){rows[i]=costs[c];x[i]=c;go(i+1);}return;}
    count++;let best=1<<30;
    for(const h of H){
      const a=rows[0][h[0]]+rows[1][h[1]]+rows[2][h[2]]+rows[3][h[3]]+rows[4][h[4]]+rows[5][h[5]];
      if(a<best)best=a;
    }
    if(best>current)current=best;
    if(best>max){max=best;witness=x.map(c=>[M(c,8),Math.floor(c/8)]);}
  }
  go(0); allMax.push(current); total+=count;
}
assert.equal(max,34);
const c=17n, d=8n, p=4n;
// b=(p^2/(p^2-r))*c/d; F(b,b*r/(2p^2)) = b^6 * p^4 /( 54 * 4^8) *4/r^2.
// Compare (4/3)F<1/16 even if only continuous projected normal bound is used.
for(const r of [1n,2n]){
 const n=p*p*c,e=(p*p-r)*d;
 assert(16n*16n*p**4n*n**6n<3n*54n*4n**8n*r*r*e**6n);
}
assert(16n*16n*4n**8n*d*d<3n*54n*4n**8n*(4n*d-c)**2n);
for(const m of [3n,4n]){
 const a=5n*m**11n*32n**5n,b=1458n*4n**11n*(8n*m-17n)**5n;
 assert(4n*a*a*12n**22n<3n*b*b*11n**22n);
}
console.log('PASS: K12/4K orbit representatives',total,allMax,'depth 34/16; witness',witness);
