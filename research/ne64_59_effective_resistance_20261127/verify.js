'use strict';
const assert = require('node:assert/strict');

function parameters(A) {
  const n=5,H=A.slice(1).map((row,i)=>row.slice(1).map((v,j)=>i===j?12:-v));
  const I=H.map((row,i)=>row.map((_,j)=>i===j?1:0));
  const mul=(B,C)=>B.map(row=>C[0].map((_,j)=>row.reduce((s,t,k)=>s+t*C[k][j],0)));
  let U=I,adj,det;
  // Faddeev--LeVerrier, with all small matrix arithmetic integral.
  // Spectrum is in [0,24], so these degree-at-most-five polynomial
  // intermediates are well within Number's exact integer range.
  for(let k=1;k<=n;k++){
    const V=mul(H,U),s=-V.reduce((acc,x,i)=>acc+x[i],0);
    assert(s%k===0);
    if(k===n){det=-s/k;adj=U;break;}
    U=V.map((x,i)=>x.map((v,j)=>v+(i===j?s/k:0)));
  }
  assert(det>0);
  assert.deepEqual(mul(H,adj),I.map(row=>row.map(x=>x*det)));
  const r=[];
  for(let i=0;i<6;i++)for(let j=0;j<i;j++){
    const a=3*((i?adj[i-1][i-1]:0)+(j?adj[j-1][j-1]:0)-(i&&j?2*adj[i-1][j-1]:0));
    r.push({i,j,a});
  }
  const top=Math.max(...r.map(x=>x.a));
  assert(top>0&&top<12*det);
  return {det,r,top};
}

const cmp=(d,a)=>{const p=98304n*(16n*BigInt(12*d-a))**59n,
                         q=BigInt(d)*BigInt(185*d)**59n;
                  return p<q?-1:p===q?0:1;};
const matching=A=>A.every((r,i)=>r.filter(x=>x===4).length===1&&
    r.every((x,j)=>i===j?x===0:x===2||x===4));
const wg=f=>Array.from({length:6},(_,i)=>Array.from({length:6},(_,j)=>i===j?0:f(i,j)));
const E=parameters(wg((i,j)=>2+2*Number((i>>1)===(j>>1))));
assert.equal(E.det,98304);assert.equal(E.top*16,7*E.det);
assert.equal(cmp(E.det,E.top),0);

const C=parameters(wg((i,j)=>2+Number((i+1)%6===j||(j+1)%6===i)));
assert.equal(C.det,101400);
for(const {i,j,a} of C.r){
  const distance=Math.min(i-j,6-i+j);
  const [u,v]=[[0,1],[209,520],[28,65],[45,104]][distance];
  assert.equal(a*v,u*C.det);
}
assert.equal(cmp(C.det,C.top),-1);
const lowerTop=Math.max(...C.r.filter(({i,j})=>i-j!==3).map(x=>x.a));
assert.equal(cmp(C.det,lowerTop),1);

const A=wg(()=>0),left=Array(6).fill(12);
let nodes=0,valid=0,eq=0;
function extend(i,j){
  nodes++;
  if(i===5){
    if(left[i]!==0)return;
    for(let a=0;a<6;a++)for(let b=a+1;b<6;b++)for(let c=b+1;c<6;c++)
      if(A[a][b]+A[a][c]+A[b][c]>12)return;
    for(const row of A)assert.equal(row.reduce((s,x)=>s+x,0),12);
    valid++;
    const {det,top}=parameters(A),side=cmp(det,top);
    assert(side<=0);
    if(side===0){eq++;assert(matching(A));}
    return;
  }
  if(j===6){if(left[i]===0)extend(i+1,i+2);return;}
  const upper=Math.min(6,left[i],left[j]);
  let lo=Math.max(0,left[i]-6*(5-j),left[j]-6*(4-i));
  if(i===0&&j>1)lo=Math.max(lo,A[0][j-1]);
  for(let w=lo;w<=upper;w++){
    if(i===0&&w*(5-j)>left[i]-w)continue;
    A[i][j]=A[j][i]=w;
    left[i]-=w;left[j]-=w;
    extend(i,j+1);
    left[i]+=w;left[j]+=w;
  }
}
extend(0,1);
assert.equal(nodes,1258711);
assert.equal(valid,55176);
assert.equal(eq,3);
console.log({nodes,valid,eq,classification:'only the perfect matching attains equality'});
