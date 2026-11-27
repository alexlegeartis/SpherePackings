'use strict';
const assert = require('node:assert/strict');
const pc=x=>{let b=0;while(x){x&=x-1;b++;}return b;};
const palettes=[1,2,3].map(k=>Array.from({length:64},(_,i)=>i).filter(i=>pc(i)===k));
const A=Array.from({length:6},()=>Array(6).fill(0)), left=Array(6).fill(12);
const dBase=98304n, betaBaseMul=185n;
function better(d,a,s) {
  return dBase*BigInt(s*s)*(16n*BigInt(12*d-a))**59n >
         36n*BigInt(d)*(betaBaseMul*BigInt(d))**59n;
}
function normal(A) {
  const H=A.slice(1).map((r,i)=>r.slice(1).map((v,j)=>i===j?12:-v));
  let U=H.map((r,i)=>r.map((v,j)=>Number(i===j)));
  for(let k=1;k<=5;k++){
    const V=H.map((r,i)=>r.map((v,j)=>r.reduce((s,a,t)=>s+a*U[t][j],0)));
    const c=-V.reduce((s,r,i)=>s+r[i],0)/k;
    assert(Number.isInteger(c));
    if(k===5){assert(c<0);return [-c,U];}
    U=V.map((r,i)=>r.map((v,j)=>i===j?v+c:v));
  }
  throw Error('characteristic polynomial');
}
let graphs=0, examined=[0,0,0];
function check() {
  graphs++;
  const [d,U]=normal(A), limit=Math.floor((3*d-1)/4);
  const cut=Array(16).fill(-1);
  for(let s=1;s<=15;s++){
    let a=Math.max(-1,Math.min(limit,Math.floor(
      d*(12-(185/16)*(36*d/(98304*s*s))**(1/59)))));
    while(a<limit&&better(d,a+1,s))a++;
    while(a>=0&&!better(d,a,s))a--;
    assert(a<0||better(d,a,s));
    assert(a===limit||!better(d,a+1,s));
    cut[s]=a;
  }
  for(let mode=0;mode<3;mode++){
    const V=palettes[mode], n=V.length, max=[6,15,10][mode];
    const D=V.map((p,i)=>V.map((q,j)=>{
      if(j<=i)return 0;
      const v=Array.from({length:5},(_,t)=>(p>>(t+1)&1)-(q>>(t+1)&1));
      let x=0;
      for(let t=0;t<5;t++)if(v[t]){
        x+=v[t]*v[t]*U[t][t];
        for(let u=0;u<t;u++)x+=2*v[t]*v[u]*U[t][u];
      }
      return 3*x;
    }));
    for(let s=1;s<=max;s++){
      const a=cut[s];
      if(a<0 || (mode>0 && 2*a<=d))continue;
      examined[mode]++;
      const N=V.map((_,i)=>V.reduce((z,v,j)=>j>i&&D[i][j]<=a?z|1<<j:z,0));
      function has(P,k){
        if(k===0)return true;
        while(pc(P)>=k){
          const i=31-Math.clz32(P&-P);P&=P-1;
          if(has(P&N[i],k-1))return true;
        }
        return false;
      }
      assert(!has((1<<n)-1,s),JSON.stringify({A,s,mode,d,a}));
    }
  }
}
function extend(i,j){
  if(i===5){
    if(left[i]!==0)return;
    for(let a=0;a<6;a++)for(let b=a+1;b<6;b++)for(let c=b+1;c<6;c++)
      if(A[a][b]+A[a][c]+A[b][c]>12)return;
    check();return;
  }
  if(j===6){if(left[i]===0)extend(i+1,i+2);return;}
  let lo=Math.max(0,left[i]-6*(5-j),left[j]-6*(4-i));
  const hi=Math.min(6,left[i],left[j]);
  if(i===0&&j>1)lo=Math.max(lo,A[0][j-1]);
  for(let w=lo;w<=hi;w++){
    if(i===0&&w*(5-j)>left[i]-w)continue;
    A[i][j]=A[j][i]=w;left[i]-=w;left[j]-=w;
    extend(i,j+1);left[i]+=w;left[j]+=w;
  }
}
extend(0,1);
assert.equal(graphs,55176);
console.log({graphs,examined,scope:'all clusters with diameter squared below 3/4'});
