"use strict";
// Exact sphere enumeration, no floating-point bounds; all input is from the catalogue.
const lower=[
[4],[2,4],[0,-2,4],[0,-2,0,4],[0,0,-2,0,4],[-2,-2,0,0,0,4],[0,0,0,0,0,-2,4],[0,0,0,0,0,0,-2,4],[0,0,0,0,1,-1,0,0,4],[0,0,0,0,-1,0,0,0,0,4],[0,0,0,0,0,0,0,0,-2,0,4],[0,1,-2,0,1,0,0,-1,1,0,0,4],[-1,0,-1,-1,1,0,1,0,1,1,0,0,4],[0,0,0,0,0,0,0,0,0,1,-2,0,-1,4],[0,0,0,0,0,0,-1,1,0,1,1,0,1,0,4],[0,0,0,0,0,0,0,0,0,-1,0,1,-1,0,1,4]
].map(r=>r.map(BigInt));
const N=16,numer=[6,4,6,4,4,4,4,4,4,0,4,0,0,0,0,4].map(BigInt);
const A=(x,m)=>{if(!x)throw Error(m);};
const abs=x=>x<0n?-x:x;
function gcd(a,b){a=abs(a);b=abs(b);while(b)[a,b]=[b,a%b];return a;}
function Q(a,b=1n){if(b<0n){a=-a;b=-b;}const d=gcd(a,b);return [a/d,b/d];}
const add=(x,y)=>Q(x[0]*y[1]+y[0]*x[1],x[1]*y[1]);
const sub=(x,y)=>Q(x[0]*y[1]-y[0]*x[1],x[1]*y[1]);
const mul=(x,y)=>Q(x[0]*y[0],x[1]*y[1]);
const div=(x,y)=>Q(x[0]*y[1],x[1]*y[0]);
const lcm=(x,y)=>x/gcd(x,y)*y;
const L=Array.from({length:N},()=>Array.from({length:N},()=>Q(0n))),D=[];
for(let i=0;i<N;i++){
 for(let j=0;j<=i;j++){
  let x=Q(lower[i][j]);
  for(let k=0;k<j;k++)x=sub(x,mul(mul(L[i][k],D[k]),L[j][k]));
  if(i===j){D[i]=x;A(x[0]>0n,"PD");L[i][i]=Q(1n);}
  else L[i][j]=div(x,D[j]);
 }
}
A(String(D.reduce(mul,Q(1n)))==="512,1","det");
const S=L.flat().map(x=>x[1]).reduce(lcm,1n);
const E=D.map(x=>x[1]).reduce(lcm,1n);
const M=L.map(r=>r.map(x=>x[0]*(S/x[1])));
const dia=D.map(x=>x[0]*(E/x[1]));
function floor(a,b){return a>=0n?a/b:-((-a+b-1n)/b);}
function isqrt(n){if(n<2n)return n;let x=1n<<BigInt((n.toString(2).length+1)>>1);for(;;){const y=(x+n/x)>>1n;if(y>=x)return x;x=y;}}
function countAtOrBelow(n,q,bn,bd){
 const target=E*(S*q)**2n*bn;A(target%bd===0n,"threshold");
 const max=target/bd, y=Array(N).fill(0n);
 let equal=0,below=0;const shell=[];
 function walk(i,v){
  if(i<0){if(v===max){equal++;shell.push([...y]);}else below++;return;}
  let offset=-S*n[i];
  for(let j=i+1;j<N;j++)offset+=M[j][i]*y[j];
  const rad=isqrt((max-v)/dia[i]);
  const min=-floor(rad+offset,S*q), maxZ=floor(rad-offset,S*q);
  for(let z=min;z<=maxZ;z++){
   y[i]=q*z-n[i];
   const w=S*q*z+offset;
   walk(i-1,v+dia[i]*w*w);
  }
 }
 walk(N-1,0n);return {equal,below,shell};
}
const a=countAtOrBelow(numer,8n,7n,2n);
const b=countAtOrBelow(numer.map(x=>2n*x),8n,2n,1n);
const c=countAtOrBelow(Array(N).fill(0n),1n,4n,1n);
A(a.below===0&&a.equal===1152,"first layer");
A(b.below===0&&b.equal===30,"second layer");
A(c.below===1&&c.equal===2982,"horizontal minimum");
const P=101n;
function rank101(input){
 const B=input.map(row=>row.map(x=>(x%P+P)%P));
 let k=0;
 for(let j=0;j<N;j++){
  const i=B.findIndex((r,i)=>i>=k&&r[j]!==0n);
  if(i<0)continue;
  [B[i],B[k]]=[B[k],B[i]];
  for(let t=k+1;t<B.length;t++){
   const q=B[t][j], p=B[k][j];
   if(q===0n)continue;
   for(let u=j;u<N;u++)B[t][u]=(p*B[t][u]-q*B[k][u]+P*P)%P;
  }
  k++;
  if(k===N)break;
 }
 return k;
}
A(a.shell[0].every((_,j)=>a.shell.reduce((sum,w)=>sum+w[j],0n)===0n),"centered");
A(rank101(a.shell)===16,"spanning");
console.log("PASS: catalogue det 512, min 4; exact first/second layer squared minima 7/2 and 2, with "+a.equal+" and "+b.equal+" witnesses.");

