'use strict';const assert=require('node:assert/strict');
const E=[96,1,4,9,16,25,36,17,32,17,36,25,16,9,4,1],cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
let verts=0,good=[];
for(let a=0;a<16;a++)for(let b=0;b<16;b++){
 if(a!==1)continue;
 const C=[[1,0,0,E[a]],[0,0,1,E[b]>=16?E[b]:0],[0,1,0,0],[1,-2,0,0],[-1,0,1,0],[-1,0,0,-5],[0,0,-1,-40]];
 for(let j=0;j<=1;j++)for(let i=-8;i<=8;i++){if(j===0&&i<=1||j===1&&i===0)continue;
 const e=E[((a*i+b*j)%16+16)%16];if(j&&e<16)continue;if((i*i+j*j+Math.min(0,i*j))*E[a]<e)C.push([i*i,2*i*j,j*j,e]);}
 for(let i=0;i<C.length;i++)for(let j=i+1;j<C.length;j++){const x=cross(C[i],C[j]);for(let k=j+1;k<C.length;k++){
 let w=dot(x,C[k]);if(w===0)continue;const y=cross(C[j],C[k]),z=cross(C[k],C[i]);let g=x.map((s,l)=>s*C[k][3]+y[l]*C[i][3]+z[l]*C[j][3]);if(w<0){w=-w;g=g.map(t=>-t);}
 if(C.some(v=>dot(v,g)<v[3]*w))continue;
 verts++;const [A,B,Cc]=g.map(BigInt),W=BigInt(w),num=A*Cc-B*B-32n*W*W;assert(num>=0n);
 if(num===0n)good.push([a,b,g.map(v=>v/w)])}}}
console.log(verts,JSON.stringify([...new Set(good.map(JSON.stringify))]));

