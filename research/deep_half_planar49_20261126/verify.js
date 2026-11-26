'use strict';
const fs=require('node:fs'), path=require('node:path'), assert=require('node:assert/strict');
const at=(p)=>fs.readFileSync(path.join(__dirname,p),'utf8');
const w=JSON.parse(at('witness.json'));
const html=at('../density_breakthroughs_20261113/seed_deformations/sources/CQ48a.html');
const tri=html.split('<a NAME="GRAM">')[1].split('<a NAME="DIVISORS">')[0].split('\n').slice(2).map(l=>l.replace(/<br\s*\/?>/ig,'').trim()).filter(l=>/^[-+\d ]+$/.test(l)).map(l=>l.split(/\s+/).map(Number));
assert.equal(tri.length,48); assert(tri.every((r,i)=>r.length===i+1));
const parent=tri.map((r,i)=>tri.map((s,j)=>i>=j?r[j]:s[i]));
const bil=(x,y)=>x.reduce((z,c,i)=>z+c*parent[i].reduce((s,d,j)=>s+d*y[j],0),0);
const old=JSON.parse(at('../p48n_hurwitz_20261122/results.json')).rows[0].half_coset_numerators_in_catalogue_basis[0];
assert(w.v.every((a,i)=>(a-old[i])%2===0)); assert.equal(bil(w.v,w.v),16);
assert.deepEqual(w.vs.map(x=>[bil(x,w.v),bil(x,x)]),[[0,6],[1,6],[2,6],[3,6],[4,6],[5,6],[6,6],[7,8],[8,8]]);
const C0=[14550592,23003136,11169024,3178496,481216,32768,768]; let no=52416000n,de=1n;
for(let p=0;p<=5;p++){const z=C0.reduce((u,x,j)=>u+BigInt(x)*BigInt(j*j)**BigInt(p),0n);assert.equal(z*de,no);no*=96n*BigInt(2*p+1);de*=BigInt(48+2*p);}
assert(32n*95n**47n < 27n*96n**47n);
const E=[96,1,4,9,16,25,36,0,32,0,36,25,16,9,4,1],cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
let verts=0,good=[];
for(let a=0;a<16;a++)for(let b=0;b<16;b++){
 if(E[a]>7)continue;
 const s=Math.max(1,E[a]), C=[[1,0,0,s],[0,0,1,E[b]],[0,1,0,0],[1,-2,0,0],[-1,0,1,0],[-1,0,0,-7],[0,0,-1,-40]];
 for(let j=0;j<=8;j++)for(let i=-8;i<=8;i++){if(j===0&&i<=1||j===1&&i===0)continue;
 const e=E[((a*i+b*j)%16+16)%16];if((i*i+j*j+Math.min(0,i*j))*s<e)C.push([i*i,2*i*j,j*j,e]);}
 for(let i=0;i<C.length;i++)for(let j=i+1;j<C.length;j++){const x=cross(C[i],C[j]);for(let k=j+1;k<C.length;k++){
 let w=dot(x,C[k]);if(w===0)continue;const y=cross(C[j],C[k]),z=cross(C[k],C[i]);let g=x.map((s,l)=>s*C[k][3]+y[l]*C[i][3]+z[l]*C[j][3]);if(w<0){w=-w;g=g.map(t=>-t);}
 if(C.some(v=>dot(v,g)<v[3]*w))continue;
 verts++;const [A,B,Cc]=g.map(BigInt),W=BigInt(w),num=A*Cc-B*B-32n*W*W;assert(num>=0n);
 if(num===0n)good.push([a,b,g.map(v=>v/w)])}}}
console.log(verts,JSON.stringify([...new Set(good.map(JSON.stringify))]));

