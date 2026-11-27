"use strict";
const ok=(z,m)=>{if(!z)throw Error(m);};
const es=(a,b)=>a*a-a*b+b*b, md=(x,n)=>(x%n+n)%n;
function mul4(a,b){let x=(b&1?a:0)^(b&2?a<<1:0);return x&4?x^7:x;}
const hex=(a,b,c)=>[a,b,c,a^b^c,a^mul4(2,b)^mul4(3,c),a^mul4(3,b)^mul4(2,c)];
const Hex=[];for(let a=0;a<4;a++)for(let b=0;b<4;b++)for(let c=0;c<4;c++)Hex.push(hex(a,b,c));
const wh={};
for(const z of Hex){const h=z.filter(Boolean).length;wh[h]=(wh[h]??0)+1;}
ok(wh[0]===1&&wh[4]===45&&wh[6]===18,"hex weight");
for(const x of Hex)for(const y of Hex)ok(x.reduce((z,a,i)=>z^mul4(a,mul4(y[i],y[i])),0)===0,"Hermitian dual");
const cost=Array.from({length:9},(_,r)=>Array.from({length:4},(_,h)=>{
 let v=99;for(let a=-6;a<=6;a++)for(let b=-6;b<=6;b++)if(md(a,3)===r%3&&md(b,3)===Math.floor(r/3)
   &&(md(a,2)+2*md(b,2))===h)v=Math.min(v,es(a,b));ok(v<=27,"coordinate cost");return v;
}));
const want={0:1,4:756,6:4032,8:20412,10:60480,12:119196,14:163296,16:122472,18:40544,24:252};
let seen={},labels=[],ho={};
for(let i=0;i<9**6;i++){
 let z=i,rows=[],dual=true,w=0,coord=[];
 for(let j=0;j<6;j++){let r=z%9,a=r%3,b=Math.floor(r/3);z=Math.floor(z/9);
   rows.push(cost[r]);dual&&=((a+b)%3===0);w+=+(a!==0);coord.push(a);}
 let d=99;
 for(const x of Hex){const q=rows[0][x[0]]+rows[1][x[1]]+rows[2][x[2]]+rows[3][x[3]]+rows[4][x[4]]+rows[5][x[5]];if(q<d)d=q;}
 seen[d]=(seen[d]??0)+1;
 ok((d===24)===(dual&&w%3===2),"canonical deep identification");
 if(d===24)labels.push(coord);
 if(dual)ho[d]=(ho[d]??0)+1;
}
ok(JSON.stringify(seen)===JSON.stringify(want)&&JSON.stringify(ho)===JSON.stringify({0:1,12:252,18:224,24:252}),"depth histograms");

const spectrum={};
for(let i=0;i<729;i++){
 let w=[],t=i;for(let j=0;j<6;j++){w.push(t%3);t=Math.floor(t/3);}
 let z=0;for(const a of labels)if(a.reduce((n,u,j)=>n+u*w[j],0)%3===0)z++;
 const l=(3*z-252)/2;ok(Number.isInteger(l),"integer Fourier coefficient");spectrum[l]=(spectrum[l]??0)+1;
}
ok(spectrum[252]===1&&spectrum[9]===476&&spectrum[-18]===252,"CT spectrum");

const E=[
 [1,-1,-1,-1,-1,-1,-1,1],
 [2,2,0,0,0,0,0,0],
];
for(let i=0;i<6;i++){let x=Array(8).fill(0);x[i]=-2;x[i+1]=2;E.push(x);}
function det(a){if(a.length===1)return a[0][0];return a[0].reduce((s,x,j)=>
 s+(x?(-1)**j*x*det(a.slice(1).map(r=>r.filter((_,k)=>k!==j))):0),0);}
ok(Math.abs(det(E))===256,"E8 basis determinant");
for(const a of E)ok(a.every(x=>md(x,2)===md(a[0],2))&&a.reduce((u,v)=>u+v,0)%4===0,"E8 membership");
const roots=[];
for(let x=0;x<256;x++){
 const v=Array.from({length:8},(_,j)=>E.reduce((s,a,i)=>s+((x>>i&1)?a[j]:0),0));
 const q=v.reduce((s,x)=>s+x*x,0)/8;ok(Number.isInteger(q),"even E8");
 if(q%2)roots.push(x);
}
ok(roots.length===120,"120 nonsingular E8 classes");
const pw=x=>{let y=0;while(x){x&=x-1;y^=1;}return y},eig={};
for(let c=0;c<256;c++){let a=roots.reduce((s,x)=>s+(pw(c&x)?-1:1),0);eig[a]=(eig[a]??0)+1;}
ok(eig[120]===1&&eig[8]===120&&eig[-8]===135,"E8 spectrum");

const A=[Array(8).fill(0)];
for(let j=0;j<6;j++){let v=Array(8).fill(0);v[j]=2;v[7]=-2;A.push(v);}
A.push([0,0,0,0,0,0,-2,-2],[1,1,1,1,1,1,-1,-1]);
for(let i=0;i<A.length;i++){
 let v=A[i];ok(v.every(x=>md(x,2)===md(v[0],2))&&v.reduce((x,y)=>x+y,0)%4===0,"simplex membership");
 for(let j=0;j<i;j++)ok(v.reduce((f,t,k)=>f+(t-A[j][k])**2,0)===8,"simplex side");
}
ok(A.length===9,"9 simplex points");

function mul(a,b){return [a[0]*b[0],a[1]*b[1]];}
function power(a,n){return [a[0]**BigInt(n),a[1]**BigInt(n)];}
function less(a,b){return a[0]*b[1]<a[1]*b[0];}
function flip(a){return [a[1],a[0]];}
function D(a,b){return mul(power(mul(a,[1n,4n]),10),mul([1n,27n],power(flip(b),4)));}
const U=[3n,50n],a=[5n,12n],b=[5n,8n],mu=[8n,3n];
const lowA=mul(U,power([6n,5n],4)),lowB=mul(U,power([4n,5n],4));
const middleA=mul([7n,1n],D(mu,a)),middleB=mul([7n,1n],D(mu,b));
const cHi=mul([102n,432n],power([8n,9n],6));
const cA=middleA;
const pHi=mul([7n,9n],D([47n,12n],b));
const target2=[3n**40n,(2n**20n)*(5n**21n)];
for(const t of [lowA,cHi,D(mu,[1n,3n]),[1n,27n]])ok(less(t,cA),"compact endpoint");
for(const t of [lowB,middleB,mul([1n,16n],D([4n,1n],b)),D([13n,4n],b),[1n,27n]])ok(less(t,pHi),"period-free endpoint");
for(const t of [cHi,middleA,pHi,mul([1n,16n],D([4n,1n],[7n,20n])),[7n,48n]])ok(less(power(t,2),target2),"benchmark comparison");
ok(less(mul([798n,10000n],mul([720n,1n],power([100n,314n],6))),U),"12D upper input conversion");
ok(!less(power([3n,16n],2),target2),"full occupancy positive trial");
ok(less(flip(mul([729n,1n],target2)),power([2n,3n],7)),"all PSD metrics gain gate");
console.log("PASS: 531441 CT ternary depths, 252 all-canonical deep classes, E8/CT spectra and all exact 20D endpoints.");
