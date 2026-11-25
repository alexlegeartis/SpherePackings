'use strict';
function verify(d){
const G=d.E7_gram,M=d.phase_isometry,B=d.phase_syndromes,A=d.alpha.slice(1),field=d.field_polynomial;
const F=[[2,-1,-1,0,-1,0,0],[0,3,-1,0,-1,0,0],[0,0,4,-3,-2,0,0],[0,0,0,5,-2,0,0],[0,0,0,0,4,-5,0],[0,0,0,0,0,3,-4],[0,0,0,0,0,0,1]],D=[30,10,5,3,3,5,40];
const range=n=>Array.from({length:n},(_,i)=>i),bits=x=>range(7).map(i=>x>>i&1);
function ok(a){if(!a)throw Error("certificate");}
for(const i of range(7))for(const j of range(7))ok(F.reduce((v,r,k)=>v+D[k]*r[i]*r[j],0)===60*G[i][j]);
const det=a=>a.length===1?a[0][0]:a[0].reduce((s,t,i)=>s+(i%2?-1:1)*t*det(a.slice(1).map(r=>r.filter((_,j)=>i!==j))),0);
ok(det(G)===2&&M.length===128&&new Set(M).size===128&&M[64]===82&&M.every(q=>q>=0&&q<128));
for(const i of range(128))for(const j of range(128))ok(M[i^j]===(M[i]^M[j]));
function gf(x,y){let r=0;while(y){if(y&1)r^=x;x<<=1;if(x&64)x^=field;y>>=1;}return r;}
for(let i=1;i<64;i++)ok(new Set(range(64).map(j=>gf(i,j))).size===64);
ok(new Set(A).size===4&&A.every(x=>x>=0&&x<64));
const C=[range(128).map(()=>99),range(128).map(()=>99)],x=range(7).map(()=>0),r=bits(82);
function enumerate(k,b){if(k<0){ok((1440-b)%60===0);const norm=(1440-b)/60;for(const m of [0,1]){let mask=0;for(const i of range(7)){const v=x[i]-m*r[i];if(v&1){mask=-1;break;}mask|=(v>>1&1)<<i;}if(mask>=0)C[m][mask]=Math.min(C[m][mask],norm);}return;}
let s=0;for(let i=k+1;i<7;i++)s+=F[k][i]*x[i];let h=0;while(D[k]*(h+1)*(h+1)<=b)h++;for(let y=Math.ceil((-h-s)/F[k][k]);y<=Math.floor((h-s)/F[k][k]);y++){x[k]=y;enumerate(k-1,b-D[k]*(F[k][k]*y+s)**2);}}
enumerate(6,1440);ok(C.every(a=>a.every(z=>z<=24)));
function merge(types){let a=Array(8192).fill(999),p=Array(8192).fill(999);a[0]=0;
for(let j=0;j<4;j++){const n=Array(8192).fill(999),q=Array(8192).fill(999);for(let k=0;k<128;k++){const s=(k&63)|(gf(A[j],k&63)<<6)|((k>>6)<<12), w=C[types[j]][M[k]];for(let l=0;l<8192;l++)if(a[l]<999){const z=l^s;n[z]=Math.min(n[z],a[l]+w);q[z]=Math.min(q[z],k?a[l]+w:p[l]+w);}}a=n;p=q;}return [a,p];}
const [ev,no0]=merge([0,0,0,0]),[od]=merge(d.half_shift_blocks);ok(no0[0]>=32&&ev[0]===0&&d.half_shift_blocks.join(",")==="0,1,1,1"&&d.half_radical_mask===82);
const tag=m=>B.reduce((s,x,k)=>m>>k&1?s^x:s,0);
ok(B.length===6&&d.local_lattice==="D6"&&tag(63)===4096);
let min=999;const out={};
for(let x=0;x<64;x++){
 let p=bits(x&31).reduce((s,v)=>s+v,0);p+=p%2;const k=x>>5,s=tag(x);
 for(const o of [0,1]){
 const local=o?(k===(p/2%2)?6:14):(p?4*p:16*k);
 const ext=o?od[s^d.offset_syndrome]:ev[s];
 if(x||o){ok(local+ext>=32);min=Math.min(min,local+ext);}
 const key=[o,local,ext].join(",");out[key]=(out[key]??0)+1;
 }
}
return {minimumTimes4:min,kernelTimes4:no0[0],doubledTag:tag(63),histogram:out};
}
console.log(JSON.stringify(verify(require('./cert.json')),null,2));
