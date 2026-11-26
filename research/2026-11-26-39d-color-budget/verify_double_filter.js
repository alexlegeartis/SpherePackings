// Exact integer necessary filters. Run with Node.js from any directory.
const fs=require("node:fs"), path=require("node:path");
const read=n=>JSON.parse(fs.readFileSync(path.join(__dirname,"..",n),"utf8"));
const f=read("e7_hermitian35_20261123/certificate.json"),g=read("e7_lattice38_20261123/certificate38.json");
const C=f.E7_gram,M=f.phase_isometry,A=f.alpha,B=g.horizontal_generators_times_2[1],I=[0,911,4096,5007];
function mul(a,b){let c=0;for(;b;b>>=1,a<<=1){if(a&64)a^=67;if(b&1)c^=a;}return c}
const S=A.map(a=>M.map((_,i)=>(i&63)|(mul(a,i&63)<<6)|(i&64?4096:0)));
function norm(x){let k=0;for(let i=0;i<7;i++)for(let j=0;j<7;j++)k+=x[i]*C[i][j]*x[j];return k}
const roots=M.map((b,i)=>b!==82&&norm(Array.from({length:7},(_,j)=>(b>>j)&1))%4===2?i:-1).filter(i=>i>=0);
const low=new Set([0]);for(let i=0;i<5;i++){for(let l=0;l<128;l++)if(M[l]!==0&&M[l]!==82)low.add(S[i][l]);for(let j=0;j<i;j++)for(const x of roots)for(const y of roots)low.add(S[i][x]^S[j][y])}
const reps=[];for(let x=1;x<8192;x++){const c=I.map(k=>x^k);if(x!==Math.min(...c))continue;let t=new Set([0,5007]);for(const y of c)if(low.has(y))for(const a of [...t])t.add(a^y);if(t.size<=4)reps.push(x)}
let F=new Int32Array(8192).fill(100000);F[0]=0;
for(let b=0;b<5;b++){
 const D=new Int32Array(128).fill(100000),x=[];
 function rec(i,m){if(i===7){D[m]=Math.min(D[m],norm(x));return}const start=-4+(B[b][i]&1);for(let z=start;z<=4;z+=2){x[i]=z;rec(i+1,m+(((z-B[b][i])/2)&1)*(1<<i))}}
 rec(0,0);const next=new Int32Array(8192).fill(100000);for(let s=0;s<8192;s++)if(F[s]<100000)for(let l=0;l<128;l++){const j=s^S[b][l];next[j]=Math.min(next[j],F[s]+D[M[l]])}F=next;
}
const assert=require("node:assert/strict");assert.equal(low.size,7939);assert.equal(roots.length,63);assert.equal(reps.length,18);
assert(I.every(i=>low.has(117^i)));
for(const d of reps){const cheap=I.map(i=>d^i).filter(x=>low.has(x));assert.equal(cheap.length,2);assert.equal(cheap[0]^cheap[1],5007);}
assert.deepEqual(reps.filter(x=>F[x]>16),[1496]);
console.log({representatives:reps, witnessCosts:reps.map(x=>F[x]), onlySurvivingCoset:I.map(i=>i^1496)});
