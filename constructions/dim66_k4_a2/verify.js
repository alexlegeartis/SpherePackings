"use strict";

// Replay with Node (no dependencies):
//   node verify.js ../../research/novel_routes_20261113/higher_parent/neb72.html
// This uses the published minimum 8 of Nebe's Gamma72 as a source theorem.

const SHA256 = "92daeaf0d93672e82ca5e5bc4479e5775daf59029aa4e983fb464e36ef002044";
const A = [
  [[71,1]],
  [[59,1],[66,-1]],
  [[52,1],[54,-1]],
  [[51,1],[65,1]],
  [[1,1],[7,-1],[13,1]],
  [[5,-1],[9,1],[13,-1]],
];
const PIVOTS = [71,59,52,51,1,5];
const simplex = n => Array.from({length:n+1},(_,i)=>Array.from({length:n},(_,j)=>i===j+1?1n:0n));
const LABELS = simplex(4).flatMap(a => simplex(2).map(b => [...a,...b]));
const fail = (condition, description) => { if (!condition) throw Error(description); };
const q = (M,a,b) => a.reduce((s,u,i) => s+u*b.reduce((t,v,j) => t+v*M[i][j],0n),0n);
const eq = (u,v) => u.length === v.length && u.every((x,i) => x === v[i]);
const sub = (u,v) => u.map((x,i) => x-v[i]);
const expected = (i,j) => i===j ? 8n : i<4 && j<4 ? -2n : i>=4 && j>=4 ? -4n : 0n;
const dot = (G,a,b) => a.reduce((s,[i,c]) => s+b.reduce((t,[j,d]) => t+BigInt(c*d)*G[i][j],0n),0n);

// Fraction-free elimination. For the symmetric G, each successive pivot is
// its corresponding positive leading principal minor (no row exchanges).
function determinant(input, positive=false) {
  const M = input.map(r => r.slice());
  const n = M.length;
  if (!n) return 1n;
  let prev=1n, sign=1n;
  for(let k=0;k<n-1;k++){
    if(positive) fail(M[k][k]>0n,"Gram matrix is not positive definite");
    if(M[k][k]===0n){
      const i=M.findIndex((r,j)=>j>k&&r[k]!==0n);
      if(i<0) return 0n;
      [M[k],M[i]]=[M[i],M[k]];
      sign=-sign;
    }
    const p=M[k][k];
    for(let i=k+1;i<n;i++){
      for(let j=k+1;j<n;j++){
        let v=M[i][j]*p-M[i][k]*M[k][j];
        fail(v%prev===0n,"Nonexact Bareiss division");
        M[i][j]=v/prev;
      }
      M[i][k]=0n;
    }
    prev=p;
  }
  if(positive) fail(M[n-1][n-1]>0n,"Nonpositive full determinant");
  return sign*M[n-1][n-1];
}

function verify(html) {
  const i = html.indexOf('<a NAME="GRAM">');
  const j = html.indexOf('<a NAME="DET">',i+1);
  fail(0<=i&&i<j,"Primary Gram section missing");
  const seq = html.slice(i,j).match(/-?\d+/g).map(BigInt);
  fail(seq.length===1+72*72 && seq[0]===72n,"Expected full 72 by 72 row-major Gram matrix");
  const G = Array.from({length:72},(_,r)=>seq.slice(1+72*r,73+72*r));
  fail(G.every((r,i)=>r[i]%2n===0n && r.every((x,j)=>x===G[j][i])),
       "Not symmetric even integral");
  fail(determinant(G,true)===1n,"Not unimodular");
  const normal = A.map(u=>A.map(v=>dot(G,u,v)));
  fail(normal.every((r,i)=>r.every((x,j)=>x===expected(i,j))),"Normal Gram not K4 direct sum 2A2");
  const minor = A.map(r=>PIVOTS.map(i=>BigInt(r.find(([j])=>i===j)?.[1]??0)));
  fail(determinant(minor)===-1n,"Coefficient rows not certified primitive");
  const d = determinant(normal,true);
  fail(d===96000n,"Wrong rank-six determinant");
  const adj = normal.map((r,i)=>r.map((_,j)=>
    (-1n)**BigInt(i+j)*determinant(normal.filter((_,s)=>s!==j).map(r=>r.filter((_,t)=>t!==i)))));
  fail(normal.every((r,i)=>r.every((_,j)=>r.reduce((s,x,k)=>s+x*adj[k][j],0n)===(i===j?d:0n))),
       "Incorrect adjugate");
  const hs={};
  for(let i=0;i<LABELS.length;i++)for(let j=i+1;j<LABELS.length;j++){
    let a=sub(LABELS[i],LABELS[j]), v=q(adj,a,a);
    fail(0n<v && v<=35200n,"Bad pairwise normal diameter");
    let k=String(v);
    hs[k]=(hs[k]??0)+1;
  }
  fail(JSON.stringify(hs)===JSON.stringify({16000:15,19200:30,35200:60}),
       "Unexpected exact distance histogram");
  const root=[[24,1],[25,1],[40,1]];
  const b=A.map(u=>dot(G,u,root));
  fail(dot(G,root,root)===8n&&eq(b,[0n,1n,-1n,0n,1n,-1n])&&q(adj,b,b)===35200n,
       "Wrong equality witness");
  fail(LABELS.some(x=>LABELS.some(y=>eq(b,sub(x,y)))),"Wrong equality witness label");
  const doubled=[2n,0n,0n,0n,0n,0n];
  fail(LABELS.every(y=>adj.map(r=>r.reduce((s,x,i)=>s+x*(doubled[i]-y[i]),0n)).some(x=>x%d!==0n)),
       "This witness must not be an additive subgroup modulo the period");
  const numerator=225n*90000n*229n**66n, denominator=256n*96000n*228n**66n;
  fail(numerator>denominator,"No strict gain over 16-site");
  fail(225n*12288n*229n**66n>81n*96000n*225n**66n,"No strict gain over 9-site");
  return { dimension:66, cosets:15, normalDeterminant:String(d), covolume:"sqrt(96000)",
    minimumSquared:"229/30 (exact)", normalSquaredDiameter:"11/30",
    exactCenterDensity:"(15/sqrt(96000)) * (229/120)^33",
    normalPairNumeratorHistogram:hs,
    ratioTo16SiteApprox:(15/16)*Math.sqrt(90000/96000)*(229/228)**33,
    ratioTo9SiteApprox:(15/9)*Math.sqrt(12288/96000)*(229/225)**33 };
}
if (typeof module !== "undefined" && module.exports) {
  module.exports={verify};
  if(require.main===module){
    const fs=require("node:fs"), crypto=require("node:crypto");
    const filename=process.argv[2];
    fail(filename,"Pass the archived primary neb72.html (see README).");
    const bytes=fs.readFileSync(filename);
    fail(crypto.createHash("sha256").update(bytes).digest("hex")===SHA256,"Primary input SHA-256 mismatch");
    console.log(JSON.stringify(verify(bytes.toString("utf8")),null,2));
  }
}
