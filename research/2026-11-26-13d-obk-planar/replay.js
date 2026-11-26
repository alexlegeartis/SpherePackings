"use strict";
const fs=require("node:fs"), crypto=require("node:crypto");
function verify(data) {
  const ok=(b) => { if(!b) throw Error("verification failed"); };
  const z=data.trim().split(/\s+/).map(Number);
  ok(z.length===7010 && z[0]===96 && z[1]===11);
  const pc=new Uint8Array(4096);
  for(let i=1;i<4096;i++) pc[i]=1+pc[i&(i-1)];
  const short=Array.from({length:4095},(_,j)=>j+1)
      .filter(x=>pc[x]<=3).sort((a,b)=>pc[a]-pc[b]||a-b);
  const banks=[], cornerCounts=[0,0,0], extEq=[0,0,0,0];
  for(let j=0;j<96;j++){
    const C=z.slice(3+73*j,75+73*j);
    ok(z[2+73*j]===72 && C.every(x=>Number.isInteger(x)&&x>=0&&x<2048)
       && new Set(C).size===72 && C.every((x,i)=>C.slice(0,i).every(y=>pc[x^y]>=4)));
    const D=[...new Set(C.flatMap(x=>C.map(y=>x^y)))];
    const deep=[];
    for(let v=0;v<2048;v++) if(D.every(e=>pc[e^v]>=3)) deep.push(v);
    if(!deep.length) continue;
    ok(deep.length===1 && deep[0]===2047);
    banks.push(j+1);
    const E=[...D,...D.map(e=>e^4095)];
    const on=new Uint8Array(4096); for(const e of E)on[e]=1;
    for(let u=0;u<4096;u++){
      if(on[u]){cornerCounts[0]++;continue;}
      const A=short.filter(m=>on[u^m]);
      function four(start,needed,used){
        if(!needed)return true;
        for(let k=start;k<A.length;k++)
          if(!(A[k]&used)&&four(k+1,needed-1,used|A[k]))return true;
        return false;
      }
      ok(A.length && pc[A[0]]<=2 && four(0,4,0));
      cornerCounts[pc[A[0]]]++;
    }
    const local=[0,0,0,0];
    for(let F=0;F<4096;F++){
      const m=pc[F];if(m%4)continue;
      const k=4095^F, p=new Uint8Array(4096);
      for(const e of E)p[e&k]=1;
      const err=3-m/4;
      const flips=[0,...short].filter(t=>(t&F)===0 && pc[t]<err);
      for(let u=k;;u=(u-1)&k){
        if(flips.every(t=>!p[t^u]))local[m/4]++;
        if(u===0)break;
      }
    }
    ok(local.join()==='0,0,0,1');for(let i=0;i<4;i++)extEq[i]+=local[i];
  }
  ok(banks.join()==='6,8,13,15,20,21,22,24,30,32,34,37,39,43,44,47,59,61,62,63,64,65,66,67,68,69,70');
  ok(cornerCounts.reduce((a,b)=>a+b,0)===110592);
  return {viable_codes:banks,cornerCounts,folded_equality_halfgrid:extEq};
}
(async()=>{
 const input=process.argv[2]?fs.readFileSync(process.argv[2],"utf8"):await (async()=>{ const r=await fetch("https://users.aalto.fi/~pat/opt/72_4.txt"); if(!r.ok)throw Error(r.status);return r.text();})();
 const norm=input.trim().split(/\s+/).map(Number).join("\n")+"\n";
 if(crypto.createHash("sha256").update(norm).digest("hex")!=="0f3fe74e631f11d3d3263c1c476d86ad142298ed6e2184f2b06370c771519d23") throw Error("catalogue hash mismatch");
 console.log(JSON.stringify(verify(norm)));
})().catch(e=>{console.error(e);process.exitCode=1;});
