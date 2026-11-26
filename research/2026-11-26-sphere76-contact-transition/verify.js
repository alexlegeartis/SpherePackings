"use strict";
// Standalone Node.js, exact integers. Leech coordinates divided by sqrt(8).
let W=[[1,-1,-1,-1,-1,1,1,1,1,1,1,1,-1,-1,1,-1,1,1,3,-1,-1,-1,1,-1],[-1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,-1,1,-1,-1,-1,3,-1,1,-1,-1,1,-1],[-1,1,-3,-1,-1,1,1,-1,-1,1,-1,-1,1,-1,-1,-1,1,-1,-1,1,-1,-1,-1,-1],[0,2,0,2,2,0,0,0,0,0,0,0,0,0,2,0,0,-2,0,-2,0,2,0,2],[1,-3,1,1,-1,-1,-1,1,1,-1,1,-1,-1,1,-1,1,-1,-1,-1,-1,1,1,-1,1]],pc=n=>{n-=(n>>>1)&1431655765;n=(n&858993459)+((n>>>2)&858993459);return ((n+(n>>>4)&252645135)*16843009)>>>24;};let roots=[],xx;for(let i=0;i<24;i++)for(let j=0;j<i;j++)for(let a of [-4,4])for(let b of [-4,4]){let x=new Int8Array(24);x[i]=a;x[j]=b;roots.push(x)}
for(let n=0;n<4096;n++){let m=0;for(let j=0;j<12;j++)if(n>>j&1)m^=0xAE3<<j;m|=(pc(m)&1)<<23;if(pc(m)==8){let S=[];for(let j=0;j<24;j++)if(m>>j&1)S.push(j);for(let s=0;s<256;s++)if(!(pc(s)&1)){let x=new Int8Array(24);for(let j=0;j<8;j++)x[S[j]]=2-4*(s>>j&1);roots.push(x)}}for(let k=0;k<24;k++){let x=new Int8Array(24);for(let j=0;j<24;j++)x[j]=(j==k?-3:1)*(1-2*(m>>j&1));roots.push(x)}}
let dot=(a,b)=>{let s=0;for(let j=0;j<24;j++)s+=a[j]*b[j];return s/8};let B=[[21,9,9,12,12],[9,29,13,20,16],[9,13,29,16,20],[12,20,16,32,16],[12,16,20,16,32]],labs=Array.from({length:6},(_,i)=>B.map((_,j)=>Number(i===j+1))),E={};let keys=[];for(let i of [0,4,5])for(let j of [1,2,3]){let q=labs[j].map((x,k)=>x-labs[i][k]);let beta=q.reduce((s,x,h)=>s+x*B[h].reduce((t,v,g)=>t+v*q[g],0),0);if(beta===29){const k=q.join();E[k]=[];keys.push(k)}}

function det(A){let a=A.map(r=>Array.from(r,z=>BigInt(z))),n=a.length,den=1n,sg=1n;for(let k=0;k<n-1;k++){if(a[k][k]===0n){let i=k+1;while(i<n&&!a[i][k])i++;if(i==n)return 0n;[a[k],a[i]]=[a[i],a[k]];sg=-sg}let piv=a[k][k];for(let i=k+1;i<n;i++)for(let j=k+1;j<n;j++)a[i][j]=(piv*a[i][j]-a[i][k]*a[k][j])/den;den=piv}return a[n-1][n-1]*sg}


function assert(ok){if(!ok)throw Error("failed");}assert.equal=(a,b)=>assert(a===b);assert.deepEqual=(a,b)=>assert(JSON.stringify(a)===JSON.stringify(b));
assert.equal(roots.length,196560);
let Z=[0,0,0,0,0,2,2,0,0,0,0,0,2,0,-2,2,0,2,0,-2,0,2,0,0],Wgram=W.map(a=>W.map(b=>dot(a,b)));assert.equal(det(Wgram),300n);assert.equal(dot(Z,Z),4);assert(W.every(w=>dot(w,Z)===0));
for(let i=0;i<5;i++)for(let j=0;j<5;j++)assert.equal(Wgram[i].reduce((s,x,k)=>s+x*B[k][j],0),i===j?60:0)
assert(W.concat([Z]).every(w=>roots.some(x=>x.every((t,j)=>t===w[j]))));
let P={};for(let x of roots){let k=W.map(w=>dot(w,x)).join();(P[k]??=[]).push(x)}
for(let i=1;i<6;i++)assert(P[labs[i].join()].length>0);
let cm=0,old=0;let got=[];for(let i=0;i<6;i++)for(let j=i+1;j<6;j++){let q=labs[j].map((z,k)=>z-labs[i][k]);let beta=q.reduce((ss,y,k)=>ss+y*B[k].reduce((s1,t,h)=>s1+t*q[h],0),0);if(beta===32){assert.equal(P[q.join()].length,2011);old+=2011;}}
let expected=[78,62,78,78,62,78];
for(let k of keys){let xs=P[k];assert.equal(xs.length,2266);let m=xs.filter(x=>dot(x,Z)===2).length;got.push(m);cm+=m}assert.deepEqual(got,expected);assert.equal(cm,436);assert.equal(old,12066);
let ds=[{"key":"0,0,0,1,0","points":2011,"gcd":"1","hit":[[[1260,58,1355,537,201,1682,990,371,1477,1086,1410,597,1596,985,681,855,1677,1767,478],"14"],[[90,710,439,1989,1858,920,1833,1495,598,175,106,1622,1806,276,40,1483,483,1065,884],"6"],[[1908,1721,599,1254,1194,1219,1770,1421,863,41,1025,1233,1650,1184,88,789,919,1822,236],"11"]]},{"key":"-1,1,0,0,0","points":2011,"gcd":"1","hit":[[[1157,1285,1704,1751,1008,560,1078,458,1996,1961,381,1843,1084,774,1359,419,486,1516,209],"1"]]}];let gcd=(a,b)=>b?gcd(b,a%b):a;for(let z of ds){let n=0n;let xs=P[z.key];for(let [idx,idxvalue] of z.hit){let dw=idx.map(i=>Array.from(xs[i],(x,k)=>x-xs[0][k]));let vv=det(W.concat(dw));if(vv<0)vv=-vv;assert.equal(vv,BigInt(idxvalue)*300n*8n**12n);n=gcd(n,BigInt(idxvalue));}assert.equal(n,1n)}
console.log("Verified two connected 3-coset components, 12066 internal unoriented edges per cell, 436 new cross contacts at endpoint; profiles",got.join(","));
