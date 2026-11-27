"use strict";
const ok = (a, msg="assertion") => { if (!a) throw new Error(msg); };
const range = n => Array.from({length:n}, (_,i)=>i);
const weight = n => { let z=0; for (;n;n &= n-1) ++z; return z; };
const shift = (w,i,n) => ((w<<i)|(w>>(n-i)))&((1<<n)-1);
const parse = s => parseInt(s,2);
const C10 = ["0100000011","0011111101","1100101100","0001010111"].flatMap(s =>
  range(10).map(i=>shift(parse(s),i,10)));
const HR = [0,...range(11).map(i=>shift(parse("11011100010"),i,11))];
const C12 = [];
for (let i=0;i<12;i++) for (let j=0;j<i;j++) {
  const s=HR[i]^HR[j]; C12.push(s,4095^s);
}
for (let i=0;i<6;i++) {const s=3<<(2*i); C12.push(s,4095^s);}
const C11 = C12.filter(c=>c<2048), core=C11.filter(c=>c!==1023);
const U=2047;
function checkCode(C,n,size) {
  ok(C.length===size && new Set(C).size===size, "code cardinality");
  let dist=Infinity; const cover=new Set();
  for (let i=0;i<C.length;i++) for (let j=0;j<i;j++) {
    const d=C[i]^C[j];
    dist=Math.min(dist,weight(d));
    if(weight(d)!==4) continue;
    for(const a of range(n))if((d>>a)&1)
      for(let b=0;b<a;b++)if((d>>b)&1) cover.add(a*n+b);
  }
  ok(dist===4 && cover.size===n*(n-1)/2, "minimum and pair support");
  return range(n).map(i=>{
    const T=new Set(C), comp=[];
    while(T.size){
      const seed=T.values().next().value, Q=[seed];T.delete(seed);
      for(let j=0;j<Q.length;j++)for(const c of T) {
        const d=Q[j]^c;
        if (((d>>i)&1)&&weight(d)===4){Q.push(c);T.delete(c);}
      }
      comp.push(Q.length);
    }
    return comp.sort((a,b)=>b-a);
  });
}
ok(checkCode(C10,10,40).every(s=>s.length===1), "P10c connectivity");
ok(checkCode(C12,12,144).every(s=>s.length===1), "Witt extension connectivity");
ok(checkCode(core,11,71).every(s=>s.length===1), "71-core connectivity");
const comp=checkCode(C11,11,72);
ok(comp.slice(0,10).every(s=>s.length===1)&&String(comp[10])==="71,1", "rattler");
ok(!core.some(c=>weight(c^1023)===4 && ((c^1023)&1024)), "rattler index");

const S5=core.map(c=>c^U).filter(s=>weight(s)===5);
ok(core.every(c=>weight(c^U)>=5)&&S5.length===66,"Witt supports");
for(let t=0;t<=U;t++)if(weight(t)===4)
  ok(S5.filter(s=>(t&s)===t).length===1,"Steiner 4-design");
for(const i of range(11)) {
  const inc=S5.filter(s=>(s>>i)&1).map(s=>s^(1<<i));
  ok(inc.length===30,"derived design count");
  for(const j of range(11))if(i!==j)
    ok(inc.filter(s=>(s>>j)&1).length===12,"derived 1-design");
  for(let t=0;t<=U;t++)if(weight(t)===3&&!((t>>i)&1))
    ok(inc.filter(s=>(t&s)===t).length===1,"derived 3-design");
}

let strict=0, ties=[];
for(let b=0;b<=U;b++) {
  if(weight(U^b)<=1)continue;
  // ES2019 stable sort preserves the original core order for tied weights.
  const S=core.map(c=>b^c).sort((a,b)=>weight(a)-weight(b));
  const counts=range(11).map(()=>0);
  let hasStrict=false, hasEqual=false;
  for(let k=1;k<=S.length;k++){
    for(const i of range(11))counts[i]+=(S[k-1]>>i)&1;
    const N=counts.reduce((s,p)=>s+Math.max(4*p,k),0);
    if(N<16*k){hasStrict=true;break;}
    if(N===16*k)hasEqual=true;
  }
  if(hasStrict){strict++;continue;}
  ok(hasEqual,"corner prefix");
  ties.push(b);
  let worst=0;
  for(let m=0;m<=U;m++){
    let v=Infinity;
    for(const s of S)v=Math.min(v,weight(m)+4*weight(s&(U^m)));
    worst=Math.max(v,worst);
  }
  ok(worst===12,"half-vertex torus coverage");
}
ok(strict===2032 && String(ties)==="1025,1028,1088,1536","rounding exclusions");

const K10=[
 [0,0,0,0,0,0,0,0,0,0],
 [-2,0,0,0,0,0,0,0,0,0],
 [-1,0,-1,0,0,0,0,0,-1,-1],
 [-1,0,1,0,0,0,0,0,-1,-1],
 [-1,-1,0,0,0,0,-1,0,0,-1],
 [-1,1,0,0,0,0,-1,0,0,-1],
 [-1,0,0,0,-1,-1,0,0,0,-1],
 [-1,0,0,0,1,-1,0,0,0,-1],
 [-1,0,0,-1,0,0,0,-1,0,-1],
 [-1,0,0,1,0,0,0,-1,0,-1],
];
const h=(v,w)=>v.reduce((s,a,i)=>s+(a-w[i])**2,0);
const mod2 = v => v.reduce((a,x,i)=>a+((x&1)<<i),0);
for(let i=0;i<10;i++){
  ok(C10.includes(mod2(K10[i])^C10[0]),"Best clique membership");
  for(let j=0;j<i;j++)ok(h(K10[i],K10[j])===4,"Best clique distance");
}
const K9=[range(8).map(()=>0),...range(7).map(i=>range(8).map(j=>2*(j===0 || j===1+i))),range(8).map(()=>1)];
for(let i=0;i<9;i++){
  const v=K9[i];
  ok(v.every(c=>c%2===v[0]%2) && v.reduce((a,b)=>a+b,0)%4===0,"double E8 membership");
  for(let j=0;j<i;j++)ok(h(v,K9[j])===8,"double E8 clique distance");
}
function rank3(A) {
  const n=A.length;
  let r=0;
  for(let j=0;j<n;j++){
    const i0=A.findIndex((x,i)=>i>=r&&x[j]%3!==0);
    if(i0<0)continue;
    [A[i0],A[r]]=[A[r],A[i0]];
    const z=A[r][j];for(let k=j;k<n;k++)A[r][k]=(A[r][k]*z)%3;
    for(let i=r+1;i<n;i++){const z0=A[i][j];for(let k=j;k<n;k++)A[i][k]=(6+A[i][k]-z0*A[r][k])%3;}
    ++r;
  }
  return r;
}
ok(rank3(range(8).map(i=>range(8).map(j=>i===j?1:2)))===7,"quadratic Gram rank");

const cmp=([a,b],[c,d])=>a*d<c*b;
const mul=([a,b],[c,d])=>[a*c,b*d];
const pow=([a,b],e)=>[a**BigInt(e),b**BigInt(e)];
const R21=[43n**21n,2n**82n*3n**23n], R22=[11n**22n,2n**46n*3n**21n];
const U21=mul(mul(pow([4n,5n],21),pow([270n,112n],9)),[1n,27n**2n]);
const U22=mul(mul(pow([4n,5n],22),pow([540n,112n],10)),pow([40n,27n*1024n],2));
ok(cmp(U21,pow([18653n,100000n],2))&&cmp(U21,R21),"21D exact comparison");
ok(cmp(U22,pow([32380n,100000n],2))&&cmp(U22,R22),"22D exact comparison");
ok(8n*75n < 44n*27n && 48n*27n < 21n*64n && 48n*27n < 22n*64n, "gap and monotonicity");
console.log("PASS: 11D fixed core, all-rounding exclusions, 10D/11D superperiod criterion, two literal contact cliques, and exact 21D/22D benchmarks");
