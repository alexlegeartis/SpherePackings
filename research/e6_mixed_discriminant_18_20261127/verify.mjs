/* Exact standalone (plain Node.js). No external packages or data. */
"use strict";
const ok=(b,m="assertion")=>{if(!b)throw new Error(m);};
const range=n=>Array.from({length:n},(_,i)=>i);
const all=range(64);
const slots=[0,1,2];
const depth=[0,3,6,8,5,2];
const perms=x=>Array.from(new Set(["012","021","102","120","201","210"].map(p=>[...p].map(i=>x[i]).join("")))).map(s=>[...s].map(Number));
const shapes=n=>range(n).flatMap(i=>range(n).slice(i).flatMap(j=>range(n).slice(j).map(k=>[i,j,k])));
const G=[[2,-1,0,0,0,0],[-1,2,-1,0,0,0],[0,-1,2,-1,0,-1],[0,0,-1,2,-1,0],[0,0,0,-1,2,0],[0,0,-1,0,0,2]];
const H=[[4,5,6,4,2,3],[5,10,12,8,4,6],[6,12,18,12,6,9],[4,8,12,10,5,6],[2,4,6,5,4,3],[3,6,9,6,3,6]];
const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
const norm=x=>G.reduce((s,g,i)=>s+x[i]*dot(g,x),0);
const mask=(x)=>x.reduce((s,y,i)=>s|((y&1)<<i),0);
const quad=all.map(x=>(norm(range(6).map(i=>x>>i&1))/2)&1);
const cls=quad.map((q,i)=>i?(q?1:2):0);
for(const i of range(6))for(const j of range(6))ok(dot(G[i],H.map(r=>r[j]))===3*(i===j),"Gram inverse");
ok(quad.reduce((a,b)=>a+b)===36,"minus-type count");
const ranges=H.map((r,i)=>Math.floor(Math.sqrt(16*r[i])));
const nearest=[];
for(let t=0;t<=1;t++) {
  const r=H.map(row=>2*t*row[0]); // 3*(2 t w)
  const xx=Array(6), nn=Array(6), ans=all.map(()=>Infinity);
  let stock=Infinity;
  const scan=i=>{
    if(i===6){
      const z=norm(xx);
      if(z<=48){
        const mm=mask(nn);
        if(z<ans[mm])ans[mm]=z;
        if(t===0&&z>0&&z<stock)stock=z;
      }
      return;
    }
    for(let n=Math.ceil((-ranges[i]-r[i])/3);n<=Math.floor((ranges[i]-r[i])/3);n++){
      nn[i]=n;xx[i]=3*n+r[i];scan(i+1);
    }
  };
  scan(0);
  for(const x of all)ok(ans[x]===6*depth[3*t+cls[x]],"torus depth");
  if(t===0)ok(stock===18,"stock minimum");
  nearest.push(ans);
}

const P=[[27,27,27],[27,-3,3],[27,3,-5]];
const six=6n;
const rows=h=>[...h].map(Number);
const ff=(h,s)=>{const f=perms(s);return six/BigInt(f.length)*f.map(x=>x.reduce((a,j,i)=>a*BigInt(P[h[i]%3][j%3]*(h[i]>=3&&j>=3?-1:2)),1n)).reduce((a,b)=>a+b,0n);};
const ffden=6n*54n**3n;
function certificate(k,cert,den,eq,ed) {
  for(const s of shapes(k))if(s.reduce((a,j)=>a+depth[j],0)>=12) {
    ok(-cert.reduce((a,[h,y])=>a+ff(rows(h),s)*BigInt(y),0n)>=BigInt(den)*ffden,"dual half/full");
  }
  ok((den+cert.reduce((a,[_,y])=>a+y,0))*ed===eq*den,"objective");
}
certificate(3,[["001",72],["002",216],["022",729]],7,1024,7);
certificate(6,[["001",314496],["003",7168],["005",1161216],["025",874800],["055",2694141],["235",35235]],24752,24576,119);
const parallel=["021","012","000"];
function check_line_system(){
  const ors=["111","112","121","122"].map(rows), code=shapes(3).flatMap(perms);
  const sep=(x,y,s)=>x.reduce((a,_,i)=>a+depth[s[i]+3*(x[i]!==y[i])],0);
  for(let i=0;i<4;i++)for(let a=0;a<3;a++)for(let b=0;b<3;b++){
    const line0=range(3).map(j=>[j,j,j]),line1=range(3).map(j=>ors[i].map((k,l)=>(j*k+[0,a,b][l])%3));
    const C=shapes(3).filter(s=>perms(s).some(x=>line0.every(v=>line1.every(w=>sep(v,w,x)>=12)))).map(x=>x.join("")).join(",");
    ok(C===(i!==0?"":a+b===0?"022,112":a*b===2?"000,001,111":""),"line relations");
  }
  for(const [z,s0] of [[0,"022"],[0,"112"],[1,"000"],[1,"001"],[1,"111"]]){
    let s=rows(s0);
    let a=-sumLine([["001",0,3024],["022",0,27],["022",1,13095]],s,z);
    ok(4n*a>=BigInt(4+3*z)*238n*2n*ffden,"line dual");
  }
  ok(119*(238+3024+27+13095)===8192*238,"line objective");
  function sumLine(Y,s,z) {return Y.reduce((a,[h,t,y])=>a+ff(rows(h),s)*BigInt(y*(z&&t?-1:2)),0n);}
}
check_line_system();

const q0=x=>Number(Boolean((x&1)&&(x&2)))^Number(Boolean((x&4)&&(x&8)))^Number(Boolean((x&16)&&(x&32)));
const T=all.filter(q0);
const forms=T.map(t=>all.map(x=>q0(x^t)^1));
const D=forms.map(r=>r.map((x,i)=>i?(x?1:2):0));

function lagrangian(g,ab,ac,bc) {
  const N=[0,0,0],B=[];
  const next=i=>1<<(2*N[i]+++6*i);
  for(let j=0;j<g;j++){
    const [a,b,c]=slots.map(next);
    B.push(a^b^c,(a^b)<<1,(a^c)<<1);
  }
  for(const [j,k,e] of [[0,1,ab],[0,2,ac],[1,2,bc]])for(let i=0;i<e;i++){
    const z=next(j)^next(k);B.push(z,z<<1);
  }
  for(let j=0;j<3;j++)while(N[j]<3)B.push(next(j));
  const U=[0];
  for(const b of B)for(let i=0,n=U.length;i<n;i++)U.push(U[i]^b);
  ok(U.length===512&&new Set(U).size===512,"stabilizer cardinality");
  return U;
}
const F=[];
for(let g=0;g<=3;g++)for(let ab=0;ab<=3-g;ab++)for(let ac=ab;ac<=3-g-ab;ac++)for(let bc=ac;bc<=Math.min(3-g-ab,3-g-ac);bc++)F.push([g,ab,ac,bc]);
ok(F.length===16,"complete shapes");
let total=0,exception=0;
for(const f of F) {
  const U=lagrangian(...f),x=U.map(i=>i&63),y=U.map(i=>i>>6&63),z=U.map(i=>i>>12);
  const bad=new Uint8Array(512),reps=new Uint32Array(512),used=[],cand=[];
  for(let a=0;a<28;a++)for(let b=0;b<28;b++)for(let c=0;c<28;c++){
    total++;
    const da=D[a],db=D[b],dc=D[c];
    used.length=0;cand.length=0;bad.fill(0);reps.fill(0);
    for(let i=1;i<512;i++){
      const j=da[x[i]], k=db[y[i]], l=dc[z[i]], d=j+k+l;
      if(!(d%2)&&d!==4){bad[i]=1;used.push(i);}
      else if(d===1||(j===1&&k===1&&l===1))cand.push(i);
    }
    let square=false, couple=cand.length<=62, stamp=0;
    search: for(let i=0;i<used.length;i++)for(let j=0;j<i;j++){
      const p=used[i],q=used[j],r=p^q;
      if(!bad[r])continue;
      square=true;
      if(couple)break search;
      let n=0;stamp++;
      for(const t of cand){
        const v=Math.min(t,t^p,t^q,t^r);
        if(reps[v]!==stamp){reps[v]=stamp;n++;}
        if(n>62)break;
      }
      if(n<=62){couple=true;break search;}
    }
    ok(square,"2-dimensional conflict subgroup");
    if(!couple){
      exception++;
      ok(f.join(",")==="1,0,0,0","exception shape");
      for(const d of [da,db,dc])for(const u of all)if(!(u&40))
        ok(d[u]===(u&3?1:u?2:0),"exception is anisotropic on GHZ pair, zero on local slots");
    }
  }
}
ok(total===351232&&exception===64,"complete stabilizer replay");
const bits=(x,i)=>(x>>i)&1;
const qm=x=>(bits(x,0)&bits(x,1))^(bits(x,2)&bits(x,3))^Number(Boolean(x&48));
const tp=x=>x?(qm(x)?1:2):0;
for(let x=1;x<64;x++){
  const [a,b,c,d,e,f]=range(6).map(i=>bits(x,i));
  const y=x^((a^c)<<1)^(c<<3),z=(e<<1)|(f<<3)|(a<<4)|(c<<5);
  ok([tp(x),tp(y),tp(z)].sort().reduce((s,j)=>s+j,0)===4,"constructive binary min");
  ok((qm(x)^qm(y)^qm(z))===0,"constructive ternary compatibility");
}

console.log("PASS: E6 depth table, full and half-only duals, complete-line dual, 351232 stabilizer cases, and 192-phase witness");

