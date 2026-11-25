'use strict';
// Exact ternary linear-code verification. All arithmetic is on small integers.

function W3(n){return Array.from({length:3**n},(_,i)=>Array.from({length:n},(_,j)=>Math.floor(i/3**j)%3));}
const chi=(u,v)=>u.reduce((s,x,i)=>s+x*v[i],0)%3? -1:2;
function groupBy(ws,fn){const A=new Map();for(const w of ws){const k=fn(w);if(!A.has(k))A.set(k,[]);A.get(k).push(w);}return [...A.values()];}

const m3=x=>(x%3+3)%3;
const prj=(r)=>r.every(x=>x===0)||r.find(Boolean)===1;
const id3=x=>x.reduce((s,y,i)=>s+y*3**i,0);
function det3(A){
 A=A.map(a=>a.slice());let d=1;
 for(let i=0;i<4;i++){let p=i;while(p<4&&!A[p][i])p++;if(p===4)return 0;
  if(p!==i){[A[i],A[p]]=[A[p],A[i]];d=-d;}
  d=m3(d*A[i][i]);const u=A[i][i];
  for(let j=i+1;j<4;j++){const c=m3(u*A[j][i]);for(let k=i;k<4;k++)A[j][k]=m3(A[j][k]-c*A[i][k]);}
 }return m3(d);
}
const U4=W3(4).filter(prj),pos4=new Map();
for(let i=0;i<U4.length;i++){pos4.set(id3(U4[i]),i);pos4.set(id3(U4[i].map(x=>m3(-x))),i);}
const qs=[];
for(const e of W3(10)){
 let l=0,A=Array.from({length:4},()=>Array(4).fill(0));
 for(let i=0;i<4;i++)for(let j=i;j<4;j++)A[i][j]=A[j][i]=e[l++];
 if(det3(A)!==1)continue;
 qs.push(U4.map((u,i)=>i===0?0:(u.reduce((s,x,i)=>s+x*u.reduce((t,y,j)=>t+A[i][j]*y,0),0)%3||3)));
}
const V=W3(6).filter(prj).slice(1);
function formsAt(t){
 const fs=[([a,b,c,d,e,f])=>[c,d,e,f],([a,b,c,d,e,f])=>[a,b,e,f],
 [([a,b,c,d,e,f])=>[a,b,c,d],([a,b,c,d,e,f])=>[m3(a-c),b,d,f],([a,b,c,d,e,f])=>[m3(a-c),m3(b-d),e,f]][t]];
 const map=fs.map(f=>V.map(v=>pos4.get(id3(f(v)))));const K=map.map(u=>u.flatMap((a,i)=>a===0?[i]:[]));
 const F=map.map((u,i)=>qs.flatMap((q,l)=>{
 let iso=0;
 for(let j=0;j<3;j++)if(j!==i){const d=K[j].map(k=>q[u[k]]);if(d.some(x=>x<=1))return [];
  if(d.every(x=>x===3))iso|=1<<j;
 }const depths=Uint8Array.from(u.map(j=>q[j]));
 return [{l,d:depths,iso}];
 }));
 return {F,K,map};
}
function rred(t,shared){
 const {F,K,map}=shared, qKeys=F[0].map(f=>f.l), qIdx=new Map(qKeys.map((l,i)=>[qs[l].join(""),i]));
 const vmap=Array(41); for(let i=0;i<map[0].length;i++)vmap[map[0][i]]=i;
 const fs=v=>[v[2],v[3],v[4],v[5]], W6=W3(6).slice(1);
 const ord=V.map((v,i)=>[i,v.reduce((s,y)=>s+(y?1:0),0)]).sort((a,b)=>a[1]-b[1]).map(x=>x[0]);
 W6.sort((x,y)=>x.filter(Boolean).length-y.filter(Boolean).length);
 const dot=(a,b)=>m3(a.reduce((s,x,i)=>s+x*b[i],0));
 const masks=Array(V.length).fill(0);for(let j=0;j<3;j++)for(const k of K[j])masks[k]|=1<<j;
 const H=new Map();
 for(const vi of ord){const u=V[vi];if(masks[vi]&1)continue;for(const f of W6){
  if(dot(u,f)===2)continue;
  if(K.some((ks,j)=>!(masks[vi]>>j&1)&&ks.some(k=>dot(V[k],f))))continue;
  const a=U4.map((_,i)=>{if(i===0)return 0;const x=V[vmap[i]],b=dot(x,f);
   return pos4.get(id3(fs(x.map((z,j)=>m3(z+b*u[j])))));});
  const key=a.join(",");
  if(!H.has(key))H.set(key,a);
 }}
 const hs=[...H.values()].filter(p=>p.some((x,i)=>x!==i));
 const gen=hs.slice(0,100);
 const visited=Array(F[0].length).fill(false),reps=[];
 for(let i=0;i<F[0].length;i++){if(visited[i])continue;
  const list=[i];visited[i]=true;
  for(let j=0;j<list.length;j++){
   const Q=qs[qKeys[list[j]]];
   for(const h of gen){const k=qIdx.get(h.map(i=>Q[i]).join(""));
    if(k===undefined)throw Error("orbit step escaped");
    if(!visited[k]){visited[k]=true;list.push(k);}
   }
  }reps.push([i,list.length]);}
 return {reps,gen,ng:hs.length};
}
function searchAdded(t,S,R){
 const {F,K}=S,{reps}=R, out={pair:0,triple:0,byRep:[]};const Vn=V.length;
 const Fs=F[2];const indexes=Array.from({length:Vn},(_,p)=>Fs.flatMap((x,i)=>x.d[p]===3?[i]:[]));
 const indexes2=Array.from({length:Vn},(_,p)=>Fs.flatMap((x,i)=>x.d[p]>=2?[i]:[]));
 const dim3=t===2;const testA=[];
 for(const [i,on] of reps){const a=F[0][i];let here=0,tries=0;
  for(const b of F[1]){
   if(K[2].some(p=>a.d[p]+b.d[p]<5))continue;
   let best=null, wrong=false;
   const req2=[],req3=[];
   for(let j=0;j<Vn;j++){const s=a.d[j]+b.d[j];
    if(s<=1){wrong=true;break;}
    if(s===2){req3.push(j);if(best===null||indexes[j].length<best.length)best=indexes[j];}
    if(s===3){req2.push(j);if(best===null||indexes2[j].length<best.length)best=indexes2[j];}
   }
   if(wrong)continue;
   out.pair++;here++;
   for(const k of best??Fs.map((_,i)=>i)){
    const c=Fs[k];out.triple++;tries++;
    if(req3.some(j=>c.d[j]<3)||req2.some(j=>c.d[j]<2))continue;
    if(Array.from({length:Vn},(_,j)=>j).every(j=>a.d[j]+b.d[j]+c.d[j]>=5))return {...out,hit:[a.l,b.l,c.l],t};
   }
  }
  out.byRep.push([a.l,on,here,tries]);
 }
 return out;
}
if(qs.length!==21060||V.length!==364)throw Error("quadratic or projective census");
const expected=[{"orbits":[48,192,192,192],"pair":1728,"triple":393984},{"orbits":[108,486,108,108],"pair":2052,"triple":554040},{"orbits":[1296,648,648,648,648],"pair":9072,"triple":13716864}];
for(let t=0;t<3;t++){
 const S=formsAt(t),R=rred(t,S),o=searchAdded(t,S,R);
 const s=R.reps.map(x=>x[1]);
 if(o.hit||JSON.stringify(s)!==JSON.stringify(expected[t].orbits)||o.pair!==expected[t].pair||o.triple!==expected[t].triple)
   throw Error("unexpected exact count, or feasible triple");
 if(s.reduce((a,b)=>a+b)!==S.F[0].length||S.K.some(x=>x.length!==4))throw Error("orbit coverage");
 console.log(JSON.stringify({kernelIntersection:t,localFormCounts:S.F.map(x=>x.length),firstFormOrbitSizes:s,pair:o.pair,triple:o.triple,hit:null}));
}
