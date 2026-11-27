'use strict';
const assert=require('node:assert/strict');
const fs=require('node:fs'), path=require('node:path'), zlib=require('node:zlib');
const source=process.argv[2]||path.join(__dirname,'../p48p_k9_actual_39_20261126');
const matrix={"rs":[[-2,-2,-2,1,1,1,1,1,1],[1,1,1,-2,-2,-2,1,1,1],[3,0,0,-3,0,0,0,0,0],[-1,2,2,-1,-1,-1,2,-1,-1],[1,1,1,-2,-2,1,1,-2,1],[0,3,0,0,-3,0,0,0,0],[1,1,1,1,-2,-2,1,-2,1],[-1,-1,-1,2,-1,-1,2,2,-1]],"G":[[6,-3,-3,-3,-3,-3,-3,3,-2],[-3,6,3,3,3,3,3,0,1],[-3,3,6,0,3,0,0,-3,3],[-3,3,0,6,3,3,3,0,-1],[-3,3,3,3,6,3,3,-3,1],[-3,3,0,3,3,6,3,0,0],[-3,3,0,3,3,3,6,0,1],[3,0,-3,0,-3,0,0,6,-1],[-2,1,3,-1,1,0,1,-1,6]],"u":[[-2,-2,0,0,0,0,0,0],[2,0,-2,0,0,0,0,0],[2,0,0,-2,0,0,0,0],[2,0,0,2,0,0,0,0],[2,0,0,0,-2,0,0,0],[1,1,-1,1,-1,1,-1,-1],[1,1,-1,1,-1,-1,1,-1],[-1,-1,-1,1,1,-1,-1,-1]],"P":[[2,2,0,-4,0,-2,0,-2],[2,2,0,2,0,4,0,-2],[2,2,0,2,0,-2,0,4],[-4,2,0,2,0,-2,0,-2],[-1,-1,3,-1,3,1,3,1],[-1,-1,3,-1,-3,1,-3,1],[2,-4,0,2,0,-2,0,-2],[-1,-1,-3,-1,3,1,-3,1],[-1,-1,-3,-1,-3,1,3,1]]};
const ks=JSON.parse(fs.readFileSync(path.join(source,'embedding.json'))).K;
assert.deepEqual(ks,matrix.G);
const dot=(x,y)=>x.reduce((s,t,i)=>s+t*y[i],0);
const {rs,G,u,P}=matrix;
for(let i=0;i<8;i++){
 assert(u[i].every(x=>(x-u[i][0])%2===0) && u[i].reduce((a,b)=>a+b,0)%4===0);
 assert(rs[i].every(x=>(x-rs[i][0])%3===0) && rs[i].reduce((a,b)=>a+b,0)===0);
 assert.equal(rs[i][0],G[i][8]);
 for(let j=0;j<8;j++){assert.equal(dot(rs[i],rs[j]),3*G[i][j]);assert.equal(3*dot(u[i],u[j]),4*G[i][j]);}
 for(let j=0;j<9;j++)assert.equal(dot(P[j],u[i]),4*rs[i][j]);
}
function det(A){
 if(A.length===1)return A[0][0];
 return A[0].reduce((s,x,j)=>s+(j%2?-1:1)*x*det(A.slice(1).map(r=>r.filter((t,k)=>k!==j))),0);
}
assert.equal(det(G.slice(0,8).map(r=>r.slice(0,8).map(x=>x/3))),1);
assert.equal(det(G),21870);
assert.equal(3**8*10/3,21870);
assert(3**7*4*3>21870 && 3**6*8*4>21870 && 3**6*3*10**3===40500*2*3**3);
assert.equal(3**5*4*10**4,30000*4*3**4);


function atlas39(cache,G,gens){
 const H=new Int16Array(21870),w=new Int16Array(39);
 if(cache.length!==21869*39)throw 'bad source length';
 for(let i=1;i<H.length;i++){
   let k=i%10,y=(i/10)|0;
   for(let l=0;l<39;l++)w[l]=k*gens[7][l];
   for(let j=0;j<7;j++){let d=y%3;y=(y/3)|0;for(let l=0;l<39;l++)w[l]+=d*gens[j][l];}
   for(let l=0;l<39;l++) w[l]=(w[l]%30)+30*cache[(i-1)*39+l];
   let Q=0;
   for(let l=0;l<39;l++){let q=0;for(let m=0;m<39;m++)q+=G[l][m]*w[m];Q+=w[l]*q;}
   if(Q%30||Q<=0)throw 'bad witness';
   H[i]=Q/30;
 }
 if(Math.max(...H)!==171)throw 'wrong max';
 return H;
}
function phases5(){
 const C=[];
 for(let p=0;p<2;p++)for(let i=0;i<3**5;i++){let j=i,a=[];for(let k=0;k<5;k++){a.push(2*(j%3)+p);j=Math.floor(j/3)}for(let b=0;b<2;b++)C.push(a.map((v,k)=>v+(k===0?6*b:0)))}
 return C;
}
const mod=(x,m)=>((x%m)+m)%m;
function kp43(C,H,P){
 const out=[];
 const A=P.map(row=>row.slice(0,5)),T=P.map(row=>row.slice(5));
 const dot=(a,b)=>a.reduce((s,t,i)=>s+t*b[i],0);
 for(const [i,u] of C.entries()){
  let best=Infinity,example; const p=u[0]&1,sum=u.reduce((s,x)=>s+x,0),YY=A.map(a=>dot(a,u));
  for(let x=-8+p;x<=8;x+=2)for(let y=-8+p;y<=8;y+=2)for(let z=-8+p;z<=8;z+=2) {
    const k=5*(x*x+y*y+z*z);if(k>360||k>=best||(sum+x+y+z)%4!==0)continue;
    const XX=T.map(([a,b,c],i)=>YY[i]+a*x+b*y+c*z);
    if(!XX.every(t=>t%4===0))throw 'bad map';const jj=mod(-XX[0]/4,3);
    let idx=0;
    for(let i=7;i>0;i--){const v=XX[i]-XX[0];if(v%12!==0)throw 'bad class';idx=3*idx+mod(v/12,3);}
    idx*=10;
    for(let j=-13+mod(jj+13,3);j<=13;j+=3){const base=k+2*j*j;if(base>360||base>=best)continue;
      const id=idx+mod(j,10),h=2*H[id]+base;
      if(h<best){best=h;example=[x,y,z,j,id]}
    }
  }
  out.push({h:best,w:example});
 }
 return out;
}
function histA(U){let h={};for(let d of U)h[d.h]=(h[d.h]||0)+1;return h;}


function phasesA5(){
 const C=[];
 for(let u=0;u<3**5;u++){let z=u,a=[];for(let i=0;i<5;i++){a.push(z%3);z=(z/3)|0;}for(let b=0;b<6;b++)C.push(a.map((v,i)=>v+(i===0?3*b:0)))}
 return C;
}
function keyA5(a){
 let u=0,m=1,s=0;for(let x of a){let v=mod(x,3);u+=m*v;m*=3;s+=(x-v)/3;}return u*6+mod(s,6);
}
function kp43A5(C,H,P){
 const out=[];
 const A=P.map(row=>row.slice(0,5)),T=P.map(row=>row.slice(6)),Q=P.map(row=>row.slice(0,6).reduce((s,x)=>s+x,0));
 const dot=(a,b)=>a.reduce((s,t,i)=>s+t*b[i],0);
 for(const [i,u] of C.entries()){
  let best=Infinity,example; const sum=u.reduce((s,x)=>s+x,0),YY=A.map(a=>2*dot(a,u));
  for(let t=Math.ceil((-10-sum)/3);t<=Math.floor((10-sum)/3);t++)for(let x=-8+mod(t,2);x<=8;x+=2)for(let y=-8+mod(t,2);y<=8;y+=2) {
    const k=10*(sum+3*t)**2+15*(x*x+y*y);if(k>1080||k>=best||(2*sum+6*t+x+y)%4!==0)continue;
    const XX=T.map(([a,b],i)=>YY[i]+Q[i]*t+a*x+b*y);
    if(!XX.every(t=>t%4===0))throw 'bad map';const jj=mod(-XX[0]/4,3);
    let idx=0;
    for(let i=7;i>0;i--){const v=XX[i]-XX[0];if(v%12!==0)throw 'bad class';idx=3*idx+mod(v/12,3);}
    idx*=10;
    for(let j=-13+mod(jj+13,3);j<=13;j+=3){const base=k+6*j*j;if(base>1080||base>=best)continue;
      const id=idx+mod(j,10),h=6*H[id]+base;
      if(h<best){best=h;example=[t,x,y,j,id]}
    }
  }
  out.push({h:best,w:example});
 }
 return out;
}


function key5(d){
 let p=mod(d[0],2),j=0,m=1,t=0;
 for(let x of d){let r=mod(x,6);if((r&1)!==p)throw 'mixed parity';j+=m*((r-p)>>1);m*=3;t+=(x-r)/6;}
 return p*486+2*j+mod(t,2);
}


function normalSmall5(d){
 let par=0,s=0,m=0;for(let x of d){let y=Math.floor(x/6+0.5),z=x-6*y;s+=z*z;par^=y&1;m=Math.max(m,Math.abs(z));}
 return s+par*(36-12*m); // 12 times squared distance
}
function normalSmallA5(a){
 let s=a.reduce((t,x)=>t+x,0),z=[...a,0].map(x=>mod(6*x-s+9,18)-9);
 return z.reduce((t,x)=>t+x,0)===0?z.reduce((t,x)=>t+x*x,0):Infinity; // 108 q; thresholds for our 6-site witness
}


const N=[[14,7,7,8,6],[7,14,8,7,6],[7,8,14,7,6],[8,7,7,14,6],[6,6,6,6,12]];
const q=(x)=>N.reduce((z,R,i)=>z+x[i]*R.reduce((y,t,j)=>y+t*x[j],0),0);
const key=(x)=>x.map(t=>((t%48)+48)%48).reduce((z,t)=>48*z+t,0);
const v=(x)=>N.map(R=>R.reduce((a,t,i)=>a+t*x[i],0));
const minus=(x,y)=>x.map((t,i)=>t-y[i]);
const gcd=(a,b)=>b?gcd(b,a%b):a;
function balls(limit){
 const b=Math.floor(Math.sqrt(limit/8)),H=new Map();let xx=[0,0,0,0,0];
 function rec(i){if(i<5){for(let j=-b;j<=b;j++){xx[i]=j;rec(i+1)}return}
 const a=q(xx);if(a>limit)return;const p=v(xx), k=key(p);
 if(!H.has(k)||H.get(k).s>a)H.set(k,{s:a,p,x:xx.slice()});
 }
 rec(0);return H;
}



function imposs(A,t,bound=1000000) {
 const n=A.length, ids=new Map(),one=Array.from({length:n},(_,i)=>1n<<BigInt(i));
 one.forEach((bit,i)=>ids.set(bit,i));
 let visits=0;
 const walk=(P,k)=>{
   if(++visits>bound)throw Error('node cap');
   if(k===0)return true;
   let U=P,order=[],color=[],j,d=0;
   while(U){
      let Q=U;++d;
      while(Q){const b=Q&-Q;j=ids.get(b);U^=b;Q&=~(b|A[j]);order.push(j);color.push(d);}
   }
   for(j=order.length-1;j>=0;j--){
      if(color[j]<k)return false;
      const i=order[j];if(walk(P&A[i],k-1))return true;P^=one[i];
   }
   return false;
 };
 return {found:walk((1n<<BigInt(n))-1n,t),visits};
}
function graph(X,f,s) {
 const A=X.map(()=>0n);
 for(let i=0;i<X.length;i++)for(let j=i+1;j<X.length;j++) if(f(X[i],X[j])<=s){A[i]|=1n<<BigInt(j);A[j]|=1n<<BigInt(i);}
 return A;
}
function target(den,D,s) {
 // delta^2 = M^2 (6den-s)^43 / ((4den)^43 D); comparator^2=137^43/(2^223 3^42).
 const a=BigInt(6*den-s)**43n*2n**223n*3n**42n;
 const b=137n**43n*BigInt(4*den)**43n*BigInt(D);
 let lo=1,hi=D+1;
 while(lo<hi){const m=(lo+hi)>>1;if(BigInt(m)**2n*a>b)hi=m;else lo=m+1;}
 return lo;
}
function run43(){
 let H=balls(85),hist=[...new Set([...H.values()].map(x=>x.s))].sort((a,b)=>a-b);
 if(target(48,3072,86)!==3073)throw Error('tail');
 return hist.map(s=>{
   const t=target(48,3072,s),X=[...H.values()].filter(x=>x.s>0&&x.s<=s);
   const f=(a,b)=>H.get(key(minus(a.p,b.p)))?.s??Infinity;
   return [s,t,X.length,t>X.length+1?{}:imposs(graph(X,f,s),t-1)];
 })
}


function verify(C,A,key,unit,D,m,s,witness,normalCheck,normalCut,expected) {
 assert(C.every((u,i)=>key(u)===i) && C.length===D);
 assert.equal(A.length,D); assert.equal(A[0].h,0);
 assert(A.slice(1).every(x=>0<x.h&&x.h<6*unit));
 assert.deepEqual(histA(A),expected);
 const H=A.map(x=>x.h);
 const dist=(u,v)=>H[key(u.map((x,i)=>x-v[i]))];
 for(let i=0;i<D;i++)assert.equal(dist(C[0],C[i]),H[i]);
 let n=0;
 for(let i=0;i<witness.length;i++)for(let j=0;j<i;j++){
   const a=C[witness[i]],b=C[witness[j]];
   assert(normalCheck(a.map((x,k)=>x-b[k]))<=normalCut);
   assert(dist(a,b)>=s);
   if(dist(a,b)===s)n++;
 }
 assert(n>0);
 const coeff=BigInt(m)**2n*BigInt(s)**43n;
 assert(coeff<BigInt(6*unit)**43n);
 // Both singletons (radius sqrt(6)/2) fall below DSW's record.
 assert(6n**43n*2n**223n*3n**42n < 137n**43n*4n**43n*BigInt(D));
 const vals=[...new Set(H.slice(1))].sort((a,b)=>b-a),report=[];
 for(let h of vals){
  const f=(a,b)=>-dist(a,b),X=C.filter(u=>dist(u,C[0])>=h),Y=BigInt(h)**43n;
  let lo=2,hi=D+1;
  while(lo<hi){const mid=(lo+hi)>>1;if(BigInt(mid)**2n*Y>coeff)hi=mid;else lo=mid+1;}
  if(lo>X.length+1)continue;
  const q=imposs(graph(X,f,-h),lo-1);
  assert(!q.found);
  report.push([h,lo,X.length,q.visits]);
 }
 console.log('det='+D,'histogram',histA(A));
 console.log('excluded [height numerator, forbidden motif size, induced graph size, search nodes]',JSON.stringify(report));
 console.log('non-singleton witness:',witness,'squared min',s+'/'+unit);
}
function run() {
 const gram=JSON.parse(fs.readFileSync(path.join(source,'gram.json')));
 const gens=JSON.parse(fs.readFileSync(path.join(source,'gens.json')));
 const c=zlib.gunzipSync(fs.readFileSync(path.join(source,'witness.gz')));
 const H=atlas39(new Int8Array(c.buffer,c.byteOffset,c.byteLength),gram,gens);
 assert(gram.some((r,i)=>r[i]===6));
 verify(phases5(),kp43(phases5(),H,P),key5,60,972,4,335,[0,971,966,325],normalSmall5,5,{"0":1,"135":2,"180":1,"200":32,"215":112,"220":80,"240":80,"255":160,"260":72,"280":90,"295":180,"300":80,"320":40,"335":32,"340":10});
 verify(phasesA5(),kp43A5(phasesA5(),H,P),keyA5,180,1458,6,990,[0,1453,1442,1407,1300,977],normalSmallA5,54,{"0":1,"270":1,"360":2,"600":60,"630":42,"640":132,"670":150,"720":60,"750":90,"760":162,"790":162,"840":90,"870":90,"880":162,"910":162,"960":30,"990":20,"1000":30,"1030":12});
 // This separate calculation uses only the normal-distance guarantee in DSW's determinant-3072 carrier.
 const normal=run43();
 assert(normal.length===23 && normal.every(r=>!r[3].found));
 console.log('det=3072 wrapped guarantee cannot beat DSW at any radius (not an actual-metric claim)');
 console.log('PASS');
}
if(require.main===module)run();
module.exports={atlas39,P,source,imposs,graph,mod,histA,normalSmall5};

