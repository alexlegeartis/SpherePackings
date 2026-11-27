const pwr=[1,3,9,27,81,243];
const digits=Array.from({length:729},(_,i)=>pwr.map(p=>Math.floor(i/p)%3));
const Fadd=new Uint16Array(729*729),Fmul=new Uint16Array(729*729),Fneg=new Uint16Array(729),inv=new Uint16Array(729);
const mod=[2,2,1,0,2,0,1],a6=mod.slice(0,6).map(a=>(3-a)%3);
for(let i=0;i<729;i++)Fneg[i]=digits[i].reduce((s,x,j)=>s+(x===0?0:3-x)*pwr[j],0);
function encode(ds){return ds.reduce((s,x,j)=>s+((x%3+3)%3)*pwr[j],0);}
function timesGen(x){const [a,b,c,d,e,f]=digits[x];return encode([f*a6[0],a+f*a6[1],b+f*a6[2],c+f*a6[3],d+f*a6[4],e+f*a6[5]]);}
for(let i=0;i<729;i++)for(let j=0;j<729;j++){let s=0;for(let k=0;k<6;k++)s+=((digits[i][k]+digits[j][k])%3)*pwr[k];Fadd[729*i+j]=s;}
const add=(x,y)=>Fadd[729*x+y],neg=x=>Fneg[x],sub=(x,y)=>add(x,neg(y));
for(let i=0;i<729;i++){
 const v=[i];for(let k=1;k<6;k++)v.push(timesGen(v[k-1]));
 const w=v.map(x=>neg(x));
 for(let j=0;j<729;j++){let s=0;for(let k=0;k<6;k++){const t=digits[j][k];if(t)s=add(s,t===1?v[k]:w[k]);}
 Fmul[729*i+j]=s;if(s===1)inv[i]=j;}
}
const mul=(x,y)=>Fmul[729*x+y];
function exp(x,k){let v=1;while(k){if(k&1)v=mul(v,x);x=mul(x,x);k>>=1;}return v;}
const cube=x=>mul(mul(x,x),x);
const field=Array.from({length:729},(_,i)=>i),rootsSq=Array.from({length:729},()=>[]),rootsAdd=Array.from({length:729},()=>[]);
for(const x of field){rootsSq[mul(x,x)].push(x);rootsAdd[add(cube(x),x)].push(x);}
function trim(x){while(x.length&&x.at(-1)===0)x.pop();return x;}
function polyAdd(a,b){const c=Array(Math.max(a.length,b.length));for(let i=0;i<c.length;i++)c[i]=add(a[i]??0,b[i]??0);return trim(c);}
function polySub(a,b){return polyAdd(a,b.map(neg));}
function polyMul(a,b){if(!a.length||!b.length)return [];const c=Array(a.length+b.length-1).fill(0);for(let i=0;i<a.length;i++)if(a[i])for(let j=0;j<b.length;j++)if(b[j])c[i+j]=add(c[i+j],mul(a[i],b[j]));return trim(c);}
function polyCube(a){if(!a.length)return [];const c=Array(3*(a.length-1)+1).fill(0);for(let i=0;i<a.length;i++)c[3*i]=cube(a[i]);return trim(c);}
function shift(a,t){let b=[];for(let i=a.length-1;i>=0;i--){b=polyMul(b,[t,1]);if(a[i])b[0]=add(b[0]??0,a[i]);}return b;}
const f=Array(29).fill(0);f[28]=1;
const omega=rootsSq[2][0];if(!omega||rootsAdd[0].length!==3||rootsAdd[0][1]!==omega&&rootsAdd[0][2]!==omega)throw Error("omega");
let verifyField=0;for(let a=1;a<729;a++)if(mul(a,inv[a])===1)verifyField++;if(verifyField!==728)throw Error("F727 fields");
function onCurve(x,y){return polySub(polySub(polySub(polyMul(y,y),polyCube(x)),x),f).length===0;}
function groups(){
 const gs=[];
 for(const a of field){
  const rr=rootsAdd[exp(a,28)];
  if(rr.length!==3)throw Error("translation constants");
  // a^27 t - a^81 t^3 + a^243 t^9
  const R=Array(10).fill(0);R[1]=exp(a,27);R[3]=neg(exp(a,81));R[9]=exp(a,243);
  R[0]=rr[0];
  if(polySub(polyAdd(polyCube(R),R),polySub(shift(f,a),f)).length)throw Error("correction");
  for(const r of rr){const S=R.slice();S[0]=r;gs.push({a,R:S});}
 }return gs;
}
function normalized(){
 const S=[],seen=new Set();
 const nine=field.filter(x=>exp(x,9)===x);if(nine.length!==9)throw Error("subfield");
 for(const u of field)for(const v of nine){
   const u3=exp(u,3),b=u?mul(exp(v,6),inv[u3]):0,c=add(exp(u,15),u);
   const wCube=u?rootsSq[sub(mul(b,b),c)].map(r=>sub(b,r)):[0];
   for(const y0 of wCube){
    if(u&&add(add(mul(y0,y0),mul(b,y0)),c)!==0)throw Error("quadratic");
    const w=u?exp(y0,243):0,z=u?neg(mul(mul(v,w),inv[u])):0;
    const y=Array(15).fill(0);y[14]=1;y[10]=u;y[7]=v;y[6]=mul(u,u);y[4]=w;y[3]=neg(mul(u,v));y[2]=neg(u3);y[1]=z;y[0]=sub(mul(v,v),mul(u,w));
    const a=polySub(polyMul(y,y),f), x=Array(10).fill(0);
    let good=true;
    for(let j=1;j<=27;j++){
      const c=sub(a[j]??0,j%3===0?cube(x[j/3]):0);
      if(j<=9)x[j]=c;else if(c){good=false;break;}
    }
    if(!good||!rootsAdd[a[0]??0].length)continue;
    if(seen.has(y.join(",")))throw Error("duplicate y");
    seen.add(y.join(","));x[0]=rootsAdd[a[0]??0][0];
    if(!onCurve(x,y))throw Error("normalized curve");
    S.push({x,y});
   }
 }
 if(S.length!==1143)throw Error("norm shell size "+S.length);return S;
}
function transform(v,g){
const x=polyAdd(shift(v.x,g.a),g.R),y=shift(v.y,g.a);return {x,y};
}
// Binary certificate: x(P)-x(V) divisor of y(P)-y(V) iff pairing -5.
// Scratch vectors reused; no floating point arithmetic.
function divides(w,v,ix,dx,dy){
for(let i=0;i<11;i++)dx[i]=sub(v.x[i]??0,w.x[i]??0);
for(let i=0;i<16;i++)dy[i]=sub(v.y[i]??0,w.y[i]??0);
for(let j=15;j>=10;j--){const c=mul(dy[j],ix);if(c)for(let k=0;k<=10;k++)dy[j-10+k]=sub(dy[j-10+k],mul(c,dx[k]));}
return dy.every(x=>x===0);
}
function probe(incoming,ids=[30,0]){
 const gs=groups(), seeds=normalized(), result=[];
 if(gs.length!==2187)throw Error("number of group words");
 for(const i of ids){
  const v=incoming.points[i];if(v.x.length!==11||v.y.length!==16||!onCurve(v.x,v.y))throw Error("bad candidate");
  let count=0,first=null;
  const dx=new Uint16Array(11),dy=new Uint16Array(16);
  const L=gs.map(g=>({g,v:transform(v,g)}));const ix=inv[v.x[10]];
  for(const [j,{g,v}] of L.entries())if(v.x.length!==11||!onCurve(v.x,v.y))throw Error("invalid transformed candidate "+j);
  for(const [j,{g,v}] of L.entries())for(const [k,w] of seeds.entries())if(divides(w,v,ix,dx,dy)){
   count++;if(first===null)first={j,k,normalSeed:w,group:g,transformedCandidate:v};}
  result.push({sourceIndex:i,candidate:v,tested:gs.length*seeds.length,exactMinusFive:count,first});
 }
 return result;
}

const checks=[{"sourceIndex":28,"candidate":{"x":[53,495,611,65,480,321,641,641,613,390,65],"y":[627,31,202,530,99,396,657,0,0,133,570,375,631,302,0,6]},"tested":2499741,"exactMinusFive":341,"first":{"j":2,"k":552,"normalSeed":{"x":[35,0,507,0,81,0,476,0,110,0],"y":[19,0,543,0,152,0,534,0,0,0,364,0,0,0,1]},"group":{"a":0,"R":[555,0,0,0,0,0,0,0,0,0]},"transformedCandidate":{"x":[491,495,611,65,480,321,641,641,613,390,65],"y":[627,31,202,530,99,396,657,0,0,133,570,375,631,302,0,6]}}},{"sourceIndex":30,"candidate":{"x":[201,571,80,684,81,247,693,693,506,507,65],"y":[608,430,227,486,37,396,646,557,0,60,430,375,353,302,0,6]},"tested":2499741,"exactMinusFive":258,"first":{"j":18,"k":283,"normalSeed":{"x":[137,0,461,0,649,0,135,0,459,0],"y":[541,0,541,0,270,0,270,0,0,0,191,0,0,0,1]},"group":{"a":6,"R":[152,329,0,702,0,0,0,0,0,418]},"transformedCandidate":{"x":[101,663,202,392,258,27,586,572,506,55,65],"y":[584,347,144,294,581,396,59,557,0,364,641,375,271,302,0,6]}}},{"sourceIndex":4,"candidate":{"x":[158,190,271,137,1,137,326,326,135,0,648],"y":[675,355,0,268,91,588,420,0,0,332,0,149,0,58,0,469]},"tested":2499741,"exactMinusFive":341,"first":{"j":13,"k":389,"normalSeed":{"x":[35,0,488,0,448,0,129,0,294,0],"y":[11,0,357,0,412,0,524,0,0,0,264,0,0,0,1]},"group":{"a":4,"R":[475,656,0,353,0,0,0,0,0,591]},"transformedCandidate":{"x":[68,0,627,242,386,603,309,110,135,167,648],"y":[430,0,604,214,633,588,139,0,0,591,505,149,567,58,0,469]}}},{"sourceIndex":0,"candidate":{"x":[0,0,1,0,0,0,0,0,1,0,1],"y":[0,1,0,0,0,2,0,2,0,0,0,2,0,1,0,2]},"tested":2499741,"exactMinusFive":260,"first":{"j":4,"k":422,"normalSeed":{"x":[0,0,1,293,1,1,557,291,293,0],"y":[0,2,293,556,556,0,555,291,0,0,292,0,0,0,1]},"group":{"a":1,"R":[293,1,0,2,0,0,0,0,0,1]},"transformedCandidate":{"x":[293,0,2,1,1,2,1,2,1,2,1],"y":[1,0,1,2,0,2,1,2,0,2,2,2,2,1,0,2]}}}];
const rootIds=[0,2187,4374,6561];
function require(x,msg){if(!x)throw Error(msg);}
function eq(a,b){return a.length===b.length&&a.every((x,i)=>x===b[i]);}
function verifyLiteral(){const gs=groups(),seen=new Set();
 for(let k=0;k<checks.length;k++){
  const C=checks[k],v=C.candidate,w=C.first.normalSeed,g=C.first.group,a=C.first.transformedCandidate;
  require(v.x.length===11&&v.y.length===16&&onCurve(v.x,v.y),'narrow root');
  require(w.x.length<=10&&w.y.length===15&&w.y[14]===1&&onCurve(w.x,w.y),'non-narrow Q contact');
  const G=gs[C.first.j];require(G.a===g.a&&eq(G.R,g.R),'translation match');
  const vv=transform(v,g);require(eq(vv.x,a.x)&&eq(vv.y,a.y)&&onCurve(a.x,a.y),'translation of root');
  const x0=polySub(shift(w.x,neg(g.a)),shift(g.R,neg(g.a))),y0=shift(w.y,neg(g.a));
  require(x0.length<=10&&y0.length===15&&y0[14]===1&&onCurve(x0,y0),'inverse Q witness');
  const dx=new Uint16Array(11),dy=new Uint16Array(16);
  require(divides(w,a,inv[a.x[10]],dx,dy)&&divides({x:x0,y:y0},v,inv[v.x[10]],dx,dy),'exact -5 contact');
  let n=0;for(const g of gs){const f=transform(v,g),key=f.x.join(',')+';'+f.y.join(',');require(!seen.has(key),'orbit repetition');seen.add(key);n++;}
  require(n===2187,'orbit size');
 }
 require(seen.size===8748,'total tested directions');
 require(14n*75n**2n>375n*14n**2n,'strict primary all-radius threshold');
 require(722n*2n*99n**53n<75n*15n*100n**53n,'lower-radius period bound');
 return {witnesses:checks.length,distinctBlockedDirections:seen.size};
}
const verified=verifyLiteral();
if(typeof process!=='undefined'&&process.argv.includes('--full')){
 const r=probe({points:checks.map(c=>c.candidate)},[0,1,2,3]);
 require(r.every((c,i)=>c.tested===checks[i].tested&&c.exactMinusFive===checks[i].exactMinusFive),'full orbit contact counts');
 verified.fullCounts=r.map(c=>c.exactMinusFive);
}
if(typeof console!=='undefined') console.log('PASS',JSON.stringify(verified));
