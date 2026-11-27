'use strict';
const P=729,N=6,m=[2,2,1,0,2,0,1],pow3=[1,3,9,27,81,243];const digits=Array.from({length:P},(_,n)=>pow3.map(k=>(n/k|0)%3));const A=new Uint16Array(P*P),M=new Uint16Array(P*P),S=new Uint16Array(P*P),I=new Uint16Array(P);
for(let a=0;a<P;a++)for(let b=0;b<P;b++){let val=0,sub=0;let x=digits[a],y=digits[b];for(let k=0;k<N;k++){val+=(x[k]+y[k])%3*pow3[k];sub+=(x[k]+3-y[k])%3*pow3[k]}A[a*P+b]=val;S[a*P+b]=sub;}
for(let a=1;a<P;a++)for(let b=1;b<=a;b++){const v=Array(11).fill(0),x=digits[a],y=digits[b];for(let i=0;i<N;i++)if(x[i])for(let j=0;j<N;j++)v[i+j]+=x[i]*y[j];for(let i=10;i>=6;i--){const z=v[i]%3;if(z)for(let k=0;k<6;k++)if(m[k])v[i-6+k]+=(3-z)*m[k];}let val=0;for(let k=0;k<N;k++)val+=v[k]%3*pow3[k];M[a*P+b]=M[b*P+a]=val;if(val===1){I[a]=b;I[b]=a}}
const ag=(a,b)=>A[a*P+b],sg=(a,b)=>S[a*P+b],mg=(a,b)=>M[a*P+b],ng=a=>sg(0,a),ig=a=>I[a];function pg(a,n){let z=1;for(;n;n>>=1,a=mg(a,a))if(n&1)z=mg(z,a);return z;}
const trim=a=>{while(a.length&&!a.at(-1))a.pop();return a};const add=(p,q)=>trim(Array.from({length:Math.max(p.length,q.length)},(_,i)=>ag(p[i]||0,q[i]||0))),sub=(p,q)=>trim(Array.from({length:Math.max(p.length,q.length)},(_,i)=>sg(p[i]||0,q[i]||0))),scale=(p,s)=>s?p.map(a=>mg(a,s)):[];function mul(p,q){if(!p.length||!q.length)return [];let a=new Array(p.length+q.length-1).fill(0);for(let i=0;i<p.length;i++)if(p[i])for(let j=0;j<q.length;j++)if(q[j])a[i+j]=ag(a[i+j],mg(p[i],q[j]));return trim(a);}
function div(p,q){if(!q.length)throw 'div0';const a=p.slice(),g=new Array(Math.max(0,p.length-q.length+1)).fill(0),qi=ig(q.at(-1));for(let i=p.length-q.length;i>=0;i--)if(a[i+q.length-1]){let z=mg(a[i+q.length-1],qi);g[i]=z;for(let j=0;j<q.length;j++)a[i+j]=sg(a[i+j],mg(z,q[j]));}return [trim(g),trim(a)];}
const poly=(i,j)=>Array.from({length:i+1},(_,k)=>i===k?j:0); const check=(x,y)=>sub(mul(y,y),add(add(mul(mul(x,x),x),x),poly(28,1))).length===0;
function sum(a,b){let d=sub(a.x,b.x),e=sub(a.y,b.y);if(!d.length)return null;let [k,rem]=div(e,d);if(rem.length)return null;let x=sub(mul(k,k),add(a.x,b.x));let y=sub(mul(k,sub(a.x,x)),a.y);if(!check(x,y))throw 'bad point';return {x,y};}

const elems=Array.from({length:P},(_,i)=>i),F9=elems.filter(a=>pg(a,9)===a),cubeRoots=Array.from({length:P},()=>[]),root2=Array.from({length:P},()=>[]);
for(let a of elems){cubeRoots[ag(pg(a,3),a)].push(a);root2[mg(a,a)].push(a)}
function normalized(){const out=[];for(let u of elems)for(let v of F9){let u2=mg(u,u),u3=mg(u2,u);let wvals;if(!u)wvals=[0];else {let k=mg(pg(v,6),ig(u3)),b=ag(pg(u,15),u),disc=sg(mg(k,k),b);wvals=root2[disc].map(a=>pg(sg(k,a),243));}for(let w of wvals){let z=u?ng(mg(mg(v,w),ig(u))):0,yy=Array(15).fill(0);yy[14]=1;yy[10]=u;yy[7]=v;yy[6]=u2;yy[4]=w;yy[3]=ng(mg(u,v));yy[2]=ng(u3);yy[1]=z;yy[0]=sg(mg(v,v),mg(u,w));let T=sub(mul(yy,yy),poly(28,1)),x=Array(10).fill(0),ok=true;for(let j=1;j<=27;j++){let q=T[j]||0;if(j%3===0)q=sg(q,pg(x[j/3],3));if(j<10)x[j]=q;else if(q){ok=false;break;}}const cc=cubeRoots[T[0]||0];if(ok&&cc.length){x[0]=cc[0];if(!check(x,yy))throw 'norm';out.push({x,y:yy,u,v,w,z});}}}if(out.length!==1143)throw ["count",out.length];return out;}
const bin=Array.from({length:29},(_,i)=>Array(i+1).fill(1));for(let i=2;i<29;i++)for(let j=1;j<i;j++)bin[i][j]=(bin[i-1][j-1]+bin[i-1][j])%3;
const COEFF=elems.map(a=>Array.from({length:29},(_,i)=>pg(a,i)));
function trans(p,a){const res=Array(p.length).fill(0),pp=COEFF[a];for(let i=0;i<p.length;i++)if(p[i])for(let j=0;j<=i;j++)if(bin[i][j])res[j]=ag(res[j],mg(mg(p[i],bin[i][j]),pp[i-j]));return res}
function trpoint(v,a){let x=trans(v.x,a),y=trans(v.y,a);let cs=cubeRoots[pg(a,28)];if(cs.length!==3)throw 'const';x[1]=ag(x[1],pg(a,27));x[3]=sg(x[3],pg(a,81));x[9]=ag(x[9],pg(a,243));return {x,y,cs}}
function remainderDivides(d,e){const a=e.slice(),qi=ig(d.at(-1));for(let i=e.length-d.length;i>=0;i--)if(a[i+d.length-1]){let z=mg(a[i+d.length-1],qi);for(let j=0;j<d.length-1;j++)a[i+j]=sg(a[i+j],mg(z,d[j]));}for(let j=0;j<d.length-1;j++)if(a[j])return false;return true;}

function gcd(a,b){while(b.length){const z=div(a,b)[1];a=b;b=z;}return scale(a,ig(a.at(-1)));}
// Raw group law, with a monic reduced denominator B.
function rational(a,b){let d=sub(a.x,b.x),e=sub(a.y,b.y);if(!d.length)return null;let g=gcd(d,e);let [B]=div(d,g),[k]=div(e,g);if(B.length!==2)return null;const sc=ig(B.at(-1));B=scale(B,sc);k=scale(k,sc);let B2=mul(B,B),B3=mul(B2,B),X=sub(mul(k,k),mul(add(a.x,b.x),B2));if(X.length!==13)return null;let Y=sub(mul(k,sub(mul(a.x,B2),X)),mul(a.y,B3));if(Y.length!==19)throw ["degree",Y.length];return {B,X,Y};}
function verify(inp,data){
const assert=(ok,why)=>{if(!ok)throw Error(JSON.stringify(why));};
assert(inp.points.length===124&&String(inp.field_modulus)===String(m),"source");
const W=normalized(),C=[new Map(),new Map()],counts=[0,0],ps=inp.points;
for(const p of ps)assert(p.x.length===11&&p.y.length===16&&check(p.x,p.y),"source point");
for(let i=0;i<ps.length;i++)for(let j=0;j<i;j++)for(let s of [1,2]){
 let l=ps[i], r={x:ps[j].x,y:scale(ps[j].y,s)}, p=sum(l,r), q=rational(l,r);
 if(p&&p.x.length===13){C[0].set([i,j,s].join(),p);counts[0]++;}
 if(q){C[1].set([i,j,s].join(),q);counts[1]++;}
}
assert(String(counts)==="24,66","raw source coverage");
let unique=[...C].map((m,k)=>new Set([...m.values()].map(x=>k?JSON.stringify([x.B,x.X,x.Y]):JSON.stringify([x.x,x.y]))));
assert(unique[0].size===21&&unique[1].size===56,"distinct directions");
for(let k of [0,1]){
 const d=data[k?"rational":"polynomial"],seen=new Set();
 for(const c of d){
  let p=C[k].get([c.i,c.j,c.s].join());
  assert(p,"pair sum replay");
  let key=k?JSON.stringify([p.B,p.X,p.Y]):JSON.stringify([p.x,p.y]);assert(!seen.has(key),"repeat");seen.add(key);
  const h=c.hit,w=W[h.j];
  assert(w&&[w.u,w.v,w.w,w.z].join()===[h.u,h.v,h.w,h.z].join()&&[1,2].includes(h.t),"normalized witness");
  const B=k?trans(p.B,h.a):[1],B2=mul(B,B),B3=mul(B2,B);
  const corrections=Array(10).fill(0);corrections[0]=h.c;corrections[1]=pg(h.a,27);corrections[3]=ng(pg(h.a,81));corrections[9]=pg(h.a,243);
  const t28=poly(28,1);
  assert(cubeRoots[pg(h.a,28)].includes(h.c),"translation constant");
  const correctionCube=poly(0,pg(corrections[0],3));
  for(let j=1;j<corrections.length;j++)correctionCube[j*3]=pg(corrections[j],3);
  for(let j=0;j<correctionCube.length;j++)correctionCube[j]??=0;
  assert(sub(add(correctionCube,corrections),sub(trans(t28,h.a),t28)).length===0,"isometry correction");
  const X=add(trans(k?p.X:p.x,h.a),mul(B2,corrections)),Y=scale(trans(k?p.Y:p.y,h.a),h.t);
  assert(X.length===13&&Y.length===19,"direction degrees");
  const rhs=add(add(mul(mul(X,X),X),mul(X,mul(B2,B2))),mul(t28,mul(B3,B3)));
  assert(sub(mul(Y,Y),rhs).length===0,"translated point equation");
  const D=sub(X,mul(w.x,B2)),E=sub(Y,mul(w.y,B3));
  assert(D.length===13&&div(E,D)[1].length===0,"exact -6 pairing");
  if(k)assert(div(X,B)[1].length>0,"reduced pole");
 }
 assert(seen.size===unique[k].size&&[...unique[k]].every(x=>seen.has(x)),"witness coverage");
}
return {sourcePairCandidates:counts,distinctWithExactBlockers:[21,56],normalizedClosestSources:W.length};
}

if(typeof module!=='undefined') module.exports=verify;
if(typeof require!=='undefined' && require.main===module){const fs=require('node:fs'),path=require('node:path'); console.log(verify(JSON.parse(fs.readFileSync(process.argv[2],'utf8')),JSON.parse(fs.readFileSync(path.join(__dirname,'witnesses.json'),'utf8'))));}
