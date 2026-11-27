
"use strict";
const H=12; // up to q^6 for Jacobi theta, q^(12/3) at cusp zero
function ck(x,m){if(!x)throw Error(m);}
const S=f=>Array.from({length:H+1},(_,n)=>f?.(n)??new Map());
const put=(m,y,x)=>{if(!x)return;const t=(m.get(y)??0n)+x;if(t)m.set(y,t);else m.delete(y);};
const M=(q=0,y=0,x=1n)=>S(n=>n===q?new Map([[y,x]]):new Map());
const add=(a,b)=>S(n=>{const t=new Map(a[n]);for(let [y,x] of b[n])put(t,y,x);return t;});
const mult=(a,b)=>S(n=>{const t=new Map;for(let j=0;j<=n;j++)for(let [u,x] of a[j])for(let [v,y] of b[n-j])put(t,u+v,x*y);return t;});
const scale=(a,x)=>S(n=>new Map([...a[n]].map(([y,v])=>[y,v*x])));
function pow(a,e){let r=M();for(let j=0;j<e;j++)r=mult(r,a);return r;}
const scalar=a=>a.map(u=>[...u.values()].reduce((r,x)=>r+x,0n));
const term=(a,n,y=0)=>a[n]?.get(y)??0n;
function div(a,d){
 const r=S();for(let n=0;n<=H;n++){
  r[n]=new Map(a[n]);for(let j=1;j<=n;j++)for(let [u,x] of r[n-j])put(r[n],u,-d[j]*x);
  for(let [u,x] of r[n]){ck(x%d[0]===0n,"nonintegral coefficient");r[n].set(u,x/d[0]);}
 }return r;
}
function prod(step,p){let r=M();for(let j=step;j<=H;j+=step)if(p)r=mult(r,pow(add(M(),M(j,0,-1n)),p));return r;}
const mk=(n,num,den,d)=>div(prod(n,num),scalar(prod(d,den)));
const U=r=>mk(1,3*r,r,3), B=U(3), C=mult(M(1),mk(3,9,3,1)), Uzero=r=>mk(3,3*r,r,1);
const up3=a=>S(n=>n%3===0?new Map(a[n/3]):new Map());
const A=S();
for(let i=-8;i<=8;i++)for(let j=-8;j<=8;j++){let n=i*i+i*j+j*j;if(n<=H)put(A[n],0,1n);}
let X=add(add(M(0,-1),M(0,1)),M(0,0,-2n));
for(let q=1;q<=H;q++)X=mult(X,mult(pow(add(M(),M(q,-1,-1n)),2),pow(add(M(),M(q,1,-1n)),2)));
X=div(X,scalar(prod(1,4)));
let T2=S(),T3=S(),T4=S();
for(let i=-8;i<=8;i++)for(let j=-8;j<=8;j++){
 let n=i*(i+1)/2+j*(j+1)/2;if(n<=H)put(T2[n],i+j+1,1n);
 n=i*i+j*j;if(n<=H){put(T3[n],i+j,1n);put(T4[n],i+j,(i+j)&1?-1n:1n);}
}
T2=div(scale(T2,4n),scalar(T2));T3=div(T3,scalar(T3));T4=div(T4,scalar(T4));
const T34=add(T3,T4),Y=add(T2,scale(S(n=>2*n<=H?new Map(T34[2*n]):new Map()),4n));
ck(term(X,0)===-2n&&term(X,0,1)===1n&&term(X,1,2)===-2n,"Jacobi generator X");
ck(term(Y,0)===10n&&term(Y,0,1)===1n&&term(Y,1,2)===10n&&term(Y,1,1)===-64n&&term(Y,1)===108n,"Jacobi generator Y");
function cols(w,r) {
 const m=w/2, cols=[];
 for(let j=0;j<=m;j++){
  const K=24-r+2*j,d=Math.floor(K/3),s=K-3*d;
  const P=mult(pow(X,j),pow(Y,m-j)), p=mult(U(r),P), z=mult(Uzero(r),up3(P));
  for(let i=0;i<=d;i++){
   const F=mult(pow(A,s),mult(pow(B,i),pow(C,d-i)));
   const S=mult(pow(A,s),mult(pow(C,i),pow(B,d-i)));
   const e=(r-s+3*i-9*(d-i))/2;ck(Number.isInteger(e),"three exponent");
   cols.push({j,i,e,inf:mult(p,F),zero:mult(z,S)});
  }
 }
 return cols;
}

const Z=[0n,1n], O=[1n,1n];
const abs=x=>x<0n?-x:x;
function gcd(a,b){a=abs(a);b=abs(b);while(b)[a,b]=[b,a%b];return a;}
function fq(n,d=1n){if(!d)throw "zero denom";if(!n)return Z;if(d<0n){n=-n;d=-d;}const z=gcd(n,d);return [n/z,d/z];}
function ad(a,b){if(!a[0])return b;if(!b[0])return a;return fq(a[0]*b[1]+b[0]*a[1],a[1]*b[1]);}
function ng(a){return [-a[0],a[1]];}
function sb(a,b){return ad(a,ng(b));}
function ml(a,b){if(!a[0]||!b[0])return Z;return fq(a[0]*b[0],a[1]*b[1]);}
function dv(a,b){return fq(a[0]*b[1],a[1]*b[0]);}
function sgn(a){return Math.sign(Number(a[0]));}
function pc3(p){if(p>=0)return [3n**BigInt(p),1n];return [1n,3n**BigInt(-p)];}
function iF(n){return [BigInt(n),1n];}
function rr(A,n){
 const r=A.map(x=>x.slice());let p=0,piv=[];
 for(let j=0;j<n;j++){
  let k=p;while(k<r.length&&!r[k][j][0])k++;
  if(k===r.length)continue;
  [r[p],r[k]]=[r[k],r[p]];
  const d=r[p][j];for(let l=j;l<r[p].length;l++)r[p][l]=dv(r[p][l],d);
  for(let i=0;i<r.length;i++)if(i!==p&&r[i][j][0]){const c=r[i][j];for(let l=j;l<r[i].length;l++)r[i][l]=sb(r[i][l],ml(c,r[p][l]));}
  piv.push(j);p++;
 }
 return {r,p,piv};
}

function matrix(w,r,lim,gate,sel,infq,maxy){
 const C=cols(w,r),n=C.length;
 const f=(c,q,y)=>C.map(v=>ml(iF(term(v[c],q,y)),c==="inf"?pc3(-v.e):v.j%2?ng(O):O));
 let E=[],I=[];
 const eq=(c,q,y,rhs)=>{const a=f(c,q,y);if(a.some(x=>x[0])||rhs)E.push(a.concat([iF(rhs)]));};
 const nonnegative=(q,y)=>{const a=f("zero",q,y).map(ng);if(a.some(x=>x[0]))I.push({tag:[q,y],row:a.concat([Z])});};
 for(let q=0;q<=infq;q++)for(let y=0;y<=maxy;y++){
  if(q<3)eq("inf",q,y,q===0&&y===0?1:0);
  else if(2*y>w+2*q-6)eq("inf",q,y,2*q===w&&y===w?1:0);
 }
 for(let q=0;q<=lim;q++)for(let y=0;y<=maxy;y++){
  if(9*y*y>2*w*(r+3*q) || y===0&&2*(r+3*q)<gate)eq("zero",q,y,0);
  else nonnegative(q,y);
 }
 const R=rr(E,n),free=Array.from({length:n},(_,i)=>i).filter(j=>!R.piv.includes(j));
 ck(R.r.slice(R.p).every(x=>!x[n][0]),"unexpected equality contradiction");
 const Is=sel.map(tag=>{const q=I.find(x=>x.tag[0]===tag[0]&&x.tag[1]===tag[1]);ck(q,'missing nonnegative tag');return q;});
 const b=Is.map(({row})=>{
  let h=row.slice();for(let i=0;i<R.p;i++){
   const j=R.piv[i],a=h[j];if(a[0])for(let k=j;k<=n;k++)h[k]=sb(h[k],ml(a,R.r[i][k]));
  }
  return free.map(j=>h[j]).concat([h[n]]);
 });
 return {b,I,E,R,free,Is};
}
function certificate(w,r,lim,gate,sel,infq,maxy){
 const {b,free,R,Is}=matrix(w,r,lim,gate,sel,infq,maxy);
 const P=free.map((_,i)=>b.map(u=>u[i]));
 const J=rr(P,sel.length),vf=Array.from({length:sel.length},(_,i)=>i).filter(j=>!J.piv.includes(j));
 if(vf.length!==1)return {vf:vf.length,rank:R.p};
 const y=Array(sel.length).fill(Z);y[vf[0]]=O;for(let j=0;j<J.p;j++)y[J.piv[j]]=ng(J.r[j][vf[0]]);
 ck(P.every(row=>row.reduce((a,v,j)=>ad(a,ml(v,y[j])),Z)[0]===0n),"Farkas combination does not vanish");
 const obj=y.reduce((v,c,j)=>ad(v,ml(c,b[j][free.length])),Z);
 return {vf:1,rank:R.p,pos:y.every(x=>x[0]>=0n),obj:obj.map(String),labels:Is.map(x=>x.tag)};
}

const cases = [
  [6,0,9,52,[[1,1],[1,2],[2,2],[4,2],[4,3],[5,1],[5,2],[5,3],[6,1],[6,2],[6,3],[7,2],[9,1],[9,3]],4,10,"-712648611180047162763","5061498431012864"],
  [6,1,9,50,[[1,1],[1,2],[2,1],[2,2],[2,3],[3,1],[3,2],[4,1],[4,3],[5,1],[5,2],[5,3],[6,2],[6,3],[7,3]],4,10,"-31586835459","1866784"],
  [6,2,9,52,[[1,1],[1,2],[2,2],[3,3],[4,2],[4,3],[5,1],[5,2],[5,3],[6,2],[6,3],[7,3],[8,3],[9,3]],4,10,"-40544875212219","40472738"],
  [8,0,9,48,[[2,1],[2,2],[2,3],[3,1],[3,2],[3,3],[4,1],[4,2],[4,3],[5,1],[5,2],[5,3],[5,4],[6,1],[6,2],[6,3],[6,4],[7,3],[7,4],[8,2],[8,4]],4,10,"-84779040015","648658481"],
  [8,1,9,50,[[1,1],[1,2],[2,1],[2,2],[2,3],[3,1],[3,2],[3,3],[4,3],[4,4],[5,1],[5,2],[5,3],[5,4],[6,3],[6,4],[7,2],[7,3],[8,4],[9,4]],4,10,"-44660602996","622291"],
  [8,2,9,52,[[1,1],[1,2],[3,3],[3,4],[4,1],[4,2],[4,3],[4,4],[5,1],[5,2],[5,3],[5,4],[6,2],[6,3],[6,4],[7,4],[8,4],[9,1],[9,4]],4,10,"-47604717833301114400","211935764919671"],
  [10,0,9,48,[[2,1],[2,2],[2,3],[3,1],[3,2],[3,3],[3,4],[4,1],[4,2],[4,4],[5,1],[5,2],[5,3],[5,4],[5,5],[6,1],[6,2],[6,3],[6,4],[7,3],[7,4],[7,5],[8,2],[8,3],[8,4],[8,5]],4,10,"-38766480383412","1357767954535"],
  [10,1,9,50,[[1,1],[1,2],[2,1],[2,2],[2,3],[3,1],[3,2],[3,3],[4,1],[4,2],[4,4],[4,5],[5,1],[5,2],[5,3],[5,4],[5,5],[6,2],[6,3],[6,4],[7,2],[7,4],[8,4],[8,5],[9,5]],4,10,"-14892675175101","1150611856"],
  [10,2,9,52,[[1,1],[1,2],[1,3],[2,2],[2,3],[3,4],[4,1],[4,2],[4,3],[4,4],[4,5],[5,1],[5,2],[5,3],[5,4],[5,5],[6,1],[6,3],[6,4],[6,5],[7,4],[7,5],[8,4],[8,5],[9,1],[9,5]],4,10,"-8192294437835353093603483","1934853988942150377510"],
  [12,0,9,54,[[1,1],[1,2],[2,2],[2,3],[2,4],[3,2],[3,3],[4,1],[4,2],[4,3],[4,4],[4,5],[5,2],[5,3],[5,4],[5,5],[6,1],[6,2],[6,3],[6,4],[6,5],[7,2],[7,4],[7,5],[7,6],[8,4],[8,5],[8,6],[9,4],[9,5],[9,6]],4,10,"-21033703934365679897","41056764410159664"],
  [12,1,9,50,[[1,1],[1,2],[1,3],[2,1],[2,2],[2,3],[2,4],[3,1],[3,2],[3,3],[3,4],[4,1],[4,2],[4,4],[4,5],[5,1],[5,2],[5,3],[5,4],[5,5],[5,6],[6,2],[6,3],[6,4],[6,5],[7,2],[7,4],[7,5],[8,4],[8,5],[8,6],[9,6]],4,10,"-268852992846","1842056209"],
  [12,2,9,52,[[1,1],[1,2],[1,3],[2,3],[3,1],[3,4],[3,5],[4,1],[4,2],[4,3],[4,4],[4,5],[4,6],[5,1],[5,2],[5,3],[5,4],[5,5],[5,6],[6,2],[6,3],[6,4],[6,5],[6,6],[7,4],[7,5],[7,6],[8,4],[8,5],[8,6],[9,5],[9,6]],4,10,"-85378159185883243646877","227002729271235448858"],
  [14,0,9,54,[[1,1],[1,2],[1,3],[2,1],[2,2],[2,3],[2,4],[3,1],[3,2],[3,3],[3,4],[4,1],[4,2],[4,3],[4,4],[4,5],[4,6],[5,1],[5,2],[5,3],[5,4],[5,5],[6,2],[6,3],[6,4],[6,5],[7,2],[7,4],[7,5],[7,6],[7,7],[8,4],[8,5],[8,6],[8,7],[9,5],[9,6],[9,7]],4,10,"-17353961299554749","157427276410404"],
  [14,1,9,54,[[0,1],[1,1],[1,2],[1,3],[2,3],[3,1],[3,2],[3,4],[3,5],[4,1],[4,2],[4,3],[4,4],[4,5],[4,6],[5,1],[5,2],[5,3],[5,4],[5,5],[6,2],[6,3],[6,4],[6,5],[6,6],[6,7],[7,2],[7,5],[7,6],[7,7],[8,2],[8,4],[8,5],[8,6],[8,7],[9,6],[9,7]],4,10,"-880448417188325079721867437","521807225436600908199080"],
  [14,2,9,52,[[0,1],[0,2],[1,1],[1,2],[1,3],[2,2],[2,3],[2,4],[3,1],[3,2],[3,4],[3,5],[4,1],[4,2],[4,3],[4,4],[4,5],[4,6],[5,1],[5,2],[5,3],[5,4],[5,7],[6,2],[6,3],[6,4],[6,5],[6,6],[6,7],[7,4],[7,5],[7,6],[7,7],[8,5],[8,6],[8,7]],4,10,"-9154068453","441373016"],
  [16,1,12,54,[[0,1],[1,1],[1,2],[1,3],[2,4],[3,2],[3,3],[3,5],[4,1],[4,2],[4,3],[4,4],[4,5],[4,6],[5,1],[5,2],[5,3],[5,4],[5,5],[5,6],[6,2],[6,3],[6,4],[6,5],[6,7],[6,8],[7,2],[7,5],[7,6],[7,7],[7,8],[8,6],[8,7],[8,8],[9,6],[9,7],[9,8],[10,7],[10,8],[11,8],[12,8],[12,10]],6,20,"-8271574393708497483","1555208441960758676"],
  [16,2,9,54,[[0,1],[0,2],[1,1],[1,2],[1,3],[2,1],[2,3],[2,4],[2,5],[3,1],[3,2],[3,5],[3,6],[4,1],[4,2],[4,3],[4,4],[4,5],[4,6],[5,2],[5,3],[5,4],[5,5],[6,2],[6,3],[6,5],[6,6],[6,7],[6,8],[7,2],[7,4],[7,5],[7,6],[7,7],[8,2],[8,4],[8,5],[8,6],[8,7],[8,8],[9,4],[9,6],[9,7],[9,8]],4,10,"-792604424129012755769860129666291251","6196313240750342882095231052658160"]
];
for(const a of cases){
 const [w,r,,gate,,,,num,den]=a,c=certificate(...a);
 ck(c.vf===1&&c.pos&&c.obj[0]===num&&c.obj[1]===den&&BigInt(num)<0n,
    "exact Farkas contradiction failed: "+JSON.stringify([w,r,c]));
 let d=gate-2;while(d%6!==2*r)d-=2;
 ck(18n*BigInt(8*d)**47n<BigInt(w)*423n**47n,"exact density comparator");
 console.log("PASS normal="+w+" residue="+r+" raw upper="+d+" Farkas RHS="+num+"/"+den);
}
ck(9n*416n**47n>4n*423n**47n,"norm-eight trial threshold");
for(const [D,M,d] of [[6,16,48],[8,52,46],[10,59,46],[12,23,48],[14,9,50]]) {
 ck(2n*BigInt(M)**2n*BigInt(8*d)**47n<BigInt(D)*423n**47n,"arbitrary ternary motif cutoff");
 ck(2n*BigInt(M+1)**2n*BigInt(8*d)**47n>BigInt(D)*423n**47n,"next cardinality comparison");
}
console.log("PASS: all 17 scoped certificates; exact center-density comparison.");
