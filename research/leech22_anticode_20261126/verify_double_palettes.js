'use strict';
const fs=require('node:fs');
const zlib=require('node:zlib');
const readJSON=p=>{const b=fs.readFileSync(p);return JSON.parse((p.endsWith('.gz')?zlib.gunzipSync(b):b).toString('utf8'));};
const data=readJSON(process.argv[2]||'rows.json.gz'), w=readJSON(process.argv[3]||'witnesses.json');
if(data.matrices.length!==430||w.length!==90)throw Error('scope');
const rows=data.matrices.map((r,i)=>[i,r.normals,r.gram,r.high,r.low,r.undecided]);
console.log(((function(original,proves){
 const Gol=[];for(let i=0;i<4096;i++){let c=0;for(let j=0;j<12;j++)if((i>>j)&1)c^=0xae3<<j;let p=0;for(let h=c;h;h&=h-1)p^=1;Gol.push(c|(p<<23))}const Gset=new Set(Gol);const X=new Int8Array(196560*24), arr=[], pow=Array.from({length:24},(_,j)=>1<<j), ksf=new Uint8Array(1<<12);for(let i=1;i<ksf.length;i++)ksf[i]=1+ksf[i&(i-1)];
 function pop(b){let k=0;while(b){k++;b&=b-1}return k}
 let z=0;
 for(let i=0;i<24;i++)for(let j=0;j<i;j++)for(let u of [-4,4])for(let v of [-4,4]){X[24*z+i]=u;X[24*z+j]=v;z++}
 for(let c of Gol){if(pop(c)!==8)continue;let cols=pow.map((p,j)=>[p,j]).filter(x=>x[0]&c).map(x=>x[1]);for(let s=0;s<256;s++){if(pop(s)%2)continue;for(let j=0;j<8;j++) X[24*z+cols[j]]=2*(1-2*((s>>j)&1));z++}}
 for(let c of Gol)for(let j=0;j<24;j++){for(let i=0;i<24;i++)X[24*z+i]=(i===j?-3:1)*(1-2*((c>>i)&1));z++}
 if(z!==196560)throw Error('shell');

let by=new Map(proves);let num=0;const fail=(m)=>{throw Error(m)};
let count=0;
const cache=new Map();
function col(v){let k=v.join(',');if(cache.has(k))return cache.get(k);let a=new Int8Array(z),sp=v.map((h,i)=>[h,i]).filter(x=>x[0]);for(let i=0;i<z;i++){let d=0;for(let j=0;j<sp.length;j++)d+=sp[j][0]*X[24*i+sp[j][1]];if(d%8)fail('dot');a[i]=d/8}cache.set(k,a);return a}
const norm=(x,y)=>x.reduce((s,v,i)=>s+v*y[i],0)/8;
const leech=x=>{let m=(x[0]%2+2)%2,word=0,sum=0;for(let i=0;i<24;i++){if((x[i]-m)%2)return false;word|=(((x[i]-m)/2)&1)<<i;sum+=x[i]}return Gset.has(word)&&(sum-4*m)%8===0};
for(const [idx,V,[a,b,c,m],F0,lo,orig] of original){
 const D=a*c-b*b,n=2*D,s=x=>(x%D+D)%D;let labels=[],map=new Map();
 for(let p=0;p<D;p++)for(let q=0;q<D;q++)map.set(s(c*p-b*q)*D+s(a*q-b*p),0);
 labels=[...map.keys()].sort((a,b)=>a-b);labels.forEach((v,i)=>map.set(v,i));
 let diff=labels.map(x=>labels.map(y=>map.get(s(Math.floor(x/D)-Math.floor(y/D))*D+s(x%D-y%D))));
 if(!V.every(leech)||norm(V[0],V[0])!==a||norm(V[0],V[1])!==b||norm(V[1],V[1])!==c||norm(V[2],V[2])!==m||norm(V[0],V[2])!==0||norm(V[1],V[2])!==0||labels.length!==D)fail('geometry');
let F=Array(n).fill(16*D),uv=V.map(col),flags=0;
for(let j=0;j<z;j++){let [p,q,r]=[uv[0][j],uv[1][j],uv[2][j]],id=map.get(s(c*p-b*q)*D+s(a*q-b*p)),v=16*D-4*(c*p*p-2*b*p*q+a*q*q);
F[id]=Math.min(F[id],v);F[D+id]=Math.min(F[D+id],v+m*D-4*D*Math.abs(r));
if(p===1&&q===0&&r===0)flags|=1;if(p===0&&q===1&&r===0)flags|=2;if(p===0&&q===0&&r===1)flags|=4;
}
if(flags!==7)fail('primitivity');
F[D]=Math.min(F[D],D*m);F[0]=0;
if(F.join(',')!==F0.join(','))fail('single shell');
let req=F.map((v,i)=>1n<<BigInt(i));
 for(let [id,r,j,k,sign] of (by.get(idx)||[])){
  let v=Array.from({length:24},(_,i)=>X[24*j+i]+X[24*k+i]),p=V.map(x=>x.reduce((s,u,i)=>s+u*v[i],0)/8);
  if(!p.every(Number.isInteger))fail('integrality');
  let h=map.get(s(c*p[0]-b*p[1])*D+s(a*p[1]-b*p[0])),q=v.reduce((s,x)=>s+x*x,0)/8;
  let val=4*D*q - 4*(c*p[0]*p[0]-2*b*p[0]*p[1]+a*p[1]*p[1]);
  if(id>=D)val+=D*m+4*D*sign*p[2];else if(sign!==0)fail('sign');
  if(h!==id%D||r!==val)fail('invalid witness '+idx+' '+id);
  F[id]=Math.min(F[id],val);let inv=diff[0][h]+(id>=D?D:0);F[inv]=Math.min(F[inv],val);num++;
 }
 for(let q of new Set(F.slice(1))){
  if(q<=0)continue;
  let wins=k=>12n*BigInt(k)**2n*(3n*BigInt(q))**22n>9n*BigInt(D)*(44n*BigInt(D))**22n;
  if(!wins(n))continue;
  let k=1;while(!wins(k))k++;let good=F.map((v,i)=>i!==0&&v>=q);
  let adj=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>diff[i%D][j%D]+D*(Math.floor(i/D)^Math.floor(j/D))).reduce((r,h,j)=>good[h]?r|req[j]:r,0n));
  const dfs=(P,sz)=>{if(sz===k)return true;if(!P)return false;let c=0,R=P,ord=[],co=[];
   while(R){let Q=R;c++;while(Q){let u=Q&-Q,i=u.toString(2).length-1;R^=u;ord.push(i);co.push(c);Q &= R&~adj[i];}}
   for(let j=ord.length-1;j>=0;j--){if(sz+co[j]<k)return false;let i=ord[j];if(dfs(P&adj[i],sz+1))return true;P^=req[i];}return false;
  };
  count++;if(dfs(adj[0],1))fail('unresolved '+idx);
 }
}
return [original.length,num,count];}))(rows,w));
