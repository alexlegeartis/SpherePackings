
const vv=[];for(let a=0;a<=6;a++)for(let b=0;a+b<=6;b++)for(let c=0;a+b+c<=6;c++)vv.push([a,b,c,6-a-b-c]);
const order=[0,3,4,1]; // binary and ternary nonzero in last index
const indx=x=>x[0]+7*x[1]+49*x[2]+343*x[3],ind=new Map(vv.map((v,i)=>[indx(v),i]));
const cost=vv.map(v=>v.reduce((s,x,k)=>s+x*order[k],0));
const K0=[[1,1,1,1],[3,-1,3,-1],[2,2,-1,-1],[6,-2,-3,1]];
const K=vv.map(x=>{
 let m=new Map([[0,1]]);
 for(let i=0;i<4;i++)for(let n=0;n<x[i];n++){
  let g=new Map;
  for(let [k,v] of m)for(let q=0;q<4;q++){let t=k+7**q;g.set(t,(g.get(t)??0)+v*K0[q][i]);}
  m=g;
 }
 return vv.map(y=>m.get(indx(y)));
});
const z=vv.findIndex(x=>x[0]===6);

function verify(data){
 const A=(c,m)=>{if(!c)throw Error(m);};
 const range=n=>Array.from({length:n},(_,i)=>i);
 const omega=x=>[0,2,3,1][x];
 const H=range(64).map(i=>{let a=i&3,b=(i>>2)&3,c=i>>4;return [a,b,c,a^b^c,a^omega(b)^omega(omega(c)),a^omega(omega(b))^omega(c)];});
 A(H.length===new Set(H.map(String)).size&&Math.min(...H.slice(1).map(r=>r.filter(Boolean).length))===4,"hexacode");
 for(let b=0;b<4;b++)for(let t=0;t<3;t++){
   let X=3*(b&1)+2*t,Y=3*(b>>1)+2*t,m=Infinity;
   for(let i=-2;i<=2;i++)for(let j=-2;j<=2;j++){let x=X-6*i,y=Y-6*j;m=Math.min(m,x*x+x*y+y*y);}
   A(m===3*(t?(b?1:4):(b?3:0)),"torus metric");
   A(Math.max(X,Y)<=7,"no farther representatives can improve");
 }
 A(data.length===11&&data.every((r,i)=>r.D===i+2),"all thresholds");
 const out=[[1,12**6,Math.floor(64*12**6)+1]];
 for(const r of data){
  const m=BigInt(r.m),U=BigInt(r.U),c=r.c.map(([y,n])=>[y,BigInt(n)]);
  A(m>0n&&c.every(([y,n])=>y!==z&&y>=0&&y<84&&n>=0n),"dual nonnegative nontrivial weights");
  for(let j=0;j<84;j++)if(cost[j]>=r.D)
    A(-c.reduce((s,[y,n])=>s+n*BigInt(K[j][y]),0n)>=m,"dual transform on every allowed distance profile");
  A(U===m+c.reduce((s,[y,n])=>s+n*BigInt(K[z][y]),0n),"exact dual objective");
  const u=U/m, t=64n*12n**6n/BigInt(r.D)**6n+1n;
  A(u===BigInt(r.floor),"integer bound");
  if([2,3,4,10,11].includes(r.D))A(u<t,"closed threshold");
  out.push([r.D,Number(u),Number(t)]);
 }
 A(out.find(v=>v[0]===12)[1]===69,"69 phase ceiling");
 return out;
}
if(typeof module!=="undefined")module.exports=verify;
if(typeof require!=="undefined"&&require.main===module)console.log(verify(require("./eisenstein_duals.json")));

