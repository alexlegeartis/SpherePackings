'use strict';
// Standalone exact arithmetic replay for the refined 96D D4 tail.
// Imports no solver. The inherited MW92 phase-depth proof is separately
// certified in research/mw92_96_periodic12_20261124/ (commit 50c5bc5f).
const G=[[325000,-4771,-5229,-458],[-4771,325458,-5000,229],
         [-5229,-5000,324542,229],[-458,229,229,315000]];
const q=120000, T=639084, H=660000, Z=1920000;
function mod(n,m){const r=n%m;return r<0?r+m:r;}
const co=(i,j,k,l)=> {
  if((i+j+k+l)%2)throw Error('not D4');
  return [i,i+j,(i+j+k-l)/2,(i+j+k+l)/2];
};
const color=([a,b,c,d])=>[mod(c+d,3),mod(a,2),mod(b+c,2)];
const image=new Set();
for(let a=0;a<6;a++)for(let b=0;b<6;b++)for(let c=0;c<6;c++)for(let d=0;d<6;d++)
  image.add(color([a,b,c,d]).join());
if(image.size!==12)throw Error('not onto');
if(T!==q*(16*10000-106743)/10000||H!==q*(32-21)/2||Z!==16*q)
  throw Error('phase thresholds');
const diag=G.map((r,i)=>r[i]-r.reduce((s,z,j)=>s+(i===j?0:Math.abs(z)),0));
if(Math.min(...diag)*7<=16*q)throw Error('long-vector bound');
let rows={},checked=0,tight=0;
for(let i=-2;i<=2;i++)for(let j=-2;j<=2;j++)for(let k=-2;k<=2;k++)for(let l=-2;l<=2;l++){
  const v=[i,j,k,l], n=v.reduce((s,x)=>s+x*x,0);
  if(!n||n>6||(i+j+k+l)%2)continue;
  const [g,a,b]=color(co(...v));
  const goal=g?T:(a||b)?H:Z;
  const len=v.reduce((s,z,u)=>s+z*G[u].reduce((s,w,t)=>s+w*v[t],0),0);
  if(len<goal)throw Error('short vector '+JSON.stringify([v,len,goal]));
  rows[goal]??={count:0,min:Infinity};rows[goal].count++;
  rows[goal].min=Math.min(rows[goal].min,len);
  checked++;if(len===goal)tight++;
}
if(checked!==144||tight!==24||rows[T].count!==90||rows[H].count!==42||rows[Z].count!==12)throw Error('coverage');
function det(m){return m.length===1?BigInt(m[0][0]):m[0].reduce((s,v,j)=>s+(j%2?-1n:1n)*BigInt(v)*det(m.slice(1).map(r=>r.filter((_,i)=>i!==j))),0n);}
function gcd(a,b){while(b){let x=a%b;a=b;b=x;}return a;}
const Dn=4n*det(G), Dd=BigInt(q)**4n, gcdD=gcd(Dn,Dd);
const a=Dn/gcdD,b=Dd/gcdD;
if(a!==300153807922100058481n||b!==1440000000000000000n||256n*a>=53361n*b)
  throw Error('determinant or exact comparison');
const scale=10n**80n;
function arctanScaled(k,N){
  let low=0n,high=0n,p=BigInt(k);
  for(let j=0;j<N;j++){
    const d=p*BigInt(2*j+1),fl=scale/d,hi=(scale+d-1n)/d;
    if(j%2){low-=hi;high-=fl;}else{low+=fl;high+=hi;}
    p*=BigInt(k)**2n;
  }
  const d=p*BigInt(2*N+1),e=(scale+d-1n)/d;
  return N%2?[low-e,high]:[low,high+e];
}
const pi=16n*arctanScaled(5,70)[0]-4n*arctanScaled(239,25)[1];
const ceilRoot=n=>{let x=1n<<BigInt(Math.ceil(n.toString(2).length/2));
  while(true){let y=(x+n/x)>>1n;if(y>=x)return x;x=y;}};
let fact=1n;for(let i=2n;i<=48n;i++)fact*=i;
const floor=ceilRoot((pi**96n*2n**192n*b*10n**80n)/(scale**96n*fact**2n*47n**14n*a));
if(floor!==6367793872884866411n)throw Error('density arithmetic');
console.log(JSON.stringify({checked,tight,rows,gershgorin_numerator:Math.min(...diag),
  squared_tail_covolume:[a.toString(),b.toString()],occupied_significand:floor.toString(),exponent:-40,
  gain_percent:(Math.sqrt(Number(53361n*b)/Number(256n*a))-1)*100}));
