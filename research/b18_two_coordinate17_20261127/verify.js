"use strict";
function ok(x,s){if(!x)throw Error(s);}
function wt(x){let n=0;while(x){x&=x-1;n++;}return n;}
const v=(1<<18)-1, z=(1<<17)-1;
function raw(seed){
  const a=seed.reduce((r,i)=>r|(1<<i),0),S=new Set([0]);
  for(let i=0;i<17;i++){
    const w=((a<<i)|(a>>(17-i)))&z;
    for(const u of [...S])S.add(u^w);
  }
  const b=[...S];return b.concat(b.map(u=>u^v));
}
const C=[raw([1,2,4,8,9,13,15,16]),raw([3,5,6,7,10,11,12,14])];
let q=0;
for(const A of C){
  ok(new Set(A).size===512&&A.every(a=>!(wt(a)%2)),"code size and parity");
  const W=A.filter(x=>wt(x)===6);
  ok(W.length===102&&Math.min(...A.filter(x=>x).map(wt))===6,"minimum");
  for(let i=0;i<18;i++)for(let j=i+1;j<18;j++){
    const I=1<<i,J=1<<j,K=I|J;
    const L=W.filter(w=>(w&K)===I),R=W.filter(w=>(w&K)===J);
    ok(L.some(a=>R.some(b=>wt(a^b)===6)),"rectangle "+i+","+j);q++;
  }
}
for(const a of C[0])for(const b of C[1])ok(!(wt(a&b)&1),"orthogonal banks");
ok(q===306,"all pairs");

const H=(1<<17)|1,B=v^H, S=[];
for(let b=0;b<2;b++)for(const w of C[b]){
  let x=w&1,y=w>>>17&1,h=x+2*y,j=b?(h===3?1:0):(h===0?1:0);
  S.push({b,w,j,s:(2*(x-y)+4*j+8)%8});
}
let mu=32,count=0;
for(let i=0;i<S.length;i++)for(let j=0;j<i;j++){
  const x=S[i],y=S[j];let d=Math.abs(x.s-y.s),out=0;d=Math.min(d,8-d);
  if(x.b===y.b){let h=wt((x.w^y.w)&B);out=h?4*h:16*((x.j^y.j)&1);}
  else{
    const [e,o]=x.b?[y,x]:[x,y];
    const parity=(1+e.j+o.j+wt(o.w&B)+wt(e.w&o.w&B))&1;
    out=16+8*parity;
  }
  const value=d*d+out;
  if(value<mu){mu=value;count=0;}
  if(value===mu)count++;
}
ok(S.length===1024&&mu===24&&count===112640,"literal 1024-site packing");
ok(128n*3n**17n<4n**17n,"near benchmark and endpoint");
ok(8n*3n**8n<4n**8n,"rectangle small-metric gate");
ok(1024n*41n**16n<64n**16n&&1024n*41n**17n<64n**17n,"free phase gates");
ok(20*64+8*81<32*64,"large radius crossover");
console.log("PASS: 306 QR rectangles; literal 1,024-site packing at minimum 24; exact density checks.");
function mob(w){
  let t=0;
  for(let i=0;i<18;i++)if((w>>i)&1){
    let j=i===0?17:0;
    if(i>0&&i<17)while((i*j)%17!==16)j++;
    t|=1<<j;
  }
  return t;
}
const opt=[
  [
    [0,76002,1966,35828,107704,59694,84536,25972,116322,45672,70014,15922,105764,53426,95140,23784,131070,95242,55068,129104,23366,71376,46534,105098,14748,85398,61056,115148,25306,77644,35930,107286],
    [0,48142,73496,1966,107704,97730,57044,25972,116322,55674,96876,115148,25306,48032,71862,107286],
  ],
  [
    [0,1884,23476,57956,93904,48780,78392,120576,101352,29514,12194,51866,38514,72390,83502,107286,131070,129186,107594,82290,73114,37166,52678,10494,29718,101556,118876,79204,92556,58680,47568,23784],
    [0,76396,95876,106154,20214,64048,42712,29514,107286,54674,35194,15804,126432,67022,88358,118876],
  ],
];
function cliqueAtLeast(A,limit){
  const n=A.length,bits=A.map((_,i)=>1n<<BigInt(i));
  const ix=new Map(bits.map((b,i)=>[b,i]));
  let calls=0;
  function visit(P,k){
    calls++;
    if(k>=limit)return true;
    let remaining=P,col=0,order=[],bounds=[];
    while(remaining){
      let avail=remaining;col++;
      while(avail){
        const b=avail&-avail,i=ix.get(b);
        order.push(i);bounds.push(col);
        // One class of the greedy colouring is independent.
        avail&=~(A[i]|b);remaining^=b;
      }
    }
    for(let l=order.length-1;l>=0;l--){
      if(k+bounds[l]<limit)return false;
      const i=order[l];
      if(visit(P&A[i],k+1))return true;
      P^=bits[i];
    }
    return false;
  }
  return [visit((1n<<BigInt(n))-1n,0),calls];
}
let stats=[];
for(let b=0;b<2;b++){
  const src=new Set(C[b]);
  ok(C[b].every(w=>src.has(mob(w))),"Möbius inversion");
  const punct=new Set(C[b].map(w=>w&B));ok(punct.size===512,"puncture injective");
  for(const d of [7,8]){
    const goal=d===7?32:16;
    const E=[...punct].filter(w=>wt(w)>=d);
    const A=E.map(u=>E.reduce((a,w,i)=>wt(u^w)>=d?a|(1n<<BigInt(i)):a,0n));
    const [found,calls]=cliqueAtLeast(A,goal);ok(!found,"too large puncture clique");
    stats.push([b,d,E.length,calls]);
    const X=opt[b][d-7];
    ok(X.length===goal&&X.every(w=>punct.has(w)),"explicit optimum word membership");
    for(let i=0;i<X.length;i++)for(let j=0;j<i;j++)ok(wt(X[i]^X[j])>=d,"optimum code distance");
  }
}
const num=126n*379n**17n,den=504n**17n;
ok(num<den&&2n*num>den,"strict ceiling, dominates 1/2 endpoint");
ok(1024n*41n**16n*den<64n**16n*num,"first window endpoint");
ok(1024n*41n**17n*den<64n**17n*num,"rectangle left endpoint");
ok(3n*11n**17n*den<2n*12n**17n*num,"two window crossing");
ok(100000000n*num<9954n**2n*den,"rounded upper bound 0.9954");
console.log("PASS: punctured QR exact maxima 32 and 16; exhaustive anchor stats "+JSON.stringify(stats));
const Q=[1,2,4,8,9,13,15,16].reduce((a,i)=>a|(1<<i),0),N=B^Q;
const forms=new Set(["32,0,0","0,32,0","0,0,64"]);
for(let i=0;i<S.length;i++)for(let j=0;j<i;j++){
  const x=S[i],y=S[j];let d=Math.abs(x.s-y.s);d=Math.min(d,8-d);
  if(x.b===y.b){let w=x.w^y.w;forms.add([4*wt(w&Q),4*wt(w&N),d*d].join());}
  else {
    const [e,o]=x.b?[y,x]:[x,y],p=(1+e.j+o.j+wt(o.w&B)+wt(e.w&o.w&B))&1;
    forms.add([8+8*p,8,d*d].join());forms.add([8,8+8*p,d*d].join());
  }
}
const allForms=[...forms].map(x=>x.split(",").map(Number));
const uf=allForms.filter(u=>!allForms.some(v=>u!==v&&u.every((x,i)=>x>=v[i])));
ok(uf.map(x=>x.join()).sort().join(";")===[
 "32,0,0","0,32,0","0,0,64","8,16,0","16,8,0","4,16,4","16,4,4",
 "12,8,4","8,12,4","8,8,16","16,0,16","0,16,16"
].sort().join(";"),"exact twelve forms");
const dot=(u,v)=>u.reduce((s,a,i)=>s+a*v[i],0);
const cross=(u,v)=>u.map((_,i)=>u[(i+1)%3]*v[(i+2)%3]-u[(i+2)%3]*v[(i+1)%3]);
const gcd=(x,y)=>y?gcd(y,x%y):x;
const vs=new Set();
for(let i=0;i<uf.length;i++)for(let j=0;j<i;j++)for(let k=0;k<j;k++){
  const [u,v,w]=[uf[i],uf[j],uf[k]],a=cross(v,w),b=cross(w,u),c=cross(u,v);
  let d=dot(u,a),r=a.map((z,l)=>24*(z+b[l]+c[l]));
  if(!d)continue;
  if(d<0){d=-d;r=r.map(z=>-z);}
  if(uf.some(z=>dot(z,r)<24*d))continue;
  const g=r.reduce((p,q)=>gcd(p,q),d);
  vs.add(r.map(x=>x/g).concat(d/g).join());
  ok(BigInt(r[0])**8n*BigInt(r[1])**8n*BigInt(r[2])>=BigInt(d)**17n,"covolume");
}
ok([...vs].sort().join(";")===["6,3,6,4","1,1,1,1","9,9,3,8","3,6,6,4","3,9,3,4","9,3,3,4"].sort().join(";"),"all six vertices");
console.log("PASS: twelve metric forms, six exact vertices, no three-scale improvement.");
