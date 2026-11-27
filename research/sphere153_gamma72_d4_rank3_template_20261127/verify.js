'use strict';
const assert=require('node:assert/strict');
const E=[[0,1],[1,2],[2,3],[3,0],[0,2],[1,3]];
const dot=(a,b)=>a.reduce((v,x,i)=>v+x*b[i],0);
const len=v=>dot(v,v);
const ww=v=>E.map(([a,b])=>(v[a]-v[b])**2);
const cuts=[];
for(let m=1;m<15;m+=2) cuts.push(ww([0,1,2,3].map(j=>(m>>j)&1)));
const nat=Array.from({length:6},(_,i)=>Array.from({length:6},(_,j)=>+(i===j)));
const constraints=[...nat.map(x=>[x,0]),...cuts.map(x=>[x,10]),[[1,0,1,0,-1,-1],0],[[0,1,0,1,-1,-1],0]];
const conormCoords=([x,y,b])=>[(x+y+b)/2,b,(x-y+b)/2,0];
const repQ=[2,4,6].map(L=>{
 const a=new Map();
 for(let x=-2;x<=2;x++)for(let y=-2;y<=2;y++)for(let z=-2;z<=2;z++)
 if(len([x,y,z])===L){
  const v=conormCoords([x,y,z]);assert(v.every(Number.isInteger));
  const t=ww(v);a.set(t.join(),[t,[x,y,z]]);
 }
 return [...a.values()];
});
assert.deepEqual(repQ.map(x=>x.length),[6,3,12]);
const integer=x=>{assert(Number.isSafeInteger(x));return x;};
function exactDiv(a,b){assert(b&&a%b===0);return integer(a/b);}
function solve(L){
 // Bareiss on the integer augmented matrix. Result x=N/D in scale x=5p.
 const A=L.map(([a,b])=>[...a,b]);let prev=1;
 for(let j=0;j<5;j++){
  let k=j;while(k<6&&A[k][j]===0)k++;if(k===6)return null;
  if(j!==k)[A[j],A[k]]=[A[k],A[j]];
  const pj=A[j][j];
  for(let m=j+1;m<6;m++){
   for(let i=j+1;i<7;i++)
     A[m][i]=exactDiv(integer(integer(pj*A[m][i])-integer(A[m][j]*A[j][i])),prev);
   A[m][j]=0;
  }
  prev=pj;
 }
 if(!A[5][5])return null;
 const D=Math.abs(A[5][5]),N=Array(6);
 for(let i=5;i>=0;i--){
   let b=integer(D*A[i][6]);
   for(let j=i+1;j<6;j++)b=integer(b-integer(A[i][j]*N[j]));
   N[i]=exactDiv(b,A[i][i]);
 }
 // Independently check A N = D b against the six untouched input equalities.
 for(const [a,b] of L)assert(dot(a,N)===D*b);
 return [N,D];
}
function detNum(p){
 const a=p[0]+p[3]+p[4],b=p[0]+p[1]+p[5],c=p[1]+p[2]+p[4],
 x=-p[0],y=-p[4],z=-p[1];
 return a*b*c+2n*x*y*z-a*z*z-b*y*y-c*x*x;
}
function minimum(L){
 let best=null, count=0, w;
 const ix=[];
 function go(k) {
  if(ix.length===6){
   const ans=solve(ix);if(!ans)return;
   const [N,D]=ans;
   if(L.some(([a,b])=>dot(a,N)<b*D))return;
   count++;
   const num=detNum(N.map(BigInt)),den=(5n*BigInt(D))**3n;
   assert(num>0n);
   if(!best||num*best[1]<den*best[0]){best=[num,den];w=[N,D];}
   return;
  }
  for(let i=k;i<=L.length-6+ix.length;i++){ix.push(L[i]);go(i+1);ix.pop();}
 }
 go(0);assert(best&&count>0);
 return [best,count,w];
}
function reduce([a,b]){const g=(x,y)=>y?g(y,x%y):x;const k=g(a,b);return [a/k,b/k];}
const caps=[15,29,49], expected=[[9n,2n],[779n,150n],[1053n,200n]];
for(let i=0;i<3;i++){
 let z=null, count=0;
 for(const [v,x] of repQ[i]){
  const [q,n]=minimum([...constraints,[v,caps[i]]]);
  assert(2n*q[0]>=9n*q[1]);
  count+=n;
  if(!z||q[0]*z[1]<q[1]*z[0])z=q;
 }
 assert.deepEqual(reduce(z),expected[i]);
 console.log('FCC norm',2+2*i,'min det above cutoff',reduce(z).join('/'),'feasible basic solutions',count);
}

// Nominal D3 tail congruences, all seven planes, and literal configurations.
function aSub(a,b){return a.map((x,i)=>x-b[i]);}
const bits=y=>{const [a,b,c]=y;
 const x=[(a+b-c)/2,(a+c-b)/2,(b+c-a)/2];
 assert(x.every(Number.isInteger));
 return x.reduce((s,x,i)=>s+((x&1)<<i),0);
};
const ks=v=>[(v&1)+((v>>1)&1),(v&1)+((v>>2)&1),((v>>1)&1)+((v>>2)&1)];
const plane=(a,b)=>[0,a,b,a^b].sort((a,b)=>a-b).join(',');
const P=new Set();for(let a=1;a<8;a++)for(let b=1;b<a;b++)P.add(plane(a,b));
const P0=new Set(),P1=new Set();
for(let i=0;i<3;i++)for(let j=0;j<3;j++)if(i!==j){
 const p=[i,j,3-i-j];
 for(let e=0;e<8;e++){
  const f=x=>bits(p.map((t,k)=>ks(x)[t]*(1-2*((e>>k)&1))));
  P0.add(plane(f(1),f(2)));
  P1.add(plane(f(1),f(6)));
 }
}
assert.equal(P0.size,4);assert.equal(P1.size,3);
assert.equal(new Set([...P0,...P1]).size,7);
assert([...P0,...P1].every(p=>P.has(p)));
const H=[[0,0,0,0],[2,0,0,0],[1,1,1,1],[3,1,1,1]], O=[0,1,2,3], T=[0,1,6,7];
const A=[
 [[0,0,0,0],[0,0,0]],
 [[-2,0,0,0],[-1,-1,0]],
 [[0,0,0,0],[-2,0,-2]],
 [[0,0,0,0],[-2,-2,0]],
 [[-1,-1,-1,-1],[-1,0,-1]],
 [[-1,-1,1,1],[-1,0,-1]],
 [[-1,1,-1,1],[-1,0,-1]],
 [[-1,1,1,-1],[-1,0,-1]]
];
const B=[
 [[0,0,0,0],[0,0,0]],
 [[0,0,0,0],[-2,-2,0]],
 [[0,0,0,0],[-2,0,-2]],
 [[-1,-1,-1,1],[-2,0,0]],
 [[-1,-1,1,-1],[-2,0,0]]
];
// Bounds for squared normal distances, in denominator 20.
// h difference raw norm^2 multiplied by (31/20)/4, tail bounds by /4.
const assertDist=(u,v,d)=>{
 const [h,y]=u,[g,w]=v;
 const H4=len(aSub(h,g)),Y4=len(aSub(y,w));
 // Bounds for norms of y: for y=2 root-vector, square scales by 4.
 const yCap=new Map([[0,0],[2,60],[4,116],[6,196],[8,240]]);
 assert(yCap.has(Y4));if(Y4===8){const z=aSub(y,w);assert(z.every(x=>x%2===0));assert(len(z.map(x=>x/2))===2);}
 assert(31*H4+yCap.get(Y4)<=80*d);
 assert(H4+Y4>0);
};
function test(X,U){
 for(const [h,y] of X){
   const i=U.indexOf(bits(y));assert(i>=0);
   const x=aSub(h,H[i]).map(t=>t/2);
   assert(x.every(Number.isInteger)&&x.reduce((s,t)=>s+t,0)%2===0);
 }
}
test(A,O);test(B,T);
for(let i=1;i<A.length;i++)for(let j=0;j<i;j++){
 // Edges to anchor, and edges from the pair {1,2} to common neighbors {3,...,7}.
 const short=j===0||(i>=3&&(j===1||j===2));
 assertDist(A[i],A[j],short?3:4);
}
for(let i=1;i<B.length;i++)for(let j=0;j<i;j++)
 assertDist(B[i],B[j],i===4&&j===3?4:3);

// Determinant budget of the comparator, and conservative simple rational caps.
const N=2n**79n*41n**11n*234n*936n**79n,D=4n*3235n**79n;
assert(81n*D<16n*N&&128n*N<729n*D); // 81/16 < D_* < 729/128
assert(729n*20n**4n<128n*31n**4n); // r < 31/20
assert(8n*729n<32n*6n**3n); // 2P short vector squared norm <6 in the low-scale branch
assert(324n*128n>729n*49n); // deleted support at r>=3/2 also loses
console.log('PASS: all seven tail planes, 8- and 5-point witnesses, determinant gate for 79D');
