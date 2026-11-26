'use strict';
// Standalone exact verified covers; simplex below only proposes distributions. No external libraries or data.
const fs=require('fs'); const C=(function(){
// Conway–Sloane, Quaternary constructions (1994), pentacode cyclic-negation form.
const Z=['01112','03110','21310','21132'],G=[0,1,3,2],C=new Set();
for(const p of Z)for(let k=0;k<5;k++)for(const s of [1,3]){
 let u=0;for(let i=0;i<5;i++)u|=G[(+p[(i+k)%5]*s)%4]<<(2*i);C.add(u);
}
const pc=x=>{let n=0;for(;x;x&=x-1)n++;return n;};
return [...C]})().sort((a,b)=>a-b);
const pc=x=>{let n=0;while(x){n++;x&=x-1}return n;};

function simplex(A,b,c){
 let m=b.length,n=c.length,T=Array.from({length:m+1},()=>new Float64Array(n+m+1)),N=n+m,R=N;
 for(let j=0;j<n;j++)T[m][j]=-c[j];
 for(let i=0;i<m;i++){if(b[i]<-1e-12)throw Error("negative rhs");T[i].set(A[i]);T[i][n+i]=1;T[i][R]=b[i];}
 let B=Array.from({length:m},(_,i)=>n+i);let steps=0;
 for(;;){if(steps++>1000)throw Error('many');let j=-1,lo=-1e-8;for(let k=0;k<N;k++)if(T[m][k]<lo){j=k;lo=T[m][k];}if(j<0)break;
 let i=-1,v=Infinity;for(let k=0;k<m;k++)if(T[k][j]>1e-8&&T[k][R]/T[k][j]<v-1e-9){v=T[k][R]/T[k][j];i=k;}if(i<0)throw Error('unbounded');
 let p=T[i][j];for(let k=0;k<=N;k++)T[i][k]/=p;for(let h=0;h<=m;h++)if(h!==i){let v=T[h][j];if(v){for(let k=0;k<=N;k++)T[h][k]-=v*T[i][k]}}B[i]=j}
 let x=Array(n).fill(0);for(let i=0;i<m;i++)if(B[i]<n)x[B[i]]=T[i][R];return {x,dual:Array.from(T[m].slice(n,n+m)),v:T[m][R],steps};}
function box(a,L,U,C){let A=[],b=[];for(let c of C){let s=a^c,t=0,r=[];for(let i=0;i<10;i++){let q=(s>>i)&1,l=L[i],u=U[i];t+=(l-q)**2;r.push((2*q-l-u)*(u-l))};r.push(1);A.push(r);b.push(t)}
for(let i=0;i<10;i++){let r=Array(11).fill(0);r[i]=1;A.push(r);b.push(1)}let R=simplex(A,b,[0,0,0,0,0,0,0,0,0,0,1]);return R}

const SCALE=4096;
const gcd=(a,b)=>b===0n?(a<0n?-a:a):gcd(b,a%b);
function approx(x){ if(x<1e-8)return [0n,1n];for(let q=1;q<=12000;q++){let p=Math.round(q*x);if(Math.abs(x-p/q)<1.e-9)return [BigInt(p),BigInt(q)]}return [BigInt(Math.round(1e8*x)),100000000n]}
function getWeights(d){let fr=d.slice(0,40).map(approx),D=1n;for(let [p,q] of fr)D=D*q/gcd(D,q);let w=fr.map(([p,q])=>p*(D/q));let g=w.reduce(gcd,0n);return w.map(a=>a/g);}
function check(a,l,u,w){
 let sum=w.reduce((s,t)=>s+t,0n);if(sum<=0n||w.some(x=>x<0n))throw Error("bad weights");
 let bound=0n;const T=BigInt(SCALE);
 for(let i=0;i<10;i++){let n=0n;for(let j=0;j<40;j++)if((C[j]^a)&(1<<i))n+=w[j];
  const evalQ=y=>sum*y*y+n*T*(T-2n*y);
  const x=evalQ(BigInt(l[i])),y=evalQ(BigInt(u[i]));bound+=x>y?x:y;
 }return 8n*bound<=20n*sum*T*T;
}
const BAD=new Set(Array.from({length:1024},(_,a)=>a).filter(a=>Math.min(...C.map(c=>pc(c^a)))===3));
const stats={boxes:0,leaves:0,maxDepth:0,maxCalls:0},glob={}, D=21/8;
function prove(a,l,u,dep){stats.boxes++;if(u.reduce((x,z)=>x+z*z,0)<SCALE*SCALE && BAD.has(a))return 'interior';stats.maxDepth=Math.max(stats.maxDepth,dep);
 let R=box(a,l.map(x=>x/SCALE),u.map(x=>x/SCALE),C),w=getWeights(R.dual);
 if(check(a,l,u,w)){stats.leaves++;return w.map((z,j)=>z?[j,String(z)]:null).filter(Boolean)}
 let j=-1,sc=-1;for(let i=0;i<10;i++){let step=(u[i]-l[i])/SCALE,t=Math.max(0,step*step*R.x[i]*(1-R.x[i]));if(t>sc&&(u[i]-l[i])>=2){sc=t;j=i}}
 if(sc<1.e-10){let ranks=Array.from({length:10},(_,i)=>i).filter(i=>u[i]-l[i]>=2).sort((i,j)=>((u[j]-l[j])*u[j]-(u[i]-l[i])*u[i]));j=ranks[0]??-1;}
 if(j<0||dep>30)throw Error('failed to cover '+[a,dep,l,u,R.v,sc]);
 let m=(l[j]+u[j])/2;if(!Number.isInteger(m))throw Error("nondyadic");
 let u1=u.slice();u1[j]=m;let l1=l.slice();l1[j]=m;
 return {j,left:prove(a,l,u1,dep+1),right:prove(a,l1,u,dep+1)};
}
function solveAll(){
 let census=[0,0,0,0],trees={};for(let a=0;a<1024;a++){
 let k=Math.min(...C.map(c=>pc(c^a)));census[k]++;trees[a]=[];
 let F=[-1];
 for(let i of F){let l=Array(10).fill(0),u=Array(10).fill(SCALE/2);if(i>=0)l[i]=SCALE/4;trees[a].push(prove(a,l,u,0))}
 }
 console.log(JSON.stringify({census,stats,code:C}));
 if(process.argv[2])fs.writeFileSync(process.argv[2],JSON.stringify({scale:SCALE,code:C,trees:trees}));
}
solveAll();

