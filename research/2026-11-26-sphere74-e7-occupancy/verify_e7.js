// Sphere74: exact common colors and mean-occupancy certificate.
const V=Array.from({length:729},(_,i)=>Array.from({length:6},(_,j)=>Math.floor(i/3**j)%3));
const D=(a,b)=>a.reduce((s,x,i)=>s+(x!=b[i]),0);
const P=[[1,1,1,0,0,0],[1,2,2,0,0,0],[2,1,2,0,0,0],[2,2,1,0,0,0]];
const A=Array(8).fill(0);
function solve(S,recurse=false){
const d=V.map(v=>S.map(u=>D(u,v)));
const ok=d.map(a=>a.every(n=>n>0&&n%3!=1));
if(S.length>7)throw Error("8 colors");
A[S.length]=Math.max(A[S.length],ok.filter(Boolean).length);
if(recurse)for(let i=0;i<729;i++)if(d[i].every(n=>n%3==2))solve([...S,V[i]],true);
} 
for(let n=1;n<4;n++)solve(P.slice(0,n));
const third=P[0].map((a,i)=>(6-a-P[1][i])%3);
solve([P[0],P[1],third]);solve(P,true);
const expected=[0,476,313,240,155,122,100,78];
if(A.join()!=expected.join())throw Error('counts');
for(let j=0;j<=7;j++)if(2*j*A[j]+18*j*j>1323+93*j)throw Error('majorant');
console.log(A,(93+Math.sqrt(770697))/288,512/19683*(93+Math.sqrt(770697))/288);
