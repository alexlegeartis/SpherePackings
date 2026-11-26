"use strict";
function assert(c){if(!c)throw Error("failed assertion")}
function det(A){if(A.length==1)return BigInt(A[0][0]);if(A.length==0)return 1n;let u=0n;for(let j=0;j<A.length;j++)if(A[0][j])u+=BigInt(A[0][j])*(j%2?-1n:1n)*det(A.slice(1).map(r=>r.filter((_,i)=>i!=j)));return u}
const mod=x=>(x%3+3)%3;
function dd(a,b){let s=0;for(let i=0;i<a.length;i++)s+=a[i]*b[i];return s}
function row3(A,rhs,cols){let n=A.length,M=A.map((a,i)=>a.map(mod).concat(mod(rhs[i]))),p=[],f=[];let r=0;for(let j=0;j<cols;j++){let k=r;while(k<n&&!M[k][j])k++;if(k==n){f.push(j);continue}[M[r],M[k]]=[M[k],M[r]];if(M[r][j]==2)M[r]=M[r].map(x=>mod(2*x));for(let i=0;i<n;i++)if(i!=r&&M[i][j]){let t=M[i][j];for(let h=j;h<=cols;h++)M[i][h]=mod(M[i][h]-t*M[r][h])}p.push(j);r++}if(M.slice(r).some(a=>a[cols]))return null;let x=Array(cols).fill(0);for(let k=0;k<p.length;k++)x[p[k]]=M[k][cols];return [x,f.map(j=>{let y=Array(cols).fill(0);y[j]=1;for(let k=0;k<p.length;k++)y[p[k]]=mod(-M[k][j]);return y})]}

function low(maxTotal=20,lim=8){let pairs=[],E=Array.from({length:8},()=>Array(8).fill(-1));for (let i=0;i<8;i++)for(let j=i+1;j<8;j++){E[i][j]=E[j][i]=pairs.length;pairs.push([i,j])}let cuts=Array.from({length:128},(_,s)=>s).filter(s=>s.toString(2).split('1').length-1==4).map(s=>pairs.map(([i,j])=>(s>>i^s>>j)&1?1:-1)),x=Array(28).fill(0),deg=Array(8).fill(0),list=[],exam=0,sort=0,pure={},rr={};
 function check(t){exam++;for(let i=1;i<8;i++)if(deg[i-1]<deg[i])return;sort++;
 let actual=t,cost=[];for(let a of cuts){let h=mod(dd(a,x));cost.push(h);actual+=h;if(actual>maxTotal)return}
 let Q=Array.from({length:7},(_,i)=>Array.from({length:7},(_,j)=>i==j?mod(2+x[E[i][7]]):mod(2*(2+x[E[i][7]]+x[E[j][7]]-x[E[i][j]]))));let rad=row3(Q,Array(7).fill(0),7)[1],r=7-rad.length;pure[actual+','+r]=(pure[actual+','+r]||0)+1;if(r>6)return;
 if(r==6){let j=rad[0].findIndex(z=>z!=0),res=Q.filter((_,i)=>i!=j).map(r=>r.filter((_,i)=>i!=j)),d=det(res);if(mod(Number(d%3n))!=1)return}
 function isRad(v){return Q.every(a=>mod(dd(a,v))==0)}let nv=0;
 for(let [i,j] of pairs){let z=Array(7).fill(0);z[i]=1;if(j<7)z[j]=-1;if(isRad(z))nv++}
 for(let s=0;s<128;s++)if(s.toString(2).split('1').length-1==4){let z=Array.from({length:7},(_,i)=>s>>i&1?1:-1);if(isRad(z))nv++}
 if(r==6)actual+=3*nv;
 if(actual>maxTotal)return;rr[r+','+actual]=(rr[r+','+actual]||0)+1;list.push([r,actual,t,nv,x.join(''),cost.join('')])}
 function recur(j,t){check(t);if(t>=lim)return;for(let k=j;k<28;k++){let [i,y]=pairs[k];x[k]=1;deg[i]++;deg[y]++;recur(k+1,t+1);if(t+1<lim){x[k]=2;deg[i]++;deg[y]++;recur(k+1,t+2);deg[i]--;deg[y]--}deg[i]--;deg[y]--;x[k]=0}}
 recur(0,0);return {exam,sort,pure,rr,list}}
const R=low(20,8);assert(R.exam==23779362 && R.sort==43164);
const pairs=[];for(let i=0;i<8;i++)for(let j=i+1;j<8;j++)pairs.push([i,j]);
const sets=Array.from({length:128},(_,i)=>i).filter(s=>s.toString(2).split("1").length-1==4);
let good=0;for(let [rank,a,t,nv,S,T] of R.list){
 let x=[...S].map(Number),y=[...T].map(Number),vals=[t];
 for(let s of sets){let c=0;for(let k=0;k<28;k++){let [i,j]=pairs[k];if(((s>>i^s>>j)&1)==0)c+=x[k];else{let q=s^(1<<i)^(1<<j);if(q&128)q^=255;c+=y[sets.indexOf(q)]}}vals.push(c)}
 if(Math.max(...vals)>9)continue;good++;assert(rank==6&&nv==0);
 let N=Array.from({length:8},()=>Array(8).fill(0));for(let i=0;i<28;i++){let [u,v]=pairs[i];N[u][v]=N[v][u]=4+2*x[i]}
 let g=Array.from({length:7},(_,i)=>Array.from({length:7},(_,j)=>(N[i][7]+N[j][7]-N[i][j])/2));for(let i=1;i<=7;i++)assert(det(g.slice(0,i).map(r=>r.slice(0,i)))>0n);assert(det(g)==1536n);
 for(let h=0;h<35;h++){let s=sets[h],sign=Array.from({length:7},(_,i)=>(s>>i&1)?1:-1);assert(dd(sign,g.map(r=>dd(r,sign)))==4*(4+2*y[h]))}
}
assert(good==7);
for(let k=1;k<=6;k++){let a=BigInt([275,189,159,143,131,121][k-1]),b=100n,K=BigInt(k),J=7n-K,c=J*a,d=7n*a-K*b;assert(2n*a**K*c**J>3n*b**K*d**J);assert(4n*(K*a*d+J*c*b)<33n*b*d)}
const [G,F]=[[[4,-2,0,0,2,2,0],[-2,4,-2,-2,0,0,0],[0,-2,4,0,0,0,0],[0,-2,0,4,0,0,2],[2,0,0,0,4,2,2],[2,0,0,0,2,4,1],[0,0,0,2,2,1,4]],[[1,-1,0,0,0,0],[0,1,-1,0,0,0],[0,0,1,-1,0,0],[0,0,1,1,0,0],[1,0,0,0,1,0],[1,0,0,0,0,1],[2,2,2,2,-1,0]]];
assert(det(G)==384n);
const Z=Array(7).fill(0),L=G.map((_,i)=>{let co=det(G.filter((_,j)=>i!=j).map(r=>r.filter((_,j)=>i!=j)));let b=0;while(BigInt((b+1)**2)*384n<=12n*co)b++;return b});
let counts={};function vec(k){if(k<7){for(let v=-L[k];v<=L[k];v++){Z[k]=v;vec(k+1)}return}if(!Z.some(a=>a))return;let value=dd(Z,G.map(r=>dd(Z,r)));if(value>12)return;let wt=Array.from({length:6},(_,j)=>mod(Z.reduce((s,b,i)=>s+b*F[i][j],0))).filter(v=>v).length;let req=!wt?12:wt%3==0?6:wt%3==1?8:4;assert(value>=req);if(value==req)counts[req]=(counts[req]||0)+1}
vec(0);assert(counts[4]==92&&counts[6]==64&&counts[8]==436&&!counts[12]);
console.log("Checked all 23779362 graph inputs, 43164 degree-sorted, 7 certified forms, sharp determinant 384 and explicit attaining seed");
