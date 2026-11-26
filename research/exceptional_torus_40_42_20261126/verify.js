'use strict';
const assert=require('node:assert/strict');
const rows={"6-10":[[2,0,1,2,1,0],[1,2,1,2,0,2],[1,0,0,1,0,0],[1,1,2,2,0,2]],"6-16":[[0,1,2,2,2,2],[2,1,2,1,2,0],[0,1,0,0,1,2]],"7-7":[[1,0,2,1,2,1,2],[0,0,2,1,1,0,2],[1,2,2,1,2,2,0],[1,1,1,2,0,2,1]],"7-11":[[1,0,2,1,2,2,0],[2,0,1,1,1,2,2],[2,0,2,1,0,0,1]],"8-4":[[2,2,2,0,1,0,2,2],[0,1,0,2,2,0,2,0],[2,1,0,0,1,2,0,0],[1,1,1,1,1,1,2,0]],"8-6":[[1,0,0,1,0,0,2,2]]};
const motifs={"6":[[0,0,0,0,0,0],[5,0,1,2,1,2],[0,3,1,2,1,2],[1,2,1,2,1,2],[0,3,1,2,0,2],[3,2,1,0,0,2],[0,2,1,2,1,0],[0,3,0,2,1,0],[0,0,0,0,1,2]],"7":[[0,0,0,0,0,0,0],[0,0,1,5,1,2,2],[1,2,1,5,0,0,0],[1,2,0,3,1,2,0],[1,2,1,2,1,0,2],[0,2,1,5,1,0,0],[1,2,1,2,1,2,0],[1,2,1,3,0,0,2],[1,0,0,5,1,0,0]],"8":[[0,0,0,0,0,0,0,0],[0,0,2,1,2,1,2,1],[1,2,0,1,2,1,2,0],[0,0,0,1,2,1,2,2],[0,1,2,1,2,1,2,0],[2,0,0,1,2,1,2,0],[0,0,0,0,2,2,2,0],[0,0,2,1,0,1,2,0],[0,0,0,0,0,1,1,0]]};
function det(a){return a.length===1?a[0][0]:a[0].reduce((s,z,j)=>z?s+z*(-1)**j*det(a.slice(1).map(w=>w.filter((_,k)=>k!==j))):s,0);}
for(const r of [6,7,8]){
 const v=Array.from({length:r},(_,i)=>Array.from({length:r},(_,j)=>i===j?2:0));
 for(let i=0;i<r-2;i++)v[i][i+1]=v[i+1][i]=-1;v[2][r-1]=v[r-1][2]=-1;
 const d=det(v),m=3*d, A=v.map((_,i)=>v.map((_,j)=>(-1)**(i+j)*det(v.filter((_,k)=>k!==i).map(w=>w.filter((_,k)=>k!==j)))));
 const K=t=>t.reduce((s,x)=>m*s+((x%m+m)%m),0),lin=c=>A.map(w=>w.reduce((s,v,i)=>s+v*c[i],0)),Q=new Map(),c=Array(r).fill(0);let seen=0;
 // q<=2 implies c_i^2 <= 2*3*2=12 (Cauchy against each simple root); no numerical enumeration.
 function scan(i){if(i<r){for(let x=-3;x<=3;x++){c[i]=x;scan(i+1)}return}const h=lin(c),n=c.reduce((s,v,j)=>s+v*h[j],0);if(n>2*m)return; const id=K(h);Q.set(id,Math.min(n,Q.get(id)??1e7));seen++}scan(0);
 const expected={6:[0,4,6,10,12,16,18],7:[0,3,4,7,8,11,12],8:[0,2,4,6]}[r];
 assert.deepEqual([...new Set(Q.values())].sort((a,b)=>a-b),expected);
 // Explicit feasible motifs, and test that their torus diameter is exactly 2/3.
 let diam=0;for(let i=0;i<9;i++)for(let j=0;j<i;j++){let n=Q.get(K(lin(motifs[r][i].map((v,k)=>v-motifs[r][j][k]))));assert(n>0);diam=Math.max(diam,n)}assert.equal(3*diam,2*m);
 function verifyH(b,B){let H=new Set(),u=Array(r).fill(0);function rec(i){if(i===B.length){let t=K(u.map(v=>d*v));assert(!H.has(t));H.add(t);if(t)assert((Q.get(t)??1e7)>b);return}for(let a=0;a<3;a++){rec(i+1);for(let j=0;j<r;j++)u[j]=(u[j]+B[i][j])%3}}rec(0);return (d*3**r)/H.size}
 const caps={};for(const [key,B]of Object.entries(rows)){const [rr,b]=key.split('-').map(Number);if(rr===r)caps[b]=verifyH(b,B)}
 if(r===7){const h=[0,0,0,3,0,3,3];assert(!Q.has(K(h))&&K(h.map(z=>2*z))===0);caps[12]=d*3**r/2;}
 const gates={6:[[4,3],[6,9],[10,caps[10]],[12,caps[16]],[16,caps[16]],[18,d*3**r]],7:[[3,2],[4,9],[7,caps[7]],[8,caps[11]],[11,caps[11]],[12,caps[12]],[13,d*3**r]],8:[[2,9],[4,caps[4]],[6,caps[6]],[7,d*3**r]]}[r];
 for(const [b,M] of gates){const M2=BigInt(M)**2n,n=BigInt(48-r);assert(M2*BigInt(6*m-b)**n<BigInt(6*m)**n);}
 console.log(r,seen,diam,caps,gates);
}
