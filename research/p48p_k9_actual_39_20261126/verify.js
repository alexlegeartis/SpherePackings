const fs=require('fs'),zlib=require('zlib'),assert=require('assert/strict');
const G=require('./gram.json'),gens=require('./gens.json');
let n=39,N=21870;const pos=[],reps=[],m3=3**7;for(let k=0;k<N;k++){let j=k%10,y=k/10|0,h=gens[7].map(x=>j*x);for(let l=0;l<7;l++){let d=y%3;y=y/3|0;for(let i=0;i<n;i++)h[i]+=d*gens[l][i];}h=h.map(x=>x%30);reps.push(h);pos.push(String.fromCharCode(...h));}assert.equal(new Set(pos).size,N);
for(let [i,h] of gens.entries()){assert(h.every(x=>x*(i===7?10:3)%30===0));for(let row of G)assert(row.reduce((s,z,j)=>s+z*h[j],0)%30===0);}
// The Gram determinant is 21870: fraction-free Bareiss, positive leading minors.
function det(M){M=M.map(r=>r.map(BigInt));let d=1n;for(let i=0;i<M.length-1;i++){if(M[i][i]<=0n)throw Error('not positive definite');for(let j=i+1;j<M.length;j++)for(let k=i+1;k<M.length;k++){let a=M[i][i]*M[j][k]-M[j][i]*M[i][k];assert(a%d===0n);M[j][k]=a/d;}d=M[i][i];}return M.at(-1).at(-1);}assert.equal(det(G),21870n);
const upper=new Int16Array(N),cache=zlib.gunzipSync(fs.readFileSync(__dirname+'/witness.gz'));assert.equal(cache.length,(N-1)*n);let max=0;for(let i=1;i<N;i++){let w=reps[i].map((x,j)=>x+30*cache.readInt8((i-1)*n+j)),v=0;for(let j=0;j<n;j++){let t=0;for(let k=0;k<n;k++)t+=G[j][k]*w[k];v+=t*w[j];}assert(v%30===0&&v>0&&v<=173*30);max=Math.max(max,v/30);upper[i]=180-v/30;}assert.equal(max,171);
const n0=2187;let differenceTable=new Int16Array(n0*n0);for(let i=0;i<n0;i++)for(let j=0;j<n0;j++)differenceTable[i*n0+j]=3*differenceTable[(i/3|0)*n0+(j/3|0)]+(3+i%3-j%3)%3;
function actualClique(D,b,target,limit=2e6){let n=2187,DL=differenceTable;let ids=[];for(let i=1;i<D.length;i++)if(D[i]<=b)ids.push(i);let N=ids.length,B=Math.ceil(N/32),A=Array.from({length:N},()=>new Uint32Array(B));for(let i=0;i<N;i++){let a=ids[i],ay=a/10|0,az=a%10;for(let j=0;j<i;j++){let v=ids[j],u=10*DL[ay*n+(v/10|0)]+(10+az-v%10)%10;if(D[u]<=b){A[i][j>>5]|=1<<(j%32);A[j][i>>5]|=1<<(i%32);}}}
let nodes=0,max=1,wit=null,mainColors=0;function pop(P,p=0){for(let k=p;k<B;k++)if(P[k])return(k<<5)+31-Math.clz32((P[k]&-P[k])>>>0);return -1;}function rec(P,depth,path){if(++nodes>limit)throw 'cap';if(depth>max){max=depth;wit=[...path];}if(depth>=target)throw 'found';let S=[],cols=[],U=P.slice(),col=0;while(pop(U)>=0){let Q=U.slice();col++;let v,i=0;while((v=pop(Q,i))>=0){i=v>>5;S.push(v);cols.push(col);U[i]&=~(1<<(v%32));Q[i]&=~(1<<(v%32));for(let j=i;j<B;j++)Q[j]&=~A[v][j];}}if(depth==1)mainColors=col;for(let i=S.length-1;i>=0;i--){if(depth+cols[i]<target)return;let v=S[i],T=P.map((z,k)=>z&A[v][k]);rec(T,depth+1,[...path,ids[v]]);P[v>>5]&=~(1<<(v%32));}}
let P=new Uint32Array(B);for(let i=0;i<N;i++)P[i>>5]|=1<<(i%32);let status='excluded';try{rec(P,1,[0])}catch(e){status=e;}return {status,nodes,max,wit,mainColors,N};}

// q=180−30h² is a necessary deficit obtained from the shorter witnessed vector. The smallest encountered q of a difference constrains nothing: a motif's *largest* q gives the guaranteed short pair and its diameter in this search graph.
const ds=[...new Set(upper)].filter(x=>x>0).sort((a,b)=>a-b);
let opt={};for(let b of ds){let goal=1;let lo=0,hi=N;while(lo<hi){let mid=(lo+hi+1)>>1;if(BigInt(mid)**2n*BigInt(180-b)**39n<=81n*160n**39n)lo=mid;else hi=mid-1;}let ub=lo;if(ub>=N)continue;
// If a larger cutoff is cheaper it is still valid, but nothing here assumes unproved inequalities or heuristic clique output. No node cap is accepted as an exclusion.
let x=actualClique(upper,b,ub+1,10000000);assert.equal(x.status,'excluded');opt[b]=[ub,x.nodes,x.N];console.log('deficit',b,'cap',ub,'nodes',x.nodes,'neighbors',x.N);}
console.log('True-metric bound: every non-singleton is <= 9*(8/9)^(39/2) times the singleton.');

