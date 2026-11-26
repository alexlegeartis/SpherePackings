// Standalone ES2020 / Node, no dependencies. Sphere100 punctured Nordstrom–Robinson metric certificate.

const pop=(n)=>{let x=n-((n>>1)&0x55555555);x=(x&0x33333333)+((x>>2)&0x33333333);return ((x+(x>>4)&0xf0f0f0f)*0x1010101)>>>24};
const pf=s=>(((s>>0)&(s>>5))^((s>>1)&(s>>4))^((s>>2)&(s>>3)))&1;
const find=(A,start)=>{if(A.length==8)return A;for(let t=start;t<64;t++)if(A.every(v=>pf(v^t))){let r=find([...A,t],t+1);if(r)return r;}};const Q=find([0],1);
const pairs=[[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]];
const evalq=(q,a)=>{let w=0;for(let k=1;k<16;k++){let c=pop(a&k)&1;for(let j=0;j<6;j++)if((q>>j&1)&&(k>>pairs[j][0]&1)&&(k>>pairs[j][1]&1))c^=1;w|=c<<(k-1);}return w;};


const W=Q.flatMap(q=>Array.from({length:16},(_,a)=>evalq(q,a))), Cset=new Set(W);
const reduced=u=>Math.min(...W.slice(0,16).map(v=>v^u)), setWS=new Set(Q.map(q=>reduced(evalq(q,0))));
const moves=[];let nAut=0;
for(let a=1;a<16;a++)for(let b=1;b<16;b++)if(a!==b)for(let c=1;c<16;c++)if(![a,b,a^b].includes(c))for(let d=1;d<16;d++)if(![a,b,a^b,c,c^a,c^b,c^a^b].includes(d)){
 const perm=Array.from({length:15},(_,i)=>{let v=i+1;return (((v&1)?a:0)^((v&2)?b:0)^((v&4)?c:0)^((v&8)?d:0))-1;});
 const apply=w=>{let z=0;while(w){let u=w&-w;w^=u;z|=1<<perm[31-Math.clz32(u)];}return z;};
 const W2=Q.map(q=>reduced(apply(evalq(q,0))));
 if([...setWS].some(v=>W2.every(w=>setWS.has(w^v)))){nAut++;if(nAut%91===31)moves.push(perm);}
}
const maskreps=[],seenM=new Uint8Array(32768);
for(let m=0;m<32768;m++)if(!seenM[m]){let data=[m];seenM[m]=1;for(let k=0;k<data.length;k++)for(let perm of moves){let z=0,w=data[k];while(w){let u=w&-w;w^=u;z|=1<<perm[31-Math.clz32(u)];}if(!seenM[z]){seenM[z]=1;data.push(z);}}maskreps.push([m,data.length]);} 

const tab={6:[4,2,1],7:[8,2,1],8:[16,2,2],9:[20,4,2],10:[40,6,2],11:[72,12,2],12:[144,24,4],13:[256,32,4],14:[512,64,8],15:[1024,128,16],16:[2048,256,32]};
const A=(n,d)=>{if(d>n)return 1;if(d<=1)return 2**n;if(d==2)return 2**(n-1);if(d%2)return A(n+1,d+1);if(3*d>2*n)return 2;return tab[n][d/2-2];};
const combinBound=(c)=>{
let [r,H,t,N,need]=c,u=Math.min(16*(N+1),32768);
for(let h=0;h<=r;h++){let j=0;while(j+t*Math.min(j,r-h)<8-1e-10)j++;u=Math.min(u,2**h*A(15-h,j));}return Math.floor(u/16);};

const WW=Q.flatMap(q=>Array.from({length:32},(_,a)=>evalq(q,a%16)^(a>=16?32767:0)));
const gcd=(a,b)=>b?gcd(b,a%b):a,params=[]; const pq=[[0,1]];for(let p=1;p<=3;p++)for(let q=1;q<=7;q++)if(gcd(p,q)==1)pq.push([p,q]);
for(let [H] of maskreps){let r=pop(H);if(r<2||r>13)continue;
for(let [p,q] of pq){let D=new Set([8*q]);for(let d=5;d<=7;d++)for(let e=0;e<=Math.min(r,d);e++)if(d*q+e*p<=8*q)D.add(d*q+e*p);
for(let m of D){let cap=256,need=Math.floor(16*(8*q/m)**7.5*(1+p/q)**(r/2)+1e-9)+1;
const left=256n*BigInt(8*q)**15n*BigInt(q+p)**BigInt(r),base=BigInt(m)**15n*BigInt(q)**BigInt(r);while(need>1&&BigInt(need-1)**2n*base>left)need--;while(BigInt(need)**2n*base<=left)need++;
if(need>cap)continue;
for(let h=0;h<=r;h++){let j=0;while(q*j+p*Math.min(j,r-h)<m)j++;cap=Math.min(cap,2**h*A(15-h,j));}
if(need<=cap)params.push([H,p,q,m,need,cap]);}}}

function certify(H,p,q,m,need){const can=Array.from({length:256},(_,i)=>Array.from({length:256},(_,j)=>{let w=WW[i]^WW[j];return q*pop(w)+p*pop(w&H)>=m;}));
let verts=Array.from({length:256},(_,i)=>i).sort((a,b)=>can[b].reduce((a,v)=>a+v,0)-can[a].reduce((a,v)=>a+v,0));
let nodes=0,cap=0,best=[],ok=true;const rec=(R,C)=>{if(!ok||best.length>=need)return;if(++nodes>50000){ok=false;return;}if(R.length>best.length)best=R.slice();let colors=[];for(let u of C){let i=0;while(colors[i]&&colors[i].some(v=>can[u][v]))i++;(colors[i]??=[]).push(u);}let cs=[],bs=[];for(let k=0;k<colors.length;k++)for(let u of colors[k]){cs.push(u);bs.push(k+1);}if(!R.length)cap=colors.length;for(let i=cs.length-1;i>=0;i--){if(bs[i]+R.length<need)return;let v=cs[i],cs2=[];for(let j=0;j<i;j++)if(can[v][cs[j]])cs2.push(cs[j]);rec([...R,v],cs2);}};
rec([],verts);return {H,p,q,m,need,ok,nodes,best:best.length,initialColor:cap};};

const results=params.map(([H,p,q,m,need])=>certify(H,p,q,m,need));if(nAut!==20160||maskreps.length!==46||params.length!==307||results.some(v=>!v.ok||v.best>=v.need))throw Error("Verification failed.");console.log({nAut,orbitCount:maskreps.length,cases:params.length,maxNodes:Math.max(...results.map(v=>v.nodes))});