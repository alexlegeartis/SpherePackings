// Sphere100, November 26, 2026. Standalone ES2020 / Node; no dependencies.
// Exact no strict improvement for the two-weight simplex-coset family described in README.md.

const pop=(n)=>{let x=n-((n>>1)&0x55555555);x=(x&0x33333333)+((x>>2)&0x33333333);return ((x+(x>>4)&0xf0f0f0f)*0x1010101)>>>24};
const pairs=[[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]];
const evalq=(q,a)=>{let w=0;for(let k=1;k<16;k++){let c=pop(a&k)&1;for(let j=0;j<6;j++)if((q>>j&1)&&(k>>pairs[j][0]&1)&&(k>>pairs[j][1]&1))c^=1;w|=c<<(k-1);}return w;};

const simplex=Array.from({length:16},(_,a)=>evalq(0,a)), seen=new Uint8Array(32768), reps=[],rowId=new Uint16Array(32768);
for(let w=0;w<32768;w++)if(!seen[w]){let id=reps.length;reps.push(w);for(let s of simplex){seen[w^s]=1;rowId[w^s]=id;}} 
for(const u of reps)for(const v of reps)if(reps[rowId[u^v]]!==(u^v))throw Error('representative set not closed');
const cols=Array.from({length:15},(_,i)=>i+1);const perm=(M,i,j)=>{let p=0;for(let c=1;c<=15;c++){let d=c^((c>>i&1)<<j);if(M&(1<<(c-1)))p|=1<<(d-1);}return p;};
let orbitReps=[],visited=new Uint8Array(32768);for(let M=0;M<32768;M++)if(!visited[M]){let all=[M];visited[M]=1;for(let o=0;o<all.length;o++){let U=all[o];for(let i=0;i<4;i++)for(let j=0;j<4;j++)if(i!==j){let m=perm(U,i,j);if(!visited[m]){visited[m]=1;all.push(m);}}}orbitReps.push([M,all.length]);}

const tab={6:[4,2,1],7:[8,2,1],8:[16,2,2],9:[20,4,2],10:[40,6,2],11:[72,12,2],12:[144,24,4],13:[256,32,4],14:[512,64,8],15:[1024,128,16],16:[2048,256,32]};
const A=(n,d)=>{if(d>n)return 1;if(d<=1)return 2**n;if(d==2)return 2**(n-1);if(d%2)return A(n+1,d+1);if(3*d>2*n)return 2;return tab[n][d/2-2];};


const gcd=(a,b)=>b?gcd(b,a%b):a;
const rats=[[0,1]];for(let p=1;p<=7;p++)for(let q=1;q<=7;q++)if(gcd(p,q)==1)rats.push([p,q]); // j-diffs relevant merely up to7 when m≤8 unless j=8 for d=8 then z=0 only
const collect=[];
for(let [H] of orbitReps){let r=pop(H);if(r<2||r>13)continue;
 const ds=[],js=[];for(let u of reps){for(let v of simplex){let w=u^v;ds.push(pop(w));js.push(pop(w&H));}}
 for(let [p,q] of rats){let hist={},Z=p/q;
 for(let id=1;id<2048;id++){let mn=8*q;for(let k=id*16;k<16+id*16;k++){let j=q*ds[k]+p*js[k];if(j<mn)mn=j;}hist[mn]=(hist[mn]||0)+1;}
 let accum=0;for(let [mstr,ct] of Object.entries(hist).sort((a,b)=>b[0]-a[0])){let m=+mstr;accum+=ct;let q0=Math.floor((8*q/m)**7.5*(1+Z)**(r/2)+1e-10)+1;const left=BigInt(8*q)**15n*BigInt(q+p)**BigInt(r),base=BigInt(m)**15n*BigInt(q)**BigInt(r); while(q0>1 && BigInt(q0-1)**2n*base>left)q0--;while(BigInt(q0)**2n*base<=left)q0++;if(accum+1<q0)continue;
 let U=32768;for(let h=0;h<=r;h++){let j=0;while(q*j+p*Math.min(j,r-h)<m)j++;U=Math.min(U,2**h*A(15-h,j));}if(U<16*q0)continue;
 collect.push([r,H,p,q,m,accum,q0,Math.floor(U/16)]);}
 }}


const scan=(case_,budget=5e4)=>{
const [r,H,p,q,m,N,need]=case_,t=p/q,allowed=[],ok=new Uint8Array(32768);
for(let j=1;j<2048;j++){let u=reps[j],vok=simplex.every(s=>{let v=u^s;return q*pop(v)+p*pop(v&H)>=m;});if(vok){allowed.push(u);ok[u]=1;}}
let adj=Object.fromEntries(allowed.map(u=>[u,allowed.filter(v=>ok[u^v])])),deg=Object.fromEntries(allowed.map(u=>[u,adj[u].length]));
let best=[0],nodes=0,exact=true;
const color=(list)=>{let a=[],colors=[];for(let u of list){let i=0;while(colors[i]?.some(v=>ok[u^v]))i++;(colors[i]??=[]).push(u);}let b=[];for(let k=0;k<colors.length;k++)for(let j of colors[k]){a.push(j);b.push(k+1);}return [a,b];};
const rec=(R,C)=>{if(!exact||best.length>=need)return;if(++nodes>budget){exact=false;return;}if(R.length>best.length)best=R.slice();if(C.length==0)return;let [cs,bs]=color(C);for(let j=cs.length-1;j>=0;j--){if(R.length+bs[j]<=best.length)return;let v=cs[j],cn=[];for(let i=0;i<j;i++)if(ok[v^cs[i]])cn.push(cs[i]);rec([...R,v],cn);}};
rec([0],allowed.sort((a,b)=>deg[b]-deg[a]));return {r,H,t,N,need,nodes,exact,k:best.length,ratio:best.length*(m/(8*q))**7.5/(1+t)**(r/2),best};};

const outcomes=collect.map(c=>scan(c));
if(orbitReps.length!==46 || collect.length!==58 || outcomes.some(e=>!e.exact||e.k>=e.need))throw Error("Verification failed.");
console.log({subsetOrbits:orbitReps.length, survivors:collect.length, largestRootNeighborhood:Math.max(...collect.map(c=>c[5])),largestSearchNodes:Math.max(...outcomes.map(e=>e.nodes))});
