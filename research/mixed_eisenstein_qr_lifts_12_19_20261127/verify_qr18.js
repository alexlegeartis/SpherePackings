
const M=1<<18, range=n=>Array.from({length:n},(_,i)=>i),pc=x=>{let a=0;for(;x;x&=x-1)a++;return a;};
const R=new Set(range(8).map(i=>(i+1)**2%17)),P=new Set(range(17).filter(x=>x&&!R.has(x))); 
function generate(S){
 const gens=range(17).map(i=>range(17).reduce((s,j)=>s|(S.has((j-i+17)%17)?1<<j:0),0)).map(x=>x|((pc(x)%2)<<17));
 let C=[0],E=new Set(C);
 for(const x of gens)if(!E.has(x)){C=[...C,...C.map(y=>y^x)];E=new Set(C);}
 return C;
}
R.add(0);P.add(0);
const C=generate(R),D=generate(P),eq=C.filter(x=>D.includes(x)),S=[...new Set([...C,...D])].filter(x=>x&&x!==(M-1));
function distances(){
 const dis=new Uint8Array(M);dis.fill(99);const que=new Int32Array(M);
 let top=0;for(const x of S){dis[x]=0;que[top++]=x;}
 for(let i=0;i<top;i++)for(let j=0;j<18;j++){let x=que[i]^(1<<j);if(dis[x]===99){dis[x]=dis[que[i]]+1;que[top++]=x;}}
 return dis;
}
const dist=distances();

const A=(c,m)=>{if(!c)throw Error(m);};
A(C.length===512&&D.length===512&&[C,D].every(c=>Math.min(...c.slice(1).map(pc))===6),"QR parameters");
A(String(eq)===String([0,M-1]),"intersection");
A(C.every(x=>D.every(y=>pc(x&y)%2===0)),"dual pair");
A(C.every(x=>!(pc(x)%2))&&D.every(x=>!(pc(x)%2)),"even");
const freq=Array(7).fill(0);
for(let r=0;r<M;r++){
 const d=dist[r];freq[d]++;
 if(d<=3||pc(r)<=2||pc(r)>=16)continue;
 A(d===4,"all nonexceptional far classes have radius four");
 let n=0;const a=Array(18).fill(0);
 for(const x of S){const u=x^r;if(pc(u)!==4)continue;n++;for(let j=0;j<18;j++)a[j]+=(u>>>j)&1;}
 A(n>0&&a.every(x=>28*x>=5*n),"uniform marginal at least 5/28");
}
A(String(freq)==="1020,17748,108936,113288,21114,36,2","coverage tally");
// All density comparisons are squared. Comparator δ0² = 13^19/(3^18 · 5^21).
const Q=13n**19n, den=3n**18n*5n**21n;
A(3n**18n*48n*den < Q*4n**18n*51n,"19D parallel-row density below incumbent");
const rows=[
 [1n,1n,2n**20n],[2n,1n,2n**20n],[3n,1n,2n**20n],
 [4n,1n,27594n],[19n,4n,4558n],[5n,1n,4558n],
 [6n,1n,1749n],[27n,4n,512n],[7n,1n,256n],[8n,1n,128n]
];
for(const [a,b,m] of rows)
 A(a**19n*m*m*den < Q*b**19n*2n**78n,"full-row family center density");
console.log("PASS: QR pair, exact Hamming rounding certificate for all 262144 binary words, and all 19D density comparisons.");

