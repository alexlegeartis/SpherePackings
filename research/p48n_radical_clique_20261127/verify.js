"use strict";
function verify(d) {
 const {G,stocks,norms,vertices,pairs,pivots}=d;
 const A=(t,m)=>{if(!t)throw Error(m)};
 const tr=M=>M[0].map((_,i)=>M.map(row=>row[i]));
 const range=n=>Array.from({length:n},(_,i)=>i);
 const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
 const mul=(M,N)=>{const Nt=tr(N);return M.map(r=>Nt.map(s=>dot(r,s)));};
 const gram=(S,V)=>mul(mul(S,G),tr(V));
 const bil=(u,M,v)=>dot(u,M.map(r=>dot(r,v)));
 function det(M) {
   const n=M.length,a=M.map(r=>r.map(BigInt));let b=1n,sgn=1n;
   for(let i=0;i<n-1;i++){
    let p=a.findIndex((r,j)=>j>=i&&r[i]!==0n);if(p<0)return 0n;
    if(p!==i){[a[p],a[i]]=[a[i],a[p]];sgn=-sgn;}
    let c=a[i][i];
    for(let j=i+1;j<n;j++)for(let k=i+1;k<n;k++){let q=c*a[j][k]-a[j][i]*a[i][k];A(q%b===0n,"Bareiss exact division");a[j][k]=q/b;}
    for(let j=i+1;j<n;j++)a[j][i]=0n;
    b=c;
   }
   return sgn*a[n-1][n-1];
 }
 A(G.length===48&&G.every((r,i)=>r.length===48&&r.every((x,j)=>Number.isInteger(x)&&x===G[j][i])),"Gram shape");
 // Sylvester also permits a source-independent check of positive definiteness.
 for(let j=1;j<=48;j++)A(det(G.slice(0,j).map(r=>r.slice(0,j)))>0n,"Gram positive");
 A(det(G)===1n,"unimodular source");
 const H=stocks.map(S=>gram(S,S));
 A(stocks.every(S=>S.length===24&&S.every(r=>r.length===48&&r.every(Number.isInteger))),"stock dimensions");
 A(H.every(M=>det(M)===1n<<24n&&M.every((r,i)=>r.every(x=>x%2===0)&&r[i]%4===0)),"scaled unimodular even stocks");
 A(gram(stocks[0],stocks[1]).flat().every(x=>x===0),"orthogonality");
 A(pivots.length===24&&new Set(pivots).size===24&&pivots.every(j=>Number.isInteger(j)&&j>=0&&j<48),"pivots");
 const bit=(u,j)=>(u>>>j)&1;
 const label=r=>pivots.reduce((s,v,i)=>s|((r[v]&1)<<i),0);
 const B=stocks.map(S=>S.map(label));
 function span(rows){
  let V=rows.map((x,i)=>[x,1<<i]),out=[];
  for(let j=0;j<24;j++){
   let k=out.length,i=V.findIndex((p,i)=>i>=k&&bit(p[0],j));if(i<0)continue;
   [V[k],V[i]]=[V[i],V[k]];
   for(let t=0;t<V.length;t++)if(t!==k&&bit(V[t][0],j)){V[t][0]^=V[k][0];V[t][1]^=V[k][1];}
   out.push(j);
  }
  let inv=Array(24).fill(0);
  V.slice(0,out.length).forEach((r,i)=>{inv[out[i]]=r[1]});
  return [out.length,inv];
 }
 let [r0,I0]=span(B[0]),[r1,I1]=span(B[1]);
 A(r0===24&&r1===24,"primitive ranks");
 function lift(l,I){let z=0;for(let i=0;i<24;i++)if(bit(l,i))z^=I[i];return range(24).map(i=>bit(z,i));}
 const ambient=P=>tr(P);
 const coord=(u,S)=>ambient(S).map(r=>dot(r,u));
 for(let i=0;i<24;i++){
   const v=coord(lift(B[1][i],I0),stocks[0]);
   A(v.every((x,j)=>(x-stocks[1][i][j])%2===0),"common full parity space");
 }
 const xorb=(x,L)=>x.reduce((y,k,i)=>y^(k&1?L[i]:0),0);
 const lab=(x,i)=>xorb(x,B[i]);
 const rad=norms.map((P,i)=>{
  const U=H[i];
  A(P.length===2&&P.every(z=>bil(z,U,z)===8)&&Math.abs(bil(P[0],U,P[1]))===4,"normal sqrt(2)A2 in Leech normalization");
  const a=lab(P[0],i),b=lab(P[1],i);A(a&&b&&a!==b,"normal mod 2");
  return [0,a,b,a^b];
 });
 A(rad.flat().slice(1).filter(m=>m!==0).every((m,i,a)=>a.indexOf(m)===i),"disjoint radical planes");
 const field=(u,v)=>((bil(lift(u,I0),H[0],lift(v,I0))/2)&1);
 const C=rad[0].slice(1,3).map(a=>rad[1].slice(1,3).map(b=>field(a,b)));
 A(((C[0][0]*C[1][1])^(C[0][1]*C[1][0]))===1,"nondegenerate cross pairing");
 const rootNorm=(x,i)=>bil(x,H[i],x)===8&&norms[i].every(y=>bil(x,H[i],y)===0);
 A(vertices.length===6&&new Set(vertices.map(String)).size===6&&String(vertices[0])==="0,0","six vertices");
 const radios=new Set(rad[0].flatMap(a=>rad[1].map(b=>a^b)));
 for(const [u,v] of vertices){
  A(rad[0].every(r=>field(u,r)===0)&&rad[1].every(r=>field(v,r)===0),"shortened quotient membership");
  A(radios.has(u^v),"radical enlargement membership");
 }
 const done=new Set;
 A(pairs.length===15,"all pairs");
 for(const {i,j,x,y} of pairs){
  A(i>j&&i<6&&j>=0&&!done.has(6*i+j),"edge index");done.add(6*i+j);
  A(rootNorm(x,0)&&rootNorm(y,1),"actual norm-8 representatives in both shortened factors");
  A(lab(x,0)===(vertices[i][0]^vertices[j][0])&&lab(y,1)===(vertices[i][1]^vertices[j][1]),"root class of edge");
 }
 // Cap = (3/2)^22 / 18; CHLWW center density squared = 25 (157/110)^44 / 605.
 A(3n**44n*605n*110n**44n < 18n**2n*2n**44n*25n*157n**44n,"ceiling below comparator");
 A(3n**44n*605n*110n**44n > 12n**2n*2n**44n*25n*157n**44n,"four radical banks would improve");
 return "PASS: exact primitive eigensubstock bases, nondegenerate radical pair, and all 15 forbidden edges of one six-clique";
}
if(typeof module!=="undefined") module.exports=verify;
if(typeof require!=="undefined"&&require.main===module) console.log(verify(require("./witness.json")));

