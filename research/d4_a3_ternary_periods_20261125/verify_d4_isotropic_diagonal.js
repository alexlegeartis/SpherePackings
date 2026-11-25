const V=Array.from({length:81},(_,i)=>Array.from({length:4},(_,j)=>Math.floor(i/3**j)%3));
const d=(x,y)=>[0,2,1,3,2][x.reduce((z,t,i)=>z+ +(t!==y[i]),0)];
const stats={};let tests=0;
for(let i=0;i<81;i++)for(let j=i+1;j<81;j++){const x=V[i],y=V[j],t=d(x,y);if(t===1)continue;
 const C=V.filter(z=>d(x,z)===1&&d(y,z)===1);
 if(C.length!==6)throw Error("common count");
 const hist=C.map(a=>C.reduce((s,b)=>s+ +(d(a,b)===1),0));
 let triangles=0;for(let a=0;a<6;a++)for(let b=a+1;b<6;b++)for(let c=b+1;c<6;c++)
  if(d(C[a],C[b])===1&&d(C[a],C[c])===1&&d(C[b],C[c])===1)triangles++;
 if(t===2){
  if(hist.some(x=>x!==4)||triangles!==8)throw Error("not K_222");
  for(const a of C){const opp=C.filter(b=>a!==b&&d(a,b)!==1);
   if(opp.length!==1||d(a,opp[0])!==2||a.some((t,j)=>(t+opp[0][j]-x[j]-y[j])%3!==0))throw Error("not a symmetric opposite pair");
  }
 }else if(hist.some(x=>x!==3)||triangles!==0)throw Error("not K_33");
 stats[t]=(stats[t]??0)+1;tests++;
}
console.log(JSON.stringify({tests,stats}));
