const result=(function(){
 const Gol=(()=>{const a=[];for(let m=0;m<4096;m++){let h=0;for(let j=0;j<12;j++)if((m>>j)&1)h^=0xae3<<j;let b=0;for(let s=h;s;s&=s-1)b^=1;h|=b<<23;a.push(h)}return a})();const X=new Int8Array(196560*24), arr=[], pow=Array.from({length:24},(_,j)=>1<<j), ksf=new Uint8Array(1<<12);for(let i=1;i<ksf.length;i++)ksf[i]=1+ksf[i&(i-1)];
 function pop(b){let k=0;while(b){k++;b&=b-1}return k}
 let z=0;
 for(let i=0;i<24;i++)for(let j=0;j<i;j++)for(let u of [-4,4])for(let v of [-4,4]){X[24*z+i]=u;X[24*z+j]=v;z++}
 for(let c of Gol){if(pop(c)!==8)continue;let cols=pow.map((p,j)=>[p,j]).filter(x=>x[0]&c).map(x=>x[1]);for(let s=0;s<256;s++){if(pop(s)%2)continue;for(let j=0;j<8;j++) X[24*z+cols[j]]=2*(1-2*((s>>j)&1));z++}}
 for(let c of Gol)for(let j=0;j<24;j++){for(let i=0;i<24;i++)X[24*z+i]=(i===j?-3:1)*(1-2*((c>>i)&1));z++}
 if(z!==196560)throw Error('shell');
 const prime=101, n=275, B=new Array(n),indices=[];
function powmod(a,b){let s=1;while(b){if(b%2)s=s*a%prime;a=a*a%prime;b>>=1}return s}
let rank=0,ct=0; const patterns=[[2,0,0],[0,2,0],[0,0,2]];
let stats=[];
for(let pp=0;pp<3;pp++){
 const M=Array.from({length:22},()=>new Int32Array(22)),sum=new Int32Array(22);let contacts=0;
 let p=patterns[pp];
 for(let x=0;x<z;x++){let j=(x*12269)%z,i=24*j;
 if(X[i]-X[i+1]!==p[0]-p[1]||X[i+1]-X[i+2]!==p[1]-p[2])continue;contacts++;
 let u=[X[i]+X[i+1]+X[i+2],...X.slice(i+3,i+24)];
 for(let k=0;k<22;k++){sum[k]+=u[k];for(let l=0;l<22;l++)M[k][l]+=u[k]*u[l]}
 if(pp||rank===n)continue;
 let f=Array.from(u);for(let k=0;k<22;k++)for(let l=0;l<=k;l++)f.push(u[k]*u[l]);f=Int32Array.from(f.map(v=>(v%prime+prime)%prime));
 let inserted=false;
 for(let k=0;k<n;k++)if(f[k]){
  if(!B[k]){const a=powmod(f[k],prime-2);for(let l=k;l<n;l++)f[l]=f[l]*a%prime;B[k]=f;rank++;indices.push(j);inserted=true;break}
  const t=f[k],b=B[k];for(let l=k;l<n;l++)f[l]=(f[l]-t*b[l]+prime*prime)%prime
 }
 }
 for(let k=0;k<22;k++)for(let l=0;l<22;l++)if(M[k][l]!==((k===l)?(k===0?82944:27648):0))throw Error('covariance');
 if(sum.some(v=>v!==0)||contacts!==20736)throw Error('centering');
 stats.push(contacts);
}
if(rank!==n)throw Error('rank '+rank);
return [rank,stats,indices];})(); console.log(result[0],result[1]);
