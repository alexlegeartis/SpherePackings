'use strict';
function verify(c){
const F=[[2,-1,-1,0,-1,0,0],[0,3,-1,0,-1,0,0],[0,0,4,-3,-2,0,0],[0,0,0,5,-2,0,0],[0,0,0,0,4,-5,0],[0,0,0,0,0,3,-4],[0,0,0,0,0,0,1]],D=[30,10,5,3,3,5,40],G=c.E7_gram,M=c.phase_isometry,A=c.alpha.slice(1);
const range=n=>Array.from({length:n},(_,i)=>i),bits=b=>range(7).map(i=>b>>i&1),chk=v=>{if(!v)throw Error("certificate verification failed")};
for(let i=0;i<7;i++)for(let j=0;j<7;j++)chk(F.reduce((s,r,k)=>s+r[i]*D[k]*r[j],0)===60*G[i][j]);
const r4=new Uint8Array(16384).fill(99),v=Array(7).fill(0);
function en(i,b){if(i<0){let k=0;for(let j=0;j<7;j++)k|=(v[j]&3)<<(2*j);r4[k]=Math.min(r4[k],(1440-b)/60);return;}let off=0;for(let j=i+1;j<7;j++)off+=F[i][j]*v[j];let q=0;while((q+1)**2*D[i]<=b)q++;for(let t=Math.ceil((-q-off)/F[i][i]);t<=Math.floor((q-off)/F[i][i]);t++){v[i]=t;en(i-1,b-D[i]*(F[i][i]*t+off)**2);}}
en(6,1440);chk(r4.every(a=>a<=24));chk(new Set(M).size===128&&M[64]===82);for(let x=0;x<128;x++)for(let y=0;y<128;y++)chk(M[x^y]===(M[x]^M[y]));
const gf=(x,y)=>{let z=0;while(y){if(y&1)z^=x;x*=2;if(x&64)x^=c.field_polynomial;y>>=1;}return z};
const sy=(j,t)=>(t&63)+(gf(A[j],t&63)<<6)+4096*(t>>6), tag=x=>c.phase_syndromes.reduce((s,t,i)=>x>>i&1?s^t:s,0);
const T=c.half_shift_blocks.map(k=>82*k);
function ext(b){let d=new Uint8Array(8192).fill(99);d[0]=0;
for(let j=0;j<4;j++){const a=new Uint8Array(8192).fill(99);for(let t=0;t<128;t++){let k=0;for(let i=0;i<7;i++)k|=((b[j]>>i&1)+2*(M[t]>>i&1))<<(2*i);let z=r4[k],s=sy(j,t);for(let u=0;u<8192;u++)a[u]=Math.min(a[u],d[s^u]+z);}d=a;}return d;}
function bank(q,b,u,m){let of=u^(m?c.offset_syndrome:0);let bs=b.map((v,j)=>{if(m)of^=sy(j,M.indexOf(T[j]&~v));return v^(m?T[j]:0)}),e=ext(bs),best=1000;
for(let k=0;k<32;k++){const wt=n=>bits(n).reduce((s,a)=>s+a,0);const p=k|((wt(k)&1)<<5),g=wt(p)/2%2;let a=0,b=10000;
for(let i=0;i<6;i++){let s=q[i]+m-2*(p>>i&1);let aa=((s+100)&7)-4,bb=((s+96)&7)-4;let x=aa*aa,y=bb*bb;[a,b]=[Math.min(a+x,b+y),Math.min(b+x,a+y)];}
best=Math.min(best,(g?b:a)+e[of^tag(k)],(g?a:b)+e[of^tag(k+32)]);}
return best;}
let q=c.new_local_shift_times_2,b=c.new_half_masks,u=c.new_integral_syndrome;
const d1=[0,1].map(m=>bank(q,b,u,m));const p=q.slice(0,5).reduce((s,v,i)=>s+((v&1)<<i),0)+32*(q.reduce((s,v)=>s+v,0)/2%2);const s=b.reduce((v,t,j)=>v^sy(j,M.indexOf(t)),0),u2=s^tag(p),d2=[0,1].map(m=>bank([0,0,0,0,0,0],[0,0,0,0],u2,m));
chk(d1[0]===26&&d1[1]===26);chk(s===919&&tag(p)===88&&u2===975);chk(d2[0]===8&&d2[1]===24);chk(c.new_layer_height_squared+''==='3,2');return {first_squared_times_four:d1,second_squared_times_four:d2,second_syndrome:u2,minimum_squared:8};}
if(typeof module!=="undefined"){module.exports=verify;if(require.main===module)console.log(verify(require("./cert.json")));}
