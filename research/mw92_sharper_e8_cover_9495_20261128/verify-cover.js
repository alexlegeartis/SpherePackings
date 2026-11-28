'use strict';
const assert=(x,m)=>{if(!x)throw Error(m||'assertion');}
const n=8,r=Array.from({length:n},(_,i)=>i),b=(x)=>BigInt(x),F=b(1000000000);
const K=r.map(i=>r.map(j=>b(i===0&&j===0?797943325:i===0||j===0?1030483:i===j?1120006184:-52646518)));
const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0n);
const fx=(v)=>dot(v,K.map(w=>dot(w,v)));
const abs=(a)=>a<0n?-a:a, gcd=(a,b)=>b===0n?abs(a):gcd(b,a%b);
const norm=(y,d)=>{const g=y.reduce((s,t)=>gcd(s,t),d);return [y.map(t=>t/g),d/g];};
const mask=(h)=>1n<<b(h);
const pop=(s)=>{let k=0;for(;s;s&=s-1n)k++;return k;};
const A=K[0][0],B=K[1][1];
let P=[];
for(let a of [-1n,1n])for(let i=0;i<=7;i++) {
 let y=r.map(j=>j===0?4n*a*A:4n*(j<=i?-B:B));
 let m=mask(a===1n?0:1);
 for(let h=1;h<=6;h++)if(y[h]===y[h+1])m|=mask(h+1);
 if(y[1]===-4n*B)m|=mask(8);
 if(y[7]===4n*B)m|=mask(9);
 P.push([y,1n,m]);
}
const V=[];
function vectors(j,arr,n2,s){
 if(j===8){if(n2&&s%4===0){let a=arr.map(b);V.push([a,fx(a),n2]);}return;}
 for(let i=-5;i<=5;i++)if((j===0||(i-arr[0])%2===0)&&(j<2||i>=arr[j-1])&&n2+i*i<=32)vectors(j+1,arr.concat(i),n2+i*i,s+i);
}
vectors(0,[],0,0);
V.sort((a,b)=>Number(a[1]-b[1]));
let mx=P.length;
for(let i=0;i<V.length;i++){
 const [v,c]=V[i],k=mask(i+10),ps=P.map(([y,d,m])=>[y,d,m,dot(v,y)-c*d]), pos=ps.filter(x=>x[3]>0n);
 let neg=ps.filter(x=>x[3]<0n),good=ps.filter(x=>x[3]<=0n).map(([y,d,m,l])=>[y,d,m|(l===0n?k:0n)]);
 let M=new Map(good.map(x=>[x[0].join(',')+'/'+x[1],x]));
 for(const [y,d,m,h] of pos)for(const [w,e,a,l] of neg){
  if(pop(m&a)<7)continue;
  const [p,q]=norm(y.map((z,j)=>h*w[j]-l*z),h*e-l*d);
  let id=p.join(',')+'/'+q, t=M.get(id),mm=(m&a)|k;
  if(t)t[2]|=mm;else M.set(id,[p,q,mm]);
 }
 P=[...M.values()];mx=Math.max(mx,P.length);
 assert(P.length>0,'nonempty; origin included');
}
const U=K[1][1]-K[1][2],V7=K[1][1]+6n*K[1][2],C=K[0][1],T=A*V7-7n*C*C;
assert(U>0n&&T>0n&&A>0n,'positive definite');
 // y^T K^-1 y = (7 sum_tail y_i^2-(sum_tail y_i)^2)/(7 U)
 //                  +(7 V7*y0^2-14 C*y0*sum_tail+ A*(sum_tail)^2)/(7 T).
const den0=112n*F*U*T;
let high=[0n,1n],top=null;
for(let [y,d] of P){
 let s=y.slice(1).reduce((a,b)=>a+b,0n),q=dot(y.slice(1),y.slice(1));
 let nom=T*(7n*q-s*s)+U*(7n*V7*y[0]**2n-14n*C*y[0]*s+A*s*s),den=den0*d*d;
 if(high[0]*den<nom*high[1]) {high=[nom,den];top=[y,d];}
}
const rr=norm([high[0]],high[1]);high=[rr[0][0],rr[1]];
const det=[U**6n*T,F**8n];
console.log({vectors:V.length,vertices:P.length,maxVertices:mx,r2:high.map(String).join('/'),r2decimal:Number(high[0])/Number(high[1]),det:det.map(String).join('/'),theta: (Number(high[0])/Number(high[1]))**4/Math.sqrt(Number(det[0])/Number(det[1]))*Math.PI**4/24});
const cover=3143n, piu=[314159265359n,100000000000n];
assert(1000n**2n*high[0]**8n*det[1]*piu[0]**8n<cover**2n*24n**2n*high[1]**8n*det[0]*piu[1]**8n,'theta<3.143');
console.log('PASS exact 8D Voronoi outer radius and covering density < 3.143');
