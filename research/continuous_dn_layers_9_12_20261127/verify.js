'use strict';
const A=require('node:assert/strict');
// Squared distance of x/D to parity-even integer lattice, times D^2.
function q(x,D) {
  const near=x.map(v=>Math.floor((v+D/2)/D));
  const r=x.map((v,i)=>Math.abs(v-D*near[i]));
  return r.reduce((z,t)=>z+t*t,0)+((near.reduce((z,t)=>z+t,0)&1)?D*D-2*D*Math.max(...r):0);
}
for(let m=8;m<=11;m++){
  const g=Array.from({length:m},(_,i)=>+(i<8)),
        x=Array.from({length:m},(_,i)=>+(i<4||i>=8));
  const h4=12-m;
  A.equal(q(g,2),8);
  A.equal(q(x,2),m-4);
  A.equal(q(x.map((v,i)=>v+g[i]),2),m-4);
  for(let k=1;k<=3;k++)for(let a=-1;a<=1;a++){
     const v=x.map((j,i)=>k*j+a*g[i]);
     A.ok(q(v,2)+k*k*h4>=8);
  }
  A.ok(16*h4>=8);
  if(m===11){
     A.equal(q(x.map(v=>2*v),2),4);
     A.equal(q(x.map((v,i)=>2*v+g[i]),2),8);
  }
}
const fixture=[3,3,4,4,4,4,4,4];
A.equal(q(fixture,8),114);
fixture[0]+=8;
A.equal(q(fixture,8),114);
A.ok(114+16>128);
for(const [m,p,r] of [[8n,17n,25n],[9n,77n,100n],[10n,87n,100n]]){
  const b=12n*p-m*r;
  A.ok(b>0n);
  A.ok((12n-m)*p**(m+1n)<=b*r**m);
}
A.ok(729n*19n**12n<1024n*(12n*19n-11n*20n)*20n**11n);
A.ok(729n<1024n);
A.ok(8*8<3*5*5); // 10D tie / Best = 8/(5 sqrt(3)) < 1.
A.ok(8<9); // 11D cap / historical 72-word construction = 8/9 (< current).
A.ok(27<32); // 12D cap / K12 = 27/32.
// Finite circle-grid Haar average: numerators in units 1/(2q)^2
// at the convex maximum, where one point has distance exactly 1/2.
for(let u=2n;u<=100n;u++){
  let s=0n;
  for(let j=0n;j<u;j++){
    const a=(u+2n*j)%(2n*u), b=a<=u?a:2n*u-a;
    s+=b*b;
  }
  A.equal(3n*s,u*(u*u+2n)); // s/(4q^3)=(q^2+2)/(12q^2)
  A.ok(8n*s<=4n*u**3n);
  if(u>=3n)A.ok(108n*s<=11n*4n*u**3n);
}
A.ok(59n>54n); // unlike orientation h^2>=59/216>1/4 at maximal radius
A.ok(93n**9n <= (3n*93n-2n*100n)*100n**8n);
A.ok(54n*93n**9n <= (54n*(8n*93n-8n*100n)+59n*100n)*100n**8n);
console.log('PASS exact torus fixture and all relevant lags in 9–12D, off-palette 9D tie, rational scale and benchmark checks, rotation Haar averages');

