"use strict";
const ok=(p,s="assertion")=>{if(!p)throw new Error(s);};
const hw=x=>{let z=0;for(;x;x&=x-1)++z;return z;};
const Q=parseInt("11011100010",2);
const H=[0,...Array.from({length:11},(_,i)=>((Q<<i)|(Q>>(11-i)))&2047)];
const A=[];
for(let i=0;i<12;i++)for(let j=0;j<i;j++){const d=H[i]^H[j];A.push(d,4095^d);}
for(let i=0;i<6;i++){const d=3<<(2*i);A.push(d,4095^d);}
const C=A.filter(u=>u<2048), V=C.filter(u=>u!==1023);
ok(C.length===72&&new Set(C).size===72&&V.length===71,"shortened stock");
for(let i=0;i<V.length;i++){
  ok(hw(2047^V[i])>=5,"core Hamming clearance");
  for(let j=0;j<i;j++)ok(hw(V[i]^V[j])>=4,"core separation");
}
ok(V.filter(x=>x<1024).length>1,"at least two transverse core residues at even axial height");
function check(k,gammaDen,stepDen,expectedN,gainNumer,gainDen){
  const B=2n**BigInt(k), E=B**4n;
  // Segment phase lengths use denominator B; all step lengths use 16*E.
  const D=16n*E, hNum=D/B;
  let pos=0n,m=0n,step=4n;
  // Allow the two tiny blocks using just axial clearance.
  ok(4n*step*E<=D,"tiny segment clearance");
  ok(hNum%step===0n,"tiny segment integrality");
  pos+=hNum;m+=hNum/step;
  for(let a=1n;a<B;a*=2n){
    const len=a*hNum;step=16n*a**4n/stepDen;
    ok(pos===len && len%step===0n,"dyadic incidence and integrality");
    // Need 4*d <= (2*a_actual^2/gammaDen)^2, cross multiplied.
    ok(step*gammaDen**2n*E<=D*a**4n,"zigzag adjacent clearance");
    ok(0n<step && step<=D/64n,"nonadjacent gap safety");
    pos+=len;m+=len/step;
  }
  ok(pos===D && 2n*m===expectedN,"exact half end and cyclic count");
  ok(gammaDen>=2n,"transverse <=1/2 and 2t<=s^2");
  const N=expectedN,F=72n*(N-1n);
  ok(N%2n===0n && 2n*(N-1n)>0n,"period and sign seam");
  const num=8n*E-11n*(F+1n),den=8n*E*F;
  ok(num===gainNumer&&den===gainDen&&num>0n,"exact Bernoulli gain");
  ok((F+1n)**2n*(4n*E-1n)**11n>F**2n*(4n*E)**11n,"direct exact squared gain");
  return [N,2n*(N-1n),num,den].map(String);
}
const c=check(13,4n,16n,24503401990400n,16622102642567949n,63563622954203754089284531912704n);
const x=check(11,2n,4n,147256021568n,24110719274253n,1492159868678706575569846272n);
ok(BigInt(c[2])*10n**17n>26n*BigInt(c[3]),"conservative decimal strict floor");
ok(BigInt(x[2])*10n**15n>16n*BigInt(x[3]),"shorter period improvement floor");
console.log("PASS: compressed exact 11D phase-slip packing; [layers, period, relative gain numerator, denominator]");
console.log(c);console.log(x);
