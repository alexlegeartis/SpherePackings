// Exact compressed replay for the 11D panel schedule.
require('./verify_phase_slip.js');
const K=1182n,M=64n,E=K**4n,B=M*K, ck=x=>{if(!x)throw Error("check");};
const T=[0n,...Array.from({length:11},(_,i)=>2n**BigInt(i)),K];
let H=0n;
for(let i=0;i+1<T.length;i++)for(let j=0n;j<M;j++){
  const w=T[i+1]-T[i], a=M*T[i]+j*w;
  const p=4n*w*B**3n, d=M**4n+a**4n, q=(p+d-1n)/d;
  ck((q-1n)*d<p && p<=q*d && 2n*w<B*q);
  H+=q;
}
const N=2n*H,F=72n*(N-1n),G=8n*E-11n*(F+1n),D=8n*E*F;
ck(N===14784265480n && G===3906505511629n && D===16622339265797396791163904n);
ck(G*10n**16n>2350n*D);
const O=2n**44n,L=72n*(147256021568n-1n);
ck((F+1n)**2n*(4n*E-1n)**11n*L**2n*(4n*O)**11n >
   (L+1n)**2n*(4n*O-1n)**11n*F**2n*(4n*E)**11n);
console.log('PASS: exact panel schedule and strict improvement over the earlier witness',String(N),String(2n*(N-1n)));
