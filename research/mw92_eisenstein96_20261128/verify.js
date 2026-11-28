// Exact replay for the 96D coupled Eisenstein tails. Node.js, no dependencies.
// Parent MW92 and the equal-ball threshold A are linked in README.md.
'use strict';
const assert=(x,m='assertion')=>{if(!x)throw Error(m);};
function gcd(a,b){while(b)[a,b]=[b,a%b];return a<0n?-a:a;}
function F(n,d=1n){n=BigInt(n);d=BigInt(d);assert(d!==0n);const g=gcd(n,d);return [d<0n?-n/g:n/g,d<0n?-d/g:d/g];}
const add=([p,q],[r,s])=>F(p*s+q*r,q*s), sub=([p,q],[r,s])=>F(p*s-q*r,q*s);
const mul=([p,q],[r,s])=>F(p*r,q*s), div=([p,q],[r,s])=>F(p*s,q*r);
const pow=([p,q],k)=>[p**k,q**k], lt=([p,q],[r,s])=>p*s<r*q;
const eq=([p,q],[r,s])=>p===r&&q===s;
const pstr=([p,q])=>p.toString()+'/'+q.toString();
const A=F(1100098526n,100000000n),s=sub(F(16),A),t=sub(F(16),mul(F(3),s));
assert(lt(F(0),t)&&lt(F(16),mul(F(4),s)));
assert(eq(sub(A,div(t,F(3))),F(32,3)));
const norm=(i,j)=>i*i-i*j+j*j;
const B=F(537,50); // 10.74
assert(lt(mul(F(3),pow(B,46n)),pow(A,46n)));
let results=[];
for(const [name,S,e,m,k] of [
 ['generic',s,F(514,125),3,27],['fixed',F(16,3),F(387,100),5,15]
]){
 const D=sub(sub(F(16),e),div(S,F(3))),W=sub(F(16),mul(F(3),e));
 assert(eq(mul(F(3),F(1,3)),F(1))); // ring norm of (2+omega)/3 is 1/3
 assert(lt(F(0),W)&&lt(W,S)); // all non-nearest norm-3 pairs are safe
 assert(lt(F(16),mul(F(5),e))); // all norm >=5
 assert(lt(F(16),add(mul(F(4),e),div(S,F(3))))); // norm 4: parity -> 2*unit -> same hole
 assert(lt(mul(F(k),pow(D,46n)),pow(A,46n)));
 // Complete N(3x+2,3y+1) check where Q(x,1)<16. The bound N(a,b)>=(a^2+b^2)/2
 // excludes |i| or |j|>=4, because then the corresponding |3i+2| or |3j+1|>=10.
 assert(lt(F(16),add(e,mul(S,F(50,9)))));
 let counts={},max=F(0),weak=F(0),fixed=new Set();
 for(let i=-3;i<=3;i++)for(let j=-3;j<=3;j++){
   let n=norm(3*i+2,3*j+1),Q=add(e,mul(S,F(n,9)));
   assert(n>=3);
   if(!lt(Q,F(16)))continue;
   let d=sub(F(16),Q);
   counts[n]=(counts[n]??0)+1;
   if(lt(max,d))max=d;
   if(n===21&&lt(weak,d))weak=d;
   if(n===3)fixed.add(((i+j)%3+3)%3);
 }
 assert(eq(max,D));
 assert(JSON.stringify(counts)===(name==='generic'?'{"3":3,"12":3,"21":6}':'{"3":3,"12":3}'));
 if(name==='generic'){
   // Midpoint dimer for each of six weakest phases: 2+sqrt(weak)<sqrt(D).
   const h=sub(sub(D,F(4)),weak);
   assert(lt(F(0),h)&&lt(mul(F(16),weak),pow(h,2n)));
   assert(k===m*(counts[3]+counts[12]+counts[21]/2));
   assert(lt(mul(F(m*m),W),mul(F(3),A)));
 } else{
   assert(fixed.size===3&&k===m*fixed.size);
   assert(lt(mul(F(m*m),W),add(mul(F(4),F(34,3)),mul(F(6),B))));
 }
 const dt=mul(F(9,16),mul(pow(S,2n),pow(e,2n)));
 results.push({name,scale:pstr(S),eta:pstr(e),D:pstr(D),W:pstr(W),cover_count:k,choices:m,shell_counts:counts,det:dt});
}
const dg=results[0].det,df=results[1].det;
const old=mul(F(4n*47n**4n),pow(F(8,9),96n));
const periodic=F(300153807922100058481n,1440000000000000000n);
assert(lt(periodic,dg)&&lt(dg,df)&&lt(df,old));
function atan(z,k){let v=F(0);for(let i=0n;i<k;i++)v=add(v,F((-1n)**i,(2n*i+1n)*z**(2n*i+1n)));return v;}
const pi=sub(mul(F(16),atan(5n,40n)),mul(F(4),atan(239n,41n)));
const fact48=Array.from({length:48},(_,i)=>BigInt(i+1)).reduce((a,b)=>a*b,1n);
const center2=div(F(2n**192n),mul(F(47n**14n),df));
const occ2=mul(center2,pow(div(pow(pi,48n),F(fact48)),2n));
function isqrt(n){let a=0n,b=1n<<BigInt(Math.ceil(n.toString(2).length/2)+1);while(b-a>1n){let c=(a+b)/2n;if(c*c<=n)a=c;else b=c;}return a;}
const floor=(v,k)=>isqrt(v[0]*10n**(2n*k)/v[1]);
const frac=(z,d)=>[z/10n**d,(z%10n**d).toString().padStart(Number(d),'0')].join('.');
const occ=floor(occ2,37n), gain=floor(div(old,df),8n)-10n**8n;
const bcsRatio2=mul(center2,F(2n**48n,3n**96n));
const bcsGain=floor(bcsRatio2,8n)-10n**8n;
assert(occ>5937109315592000n);
assert(lt(F((10n**8n+gain)**2n,10n**16n),div(old,df)));
const out={branch:results.map(x=>({...x,det:pstr(x.det)})),previous_lattice_det:pstr(old),unrestricted_det_unchanged:pstr(periodic),worst_det:pstr(df),
 occupied_density_floor:frac(occ,15n)+'e-22',previous_occupied_lattice_comparator:'5.9371093155920008e-22',gain_percent_floor:frac(gain,6n),gain_vs_BCS96_percent_floor:frac(bcsGain,6n)};
console.log(JSON.stringify(out,null,2));
