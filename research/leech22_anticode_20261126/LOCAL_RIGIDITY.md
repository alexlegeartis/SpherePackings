# Strict local optimality of the 22D Conway–Sloane 3-translate packing in the original cell

Permit arbitrary common positive-definite 22×22 Gram deformation and arbitrary independent real shifts of its three complete translates of the same lattice. Normalize squared packing minimum to μ=11/3. Modulo rigid motions this known 3-translate packing is a strict local maximum of density in those variables. The result does not allow independent site displacements in a larger supercell and is not an assertion about remote optima. No new world lower bound. The incumbent and its three normal labels are described at https://arxiv.org/pdf/math/0207182 (p. 3).

Use Leech coordinates scaled by 1/√8 and W={x₀=x₁=x₂}. The three normal label differences, up to sign, are the projections onto W^⊥ of (2,0,0), (0,2,0) and (0,0,2), in the common raw units. Indeed one may take the three mark representatives (0,0,0), (2,0,0), (0,-2,0); on quotienting the diagonal a difference (2,2,0) becomes (0,0,-2). For each oriented edge of the 3-color motif let S_ij consist of the W projections of Leech norm-four points with one of these exact normal projections. They all give actual contacts between the indicated complete translates: their squared normal norm is 1/3, leaving squared horizontal length μ. The relative placement of the anchor chosen in each translate does not change these physical contact vectors.

The verifier below certifies **separately for all three edge fibers** that |S_ij|=20,736, mean(S_ij)=0 and mean(vvᵀ)=I/6 in physical orthonormal 22-coordinates. For the first edge the matrix of 22 linear and 253 quadratic monomials evaluated on S_01 has rank 275 over F_101, hence over the reals. Computation uses coordinates u=(x0+x1+x2,x3,…,x23) instead: the physical coordinates are u0/√24 and ui/√8 for i≥1. The invertible diagonal change preserves monomial rank. The observed second moment sum in u-coordinates is diag(82944,27648,…,27648), with all cross entries zero. Thus the stated physical moments follow. The negative of each shell has identical properties.

**Proof.** Use co-moving logical coordinates for the three coset translations d_i (physical coordinates before applying the common deformation); fix d_0=0 to remove all common translations. Represent the deformation modulo rotations by the Gram I+H, H symmetric. Center density at fixed minimum is inversely proportional to sqrt(det(I+H)). If strict local optimality failed, there would be a nonconstant sequence (H,d₁,d₂)→0 of feasible configurations with det(I+H)≤1. Divide this parameter vector by its Euclidean norm and pass to a convergent unit subsequence (h,e₁,e₂). The determinant expansion gives tr h≤0. For any v∈S_ij the original contact becomes v+d_j−d_i in the new Gram, so first-order feasibility gives
vᵀ h v + 2 v·(e_j−e_i) ≥ 0.
Averaging on any one edge gives tr h/6≥0 by centered isotropy. Thus tr h=0, and every displayed inequality on every edge is equality. On the rank-275 edge all coefficients of that linear-plus-quadratic polynomial must be zero: h=0 and e_j−e_i=0. On each other edge isotropy implies its contact vectors span all of W, forcing e_j−e_i=0 as well. The color triangle is connected and e₀=0, so all e_i=0. This contradicts the unit limit. The argument actually shows positive first-order covolume increase in every nontrivial normalized infinitesimally feasible direction. No uniform finite supercell assertion is used.

**Quantitative corollary for translation only.** For any unit direction z the random variable v·z for v uniform in one edge fiber has mean zero, variance 1/6, and absolute value ≤r=√(11/3). If its smallest value is −b then the bound E[(Y+b)(r−Y)]≥0 implies br≥1/6, so min v·z≤−1/(6r)=−r/22. Moving this single *relative* translation by distance t increases at least one restriction in the wrong direction: the new squared inter-coset minimum is at most μ−(r/11)t+t². Thus any nonzero glide of length <r/11≈0.174 under the fixed metric strictly worsens that edge. This concerns relative translation of *complete* layers, not an individual site.

### Standalone Node verifier
Generates the three standard families of minimal Leech vectors from the extended binary Golay code with length-23 cyclic generator 0xae3. The raw shapes are (±4,±4), signed even-parity octads (±2)^8, and one exceptional coordinate −3s among a sign vector s from a Golay word. All are valid Leech vectors; for the local theorem a subset of actual contacts already suffices. The permutation of shell indices is coprime to 196560, and every vector is considered. Inserting pivots by modular elimination supplies rank. No optimization or floating approximation.

```js
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
```
