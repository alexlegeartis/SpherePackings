# 22D: the exact limit of all rank-two Leech diameter certificates, including winding

Sphere106, 26 November 2026. Let Λ be the Leech lattice, integral even unimodular, of minimum squared norm 4. Let K = Λ ∩ U be any primitive planar section, L=Λ∩U^⊥, and put A=1/covol(K). By unimodularity π_U Λ=K*, both complement sections have the same covolume, and K*/K labels the complete translates of L projected from Λ. On any set of s distinct labels let β be its squared diameter **in the flat torus** U/K: distances of different pairs may use different translations. The unconditional lower bound supplied solely by the ambient minimum four is μ≥4−β for the horizontal squared minimum, since one may select a shortest normal representative in each difference class before using the ambient bound. Consequently its reported center-density guarantee is s A (1−β/4)^11 for β<4, zero otherwise. For β=0 it is a single coset. This note bounds this GUARANTEE; it does not upper-bound actual horizontal distances if a shell is absent, the best deformation of a projected packing, or other 22D packings. The existing construction of [Conway–Sloane, §2, Dimension 22](https://arxiv.org/pdf/math/0207182) realizes the optimum δ0=3/√12·(11/12)^11=0.3325496245977818…, occupied density 0.00245103404412118…; see also [the project table](https://github.com/alexlegeartis/SpherePackings/blob/a3d8d1a49b771107ccf3245e64f54d4a750a196d/LOWER_BOUNDS_1_96.md).

**Theorem.** δ0 is the maximum of this construction certificate over *all* primitive rank-two normal sections, all marked subsets, including triples that wind on the normal torus. Equality requires the √2 A2 normal and three equilateral triangle marks of squared side 1/3. In particular neither alternative short planes nor adding more labels can yield a record without actually exploiting improved horizontal clearances.

Write a=λ₁(K)≥2, b=covol(K)/a, d=√β. Orient a Gauss-reduced basis with first vector (a,0), second (t,b). Then A=1/(ab)≤1/√12 by b≥√3 a/2. The dual K* is a rotation of A K, so λ₁(K*)=aA≥2A.

**1. Unwrapped marks.** If marks have plane representatives of diameter d, then for s=1 the guarantee is at most A; for collinear s≥2 their span is ≥(s−1)2A, and sA≤d. The corresponding bound is max_β √β (1−β/4)^11 = 2/√23 (22/23)^11 < 0.256. Three noncollinear points span a triangle of area ≥A/2 and ≤√3 β/4. Thus sA≤(3√3/2)β; β(1−β/4)^11 has unique maximum at 1/3, giving exactly δ0. Equality means a fundamental equilateral triangle. The resulting A=1/√12 also forces minimum 2 and equality in planar Gauss reduction, hence the indicated normal.

For ≥4 noncollinear points, triangulate the hull using all marked points; at least s−2 triangles, each of area ≥A/2. If s≥6, the planar isodiametric inequality gives sA≤(s/(s−2))(π/2)β≤(3π/4)β. For s=4 the hull has area ≤β/2 (quadrilateral diagonal product/2, or triangle), giving sA≤2β. For s=5, perimeter ≤πd by Cauchy's width formula, and the polygon isoperimetric inequality with at most five sides gives area ≤π²β/(20 tan(π/5)). Using tan(π/5)>7/10 and π<22/7 gives sA≤(3π/4)β. Thus for every such s the density is ≤(π/(2√3)) δ0, a strict 9.31% miss. This common-lift Pick argument adapts the prior 46D analysis; the treatment of the medium-scale winding obstruction below differs.

**2. Large torus diameter.** Choose each representative nearest to a fixed anchor at zero, all in the disk of radius d. The hull has area ≤πβ, perimeter ≤2πd. Pick's formula applied to the full hull lattice points, and boundary point spacing ≥2A, give sA≤πβ+πd/2+A. In a collinear hull the line count sA≤d+A is smaller. For β≥8/9 the upper estimate (πβ+π√β/2+1/√12)(1−β/4)^11 is decreasing up to β=4, since the logarithmic derivative of the first factor is ≤1/β and 11/(4−β)>1/β. At β=8/9 it is less than (176/63+11/7+3/10)(7/9)^11<0.3<δ0.

**3. Small torus diameter; rank-one winding.** Thus d<√8/3. If a>3d, all representatives relative to the anchor have ordinary pairwise distances ≤d: comparison to a shortest difference of each pair would otherwise exhibit a nonzero period of length ≤3d. Since K is even, the remaining a² can only be 4 or 6. If b>3d, unwrap the vertical coordinates: each may be chosen within d of the anchor modulo b; every pair then has plain vertical difference at most d (any other multiple of b would have magnitude >d). No short pair can now use the second period. The horizontal coordinates, taken modulo a, are in columns j/a with j mod a². For a=2 the possible cyclic column differences are only 0 or ±1 (since ad<2); any set uses at most two adjacent columns in C4. Taking adjacent representatives gives a common lift. For a=√6, cyclic column differences are at most 2. A subset of Z/6 of this kind occupies at most three columns. It fits in three consecutive positions (a common lift) unless precisely {0,2,4} or its translate. In this exception β≥2/3; vertical differences between two distinct columns are at most e=√(β−2/3). If Li is the vertical width in column i then Li+Lj≤2e, hence their summed width ≤3e. Vertical values within a column have spacing ≥1/b from pairing with the second basis vector. Thus sA≤(3/a)(e+1/b). Here e<1/2, b>3d≥√6>2, and so δ≤(3/√6)(5/6)^11<0.17. No rank-one winding can win.

**4. Finite two-period winding.** In what remains a²≤6 and b≤3d<√8, implying D=(ab)²<48. Every possible even reduced binary normal Gram is [[u,v],[v,w]], u=4 or 6, w≥u even, 0≤v≤u/2, D=uw−v²<48; this gives 21 forms, regardless of realizability in Λ. The snippet below computes K*/K=Z²/G Z² with labels adj(G)p mod D. A dual vector of dot coordinates p=(x,y) has squared norm (w x²−2vxy+u y²)/D. All required shortest distances <8/9 are enumerated: at that threshold Cauchy's inequality gives |x|≤√(8u/9) and |y|≤√(8w/9), far inside the finite box. Every label is assigned; costs over threshold may be overestimated without any consequence. At each distance q/D below 8/9, translate a maximum clique of the compatible-difference graph to contain zero and exhaust all cliques containing zero. The exact comparison for clique size s is 12s²[3(4D−q)]²² ≤ 9D(11D)²². Across all 21 forms the only equality is G=(4,2,4), q=4, clique size 3. There cannot be a winding equality at that value: 3d=√3<2. A single label was already handled in §1.

### Standalone JS verifier
Run under Node as a script, replacing the surrounding assignment by console.log if desired. It uses only integer addition, multiplication, comparisons, and BigInt in the decisive comparison; floating point in the diagnostic display of center densities does not participate in the proof.

```js
const forms=(function main(){
let out=[];
function solve(a,b,c){
 let D=a*c-b*b; const norm=(x,y)=>c*x*x-2*b*x*y+a*y*y;
 let E=new Map();
 for(let x=-2*(a+b+c);x<=2*(a+b+c);x++) for(let y=-2*(a+b+c);y<=2*(a+b+c);y++){
   let s=(p)=>((p%D)+D)%D; let key=s(c*x-b*y)+D*s(a*y-b*x);
   if(!E.has(key)||E.get(key)>norm(x,y))E.set(key,norm(x,y))
 }
 if(E.size!==D)throw Error('group order');
 let labels=[...E.keys()].sort((a,b)=>a-b),F=labels.map(l=>E.get(l)), lev=[...new Set(F)].filter(q=>0<q&&9*q<8*D);
 let s=x=>((x%D)+D)%D,sub=(x,y)=>s(x%D-y%D)+D*s(Math.floor(x/D)-Math.floor(y/D));
 let best=0,win=null,levels=[];
 for(let q of lev){
  let adj=labels.map(x=>labels.reduce((m,y,j)=>x!==y&&E.get(sub(x,y))<=q?m|1n<<BigInt(j):m,0n));
  let max=1,idx=[];
  function dfs(P,ls){
    if(!P){if(ls.length>max){max=ls.length;idx=ls;}return}
    while(P){
      let v=P&-P,j=v.toString(2).length-1;P^=v;
      dfs(P&adj[j],[...ls,j]);
    }
  }
  dfs(adj[0],[0]);
  let mu=4-q/D, delta=max/Math.sqrt(D)*(mu/4)**11;
  let lhs=12n*BigInt(max)**2n*(3n*BigInt(4*D-q))**22n,rhs=9n*BigInt(D)*(11n*BigInt(D))**22n;
  levels.push([q,max,Math.sign(Number(lhs-rhs))]);
  if(delta>best){best=delta;win=[q,max,idx]}
 }
 return {a,b,c,D,best,win,levels};
}
for(let a of [4,6])for(let b=0;b<=a/2;b++)for(let c=a;a*c-b*b<48;c+=2)out.push(solve(a,b,c));
return out;})();
if(forms.length!==21)throw Error('forms');
const equal=[];
for(const f of forms)for(const [q,s,relation] of f.levels){
 if(relation>0)throw Error('counterexample');
 if(relation===0)equal.push([f.a,f.b,f.c,q,s]);
}
if(JSON.stringify(equal)!=='[[4,2,4,4,3]]')throw Error('equality cases');
console.log(forms.length,equal);
```
