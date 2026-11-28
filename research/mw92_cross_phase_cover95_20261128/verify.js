#!/usr/bin/env node
'use strict';
const fail=(x,msg)=>{if(!x)throw Error(msg||'assertion');};
const ga=(a,b)=>b===0n?(a<0n?-a:a):ga(b,a%b);
const F=(n,d=1n)=>{n=BigInt(n);d=BigInt(d);if(d===0n)throw Error('division by zero');if(d<0n){n=-n;d=-d;}let g=ga(n,d);return [n/g,d/g];};
const A=(a,b)=>F(a[0]*b[1]+b[0]*a[1],a[1]*b[1]), M=(a,b)=>F(a[0]*b[0],a[1]*b[1]), S=(a,b)=>A(a,F(-b[0],b[1])), V=(a,b)=>F(a[0]*b[1],a[1]*b[0]), P=(a,k)=>F(a[0]**BigInt(k),a[1]**BigInt(k));
const L=(a,b)=>a[0]*b[1]<b[0]*a[1], E=(a,b)=>a[0]*b[1]===b[0]*a[1], Z=F(0), I=F(1), T=F(2), N=F(16);
const H=(a)=>V(a,T), m=(a,k)=>M(a,F(k)), q=(a)=>P(a,2), pj=a=>a.map(String).join('/');
const ca=F(1100098526n,100000000n), cb=F(1078402070n,100000000n), eps=F(1,100000), s=S(N,ca), p=F(14,3), t=S(N,cb);
const Dg=F(1055366,100000), Df=F(1055460,100000);
const generic=[[s,H(s),m(H(s),-1)],[H(s),s,m(H(s),-1)],[m(H(s),-1),m(H(s),-1),S(N,Dg)]];
const fixed=[[p,H(p),H(p)],[H(p),t,H(t)],[H(p),H(t),S(N,Df)]];
const C=(G,i)=>{let [j,k]=[0,1,2].filter(l=>l!==i);return S(M(G[j][j],G[k][k]),q(G[j][k]));};
const det=(G)=>{let [a,b,c]=G[0],d=G[1][1],e=G[1][2],f=G[2][2];return [M(M(a,d),f),m(M(M(b,c),e),2),m(M(f,q(b)),-1),m(M(d,q(c)),-1),m(M(a,q(e)),-1)].reduce(A,Z);};
const evalQ=(G,v)=>{let r=Z;for(let i=0;i<3;i++)for(let j=0;j<3;j++)r=A(r,m(G[i][j],v[i]*v[j]));return r;};
const gStrong=[[0,0],[0,1],[1,0]],gMid=[[-1,1],[1,-1],[1,1]],gWeak=[[-1,0],[-1,2],[0,-1],[0,2],[2,-1],[2,0]];
const gMatches=[[0,0],[1,2],[2,3]];
const fStrong=[[0,-1],[0,0],[2,0]],fMid=[[2,1],[2,-1],[1,-1]],fTiny=[[0,1],[1,1],[0,-2],[1,-2]],fDim=[[1,0]];
const fMatches=[[[0,0],[2,3]],[[0,0],[1,2]],[[0,0],[1,3]]];
const n3=(a)=>((a%3)+3)%3;
const k=(a)=>a.join(',');
const gReq=new Map([...gStrong.map(v=>[k(v),Dg]),...gMid.map(v=>[k(v),S(Dg,s)]),...gWeak.map(v=>[k(v),S(Dg,m(s,2))])]);
const fReq=new Map([...fStrong.map(v=>[k(v),Df]),[k(fMid[0]),S(Df,S(m(t,2),p))],...[fMid[1],fMid[2]].map(v=>[k(v),S(Df,p)]),...fTiny.map(v=>[k(v),S(Df,m(t,2))]),[k(fDim[0]),S(Df,m(p,2))]]);
function audit(G,isFix,expected){
 let d=det(G); fail(L(Z,G[0][0])&&L(Z,C(G,2))&&L(Z,d),'PD');
 for(let i=0;i<3;i++)fail(L(m(C(G,i),16),m(d,9)),'finite coordinate bounds'); // |z_l|<3 for Q<16
 let seen=new Map(), seen0=[];
 for(let i=-2;i<=2;i++)for(let j=-2;j<=2;j++)for(let l=0;l<=2;l++){
   let Q=evalQ(G,[i,j,l]);if(!L(Q,N))continue;
   if(l===2)throw Error('unexpected second layer');
   if(l===0){if(i||j){
      if(!isFix){let r=i*i+i*j+j*j; fail(r===1||r===3,'generic Eisenstein'); fail(E(Q,m(s,r)),'scale');}
      else{
        if(j===0)fail(Math.abs(i)===1&&E(Q,p),'fixed first layers');
        else {fail(Math.abs(j)===1);let ii=n3(i*(j>0?1:-1));let low=ii===1?S(cb,m(p,2)):cb;fail(!L(A(Q,low),N),'fixed j=1 or −1');}
      }
      seen0.push([i,j]);
   }continue;}
   let a=k([isFix?n3(i):i,j]), u=S(N,Q),old=seen.get(a);
   if(!old||L(old,u))seen.set(a,u);
 }
 fail(seen.size===expected.size&&[...seen].every(([k,v])=>expected.has(k)&&E(v,expected.get(k))),'exact list of forbidden phases');
 fail(seen0.length>0);
 return d;
}
let dg=audit(generic,false,gReq),df=audit(fixed,true,fReq);
for(let [i,j] of gMatches){let a=gMid[i],b=gWeak[j],x=a[0]-b[0],y=a[1]-b[1];fail(x*x+x*y+y*y===1,'generic match');}
for(let i=0;i<fMatches.length;i++){
 const [[m0,w0],[mj,wj]]=fMatches[i];
 let d0=[n3(fTiny[w0][0]-fMid[m0][0]),fTiny[w0][1]-fMid[m0][1]];
 let z=[n3(fTiny[wj][0]-fMid[mj][0]),fTiny[wj][1]-fMid[mj][1]];
 fail(d0[1]===0&&d0[0]!==0,'u match');fail(z[1]===-1&&z[0]===i,'boundary match type');
}
function upper(r,d){fail(L(d,q(r)),'upper square');return r;}
function lower(r,d){fail(L(q(r),d),'lower square');return r;}
const F1000=x=>F(x,1000);
const ag=upper(F1000(3317),A(ca,eps)), agmid=upper(F1000(2357),S(Dg,s)), agweak=upper(F1000(746),S(Dg,m(s,2))), Rg=lower(F1000(3248),Dg);
fail(L(A(A(ag,agmid),agweak),m(Rg,2)), 'generic unlike-radius merge');
fail(L(A(upper(F1000(2450),F(6)),agweak),Rg),'generic tetrahedra');
const ut=upper(F1000(579),V(S(N,m(s,3)),F(3))), uup=A(ag,ut);
const Rf=lower(F1000(3248),Df), f0=upper(F1000(2189),fReq.get(k(fMid[0]))), f1=upper(F1000(2427),fReq.get(k(fMid[1]))), fw=upper(F1000(351),fReq.get(k(fTiny[0]))), fb=upper(F1000(3284),A(cb,eps)), fr=upper(F1000(1106),fReq.get(k(fDim[0])));
fail(L(A(A(uup,f0),fw),m(Rf,2)),'fixed same-j merge');
fail(L(A(A(fb,f1),fw),m(Rf,2)),'fixed boundary merge');
fail(L(A(T,fr),Rf),'fixed dimers');fail(L(A(F(8,3),fw),Rf),'fixed nine-simplex');
const Ng=F(27,4),Nf=F(121,18);
fail(L(M(Ng,P(Dg,46)),P(ca,46))&&L(M(Nf,P(Df,46)),P(ca,46)),'covering densities');
const old=F(14260143213204819n,200000000000000n);
fail(L(dg,df)&&L(df,old)&&L(M(q(F(100099,100000)),df),old),'strict 0.099% gain');
const oldFloor=F(2133434801047687n,10n**36n),newFloor=F(213554696184854n,10n**35n);
fail(L(M(q(newFloor),df),M(q(oldFloor),old)),'published decimal floor');
console.log('PASS: complete tail, exact unlike-radius merging, strict covering contradiction.');
console.log('Generic determinant:',pj(dg));console.log('Fixed determinant:',pj(df));
console.log('95D occupied density >= 2.13554696184854e-21, improvement > 0.099%.');
